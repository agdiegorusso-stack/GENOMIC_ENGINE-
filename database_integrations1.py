#!/usr/bin/env python3
"""
Integrazione completa con database di emoglobinopatie
HbVar, ClinVar, LOVD, gnomAD, IthaGenes
Versione corretta con gestione errori migliorata e URL HbVar aggiornato
"""

import os
import re
import json
import time
import sqlite3
import hashlib
import logging
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
from contextlib import contextmanager
import xml.etree.ElementTree as ET
from urllib.parse import quote

import pandas as pd
import numpy as np
import requests
from bs4 import BeautifulSoup
from Bio import Entrez
import configparser

from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager
import time

# Configurazione logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class Configuration:
    """Gestione centralizzata della configurazione"""
    
    def __init__(self, config_file='config.ini'):
        self.config = configparser.ConfigParser()
        if os.path.exists(config_file):
            self.config.read(config_file)
        else:
            self._create_default_config(config_file)
    
    def _create_default_config(self, config_file):
        """Crea file di configurazione di default"""
        self.config['NCBI'] = {
            'email': 'ag.diegorusso@gmail.com',
            'api_key': ''  # Opzionale ma raccomandato
        }
        self.config['DATABASE'] = {
            'cache_dir': 'cache',
            'output_dir': 'integrated_db'
        }
        self.config['API'] = {
            'timeout': '30',
            'retry_attempts': '3',
            'rate_limit_delay': '1.0'
        }
        
        with open(config_file, 'w') as f:
            self.config.write(f)
        
        logger.warning(f"Created default config file: {config_file}")
        logger.warning("Please update the email address in the config file")
    
    def get_email(self):
        return self.config.get('NCBI', 'email')
    
    def get_api_key(self):
        return self.config.get('NCBI', 'api_key', fallback=None)


class DatabaseConnection:
    """Context manager per connessioni database sicure"""
    
    @staticmethod
    @contextmanager
    def get_connection(db_path: str):
        """Context manager per connessioni SQLite"""
        conn = None
        try:
            conn = sqlite3.connect(db_path)
            conn.row_factory = sqlite3.Row
            yield conn
        except sqlite3.Error as e:
            if conn:
                conn.rollback()
            logger.error(f"Database error: {e}")
            raise
        finally:
            if conn:
                conn.close()


class HbVarDatabase:
    """
    Interfaccia per HbVar - Versione 9.0 (Soluzione a due passaggi)
    Simula l'intero flusso di navigazione:
    1. Invia il form di ricerca.
    2. Nella pagina di cronologia, clicca 'Go' per visualizzare i risultati.
    3. Analizza la tabella finale.
    """
    
    # ... Il CHAIN_MAPPING e il costruttore __init__ non cambiano ...
    CHAIN_MAPPING = {
        'α': 'HBA', 'alpha': 'HBA', 'a': 'HBA', 'β': 'HBB', 'beta': 'HBB', 'b': 'HBB',
        'γ': 'HBG', 'gamma': 'HBG', 'g': 'HBG', 'δ': 'HBD', 'delta': 'HBD', 'd': 'HBD',
        'ε': 'HBE', 'epsilon': 'HBE', 'e': 'HBE', 'ζ': 'HBZ', 'zeta': 'HBZ', 'z': 'HBZ'
    }
    
    def __init__(self, cache_dir="hbvar_cache", config=None):
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        self.query_url = "https://globin.bx.psu.edu/cgi-bin/hbvar/query_vars3"
        self.db_path = os.path.join(self.cache_dir, "hbvar.db")
        self.config = config or Configuration()
        self._initialize_database()

    def _initialize_database(self):
        with DatabaseConnection.get_connection(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS variants (
                    variant_id TEXT PRIMARY KEY, hb_name TEXT NOT NULL, gene TEXT,
                    mutation_protein TEXT, position INTEGER, wild_type_aa TEXT, 
                    mutant_aa TEXT, phenotype TEXT, hgvs TEXT, data_source TEXT DEFAULT 'HbVar',
                    last_updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
            conn.commit()

    def download_all_variants(self, force_update=False) -> List[Dict]:
        logger.info("Starting HbVar download with 2-step navigation...")
        if not force_update and self._has_recent_data():
            logger.info("Using cached HbVar data")
            return self._get_cached_variants()

        all_variants = self._fetch_with_selenium()
        if all_variants:
            self._save_variants_to_db(all_variants)
        
        logger.info(f"Total variants downloaded from HbVar: {len(all_variants)}")
        return all_variants

    def _fetch_with_selenium(self) -> List[Dict]:
        """
        Usa Selenium per automatizzare la navigazione basandosi sul flusso corretto a 2 passaggi.
        Adesso cerca il testo corretto nella pagina di risultato.
        """
        variants = []
        driver = None
        try:
            logger.info("Setting up headless browser with Selenium...")
            chrome_options = Options()
            chrome_options.add_argument("--headless")
            chrome_options.add_argument("--no-sandbox")
            chrome_options.add_argument("--disable-dev-shm-usage")
            
            service = Service(ChromeDriverManager().install())
            driver = webdriver.Chrome(service=service, options=chrome_options)
            driver.set_page_load_timeout(120)

            # PASSO 1: Naviga al form e invialo
            logger.info(f"Navigating to {self.query_url}...")
            driver.get(self.query_url)
            WebDriverWait(driver, 20).until(
                EC.element_to_be_clickable((By.XPATH, '//form[1]//input[@type="submit" and @value="Submit Query"]'))
            )
            logger.info("Submitting the initial search form...")
            submit_button = driver.find_element(By.XPATH, '//form[1]//input[@type="submit" and @value="Submit Query"]')
            submit_button.click()
            
            # PASSO 2: Dalla pagina di cronologia, clicca 'Go'
            wait = WebDriverWait(driver, 30)
            logger.info("On history page. Waiting for 'Go' button...")
            go_button = wait.until(
                EC.element_to_be_clickable((By.XPATH, '//input[@type="submit" and @value="Go"]'))
            )
            logger.info("Clicking 'Go' to display results...")
            go_button.click()

            # ***** INIZIO CORREZIONE FINALE *****
            # L'attesa ora cerca "matches" con la 'm' minuscola, come nell'HTML reale.
            wait_time = 300
            logger.info(f"Waiting for final results table to load (max {wait_time} seconds)...")
            wait = WebDriverWait(driver, wait_time)
            
            wait.until(EC.presence_of_element_located((By.XPATH, "//*[contains(text(), 'matches to your query')]")))
            # ***** FINE CORREZIONE FINALE *****
            
            logger.info("Final results page loaded. Extracting data...")
            html_content = driver.page_source
            variants = self._parse_hbvar_page(html_content)

        except Exception as e:
            logger.error(f"An error occurred during Selenium execution: {e}", exc_info=True)
            if driver:
                driver.save_screenshot('selenium_error_screenshot.png')
                with open("selenium_error_page.html", "w", encoding="utf-8") as f:
                    f.write(driver.page_source)
                logger.info("Saved a screenshot and HTML of the page state at the time of error.")
        finally:
            if driver:
                driver.quit()
        
        return variants

    # ... Il resto del codice (da _parse_hbvar_page in poi) è identico a prima e corretto ...
    def _parse_hbvar_page(self, html_content: str) -> List[Dict]:
        soup = BeautifulSoup(html_content, 'html.parser')
        main_table = None
        results_header = soup.find(string=re.compile(r'\d+ matches to your query'))
        
        if results_header:
             main_table = results_header.find_next('table')
             logger.info("Found data table using the 'matches to your query' text.")
        
        if not main_table:
            with open("hbvar_parsing_ERROR.html", "w", encoding="utf-8") as f: f.write(soup.prettify())
            logger.error("CRITICAL: No data table found on the results page. Saved the response to 'hbvar_parsing_ERROR.html'.")
            return []

        variants, header_map = [], {'name': 'hb_name', 'mutation': 'mutation', 'mutation, hgvs nomenclature': 'hgvs'}
        final_headers = [h.get_text(strip=True).lower() for h in main_table.find_all('tr')[0].find_all(['td', 'th'])]
        
        for row in main_table.find_all('tr')[1:]:
            cells = row.find_all('td')
            if len(cells) != len(final_headers): continue
            variant_data = {header_map[h]: c.get_text(strip=True) for i, (h, c) in enumerate(zip(final_headers, cells)) if h in header_map}
            if 'mutation' in variant_data:
                variant_data.update(self._parse_mutation_notation(variant_data['mutation']))
            if variant_data.get('hb_name'): variants.append(variant_data)
        
        logger.info(f"Successfully parsed {len(variants)} variants from page source.")
        return variants

    def _parse_mutation_notation(self, mutation_str: str) -> Dict:
        info = {'gene': '', 'position': None, 'wild_type_aa': '', 'mutant_aa': '', 'mutation_protein': mutation_str}
        normalized = mutation_str.replace('→', '>').replace(' to ', '>').replace('-->', '>')
        patterns = [
            r'([αβγδεζa-z\s]+)\s*(\d+)\s*\(?[\w\d]+\)?\s*(\w{3})\s*>\s*(\w{3})',
            r'([αβγδεζa-z\s]+)\s*nt\s*-?(\d+)\s*(\w)\s*>\s*(\w)'
        ]
        for pattern in patterns:
            match = re.search(pattern, normalized, re.IGNORECASE)
            if match:
                groups = match.groups()
                chain_str = groups[0].lower().strip()
                for key, val in self.CHAIN_MAPPING.items():
                    if key in chain_str:
                        info['gene'] = val; break
                info['position'] = int(groups[1])
                info['wild_type_aa'], info['mutant_aa'] = self._normalize_amino_acid(groups[2]), self._normalize_amino_acid(groups[3])
                return info
        return info

    def _normalize_amino_acid(self, aa: str) -> str:
        aa_map = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val', 'X': 'Ter', '*': 'Ter', 'Stop': 'Ter'}
        return aa_map.get(aa.upper(), aa.upper()) if len(aa) == 1 else aa.capitalize()

    def _save_variants_to_db(self, variants: List[Dict]):
        with DatabaseConnection.get_connection(self.db_path) as conn:
            cursor = conn.cursor()
            try:
                for variant in variants:
                    unique_str = f"{variant.get('hb_name', '')}_{variant.get('mutation', '')}"
                    variant_id = hashlib.sha256(unique_str.encode()).hexdigest()[:32]
                    cursor.execute("""
                        INSERT OR REPLACE INTO variants (
                            variant_id, hb_name, gene, mutation_protein, position, 
                            wild_type_aa, mutant_aa, phenotype, hgvs, last_updated
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """, (
                        variant_id, variant.get('hb_name'), variant.get('gene'),
                        variant.get('mutation_protein'), variant.get('position'),
                        variant.get('wild_type_aa'), variant.get('mutant_aa'),
                        variant.get('phenotype'), variant.get('hgvs'), datetime.now()
                    ))
                conn.commit()
                logger.info(f"Saved/Updated {len(variants)} variants in the local database.")
            except sqlite3.Error as e:
                conn.rollback(); logger.error(f"Database save error: {e}"); raise
    
    def _has_recent_data(self, days=7) -> bool:
        if not os.path.exists(self.db_path): return False
        with DatabaseConnection.get_connection(self.db_path) as conn:
            cursor = conn.cursor()
            try:
                cursor.execute("SELECT MAX(last_updated) as last_update FROM variants")
                result = cursor.fetchone()
                if result and result['last_update']:
                    last_update = datetime.fromisoformat(result['last_update'])
                    return (datetime.now() - last_update).days < days
            except sqlite3.OperationalError: return False
        return False
    
    def _get_cached_variants(self) -> List[Dict]:
        with DatabaseConnection.get_connection(self.db_path) as conn:
            return pd.read_sql_query("SELECT * FROM variants", conn).to_dict('records')

    def query_variant(self, gene: str = None, position: int = None, hb_name: str = None) -> pd.DataFrame:
        query = "SELECT * FROM variants WHERE 1=1"
        params = []
        if gene: query += " AND gene = ?"; params.append(gene)
        if position: query += " AND position = ?"; params.append(position)
        if hb_name: query += " AND hb_name LIKE ?"; params.append(f"%{hb_name}%")
        with DatabaseConnection.get_connection(self.db_path) as conn:
            return pd.read_sql_query(query, conn, params=params)

# ... (The rest of the script, including ClinVarIntegrator, GnomADIntegrator, etc., remains unchanged)
# The classes below this point are assumed to be correct and are not repeated for brevity.
class ClinVarIntegrator:
    """
    Integrazione migliorata con ClinVar per varianti di emoglobina
    """
    
    def __init__(self, config=None):
        self.config = config or Configuration()
        Entrez.email = self.config.get_email()
        
        api_key = self.config.get_api_key()
        if api_key:
            Entrez.api_key = api_key
        
        self.cache_dir = "clinvar_cache"
        os.makedirs(self.cache_dir, exist_ok=True)
        self.cache = {}
    
    def search_hemoglobin_variants(self, batch_size: int = 500, 
                                  clinical_significance: str = None) -> pd.DataFrame:
        """
        Cerca varianti di emoglobina in ClinVar con filtri avanzati
        """
        logger.info("Searching ClinVar for hemoglobin variants...")
        
        # Geni di emoglobina completi
        hemoglobin_genes = [
            'HBB', 'HBA1', 'HBA2',  # Principali
            'HBG1', 'HBG2',  # Gamma (fetale)
            'HBD',  # Delta
            'HBE1',  # Epsilon
            'HBZ', 'HBM',  # Zeta e Mu
            'HBQ1'  # Theta
        ]
        
        all_variants = []
        
        for gene in hemoglobin_genes:
            logger.info(f"  Searching {gene}...")
            
            # Costruisci query con filtri
            query_parts = [f"{gene}[gene]"]
            
            if clinical_significance:
                query_parts.append(f"{clinical_significance}[Clinical_Significance]")
            
            query = " AND ".join(query_parts)
            
            try:
                variants = self._search_clinvar_batch(query, batch_size)
                all_variants.extend(variants)
                logger.info(f"    Found {len(variants)} variants")
                
                # Rate limiting
                time.sleep(float(self.config.config.get('API', 'rate_limit_delay')))
                
            except Exception as e:
                logger.error(f"    Error searching {gene}: {str(e)}")
        
        df = pd.DataFrame(all_variants)
        
        # Post-processing
        if not df.empty:
            df = self._process_clinvar_dataframe(df)
        
        logger.info(f"Total ClinVar variants found: {len(df)}")
        return df
    
    def _search_clinvar_batch(self, query: str, batch_size: int) -> List[Dict]:
        """Ricerca batch con gestione errori"""
        variants = []
        
        try:
            # Search
            with Entrez.esearch(db="clinvar", term=query, retmax=batch_size) as handle:
                result = Entrez.read(handle)
            
            if not result.get('IdList'):
                return variants
            
            # Fetch in batch più piccoli per evitare timeout
            id_list = result['IdList']
            chunk_size = 50
            
            for i in range(0, len(id_list), chunk_size):
                chunk_ids = id_list[i:i+chunk_size]
                
                try:
                    with Entrez.efetch(
                        db="clinvar",
                        id=','.join(chunk_ids),
                        rettype="vcv",
                        retmode="xml"
                    ) as handle:
                        xml_data = handle.read()
                    
                    chunk_variants = self._parse_clinvar_xml(xml_data)
                    variants.extend(chunk_variants)
                    
                except Exception as e:
                    logger.warning(f"Error fetching chunk: {e}")
                
                time.sleep(0.5)  # Rate limiting tra chunk
        
        except Exception as e:
            logger.error(f"ClinVar search error: {e}")
        
        return variants
    
    def _parse_clinvar_xml(self, xml_data: str) -> List[Dict]:
        """
        Parse XML da ClinVar con gestione errori robusta
        """
        variants = []
        
        try:
            root = ET.fromstring(xml_data)
            
            # Namespace handling
            ns = {'cv': 'https://www.ncbi.nlm.nih.gov/clinvar/'}
            
            for record in root.findall('.//VariationArchive'):
                try:
                    variant = self._extract_variant_from_record(record, ns)
                    if variant:
                        variants.append(variant)
                except Exception as e:
                    logger.warning(f"Error parsing record: {e}")
                    continue
        
        except ET.ParseError as e:
            logger.error(f"XML parsing error: {e}")
        
        return variants
    
    def _extract_variant_from_record(self, record: ET.Element, ns: Dict) -> Optional[Dict]:
        """Estrae informazioni da un record ClinVar"""
        variant = {
            'variant_id': record.get('VariationID', ''),
            'name': record.get('VariationName', ''),
            'gene': '',
            'clinical_significance': '',
            'review_status': '',
            'conditions': [],
            'hgvs_expressions': [],
            'molecular_consequence': '',
            'protein_change': '',
            'submission_count': 0
        }
        
        # Gene
        gene_elem = record.find('.//Gene')
        if gene_elem is not None:
            variant['gene'] = gene_elem.get('Symbol', '')
        
        # Clinical significance
        clin_sig = record.find('.//ClinicalSignificance/Description')
        if clin_sig is not None:
            variant['clinical_significance'] = clin_sig.text or ''
        
        # Review status
        review = record.find('.//ClinicalSignificance/ReviewStatus')
        if review is not None:
            variant['review_status'] = review.text or ''
        
        # Conditions/Phenotypes
        for trait in record.findall('.//TraitSet/Trait'):
            name_elem = trait.find('.//Name/ElementValue')
            if name_elem is not None and name_elem.text:
                variant['conditions'].append(name_elem.text)
        
        # HGVS expressions
        for hgvs in record.findall('.//HGVSExpression'):
            expr = hgvs.get('change', '')
            if expr:
                variant['hgvs_expressions'].append(expr)
                
                # Estrai protein change se presente
                if 'p.' in expr:
                    variant['protein_change'] = expr.split('p.')[-1]
        
        # Molecular consequence
        mol_cons = record.find('.//MolecularConsequence')
        if mol_cons is not None:
            variant['molecular_consequence'] = mol_cons.get('Type', '')
        
        # Conta submissions
        submissions = record.findall('.//ClinicalAssertion')
        variant['submission_count'] = len(submissions)
        
        return variant
    
    def _process_clinvar_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """Post-processing del dataframe ClinVar"""
        # Converti liste in stringhe
        if 'conditions' in df.columns:
            df['conditions_str'] = df['conditions'].apply(
                lambda x: '; '.join(x) if isinstance(x, list) else str(x)
            )
        
        if 'hgvs_expressions' in df.columns:
            df['hgvs_str'] = df['hgvs_expressions'].apply(
                lambda x: '; '.join(x) if isinstance(x, list) else str(x)
            )
        
        # Categorizza significance
        significance_categories = {
            'Pathogenic': 'pathogenic',
            'Likely pathogenic': 'pathogenic',
            'Benign': 'benign',
            'Likely benign': 'benign',
            'Uncertain significance': 'uncertain',
            'Conflicting interpretations': 'conflicting'
        }
        
        df['significance_category'] = df['clinical_significance'].map(
            lambda x: significance_categories.get(x, 'unknown')
        )
        
        # Aggiungi score di affidabilità basato su review status
        review_scores = {
            'reviewed by expert panel': 4,
            'criteria provided, multiple submitters, no conflicts': 3,
            'criteria provided, single submitter': 2,
            'no assertion criteria provided': 1,
            'no assertion provided': 0
        }
        
        df['confidence_score'] = df['review_status'].map(
            lambda x: review_scores.get(x, 0)
        )
        
        return df


class GnomADIntegrator:
    """
    Integrazione con gnomAD v3 per frequenze popolazionali
    """
    
    def __init__(self, config=None):
        self.config = config or Configuration()
        self.base_url = "https://gnomad.broadinstitute.org/api"
        self.cache_dir = "gnomad_cache"
        os.makedirs(self.cache_dir, exist_ok=True)
    
    def get_variant_by_position(self, chrom: str, pos: int, ref: str, alt: str) -> Dict:
        """
        Query diretta per variante specifica
        """
        # Costruisci variant ID in formato gnomAD
        variant_id = f"{chrom}-{pos}-{ref}-{alt}"
        
        # Controlla cache
        cache_file = os.path.join(self.cache_dir, f"{variant_id}.json")
        if os.path.exists(cache_file):
            with open(cache_file, 'r') as f:
                return json.load(f)
        
        # Query GraphQL
        query = """
        query VariantInfo($variantId: String!) {
            variant(variantId: $variantId, dataset: gnomad_r3) {
                variant_id
                chrom
                pos
                ref
                alt
                exome {
                    ac
                    an
                    af
                    filters
                    populations {
                        id
                        ac
                        an
                        af
                    }
                }
                genome {
                    ac
                    an
                    af
                    filters
                    populations {
                        id
                        ac
                        an
                        af
                    }
                }
                transcript_consequences {
                    gene_symbol
                    transcript_id
                    consequence_terms
                    hgvsc
                    hgvsp
                    lof
                    polyphen_prediction
                    sift_prediction
                }
            }
        }
        """
        
        try:
            response = requests.post(
                self.base_url,
                json={'query': query, 'variables': {'variantId': variant_id}},
                timeout=30
            )
            
            if response.status_code == 200:
                data = response.json()
                
                if data.get('data', {}).get('variant'):
                    variant_data = data['data']['variant']
                    
                    # Cache result
                    with open(cache_file, 'w') as f:
                        json.dump(variant_data, f)
                    
                    return variant_data
        
        except Exception as e:
            logger.error(f"gnomAD query error: {e}")
        
        return {}
    
    def get_hemoglobin_gene_variants(self, gene_symbol: str) -> pd.DataFrame:
        """
        Ottiene tutte le varianti per un gene di emoglobina
        """
        logger.info(f"Fetching gnomAD data for {gene_symbol}...")
        
        # Mappatura coordinate genomiche (GRCh38)
        gene_coordinates = {
            'HBB': {'chrom': '11', 'start': 5225464, 'end': 5227071},
            'HBA1': {'chrom': '16', 'start': 176680, 'end': 177522},
            'HBA2': {'chrom': '16', 'start': 172876, 'end': 173710},
            'HBG1': {'chrom': '11', 'start': 5269755, 'end': 5271400},
            'HBG2': {'chrom': '11', 'start': 5274620, 'end': 5276265},
            'HBD': {'chrom': '11', 'start': 5233577, 'end': 5234861},
        }
        
        if gene_symbol not in gene_coordinates:
            logger.warning(f"Gene {gene_symbol} not in coordinate mapping")
            return pd.DataFrame()
        
        coords = gene_coordinates[gene_symbol]
        
        # Query per regione
        query = """
        query RegionVariants($chrom: String!, $start: Int!, $stop: Int!) {
            region(chrom: $chrom, start: $start, stop: $stop) {
                variants(dataset: gnomad_r3) {
                    variant_id
                    pos
                    ref
                    alt
                    exome {
                        af
                    }
                    genome {
                        af
                    }
                    transcript_consequences {
                        gene_symbol
                        consequence_terms
                        hgvsp
                    }
                }
            }
        }
        """
        
        try:
            response = requests.post(
                self.base_url,
                json={
                    'query': query,
                    'variables': {
                        'chrom': coords['chrom'],
                        'start': coords['start'],
                        'stop': coords['end']
                    }
                },
                timeout=60
            )
            
            if response.status_code == 200:
                data = response.json()
                
                if data.get('data', {}).get('region', {}).get('variants'):
                    variants = data['data']['region']['variants']
                    
                    # Converti in DataFrame
                    df = pd.json_normalize(variants)
                    
                    # Filtra per gene specifico
                    if 'transcript_consequences' in df.columns:
                        df = df[df['transcript_consequences'].apply(
                            lambda x: any(c.get('gene_symbol') == gene_symbol 
                                        for c in x if isinstance(c, dict))
                        )]
                    
                    logger.info(f"  Found {len(df)} variants")
                    return df
        
        except Exception as e:
            logger.error(f"gnomAD region query error: {e}")
        
        return pd.DataFrame()
    
    def annotate_with_population_frequencies(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """
        Annota un dataframe di varianti con frequenze popolazionali
        """
        if variants_df.empty:
            return variants_df
        
        logger.info("Annotating variants with population frequencies...")
        
        # Colonne richieste
        required_cols = ['chrom', 'pos', 'ref', 'alt']
        
        # Verifica colonne
        missing_cols = [col for col in required_cols if col not in variants_df.columns]
        if missing_cols:
            logger.warning(f"Missing required columns: {missing_cols}")
            return variants_df
        
        # Aggiungi colonne per frequenze
        pop_columns = ['af_global', 'af_afr', 'af_amr', 'af_eas', 
                      'af_nfe', 'af_sas', 'af_asj', 'af_fin']
        
        for col in pop_columns:
            variants_df[col] = None
        
        # Annota ogni variante
        for idx, row in variants_df.iterrows():
            variant_data = self.get_variant_by_position(
                row['chrom'], row['pos'], row['ref'], row['alt']
            )
            
            if variant_data:
                # Estrai frequenze globali
                if 'genome' in variant_data and variant_data['genome']:
                    variants_df.at[idx, 'af_global'] = variant_data['genome'].get('af')
                    
                    # Frequenze per popolazione
                    if 'populations' in variant_data['genome']:
                        for pop in variant_data['genome']['populations']:
                            pop_id = pop.get('id', '')
                            if f'af_{pop_id}' in pop_columns:
                                variants_df.at[idx, f'af_{pop_id}'] = pop.get('af')
        
        # Categorizza frequenze
        variants_df['frequency_category'] = pd.cut(
            variants_df['af_global'].fillna(0),
            bins=[0, 1e-5, 1e-4, 1e-3, 1e-2, 1],
            labels=['ultra_rare', 'very_rare', 'rare', 'low_frequency', 'common']
        )
        
        logger.info(f"Annotated {len(variants_df)} variants")
        return variants_df


class IntegratedHemoglobinDatabase:
    """
    Database integrato che combina tutte le fonti con architettura migliorata
    """
    
    def __init__(self, output_dir="integrated_db", config=None):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        self.config = config or Configuration()
        
        # Inizializza componenti
        self.hbvar = HbVarDatabase(config=self.config)
        self.clinvar = ClinVarIntegrator(config=self.config)
        self.gnomad = GnomADIntegrator(config=self.config)
        
        # Database path
        self.db_path = os.path.join(self.output_dir, "hemoglobin_integrated.db")
        self._create_integrated_database()
    
    def _create_integrated_database(self):
        """
        Crea schema database integrato ottimizzato
        """
        with DatabaseConnection.get_connection(self.db_path) as conn:
            cursor = conn.cursor()
            
            # Tabella principale varianti con chiave composita
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS variants (
                    variant_id TEXT PRIMARY KEY,
                    gene TEXT NOT NULL,
                    chromosome TEXT,
                    position_grch38 INTEGER,
                    position_grch37 INTEGER,
                    ref_allele TEXT,
                    alt_allele TEXT,
                    variant_type TEXT,
                    hgvs_genomic TEXT,
                    hgvs_coding TEXT,
                    hgvs_protein TEXT,
                    protein_position INTEGER,
                    aa_change TEXT,
                    codon_change TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
            
            # Tabella annotazioni cliniche multi-sorgente
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS clinical_annotations (
                    annotation_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    variant_id TEXT NOT NULL,
                    source TEXT NOT NULL,
                    clinical_significance TEXT,
                    pathogenicity_score REAL,
                    phenotype TEXT,
                    inheritance_mode TEXT,
                    penetrance TEXT,
                    evidence_level TEXT,
                    review_status TEXT,
                    last_evaluated DATE,
                    pmids TEXT,
                    FOREIGN KEY(variant_id) REFERENCES variants(variant_id),
                    UNIQUE(variant_id, source)
                )
            """)
            
            # Tabella frequenze alleliche dettagliate
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS allele_frequencies (
                    frequency_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    variant_id TEXT NOT NULL,
                    dataset TEXT NOT NULL,
                    population TEXT NOT NULL,
                    allele_count INTEGER,
                    allele_number INTEGER,
                    allele_frequency REAL,
                    homozygote_count INTEGER,
                    hemizygote_count INTEGER,
                    coverage_mean REAL,
                    FOREIGN KEY(variant_id) REFERENCES variants(variant_id),
                    UNIQUE(variant_id, dataset, population)
                )
            """)
            
            # Tabella predizioni in silico
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS in_silico_predictions (
                    prediction_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    variant_id TEXT NOT NULL,
                    tool TEXT NOT NULL,
                    score REAL,
                    prediction TEXT,
                    confidence REAL,
                    rankscore REAL,
                    FOREIGN KEY(variant_id) REFERENCES variants(variant_id),
                    UNIQUE(variant_id, tool)
                )
            """)
            
            # Tabella dati funzionali
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS functional_data (
                    functional_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    variant_id TEXT NOT NULL,
                    assay_type TEXT,
                    measurement TEXT,
                    value REAL,
                    unit TEXT,
                    condition TEXT,
                    pmid TEXT,
                    FOREIGN KEY(variant_id) REFERENCES variants(variant_id)
                )
            """)
            
            # Indici ottimizzati
            indices = [
                "CREATE INDEX IF NOT EXISTS idx_var_gene ON variants(gene)",
                "CREATE INDEX IF NOT EXISTS idx_var_position ON variants(chromosome, position_grch38)",
                "CREATE INDEX IF NOT EXISTS idx_var_protein ON variants(gene, protein_position)",
                "CREATE INDEX IF NOT EXISTS idx_clin_significance ON clinical_annotations(clinical_significance)",
                "CREATE INDEX IF NOT EXISTS idx_freq_pop ON allele_frequencies(population, allele_frequency)",
                "CREATE INDEX IF NOT EXISTS idx_pred_score ON in_silico_predictions(tool, score)"
            ]
            
            for idx_query in indices:
                cursor.execute(idx_query)
            
            conn.commit()
    
    def build_comprehensive_database(self, update_existing=False):
        """
        Costruisce database completo con gestione incrementale
        """
        logger.info("="*80)
        logger.info("BUILDING COMPREHENSIVE HEMOGLOBIN DATABASE")
        logger.info("="*80)
        
        try:
            # 1. HbVar
            logger.info("\n[1/4] Processing HbVar database...")
            hbvar_data = self.hbvar.download_all_variants(force_update=update_existing)
            self._process_hbvar_data(hbvar_data)
            
            # 2. ClinVar
            logger.info("\n[2/4] Processing ClinVar database...")
            clinvar_data = self.clinvar.search_hemoglobin_variants()
            self._process_clinvar_data(clinvar_data)
            
            # 3. gnomAD frequencies
            logger.info("\n[3/4] Processing gnomAD frequencies...")
            self._process_gnomad_frequencies()
            
            # 4. Generate reports
            logger.info("\n[4/4] Generating reports...")
            self._generate_comprehensive_report()
            
            logger.info("\n" + "="*80)
            logger.info("DATABASE BUILD COMPLETE!")
            logger.info(f"Location: {self.db_path}")
            logger.info("="*80)
            
        except Exception as e:
            logger.error(f"Build failed: {e}")
            raise
    
    def _process_hbvar_data(self, hbvar_data: List[Dict]):
        """Processa e integra dati HbVar"""
        if not hbvar_data:
            logger.warning("No HbVar data to process")
            return
        
        with DatabaseConnection.get_connection(self.db_path) as conn:
            cursor = conn.cursor()
            
            for variant in hbvar_data:
                # Genera ID consistente
                variant_id = self._generate_variant_id(variant)
                
                # Inserisci variante
                cursor.execute("""
                    INSERT OR IGNORE INTO variants (
                        variant_id, gene, hgvs_protein,
                        protein_position, aa_change
                    ) VALUES (?, ?, ?, ?, ?)
                """, (
                    variant_id,
                    variant.get('gene'),
                    variant.get('mutation_protein'),
                    variant.get('position'),
                    f"{variant.get('wild_type_aa')}>{variant.get('mutant_aa')}"
                    if variant.get('wild_type_aa') and variant.get('mutant_aa') else None
                ))
                
                # Annotazione clinica
                if variant.get('phenotype'):
                    cursor.execute("""
                        INSERT OR REPLACE INTO clinical_annotations (
                            variant_id, source, phenotype
                        ) VALUES (?, ?, ?)
                    """, (variant_id, 'HbVar', variant.get('phenotype')))
            
            conn.commit()
            logger.info(f"  Processed {len(hbvar_data)} HbVar variants")
    
    def _process_clinvar_data(self, clinvar_df: pd.DataFrame):
        """Processa e integra dati ClinVar"""
        if clinvar_df.empty:
            logger.warning("No ClinVar data to process")
            return
        
        with DatabaseConnection.get_connection(self.db_path) as conn:
            for _, row in clinvar_df.iterrows():
                # Estrai informazioni HGVS
                hgvs_info = self._parse_hgvs(row.get('hgvs_expressions', []))
                
                variant_id = self._generate_variant_id({
                    'gene': row['gene'],
                    'clinvar_id': row['variant_id']
                })
                
                # Inserisci variante
                conn.execute("""
                    INSERT OR IGNORE INTO variants (
                        variant_id, gene, hgvs_coding, hgvs_protein
                    ) VALUES (?, ?, ?, ?)
                """, (
                    variant_id,
                    row['gene'],
                    hgvs_info.get('coding'),
                    hgvs_info.get('protein')
                ))
                
                # Annotazione clinica
                conn.execute("""
                    INSERT OR REPLACE INTO clinical_annotations (
                        variant_id, source, clinical_significance,
                        phenotype, review_status, evidence_level
                    ) VALUES (?, ?, ?, ?, ?, ?)
                """, (
                    variant_id,
                    'ClinVar',
                    row.get('clinical_significance'),
                    row.get('conditions_str', ''),
                    row.get('review_status'),
                    str(row.get('confidence_score', 0))
                ))
            
            conn.commit()
            logger.info(f"  Processed {len(clinvar_df)} ClinVar variants")
    
    def _process_gnomad_frequencies(self):
        """Processa frequenze gnomAD per varianti esistenti"""
        genes = ['HBB', 'HBA1', 'HBA2']
        
        for gene in genes:
            logger.info(f"  Processing gnomAD for {gene}...")
            gnomad_df = self.gnomad.get_hemoglobin_gene_variants(gene)
            
            if not gnomad_df.empty:
                self._integrate_gnomad_data(gnomad_df, gene)
    
    def _integrate_gnomad_data(self, gnomad_df: pd.DataFrame, gene: str):
        """Integra dati gnomAD nel database"""
        with DatabaseConnection.get_connection(self.db_path) as conn:
            cursor = conn.cursor()
            
            for _, row in gnomad_df.iterrows():
                # Genera variant ID
                variant_id = self._generate_variant_id({
                    'gene': gene,
                    'position': row.get('pos'),
                    'variant_id': row.get('variant_id')
                })
                
                # Inserisci o aggiorna variante
                cursor.execute("""
                    INSERT OR IGNORE INTO variants (
                        variant_id, gene, position_grch38,
                        ref_allele, alt_allele
                    ) VALUES (?, ?, ?, ?, ?)
                """, (
                    variant_id,
                    gene,
                    row.get('pos'),
                    row.get('ref'),
                    row.get('alt')
                ))
                
                # Aggiungi frequenza allelica
                # Use json_normalize to handle nested data
                if 'genome' in row and row['genome'] and row['genome'].get('af') is not None:
                    cursor.execute("""
                        INSERT OR REPLACE INTO allele_frequencies (
                            variant_id, dataset, population,
                            allele_frequency
                        ) VALUES (?, ?, ?, ?)
                    """, (
                        variant_id,
                        'gnomAD_v3',
                        'global',
                        row['genome']['af']
                    ))
                
            conn.commit()
            logger.info(f"    Integrated {len(gnomad_df)} variants from gnomAD for {gene}")

    
    def _generate_variant_id(self, variant_info: Dict) -> str:
        """Genera ID univoco per variante"""
        # Usa combinazione di attributi chiave
        key_parts = []
        
        if variant_info.get('gene'):
            key_parts.append(variant_info['gene'])
        
        if variant_info.get('position'):
            key_parts.append(str(variant_info['position']))
        
        if variant_info.get('clinvar_id'):
            key_parts.append(f"cv_{variant_info['clinvar_id']}")
        
        if variant_info.get('hb_name'):
            key_parts.append(variant_info['hb_name'].replace(" ", ""))

        if variant_info.get('mutation'):
             key_parts.append(variant_info['mutation'].replace(" ", ""))

        if not key_parts:
            key_parts.append(str(variant_info))

        key_string = '_'.join(filter(None, key_parts))
        
        # Genera hash SHA256 deterministico
        return hashlib.sha256(key_string.encode()).hexdigest()[:32]
    
    def _parse_hgvs(self, hgvs_list: List[str]) -> Dict:
        """Parse delle espressioni HGVS"""
        result = {'genomic': None, 'coding': None, 'protein': None}
        
        if not isinstance(hgvs_list, list):
            return result
        
        for hgvs in hgvs_list:
            if ':g.' in hgvs:
                result['genomic'] = hgvs
            elif ':c.' in hgvs:
                result['coding'] = hgvs
            elif ':p.' in hgvs:
                result['protein'] = hgvs
        
        return result
    
    def _generate_comprehensive_report(self):
        """Genera report dettagliato del database"""
        report_path = os.path.join(self.output_dir, "database_report.md")
        
        with DatabaseConnection.get_connection(self.db_path) as conn:
            cursor = conn.cursor()
            
            with open(report_path, 'w') as f:
                f.write("# Hemoglobin Database Report\n\n")
                f.write(f"Generated: {datetime.now().isoformat()}\n\n")
                
                # Statistiche generali
                cursor.execute("SELECT COUNT(*) FROM variants")
                total_variants = cursor.fetchone()[0]
                
                f.write("## Summary Statistics\n\n")
                f.write(f"- Total unique variants: {total_variants}\n")
                
                # Per gene
                cursor.execute("""
                    SELECT gene, COUNT(*) as count
                    FROM variants
                    WHERE gene IS NOT NULL AND gene != ''
                    GROUP BY gene
                    ORDER BY count DESC
                """)
                
                f.write("\n### Variants by Gene\n\n")
                f.write("| Gene | Count |\n")
                f.write("|------|-------|\n")
                for row in cursor.fetchall():
                    f.write(f"| {row[0]} | {row[1]} |\n")
                
                # Per significato clinico
                cursor.execute("""
                    SELECT clinical_significance, COUNT(*) as count
                    FROM clinical_annotations
                    WHERE clinical_significance IS NOT NULL AND clinical_significance != ''
                    GROUP BY clinical_significance
                    ORDER BY count DESC
                    LIMIT 10
                """)
                
                f.write("\n### Top 10 Clinical Significance\n\n")
                f.write("| Significance | Count |\n")
                f.write("|--------------|-------|\n")
                for row in cursor.fetchall():
                    f.write(f"| {row[0]} | {row[1]} |\n")
                
                # Data sources
                cursor.execute("""
                    SELECT source, COUNT(DISTINCT variant_id) as count
                    FROM clinical_annotations
                    GROUP BY source
                """)
                
                f.write("\n### Variants by Data Source\n\n")
                f.write("| Source | Variants |\n")
                f.write("|--------|----------|\n")
                for row in cursor.fetchall():
                    f.write(f"| {row[0]} | {row[1]} |\n")
        
        logger.info(f"  Report saved to: {report_path}")
    
    def query_variant(self, **kwargs) -> pd.DataFrame:
        """
        Query flessibile del database
        
        Parameters:
        - gene: str, position: int, chromosome: str, clinical_significance: str
        - min_frequency: float, max_frequency: float
        """
        query = """
            SELECT DISTINCT
                v.gene, v.protein_position, v.aa_change, v.hgvs_protein,
                c.clinical_significance, c.phenotype, c.source as annotation_source,
                f.population, f.allele_frequency
            FROM variants v
            LEFT JOIN clinical_annotations c ON v.variant_id = c.variant_id
            LEFT JOIN allele_frequencies f ON v.variant_id = f.variant_id
            WHERE 1=1
        """
        params = {}
        if kwargs.get('gene'):
            query += " AND v.gene = :gene"
            params['gene'] = kwargs['gene']
        if kwargs.get('position'):
            query += " AND v.protein_position = :position"
            params['position'] = kwargs['position']
        if kwargs.get('clinical_significance'):
            query += " AND c.clinical_significance LIKE :clinical_significance"
            params['clinical_significance'] = f"%{kwargs['clinical_significance']}%"
        if kwargs.get('min_frequency'):
            query += " AND f.allele_frequency >= :min_frequency"
            params['min_frequency'] = kwargs['min_frequency']
        if kwargs.get('max_frequency'):
            query += " AND f.allele_frequency <= :max_frequency"
            params['max_frequency'] = kwargs['max_frequency']
        
        with DatabaseConnection.get_connection(self.db_path) as conn:
            return pd.read_sql_query(query, conn, params=params)


def main():
    """Esecuzione principale con gestione errori migliorata"""
    try:
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[logging.FileHandler('hemoglobin_db.log'), logging.StreamHandler()]
        )
        
        print("╔" + "═" * 78 + "╗")
        print("║" + " HEMOGLOBIN CLINICAL DATABASE INTEGRATOR v2.1 ".center(78) + "║")
        print("║" + " Corrected HbVar integration and enhanced logic ".center(78) + "║")
        print("╚" + "═" * 78 + "╝\n")
        
        config = Configuration('hemoglobin_config.ini')
        
        if 'your.email' in config.get_email():
            print("\n⚠️  WARNING: Please update your email in hemoglobin_config.ini")
            return
        
        integrator = IntegratedHemoglobinDatabase(config=config)
        
        while True:
            print("\nOptions:\n1. Build/Update complete database\n2. Query existing database\n3. Generate reports\n4. Exit")
            choice = input("\nSelect option (1-4): ").strip()
            
            if choice == '1':
                integrator.build_comprehensive_database(update_existing=True)
            elif choice == '2':
                print("\nQuery examples:\n1. All HBB variants\n2. Pathogenic variants\n3. HBB position 6 (Sickle Cell)")
                query_choice = input("Select query (1-3): ").strip()
                results = pd.DataFrame()
                if query_choice == '1': results = integrator.query_variant(gene='HBB')
                elif query_choice == '2': results = integrator.query_variant(clinical_significance='Pathogenic')
                elif query_choice == '3': results = integrator.query_variant(gene='HBB', position=6)
                
                if not results.empty:
                    print(f"\nFound {len(results)} variants:")
                    print(results.head(10).to_string())
                else:
                    print("No results found.")
            elif choice == '3':
                integrator._generate_comprehensive_report()
            elif choice == '4':
                print("\nGoodbye!")
                break
            else:
                print("Invalid option.")
    
    except KeyboardInterrupt:
        print("\n\nOperation cancelled by user.")
    except Exception as e:
        logger.error(f"A fatal error occurred: {e}", exc_info=True)
        raise

if __name__ == "__main__":
    main()