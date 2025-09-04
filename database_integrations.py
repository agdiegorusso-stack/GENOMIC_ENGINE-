#!/usr/bin/env python3
"""
Integrazione completa con database di emoglobinopatie
HbVar, ClinVar, LOVD, gnomAD, IthaGenes
Autore: Assistant
Data: 2025
"""

import os
import re
import json
import time
import sqlite3
import hashlib
from datetime import datetime
from typing import Dict, List, Optional, Tuple
import xml.etree.ElementTree as ET

import pandas as pd
import numpy as np
import requests
from bs4 import BeautifulSoup
from Bio import Entrez
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# Configurazione
Entrez.email = "agdiegorusso@gmail.com"  # MODIFICA


class HbVarDatabase:
    """
    Interfaccia per HbVar - Database delle varianti di emoglobina
    http://globin.bx.psu.edu/hbvar/
    """
    
    def __init__(self, cache_dir="hbvar_cache"):
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        self.base_url = "http://globin.bx.psu.edu/cgi-bin/hbvar"
        self.variants_db = self._initialize_database()
    
    def _initialize_database(self):
        """Inizializza database SQLite locale"""
        db_path = os.path.join(self.cache_dir, "hbvar.db")
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Crea tabella varianti
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS variants (
                variant_id TEXT PRIMARY KEY,
                hb_name TEXT,
                gene TEXT,
                mutation_dna TEXT,
                mutation_protein TEXT,
                chain_type TEXT,
                position INTEGER,
                wild_type_aa TEXT,
                mutant_aa TEXT,
                clinical_significance TEXT,
                phenotype TEXT,
                geographic_origin TEXT,
                frequency TEXT,
                reference TEXT,
                last_updated TIMESTAMP
            )
        """)
        
        # Crea tabella per dati clinici
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS clinical_data (
                variant_id TEXT,
                parameter TEXT,
                value TEXT,
                unit TEXT,
                FOREIGN KEY(variant_id) REFERENCES variants(variant_id)
            )
        """)
        
        conn.commit()
        return conn
    
    def download_all_variants(self):
        """
        Scarica tutte le varianti dal database HbVar
        """
        print("Downloading HbVar database...")
        
        # URL delle diverse categorie
        categories = [
            ('query_menu.html', 'All variants'),
            ('query_menu2.html', 'Alpha chain variants'),
            ('query_menu3.html', 'Beta chain variants'),
            ('query_menu4.html', 'Gamma chain variants'),
            ('query_menu5.html', 'Delta chain variants')
        ]
        
        all_variants = []
        
        for url_suffix, category in categories:
            print(f"  Fetching {category}...")
            url = f"{self.base_url}/{url_suffix}"
            
            try:
                response = requests.get(url, timeout=30)
                if response.status_code == 200:
                    variants = self._parse_hbvar_page(response.text, category)
                    all_variants.extend(variants)
                    print(f"    Found {len(variants)} variants")
            except Exception as e:
                print(f"    Error: {str(e)}")
            
            time.sleep(1)  # Rate limiting
        
        # Salva nel database
        self._save_variants_to_db(all_variants)
        
        print(f"Total variants downloaded: {len(all_variants)}")
        return all_variants
    
    def _parse_hbvar_page(self, html_content: str, category: str) -> List[Dict]:
        """
        Parse della pagina HTML di HbVar
        """
        soup = BeautifulSoup(html_content, 'html.parser')
        variants = []
        
        # Trova tabelle con varianti
        tables = soup.find_all('table')
        
        for table in tables:
            rows = table.find_all('tr')
            
            for row in rows[1:]:  # Skip header
                cells = row.find_all('td')
                
                if len(cells) >= 5:
                    variant = {
                        'hb_name': cells[0].get_text(strip=True),
                        'mutation': cells[1].get_text(strip=True),
                        'chain': cells[2].get_text(strip=True) if len(cells) > 2 else '',
                        'phenotype': cells[3].get_text(strip=True) if len(cells) > 3 else '',
                        'category': category
                    }
                    
                    # Parse mutation notation
                    mutation_info = self._parse_mutation_notation(variant['mutation'])
                    variant.update(mutation_info)
                    
                    variants.append(variant)
        
        return variants
    
    def _parse_mutation_notation(self, mutation_str: str) -> Dict:
        """
        Parse della notazione delle mutazioni (es. β6(A3)Glu→Val)
        """
        info = {
            'gene': '',
            'position': None,
            'wild_type_aa': '',
            'mutant_aa': '',
            'mutation_protein': mutation_str
        }
        
        # Pattern per diverse notazioni
        patterns = [
            r'([αβγδ])(\d+)\(?\w*\)?\s*(\w{3})→(\w{3})',  # β6(A3)Glu→Val
            r'([αβγδ])(\d+)\s+(\w{3})>(\w{3})',  # β6 Glu>Val
            r'([αβγδ])(\d+)\s+(\w+)\s+to\s+(\w+)',  # β6 Glu to Val
        ]
        
        for pattern in patterns:
            match = re.search(pattern, mutation_str)
            if match:
                chain_map = {'α': 'HBA', 'β': 'HBB', 'γ': 'HBG', 'δ': 'HBD'}
                info['gene'] = chain_map.get(match.group(1), match.group(1))
                info['position'] = int(match.group(2))
                info['wild_type_aa'] = match.group(3)
                info['mutant_aa'] = match.group(4)
                break
        
        return info
    
    def _save_variants_to_db(self, variants: List[Dict]):
        """Salva varianti nel database SQLite"""
        cursor = self.variants_db.cursor()
        
        for variant in variants:
            # Genera ID unico
            variant_id = hashlib.md5(
                f"{variant.get('hb_name', '')}_{variant.get('mutation', '')}".encode()
            ).hexdigest()[:16]
            
            cursor.execute("""
                INSERT OR REPLACE INTO variants (
                    variant_id, hb_name, gene, mutation_protein,
                    position, wild_type_aa, mutant_aa, phenotype,
                    last_updated
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                variant_id,
                variant.get('hb_name'),
                variant.get('gene'),
                variant.get('mutation_protein'),
                variant.get('position'),
                variant.get('wild_type_aa'),
                variant.get('mutant_aa'),
                variant.get('phenotype'),
                datetime.now()
            ))
        
        self.variants_db.commit()
    
    def query_variant(self, gene: str, position: int, aa_change: str) -> Optional[Dict]:
        """
        Query per una specifica variante
        """
        cursor = self.variants_db.cursor()
        
        cursor.execute("""
            SELECT * FROM variants
            WHERE gene = ? AND position = ? 
            AND (wild_type_aa = ? OR mutant_aa = ?)
        """, (gene, position, aa_change.split('>')[0] if '>' in aa_change else aa_change,
              aa_change.split('>')[-1] if '>' in aa_change else aa_change))
        
        result = cursor.fetchone()
        
        if result:
            columns = [desc[0] for desc in cursor.description]
            return dict(zip(columns, result))
        
        return None
    
    def get_statistics(self) -> Dict:
        """
        Statistiche del database
        """
        cursor = self.variants_db.cursor()
        
        stats = {}
        
        # Totale varianti
        cursor.execute("SELECT COUNT(*) FROM variants")
        stats['total_variants'] = cursor.fetchone()[0]
        
        # Per gene
        cursor.execute("""
            SELECT gene, COUNT(*) as count
            FROM variants
            WHERE gene IS NOT NULL
            GROUP BY gene
        """)
        stats['by_gene'] = dict(cursor.fetchall())
        
        # Per fenotipo
        cursor.execute("""
            SELECT phenotype, COUNT(*) as count
            FROM variants
            WHERE phenotype IS NOT NULL
            GROUP BY phenotype
            ORDER BY count DESC
            LIMIT 10
        """)
        stats['top_phenotypes'] = dict(cursor.fetchall())
        
        return stats


class IthaGenesDatabase:
    """
    Interfaccia per IthaGenes - Database delle talassemie
    https://www.ithanet.eu/db/ithagenes
    """
    
    def __init__(self, cache_dir="ithagenes_cache"):
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        self.base_url = "https://www.ithanet.eu/db/ithagenes"
    
    def search_mutations(self, gene: str = None) -> List[Dict]:
        """
        Cerca mutazioni nel database IthaGenes
        """
        print(f"Searching IthaGenes for {gene if gene else 'all genes'}...")
        
        mutations = []
        
        # Parametri di ricerca
        params = {
            'action': 'search',
            'gene': gene if gene else '',
            'mutation_type': 'all'
        }
        
        try:
            response = requests.get(self.base_url, params=params, timeout=30)
            
            if response.status_code == 200:
                soup = BeautifulSoup(response.text, 'html.parser')
                
                # Parse risultati
                results_table = soup.find('table', {'class': 'results'})
                if results_table:
                    rows = results_table.find_all('tr')[1:]  # Skip header
                    
                    for row in rows:
                        cells = row.find_all('td')
                        if len(cells) >= 6:
                            mutation = {
                                'gene': cells[0].get_text(strip=True),
                                'mutation_name': cells[1].get_text(strip=True),
                                'hgvs_notation': cells[2].get_text(strip=True),
                                'phenotype': cells[3].get_text(strip=True),
                                'severity': cells[4].get_text(strip=True),
                                'ethnicity': cells[5].get_text(strip=True)
                            }
                            mutations.append(mutation)
            
            print(f"  Found {len(mutations)} mutations")
            
        except Exception as e:
            print(f"  Error: {str(e)}")
        
        return mutations
    
    def get_thalassemia_modifiers(self) -> List[Dict]:
        """
        Ottiene lista dei geni modificatori della talassemia
        """
        modifiers = [
            {
                'gene': 'BCL11A',
                'effect': 'Fetal hemoglobin persistence',
                'snps': ['rs11886868', 'rs4671393', 'rs7557939']
            },
            {
                'gene': 'HBS1L-MYB',
                'effect': 'Fetal hemoglobin levels',
                'snps': ['rs9399137', 'rs4895441', 'rs9402686']
            },
            {
                'gene': 'KLF1',
                'effect': 'Fetal hemoglobin persistence',
                'snps': ['rs112631212']
            },
            {
                'gene': 'GATA1',
                'effect': 'Erythropoiesis regulation',
                'snps': []
            }
        ]
        
        return modifiers


class ClinVarIntegrator:
    """
    Integrazione avanzata con ClinVar per varianti di emoglobina
    """
    
    def __init__(self):
        Entrez.email = "your.email@example.com"  # Modifica
        self.cache = {}
    
    def search_hemoglobin_variants(self, batch_size: int = 100) -> pd.DataFrame:
        """
        Cerca tutte le varianti di emoglobina in ClinVar
        """
        print("Searching ClinVar for hemoglobin variants...")
        
        # Geni di interesse
        genes = ['HBB', 'HBA1', 'HBA2', 'HBG1', 'HBG2', 'HBD', 'HBE1', 'HBZ', 'HBM']
        
        all_variants = []
        
        for gene in genes:
            print(f"  Searching {gene}...")
            
            # Search query
            query = f"{gene}[gene] AND pathogenic[Clinical_Significance]"
            
            try:
                # Search
                handle = Entrez.esearch(
                    db="clinvar",
                    term=query,
                    retmax=batch_size
                )
                result = Entrez.read(handle)
                handle.close()
                
                if result['IdList']:
                    # Fetch details
                    handle = Entrez.efetch(
                        db="clinvar",
                        id=','.join(result['IdList']),
                        rettype="vcv",
                        retmode="xml"
                    )
                    xml_data = handle.read()
                    handle.close()
                    
                    # Parse XML
                    variants = self._parse_clinvar_xml(xml_data, gene)
                    all_variants.extend(variants)
                    
                    print(f"    Found {len(variants)} pathogenic variants")
                
                time.sleep(0.5)  # Rate limit
                
            except Exception as e:
                print(f"    Error: {str(e)}")
        
        df = pd.DataFrame(all_variants)
        print(f"\nTotal ClinVar variants found: {len(df)}")
        
        return df
    
    def _parse_clinvar_xml(self, xml_data: str, gene: str) -> List[Dict]:
        """
        Parse XML da ClinVar
        """
        variants = []
        
        try:
            root = ET.fromstring(xml_data)
            
            for record in root.findall('.//VariationArchive'):
                variant = {
                    'gene': gene,
                    'variant_id': record.get('VariationID'),
                    'name': record.get('VariationName'),
                    'clinical_significance': '',
                    'review_status': '',
                    'conditions': []
                }
                
                # Clinical significance
                clin_sig = record.find('.//ClinicalSignificance/Description')
                if clin_sig is not None:
                    variant['clinical_significance'] = clin_sig.text
                
                # Review status
                review = record.find('.//ClinicalSignificance/ReviewStatus')
                if review is not None:
                    variant['review_status'] = review.text
                
                # Conditions
                conditions = record.findall('.//TraitSet/Trait/Name/ElementValue')
                variant['conditions'] = [c.text for c in conditions]
                
                # HGVS notation
                hgvs = record.find('.//SequenceLocation[@Assembly="GRCh38"]')
                if hgvs is not None:
                    variant['hgvs'] = hgvs.get('variantLength')
                
                variants.append(variant)
                
        except Exception as e:
            print(f"    XML parsing error: {str(e)}")
        
        return variants


class GnomADIntegrator:
    """
    Integrazione con gnomAD per frequenze popolazionali
    """
    
    def __init__(self):
        self.base_url = "https://gnomad.broadinstitute.org/api"
        self.cache = {}
    
    def get_hemoglobin_frequencies(self, gene: str, transcript_id: str = None) -> pd.DataFrame:
        """
        Ottiene frequenze alleliche per varianti di un gene
        """
        print(f"Fetching gnomAD frequencies for {gene}...")
        
        query = """
        query GeneVariants($gene: String!, $datasetId: DatasetId!) {
            gene(gene_symbol: $gene, reference_genome: GRCh38) {
                variants(dataset: $datasetId) {
                    variant_id
                    pos
                    ref
                    alt
                    consequence
                    hgvsp
                    exome {
                        ac
                        an
                        af
                    }
                    genome {
                        ac
                        an
                        af
                    }
                    populations {
                        id
                        ac
                        an
                        af
                    }
                }
            }
        }
        """
        
        variables = {
            "gene": gene,
            "datasetId": "gnomad_r3"
        }
        
        try:
            response = requests.post(
                self.base_url,
                json={'query': query, 'variables': variables},
                timeout=30
            )
            
            if response.status_code == 200:
                data = response.json()
                
                if data.get('data', {}).get('gene', {}).get('variants'):
                    variants = data['data']['gene']['variants']
                    
                    # Converti in DataFrame
                    df = pd.DataFrame(variants)
                    
                    # Espandi colonne nested
                    if not df.empty:
                        df['af_exome'] = df['exome'].apply(
                            lambda x: x.get('af', 0) if x else 0
                        )
                        df['af_genome'] = df['genome'].apply(
                            lambda x: x.get('af', 0) if x else 0
                        )
                        
                        # Classifica per frequenza
                        df['frequency_class'] = pd.cut(
                            df['af_genome'].fillna(0),
                            bins=[0, 0.0001, 0.001, 0.01, 0.05, 1.0],
                            labels=['ultra_rare', 'rare', 'low_freq', 'common', 'very_common']
                        )
                    
                    print(f"  Found {len(df)} variants")
                    return df
            
            return pd.DataFrame()
            
        except Exception as e:
            print(f"  Error: {str(e)}")
            return pd.DataFrame()
    
    def get_population_specific_frequencies(self, variant_id: str) -> Dict:
        """
        Ottiene frequenze popolazione-specifiche per una variante
        """
        # Implementazione semplificata
        populations = {
            'afr': 'African/African American',
            'amr': 'Latino/Admixed American',
            'asj': 'Ashkenazi Jewish',
            'eas': 'East Asian',
            'fin': 'Finnish',
            'nfe': 'Non-Finnish European',
            'sas': 'South Asian'
        }
        
        # In produzione, faresti query reale a gnomAD
        # Qui simuliamo con dati di esempio
        
        frequencies = {}
        for pop_id, pop_name in populations.items():
            frequencies[pop_name] = {
                'af': np.random.exponential(0.001),  # Simulato
                'ac': np.random.poisson(5),
                'an': 10000
            }
        
        return frequencies


class IntegratedHemoglobinDatabase:
    """
    Database integrato che combina tutte le fonti
    """
    
    def __init__(self, output_dir="integrated_db"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Inizializza tutti i database
        self.hbvar = HbVarDatabase()
        self.ithagenes = IthaGenesDatabase()
        self.clinvar = ClinVarIntegrator()
        self.gnomad = GnomADIntegrator()
        
        # Database integrato SQLite
        self.integrated_db = self._create_integrated_database()
    
    def _create_integrated_database(self):
        """
        Crea database SQLite integrato
        """
        db_path = os.path.join(self.output_dir, "hemoglobin_mutations.db")
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Tabella principale mutazioni
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS mutations (
                mutation_id TEXT PRIMARY KEY,
                gene TEXT,
                chromosome TEXT,
                position_genomic INTEGER,
                position_cds INTEGER,
                position_protein INTEGER,
                ref_nucleotide TEXT,
                alt_nucleotide TEXT,
                ref_amino_acid TEXT,
                alt_amino_acid TEXT,
                mutation_type TEXT,
                hgvs_cdna TEXT,
                hgvs_protein TEXT,
                last_updated TIMESTAMP
            )
        """)
        
        # Tabella annotazioni cliniche
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS clinical_annotations (
                mutation_id TEXT,
                source TEXT,
                clinical_significance TEXT,
                phenotype TEXT,
                severity TEXT,
                inheritance_pattern TEXT,
                penetrance TEXT,
                evidence_level TEXT,
                pmid TEXT,
                FOREIGN KEY(mutation_id) REFERENCES mutations(mutation_id)
            )
        """)
        
        # Tabella frequenze popolazionali
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS population_frequencies (
                mutation_id TEXT,
                population TEXT,
                allele_count INTEGER,
                allele_number INTEGER,
                allele_frequency REAL,
                homozygotes INTEGER,
                source TEXT,
                FOREIGN KEY(mutation_id) REFERENCES mutations(mutation_id)
            )
        """)
        
        # Tabella predizioni funzionali
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS functional_predictions (
                mutation_id TEXT,
                predictor TEXT,
                score REAL,
                prediction TEXT,
                confidence REAL,
                FOREIGN KEY(mutation_id) REFERENCES mutations(mutation_id)
            )
        """)
        
        conn.commit()
        return conn
    
    def build_comprehensive_database(self):
        """
        Costruisce database completo integrando tutte le fonti
        """
        print("\n" + "="*80)
        print("BUILDING COMPREHENSIVE HEMOGLOBIN MUTATION DATABASE")
        print("="*80)
        
        # 1. Download HbVar
        print("\n1. Downloading HbVar database...")
        hbvar_variants = self.hbvar.download_all_variants()
        
        # 2. Query ClinVar
        print("\n2. Querying ClinVar...")
        clinvar_df = self.clinvar.search_hemoglobin_variants()
        
        # 3. Query IthaGenes
        print("\n3. Querying IthaGenes...")
        ithagenes_mutations = []
        for gene in ['HBA1', 'HBA2', 'HBB', 'HBG1', 'HBG2']:
            ithagenes_mutations.extend(self.ithagenes.search_mutations(gene))
        
        # 4. Get gnomAD frequencies
        print("\n4. Fetching gnomAD frequencies...")
        gnomad_data = {}
        for gene in ['HBB', 'HBA1', 'HBA2']:
            gnomad_data[gene] = self.gnomad.get_hemoglobin_frequencies(gene)
        
        # 5. Integra tutti i dati
        print("\n5. Integrating all data sources...")
        self._integrate_all_sources(
            hbvar_variants,
            clinvar_df,
            ithagenes_mutations,
            gnomad_data
        )
        
        # 6. Genera report
        print("\n6. Generating summary report...")
        self._generate_summary_report()
        
        print("\n" + "="*80)
        print("DATABASE BUILD COMPLETE!")
        print(f"Database location: {self.output_dir}/hemoglobin_mutations.db")
        print("="*80)
    
    def _integrate_all_sources(self, hbvar, clinvar, ithagenes, gnomad):
        """
        Integra dati da tutte le fonti
        """
        cursor = self.integrated_db.cursor()
        
        # Mappatura geni standard
        gene_mapping = {
            'α': 'HBA', 'alpha': 'HBA', 'a': 'HBA',
            'β': 'HBB', 'beta': 'HBB', 'b': 'HBB',
            'γ': 'HBG', 'gamma': 'HBG', 'g': 'HBG',
            'δ': 'HBD', 'delta': 'HBD', 'd': 'HBD'
        }
        
        mutations_added = 0
        
        # Integra HbVar
        for variant in hbvar:
            if variant.get('gene') and variant.get('position'):
                mutation_id = f"{variant['gene']}_p{variant['position']}"
                
                cursor.execute("""
                    INSERT OR IGNORE INTO mutations (
                        mutation_id, gene, position_protein,
                        ref_amino_acid, alt_amino_acid,
                        mutation_type, hgvs_protein, last_updated
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    mutation_id,
                    variant['gene'],
                    variant['position'],
                    variant.get('wild_type_aa'),
                    variant.get('mutant_aa'),
                    'missense',
                    variant.get('mutation_protein'),
                    datetime.now()
                ))
                
                # Aggiungi annotazione clinica
                if variant.get('phenotype'):
                    cursor.execute("""
                        INSERT INTO clinical_annotations (
                            mutation_id, source, phenotype
                        ) VALUES (?, ?, ?)
                    """, (mutation_id, 'HbVar', variant['phenotype']))
                
                mutations_added += 1
        
        # Integra ClinVar
        if not clinvar.empty:
            for _, row in clinvar.iterrows():
                mutation_id = f"{row['gene']}_{row['variant_id']}"
                
                cursor.execute("""
                    INSERT OR IGNORE INTO mutations (
                        mutation_id, gene, hgvs_cdna,
                        mutation_type, last_updated
                    ) VALUES (?, ?, ?, ?, ?)
                """, (
                    mutation_id,
                    row['gene'],
                    row.get('hgvs'),
                    'variant',
                    datetime.now()
                ))
                
                cursor.execute("""
                    INSERT INTO clinical_annotations (
                        mutation_id, source, clinical_significance,
                        phenotype, evidence_level
                    ) VALUES (?, ?, ?, ?, ?)
                """, (
                    mutation_id,
                    'ClinVar',
                    row['clinical_significance'],
                    ','.join(row.get('conditions', [])),
                    row.get('review_status')
                ))
        
        # Integra frequenze gnomAD
        for gene, df in gnomad.items():
            if not df.empty:
                for _, row in df.iterrows():
                    mutation_id = f"{gene}_{row['pos']}_{row['ref']}>{row['alt']}"
                    
                    cursor.execute("""
                        INSERT OR IGNORE INTO mutations (
                            mutation_id, gene, position_genomic,
                            ref_nucleotide, alt_nucleotide,
                            last_updated
                        ) VALUES (?, ?, ?, ?, ?, ?)
                    """, (
                        mutation_id,
                        gene,
                        row['pos'],
                        row['ref'],
                        row['alt'],
                        datetime.now()
                    ))
                    
                    # Aggiungi frequenze
                    if row.get('populations'):
                        for pop in row['populations']:
                            cursor.execute("""
                                INSERT INTO population_frequencies (
                                    mutation_id, population,
                                    allele_count, allele_number,
                                    allele_frequency, source
                                ) VALUES (?, ?, ?, ?, ?, ?)
                            """, (
                                mutation_id,
                                pop.get('id'),
                                pop.get('ac'),
                                pop.get('an'),
                                pop.get('af'),
                                'gnomAD'
                            ))
        
        self.integrated_db.commit()
        print(f"  Integrated {mutations_added} mutations")
    
    def _generate_summary_report(self):
        """
        Genera report riassuntivo del database
        """
        cursor = self.integrated_db.cursor()
        
        report_file = os.path.join(self.output_dir, "database_summary.txt")
        
        with open(report_file, 'w') as f:
            f.write("HEMOGLOBIN MUTATION DATABASE SUMMARY\n")
            f.write("="*60 + "\n")
            f.write(f"Generated: {datetime.now()}\n\n")
            
            # Statistiche generali
            cursor.execute("SELECT COUNT(*) FROM mutations")
            total_mutations = cursor.fetchone()[0]
            f.write(f"Total mutations: {total_mutations}\n")
            
            # Per gene
            cursor.execute("""
                SELECT gene, COUNT(*) as count
                FROM mutations
                GROUP BY gene
                ORDER BY count DESC
            """)
            f.write("\nMutations by gene:\n")
            for gene, count in cursor.fetchall():
                f.write(f"  {gene}: {count}\n")
            
            # Per significato clinico
            cursor.execute("""
                SELECT clinical_significance, COUNT(*) as count
                FROM clinical_annotations
                WHERE clinical_significance IS NOT NULL
                GROUP BY clinical_significance
                ORDER BY count DESC
            """)
            f.write("\nClinical significance:\n")
            for sig, count in cursor.fetchall():
                f.write(f"  {sig}: {count}\n")
            
            # Fenotipi più comuni
            cursor.execute("""
                SELECT phenotype, COUNT(*) as count
                FROM clinical_annotations
                WHERE phenotype IS NOT NULL
                GROUP BY phenotype
                ORDER BY count DESC
                LIMIT 20
            """)
            f.write("\nTop phenotypes:\n")
            for phenotype, count in cursor.fetchall():
                if phenotype:
                    f.write(f"  {phenotype[:50]}: {count}\n")
            
            # Fonti dati
            cursor.execute("""
                SELECT source, COUNT(DISTINCT mutation_id) as count
                FROM clinical_annotations
                GROUP BY source
            """)
            f.write("\nData sources:\n")
            for source, count in cursor.fetchall():
                f.write(f"  {source}: {count} mutations\n")
        
        print(f"  Report saved to: {report_file}")
    
    def query_mutation(self, gene: str, position: int = None, 
                       amino_acid_change: str = None) -> pd.DataFrame:
        """
        Query del database integrato
        """
        query = """
            SELECT DISTINCT
                m.mutation_id,
                m.gene,
                m.position_protein,
                m.ref_amino_acid,
                m.alt_amino_acid,
                m.hgvs_protein,
                c.clinical_significance,
                c.phenotype,
                c.source as clinical_source,
                p.population,
                p.allele_frequency
            FROM mutations m
            LEFT JOIN clinical_annotations c ON m.mutation_id = c.mutation_id
            LEFT JOIN population_frequencies p ON m.mutation_id = p.mutation_id
            WHERE m.gene = ?
        """
        
        params = [gene]
        
        if position:
            query += " AND m.position_protein = ?"
            params.append(position)
        
        if amino_acid_change:
            query += " AND (m.ref_amino_acid = ? OR m.alt_amino_acid = ?)"
            params.extend([amino_acid_change, amino_acid_change])
        
        df = pd.read_sql_query(query, self.integrated_db, params=params)
        
        return df


def main():
    """
    Esegue l'integrazione completa dei database
    """
    print("╔" + "═" * 78 + "╗")
    print("║" + " HEMOGLOBIN CLINICAL DATABASE INTEGRATOR ".center(78) + "║")
    print("║" + " Complete Integration of Clinical and Population Databases ".center(78) + "║")
    print("╚" + "═" * 78 + "╝")
    
    # Inizializza sistema integrato
    integrator = IntegratedHemoglobinDatabase()
    
    # Costruisci database completo
    integrator.build_comprehensive_database()
    
    # Esempi di query
    print("\n" + "="*80)
    print("EXAMPLE QUERIES")
    print("="*80)
    
    # Query 1: Tutte le mutazioni di HBB
    print("\n1. All HBB mutations:")
    hbb_mutations = integrator.query_mutation('HBB')
    print(f"   Found {len(hbb_mutations)} HBB mutations")
    
    # Query 2: Mutazione specifica (HbS)
    print("\n2. HbS mutation (β6 Glu>Val):")
    hbs = integrator.query_mutation('HBB', position=6)
    if not hbs.empty:
        print(hbs[['gene', 'position_protein', 'clinical_significance', 'phenotype']].head())
    
    # Statistiche database
    print("\n3. Database statistics:")
    stats = integrator.hbvar.get_statistics()
    print(f"   Total variants: {stats['total_variants']}")
    print(f"   By gene: {stats['by_gene']}")
    
    print("\n" + "="*80)
    print("INTEGRATION COMPLETE!")
    print("="*80)
    print("\nOutput files:")
    print(f"  - Integrated database: {integrator.output_dir}/hemoglobin_mutations.db")
    print(f"  - Summary report: {integrator.output_dir}/database_summary.txt")
    print("\nYou can now:")
    print("  1. Query the integrated database for any mutation")
    print("  2. Cross-reference clinical significance across sources")
    print("  3. Analyze population-specific frequencies")
    print("  4. Export data for machine learning models")


if __name__ == "__main__":
    main()