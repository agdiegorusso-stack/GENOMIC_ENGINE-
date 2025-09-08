#!/usr/bin/env python3
"""
Sistema Completo di Analisi Mutazionale per il Sistema Emoglobinico
Genera tutte le mutazioni possibili, predice effetti e integra database clinici
Autore: Assistant
Data: 2025
Requisiti: pip install biopython pandas numpy requests xmltodict scikit-learn
"""

import os
import sys
import json
import time
import pickle
import hashlib
import itertools
import warnings
from datetime import datetime
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing as mp

import pandas as pd
import numpy as np
import requests
import xmltodict
from Bio import Entrez, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from Bio.Align import substitution_matrices
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler

# GitHub integrations
try:
    import hail as hl
    HAIL_AVAILABLE = True
except ImportError:
    HAIL_AVAILABLE = False
    hl = None  # Define hl as None to prevent NameError
    print("Hail not available. Install with: pip install hail")

try:
    from SigProfilerMatrixGenerator import install as sig_install
    from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGenerator as sig
    SIGPROFILER_AVAILABLE = True
except ImportError:
    SIGPROFILER_AVAILABLE = False
    sig_install = None
    sig = None
    print("SigProfilerMatrixGenerator not available. Install with: pip install SigProfilerMatrixGenerator")

# Configurazione
Entrez.email = "agdiegorusso@gmail.com"  # MODIFICA CON LA TUA EMAIL
warnings.filterwarnings('ignore')

class TranscriptProteinMapper:
    """
    Comprehensive ENST -> ENSP mapping using Ensembl REST API with local caching
    """
    
    def __init__(self, cache_dir: str = "transcript_cache"):
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        self.cache_file = os.path.join(cache_dir, "enst_ensp_mappings.json")
        self.mappings = self._load_cache()
        self.session = requests.Session()
        self.session.headers.update({
            'Content-Type': 'application/json',
            'Accept': 'application/json'
        })
    
    def _load_cache(self) -> Dict[str, str]:
        """Load cached ENST->ENSP mappings"""
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, 'r') as f:
                    mappings = json.load(f)
                print(f"📁 Loaded {len(mappings)} cached transcript mappings")
                return mappings
            except Exception as e:
                print(f"⚠️ Error loading cache: {e}")
        return {}
    
    def _save_cache(self):
        """Save mappings to cache"""
        try:
            with open(self.cache_file, 'w') as f:
                json.dump(self.mappings, f, indent=2)
            print(f"💾 Saved {len(self.mappings)} transcript mappings to cache")
        except Exception as e:
            print(f"❌ Error saving cache: {e}")
    
    def get_protein_id(self, transcript_id: str) -> Optional[str]:
        """
        Get ENSP ID for given ENST ID
        Uses cache first, then queries Ensembl REST API
        """
        # Remove version if present (ENST00000123.4 -> ENST00000123)
        clean_transcript_id = transcript_id.split('.')[0]
        
        # Check cache first
        if clean_transcript_id in self.mappings:
            return self.mappings[clean_transcript_id]
        
        # Query Ensembl REST API
        protein_id = self._query_ensembl_api(clean_transcript_id)
        
        if protein_id:
            # Cache the result
            self.mappings[clean_transcript_id] = protein_id
            # Save cache periodically
            if len(self.mappings) % 50 == 0:
                self._save_cache()
        
        return protein_id
    
    def _query_ensembl_api(self, transcript_id: str) -> Optional[str]:
        """Query Ensembl REST API for transcript->protein mapping"""
        try:
            url = f"https://rest.ensembl.org/lookup/id/{transcript_id}"
            params = {'expand': '1'}
            
            response = self.session.get(url, params=params, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                
                # Extract protein ID from Translation object
                if 'Translation' in data:
                    protein_id = data['Translation'].get('id')
                    if protein_id:
                        print(f"🔍 {transcript_id} -> {protein_id}")
                        return protein_id
                else:
                    print(f"⚠️ No protein translation for {transcript_id}")
                    # Cache negative results to avoid repeated queries
                    self.mappings[transcript_id] = None
                    return None
            
            elif response.status_code == 429:
                # Rate limiting
                print(f"⏱️ Rate limited, waiting 2 seconds...")
                time.sleep(2)
                # Retry once
                return self._query_ensembl_api(transcript_id)
            
            else:
                print(f"❌ API error {response.status_code} for {transcript_id}")
                return None
                
        except Exception as e:
            print(f"❌ Error querying Ensembl API for {transcript_id}: {e}")
            return None
    
    def bulk_map_transcripts(self, transcript_ids: List[str]) -> Dict[str, Optional[str]]:
        """
        Map multiple transcript IDs efficiently
        Uses batch API when available, falls back to individual queries
        """
        results = {}
        unmapped_transcripts = []
        
        # Check cache first
        for transcript_id in transcript_ids:
            clean_id = transcript_id.split('.')[0]
            if clean_id in self.mappings:
                results[transcript_id] = self.mappings[clean_id]
            else:
                unmapped_transcripts.append(clean_id)
        
        print(f"📊 Found {len(results)} cached, querying {len(unmapped_transcripts)} from API")
        
        # Query unmapped transcripts
        for transcript_id in unmapped_transcripts:
            protein_id = self._query_ensembl_api(transcript_id)
            results[transcript_id] = protein_id
            # Small delay to be respectful to API
            time.sleep(0.1)
        
        # Save cache after bulk operations
        if unmapped_transcripts:
            self._save_cache()
        
        return results
    
    def get_canonical_transcript_mapping(self, gene_symbol: str) -> Dict[str, str]:
        """Get canonical transcript->protein mapping for a gene"""
        try:
            url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}"
            params = {'expand': '1'}
            
            response = self.session.get(url, params=params, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                canonical_transcript = data.get('canonical_transcript')
                
                if canonical_transcript:
                    protein_id = self.get_protein_id(canonical_transcript)
                    if protein_id:
                        return {canonical_transcript: protein_id}
            
        except Exception as e:
            print(f"❌ Error getting canonical transcript for {gene_symbol}: {e}")
        
        return {}
    
    def close(self):
        """Save cache and close session"""
        self._save_cache()

class MutationGenerator:
    """
    Genera sistematicamente tutte le possibili mutazioni per ogni gene
    """
    
    def __init__(self):
        self.codon_table = CodonTable.standard_dna_table
        self.stop_codons = self.codon_table.stop_codons
        self.aa_properties = self._load_aa_properties()
        
    def _load_aa_properties(self):
        """Proprietà fisico-chimiche degli aminoacidi"""
        return {
            'A': {'hydrophobic': 1.8, 'size': 89, 'charge': 0, 'polar': 0},
            'R': {'hydrophobic': -4.5, 'size': 174, 'charge': 1, 'polar': 1},
            'N': {'hydrophobic': -3.5, 'size': 132, 'charge': 0, 'polar': 1},
            'D': {'hydrophobic': -3.5, 'size': 133, 'charge': -1, 'polar': 1},
            'C': {'hydrophobic': 2.5, 'size': 121, 'charge': 0, 'polar': 0},
            'E': {'hydrophobic': -3.5, 'size': 147, 'charge': -1, 'polar': 1},
            'Q': {'hydrophobic': -3.5, 'size': 146, 'charge': 0, 'polar': 1},
            'G': {'hydrophobic': -0.4, 'size': 75, 'charge': 0, 'polar': 0},
            'H': {'hydrophobic': -3.2, 'size': 155, 'charge': 0.5, 'polar': 1},
            'I': {'hydrophobic': 4.5, 'size': 131, 'charge': 0, 'polar': 0},
            'L': {'hydrophobic': 3.8, 'size': 131, 'charge': 0, 'polar': 0},
            'K': {'hydrophobic': -3.9, 'size': 146, 'charge': 1, 'polar': 1},
            'M': {'hydrophobic': 1.9, 'size': 149, 'charge': 0, 'polar': 0},
            'F': {'hydrophobic': 2.8, 'size': 165, 'charge': 0, 'polar': 0},
            'P': {'hydrophobic': -1.6, 'size': 115, 'charge': 0, 'polar': 0},
            'S': {'hydrophobic': -0.8, 'size': 105, 'charge': 0, 'polar': 1},
            'T': {'hydrophobic': -0.7, 'size': 119, 'charge': 0, 'polar': 1},
            'W': {'hydrophobic': -0.9, 'size': 204, 'charge': 0, 'polar': 0},
            'Y': {'hydrophobic': -1.3, 'size': 181, 'charge': 0, 'polar': 1},
            'V': {'hydrophobic': 4.2, 'size': 117, 'charge': 0, 'polar': 0},
            '*': {'hydrophobic': 0, 'size': 0, 'charge': 0, 'polar': 0}
        }
    
    def generate_all_single_nucleotide_mutations(self, gene_name: str, cds_sequence: str) -> pd.DataFrame:
        """
        Genera tutte le possibili mutazioni puntiformi per un gene
        """
        print(f"\nGenerating mutations for {gene_name}...")
        mutations = []
        seq_obj = Seq(cds_sequence)
        protein_wt = str(seq_obj.translate())
        
        # Per ogni posizione nel CDS
        for pos in range(len(cds_sequence)):
            original_nt = cds_sequence[pos]
            codon_pos = pos // 3
            pos_in_codon = pos % 3
            
            # Testa tutte le possibili sostituzioni
            for new_nt in ['A', 'T', 'G', 'C']:
                if new_nt == original_nt:
                    continue
                
                # Crea la sequenza mutata
                mutated_seq = cds_sequence[:pos] + new_nt + cds_sequence[pos+1:]
                mutated_seq_obj = Seq(mutated_seq)
                
                try:
                    protein_mut = str(mutated_seq_obj.translate())
                except:
                    protein_mut = None
                
                # Classifica la mutazione
                mutation_info = self._classify_mutation(
                    cds_sequence, mutated_seq, 
                    protein_wt, protein_mut,
                    pos, codon_pos
                )
                
                mutation_info.update({
                    'gene': gene_name,
                    'position': pos + 1,  # 1-based
                    'ref_nt': original_nt,
                    'alt_nt': new_nt,
                    'codon_position': codon_pos + 1,
                    'position_in_codon': pos_in_codon + 1,
                    'mutation_id': f"{gene_name}:c.{pos+1}{original_nt}>{new_nt}"
                })
                
                mutations.append(mutation_info)
        
        df = pd.DataFrame(mutations)
        print(f"  Generated {len(df)} mutations")
        print(f"  - Synonymous: {len(df[df['type'] == 'synonymous'])}")
        print(f"  - Missense: {len(df[df['type'] == 'missense'])}")
        print(f"  - Nonsense: {len(df[df['type'] == 'nonsense'])}")
        print(f"  - Start lost: {len(df[df['type'] == 'start_lost'])}")
        
        return df
    
    def _classify_mutation(self, wt_cds, mut_cds, wt_prot, mut_prot, nt_pos, aa_pos):
        """
        Classifica dettagliata della mutazione
        """
        # Estrai i codoni
        codon_start = (nt_pos // 3) * 3
        codon_end = codon_start + 3
        wt_codon = wt_cds[codon_start:codon_end]
        mut_codon = mut_cds[codon_start:codon_end]
        
        # Caso speciale: mutazione nel codone di start
        if aa_pos == 0 and wt_codon == 'ATG' and mut_codon != 'ATG':
            return {
                'type': 'start_lost',
                'wt_aa': 'M',
                'mut_aa': '-',
                'aa_change': 'p.Met1?',
                'severity_score': 10
            }
        
        # Traduci i codoni
        try:
            wt_aa = str(Seq(wt_codon).translate())
            mut_aa = str(Seq(mut_codon).translate())
        except:
            return {
                'type': 'unknown',
                'wt_aa': '?',
                'mut_aa': '?',
                'aa_change': 'p.?',
                'severity_score': 0
            }
        
        # Classifica
        if wt_aa == mut_aa:
            mutation_type = 'synonymous'
            aa_change = f"p.{self._aa_three_letter(wt_aa)}{aa_pos+1}="
            severity = 0
        elif mut_aa == '*':
            mutation_type = 'nonsense'
            aa_change = f"p.{self._aa_three_letter(wt_aa)}{aa_pos+1}*"
            severity = 10
        else:
            mutation_type = 'missense'
            aa_change = f"p.{self._aa_three_letter(wt_aa)}{aa_pos+1}{self._aa_three_letter(mut_aa)}"
            severity = self._calculate_missense_severity(wt_aa, mut_aa)
        
        return {
            'type': mutation_type,
            'wt_aa': wt_aa,
            'mut_aa': mut_aa,
            'wt_codon': wt_codon,
            'mut_codon': mut_codon,
            'aa_change': aa_change,
            'severity_score': severity
        }
    
    def _calculate_missense_severity(self, wt_aa, mut_aa):
        """
        Calcola un punteggio di severità per mutazioni missense
        """
        if wt_aa == '?' or mut_aa == '?' or wt_aa == '*' or mut_aa == '*':
            return 5
        
        wt_props = self.aa_properties.get(wt_aa, self.aa_properties['A'])
        mut_props = self.aa_properties.get(mut_aa, self.aa_properties['A'])
        
        # Calcola differenze nelle proprietà
        hydro_diff = abs(wt_props['hydrophobic'] - mut_props['hydrophobic'])
        size_diff = abs(wt_props['size'] - mut_props['size']) / 100
        charge_diff = abs(wt_props['charge'] - mut_props['charge']) * 2
        polar_diff = abs(wt_props['polar'] - mut_props['polar']) * 2
        
        # Punteggio composito (0-10)
        severity = min(10, hydro_diff + size_diff + charge_diff + polar_diff)
        
        return round(severity, 2)
    
    def _aa_three_letter(self, aa):
        """Converte aminoacido da 1 a 3 lettere"""
        aa_map = {
            'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
            'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
            'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
            'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
            '*': 'Ter', '?': '?'
        }
        return aa_map.get(aa, '?')
    
    def generate_double_mutations(self, gene_name: str, cds_sequence: str, 
                                 max_combinations: int = 1000) -> pd.DataFrame:
        """
        Genera mutazioni doppie (epistasi)
        """
        print(f"\nGenerating double mutations for {gene_name}...")
        
        # Genera prima tutte le singole
        single_muts = self.generate_all_single_nucleotide_mutations(gene_name, cds_sequence)
        
        # Seleziona un sottoinsieme casuale per combinazioni
        if len(single_muts) > 100:
            sample = single_muts.sample(min(100, len(single_muts)))
        else:
            sample = single_muts
        
        double_mutations = []
        count = 0
        
        for idx1, mut1 in sample.iterrows():
            for idx2, mut2 in sample.iterrows():
                if idx1 >= idx2:  # Evita duplicati
                    continue
                if mut1['position'] == mut2['position']:  # Stessa posizione
                    continue
                if count >= max_combinations:
                    break
                
                # Applica entrambe le mutazioni
                temp_seq = list(cds_sequence)
                temp_seq[mut1['position']-1] = mut1['alt_nt']
                temp_seq[mut2['position']-1] = mut2['alt_nt']
                double_seq = ''.join(temp_seq)
                
                try:
                    double_prot = str(Seq(double_seq).translate())
                    wt_prot = str(Seq(cds_sequence).translate())
                    
                    # Conta i cambiamenti aminoacidici
                    aa_changes = sum(1 for i in range(min(len(double_prot), len(wt_prot))) 
                                   if double_prot[i] != wt_prot[i])
                    
                    double_mutations.append({
                        'gene': gene_name,
                        'mutation1': mut1['mutation_id'],
                        'mutation2': mut2['mutation_id'],
                        'type1': mut1['type'],
                        'type2': mut2['type'],
                        'combined_severity': (mut1['severity_score'] + mut2['severity_score']) / 2,
                        'aa_changes_count': aa_changes,
                        'epistasis_potential': 'high' if aa_changes > 2 else 'low'
                    })
                    count += 1
                except:
                    continue
        
        df = pd.DataFrame(double_mutations)
        print(f"  Generated {len(df)} double mutations")
        
        return df


class ClinicalDatabaseIntegrator:
    """
    Integra database di mutazioni cliniche (ClinVar, HbVar, gnomAD)
    """
    
    def __init__(self, cache_dir="cache", vep_data_path=None, reference_fasta=None):
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        self.clinvar_cache = {}
        self.hbvar_data = self._load_hbvar_database()
        
        # Initialize robust genomic mapper if VEP data available
        self.genomic_mapper = None
        if vep_data_path and os.path.exists(vep_data_path):
            try:
                # Import robust mapper
                import sys
                sys.path.append(os.path.dirname(__file__))
                from robust_genomic_mapper import RobustGenomicMapper
                
                self.genomic_mapper = RobustGenomicMapper(
                    vep_data_path=vep_data_path,
                    reference_fasta=reference_fasta
                )
                print("Initialized robust genomic mapper")
            except Exception as e:
                print(f"Could not initialize robust genomic mapper: {e}")
                self.genomic_mapper = None
        else:
            print("Warning: No VEP data provided. Using legacy gnomAD mapping.")
    
    def _load_hbvar_database(self):
        """
        Carica il database HbVar delle emoglobinopatie note
        """
        # Database semplificato delle mutazioni più comuni
        # In produzione, scaricheresti il database completo da http://globin.bx.psu.edu/hbvar/
        
        hbvar_mutations = {
            'HbS': {
                'gene': 'HBB',
                'mutation': 'c.20A>T',
                'protein': 'p.Glu6Val',
                'clinical': 'Sickle cell disease',
                'severity': 'severe',
                'frequency': 'common_in_Africa'
            },
            'HbC': {
                'gene': 'HBB',
                'mutation': 'c.19G>A',
                'protein': 'p.Glu6Lys',
                'clinical': 'HbC disease',
                'severity': 'mild',
                'frequency': 'common_in_West_Africa'
            },
            'HbE': {
                'gene': 'HBB',
                'mutation': 'c.79G>A',
                'protein': 'p.Glu26Lys',
                'clinical': 'HbE disease',
                'severity': 'mild',
                'frequency': 'common_in_Southeast_Asia'
            },
            'Hb_Constant_Spring': {
                'gene': 'HBA2',
                'mutation': 'c.427T>C',
                'protein': 'p.Ter142GlnextTer31',
                'clinical': 'Alpha thalassemia',
                'severity': 'moderate',
                'frequency': 'common_in_Southeast_Asia'
            }
        }
        
        print(f"Loaded {len(hbvar_mutations)} known hemoglobinopathy mutations from HbVar")
        return hbvar_mutations
    
    def query_clinvar(self, gene_name: str, mutation: str) -> Dict:
        """
        Query ClinVar per informazioni cliniche su una mutazione
        """
        cache_key = f"{gene_name}_{mutation}"
        
        # Controlla cache
        if cache_key in self.clinvar_cache:
            return self.clinvar_cache[cache_key]
        
        max_retries = 3
        base_delay = 2.0  # Aumentato da 0.5 a 2 secondi
        
        for attempt in range(max_retries):
            try:
                # Query ClinVar via Entrez
                search_term = f"{gene_name}[gene] AND {mutation}[variant]"
                
                handle = Entrez.esearch(db="clinvar", term=search_term, retmax=10)
                record = Entrez.read(handle)
                handle.close()
                
                if record['IdList']:
                    # Fetch dettagli
                    variant_id = record['IdList'][0]
                    handle = Entrez.efetch(db="clinvar", id=variant_id, rettype="vcv", retmode="xml")
                    xml_data = handle.read()
                    handle.close()
                    
                    # Parse XML
                    parsed = xmltodict.parse(xml_data)
                    
                    # Estrai informazioni rilevanti
                    result = {
                        'found': True,
                        'variant_id': variant_id,
                        'clinical_significance': self._extract_clinical_significance(parsed),
                        'condition': self._extract_condition(parsed),
                        'review_status': self._extract_review_status(parsed)
                    }
                else:
                    result = {
                        'found': False,
                        'variant_id': None,
                        'clinical_significance': 'unknown',
                        'condition': None,
                        'review_status': None
                    }
                
                # Salva in cache
                self.clinvar_cache[cache_key] = result
                
                # Salva cache su disco periodicamente
                if len(self.clinvar_cache) % 100 == 0:
                    self._save_cache()
                
                # Rate limiting con backoff esponenziale
                delay = base_delay * (2 ** attempt)
                time.sleep(delay)
                
                return result
                
            except Exception as e:
                error_str = str(e)
                
                # Gestione specifica per rate limiting
                if "429" in error_str or "Too Many Requests" in error_str:
                    if attempt < max_retries - 1:
                        delay = base_delay * (2 ** attempt) * 2  # Doppio delay per 429
                        print(f"Rate limit hit for {gene_name} {mutation}. Waiting {delay:.1f}s before retry {attempt + 1}/{max_retries}")
                        time.sleep(delay)
                        continue
                    else:
                        print(f"Rate limit persistent for {gene_name} {mutation}. Skipping.")
                        return {
                            'found': False,
                            'error': 'Rate limit exceeded'
                        }
                else:
                    print(f"Error querying ClinVar for {gene_name} {mutation}: {error_str}")
                    return {
                        'found': False,
                        'error': error_str
                    }
        
        # Se tutti i tentativi falliscono
        return {
            'found': False,
            'error': 'Max retries exceeded'
        }
    
    def _extract_clinical_significance(self, parsed_xml):
        """Estrae significato clinico da XML ClinVar con supporto VCV"""
        try:
            # Nuovo formato VCV (VariationClassification)
            if 'VariationClassification' in parsed_xml:
                var_class = parsed_xml['VariationClassification']
                
                # Cerca GermlineClassification dentro VariationClassification
                if 'GermlineClassification' in var_class:
                    germline = var_class['GermlineClassification']
                    if 'Description' in germline:
                        return germline['Description']
                    elif 'ReviewStatus' in germline:
                        return germline['ReviewStatus']
                
                # Cerca classifications dirette
                if 'Classifications' in var_class:
                    classifications = var_class['Classifications']
                    if 'GermlineClassification' in classifications:
                        germline = classifications['GermlineClassification']
                        if isinstance(germline, dict):
                            return germline.get('Description', germline.get('ReviewStatus', 'unknown'))
            
            # Formato legacy RCV
            if 'ClinicalSignificance' in parsed_xml:
                clinical_sig = parsed_xml['ClinicalSignificance']
                if 'Description' in clinical_sig:
                    return clinical_sig['Description']
                elif '@Description' in clinical_sig:
                    return clinical_sig['@Description']
            
            # Formato intermedio
            if 'GermlineClassification' in parsed_xml:
                germline = parsed_xml['GermlineClassification']
                if isinstance(germline, dict):
                    return germline.get('Description', germline.get('@Description', 'unknown'))
                elif isinstance(germline, str):
                    return germline
            
            # Cerca in root level per strutture piatte
            root_keys = ['Description', 'ClinicalSignificance', 'Significance']
            for key in root_keys:
                if key in parsed_xml:
                    value = parsed_xml[key]
                    if isinstance(value, str):
                        return value
                    elif isinstance(value, dict) and 'Description' in value:
                        return value['Description']
            
            return 'unknown'
            
        except Exception as e:
            print(f"Error parsing clinical significance: {e}")
            return 'unknown'
    
    def _extract_condition(self, parsed_xml):
        """Estrae condizione clinica da XML ClinVar con supporto VCV"""
        try:
            # Formato VCV - Cerca in VariationClassification
            if 'VariationClassification' in parsed_xml:
                var_class = parsed_xml['VariationClassification']
                
                # Cerca TraitSet in varie posizioni
                if 'TraitSet' in var_class:
                    trait_set = var_class['TraitSet']
                elif 'Traits' in var_class:
                    trait_set = var_class['Traits']
                else:
                    trait_set = None
                
                if trait_set:
                    if isinstance(trait_set, dict):
                        trait = trait_set.get('Trait', trait_set)
                        if isinstance(trait, dict):
                            return trait.get('Name', trait.get('@Name', 'unknown'))
                        elif isinstance(trait, list) and trait:
                            first_trait = trait[0]
                            if isinstance(first_trait, dict):
                                return first_trait.get('Name', first_trait.get('@Name', 'unknown'))
            
            # Formato legacy RCV
            if 'TraitSet' in parsed_xml:
                trait_set = parsed_xml['TraitSet']
                if isinstance(trait_set, dict):
                    trait = trait_set.get('Trait', {})
                    if isinstance(trait, dict):
                        return trait.get('Name', trait.get('@Name', 'unknown'))
                    elif isinstance(trait, list) and trait:
                        return trait[0].get('Name', trait[0].get('@Name', 'unknown'))
            
            # Cerca nomi alternativi
            condition_keys = ['Condition', 'Phenotype', 'Disease', 'DiseaseName']
            for key in condition_keys:
                if key in parsed_xml:
                    value = parsed_xml[key]
                    if isinstance(value, str):
                        return value
                    elif isinstance(value, dict):
                        name = value.get('Name', value.get('@Name', ''))
                        if name:
                            return name
            
            return 'unknown'
            
        except Exception as e:
            print(f"Error parsing condition: {e}")
            return 'unknown'
    
    def _extract_review_status(self, parsed_xml):
        """Estrae stato di review da XML ClinVar con supporto VCV"""
        try:
            # Formato VCV - Cerca in VariationClassification
            if 'VariationClassification' in parsed_xml:
                var_class = parsed_xml['VariationClassification']
                
                # Cerca ReviewStatus in GermlineClassification
                if 'GermlineClassification' in var_class:
                    germline = var_class['GermlineClassification']
                    if isinstance(germline, dict) and 'ReviewStatus' in germline:
                        return germline['ReviewStatus']
                
                # Cerca ReviewStatus diretto in VariationClassification
                if 'ReviewStatus' in var_class:
                    return var_class['ReviewStatus']
                
                # Cerca in Classifications
                if 'Classifications' in var_class:
                    classifications = var_class['Classifications']
                    if 'GermlineClassification' in classifications:
                        germline = classifications['GermlineClassification']
                        if isinstance(germline, dict) and 'ReviewStatus' in germline:
                            return germline['ReviewStatus']
            
            # Formato legacy RCV
            if 'ReviewStatus' in parsed_xml:
                review_status = parsed_xml['ReviewStatus']
                if isinstance(review_status, str):
                    return review_status
                elif isinstance(review_status, dict):
                    return review_status.get('@Description', review_status.get('Description', 'unknown'))
            
            # Cerca attributi alternativi
            status_keys = ['@ReviewStatus', 'Status', 'Review', 'Evidence']
            for key in status_keys:
                if key in parsed_xml:
                    value = parsed_xml[key]
                    if isinstance(value, str):
                        return value
                    elif isinstance(value, dict) and 'Description' in value:
                        return value['Description']
            
            return 'unknown'
            
        except Exception as e:
            print(f"Error parsing review status: {e}")
            return 'unknown'
    
    def query_gnomad(self, gene_name: str, hgvsp: str, position: int = None, ref: str = None, alt: str = None) -> Dict:
        """
        Query gnomAD using robust genomic mapping
        Uses RobustGenomicMapper if available, falls back to approximation
        """
        # Try to use robust mapper if available
        if hasattr(self, 'genomic_mapper') and self.genomic_mapper:
            try:
                return self.genomic_mapper.query_gnomad_with_mapping(gene_name, hgvsp)
            except Exception as e:
                print(f"Error with robust mapper: {e}")
        
        # Fallback to legacy approximation (deprecated)
        print(f"Warning: Using legacy gnomAD mapping for {gene_name}:{hgvsp}")
        
        # gnomAD API endpoint  
        url = "https://gnomad.broadinstitute.org/api"
        
        query = """
        query VariantInfo($variantId: String!) {
            variant(variantId: $variantId) {
                variantId
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
        """
        
        # Legacy mapping (approximate)
        gene_coords = {
            'HBA1': ('16', 176680),
            'HBA2': ('16', 172847),   
            'HBB': ('11', 5225464),
            'HBG1': ('11', 5274488),
            'HBG2': ('11', 5266047),
        }
        
        if gene_name in gene_coords and position:
            chrom, gene_start = gene_coords[gene_name]
            genomic_pos = gene_start + position
            variant_id = f"{chrom}-{genomic_pos}-{ref or 'A'}-{alt or 'T'}"
        else:
            return {'found': False, 'af_exome': 0, 'af_genome': 0, 'error': 'No mapping available'}
        
        try:
            response = requests.post(
                url,
                json={'query': query, 'variables': {'variantId': variant_id}},
                timeout=15
            )
            
            if response.status_code == 200:
                data = response.json()
                if data.get('data', {}).get('variant'):
                    variant = data['data']['variant']
                    return {
                        'found': True,
                        'variant_id': variant_id,
                        'af_exome': variant.get('exome', {}).get('af', 0),
                        'af_genome': variant.get('genome', {}).get('af', 0),
                        'populations': variant.get('populations', [])
                    }
            
            return {'found': False, 'variant_id': variant_id, 'af_exome': 0, 'af_genome': 0}
            
        except Exception as e:
            print(f"Error querying gnomAD: {str(e)}")
            return {'found': False, 'error': str(e), 'variant_id': variant_id}
    
    def _save_cache(self):
        """Salva cache su disco"""
        cache_file = os.path.join(self.cache_dir, "clinvar_cache.pkl")
        with open(cache_file, 'wb') as f:
            pickle.dump(self.clinvar_cache, f)
    
    def _load_cache(self):
        """Carica cache da disco"""
        cache_file = os.path.join(self.cache_dir, "clinvar_cache.pkl")
        if os.path.exists(cache_file):
            with open(cache_file, 'rb') as f:
                self.clinvar_cache = pickle.load(f)


class MutationEffectPredictor:
    """
    Predice l'effetto funzionale delle mutazioni usando ML e analisi strutturale
    """
    
    def __init__(self):
        self.conservation_scores = {}
        self.structure_features = {}
        self.ml_model = None
        self.scaler = StandardScaler()
    
    def calculate_conservation_score(self, gene_name: str, position: int, 
                                    alignment_file: Optional[str] = None) -> float:
        """
        Calcola score di conservazione evolutiva
        """
        if alignment_file and os.path.exists(alignment_file):
            # Carica allineamento multiplo
            alignment = AlignIO.read(alignment_file, "fasta")
            
            # Calcola conservazione per posizione
            conservation = []
            for i in range(alignment.get_alignment_length()):
                column = alignment[:, i]
                # Shannon entropy
                unique_chars = set(column)
                entropy = 0
                for char in unique_chars:
                    if char != '-':
                        freq = column.count(char) / len(column)
                        if freq > 0:
                            entropy -= freq * np.log2(freq)
                
                # Converti entropia in score di conservazione (0-1)
                max_entropy = np.log2(20)  # 20 aminoacidi
                conservation_score = 1 - (entropy / max_entropy)
                conservation.append(conservation_score)
            
            # Ritorna score per la posizione richiesta
            if position < len(conservation):
                # Convert 1-based position to 0-based array index
                if 1 <= position <= len(conservation):
                    return conservation[position - 1]
                else:
                    return 0.5  # Default conservation score for out-of-bounds
        
        # Default: score medio
        return 0.5
    
    def predict_structural_impact(self, wt_aa: str, mut_aa: str, 
                                 position: int, gene_name: str) -> Dict:
        """
        Predice l'impatto sulla struttura proteica
        """
        # Caratteristiche strutturali basate sulla posizione
        # Per emoglobina, conosciamo alcune regioni critiche
        
        critical_regions = {
            'HBB': {
                'heme_binding': list(range(58, 66)) + list(range(87, 95)),
                'alpha_beta_interface': list(range(30, 45)) + list(range(112, 125)),
                'oxygen_binding': [87, 92],
                'polymerization_site': [6, 85, 88]  # Per HbS
            },
            'HBA1': {
                'heme_binding': list(range(56, 64)) + list(range(85, 93)),
                'alpha_beta_interface': list(range(31, 42)) + list(range(110, 122)),
                'oxygen_binding': [87, 92]
            },
            'HBA2': {
                'heme_binding': list(range(56, 64)) + list(range(85, 93)),
                'alpha_beta_interface': list(range(31, 42)) + list(range(110, 122)),
                'oxygen_binding': [87, 92]
            }
        }
        
        impact_scores = {
            'position_criticality': 0,
            'chemical_change': 0,
            'structural_change': 0,
            'functional_region': 'none'
        }
        
        # Controlla se la posizione è in una regione critica
        if gene_name in critical_regions:
            for region, positions in critical_regions[gene_name].items():
                if position in positions:
                    impact_scores['position_criticality'] = 0.8
                    impact_scores['functional_region'] = region
                    break
        
        # Valuta il cambio chimico
        aa_groups = {
            'hydrophobic': ['A', 'V', 'I', 'L', 'M', 'F', 'W', 'P', 'G'],
            'polar': ['S', 'T', 'Y', 'N', 'Q', 'C'],
            'positive': ['K', 'R', 'H'],
            'negative': ['D', 'E']
        }
        
        wt_group = [k for k, v in aa_groups.items() if wt_aa in v]
        mut_group = [k for k, v in aa_groups.items() if mut_aa in v]
        
        if wt_group != mut_group:
            impact_scores['chemical_change'] = 0.6
            if ('positive' in wt_group and 'negative' in mut_group) or \
               ('negative' in wt_group and 'positive' in mut_group):
                impact_scores['chemical_change'] = 0.9
        
        # Valuta cambio strutturale
        structure_breaking = {
            'P': 0.9,  # Prolina rompe alpha-elica
            'G': 0.7   # Glicina aumenta flessibilità
        }
        
        if mut_aa in structure_breaking:
            impact_scores['structural_change'] = structure_breaking[mut_aa]
        elif wt_aa in structure_breaking and mut_aa not in structure_breaking:
            impact_scores['structural_change'] = 0.8
        
        # Score complessivo
        impact_scores['overall_impact'] = np.mean([
            impact_scores['position_criticality'],
            impact_scores['chemical_change'],
            impact_scores['structural_change']
        ])
        
        return impact_scores
    
  

    def train_ml_model(self, all_generated_mutations: pd.DataFrame):
        """
        Addestra il modello ML usando dati reali da ClinVar come ground truth.
        """
        print("\nTraining machine learning model on REAL data...")
        
        # Carica il nostro ground truth da ClinVar, creato dalla Fase 1
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)  # Up from scripts/ to hemoglobin_analysis/
        training_file = os.path.join(project_root, 'data', 'processed', 'clinvar_training_set.csv')
        
        print(f"Looking for training file at: {training_file}")
        
        if not os.path.exists(training_file):
            print(f"ERRORE: File di addestramento 'clinvar_training_set.csv' non trovato in {training_file}")
            print("Il file dovrebbe essere in: hemoglobin_analysis/data/processed/")
            print("Disponibili alternative:")
            for alt_path in [
                os.path.join(os.getcwd(), 'hemoglobin_analysis', 'data', 'processed', 'clinvar_training_set.csv'),
                os.path.join(os.path.dirname(os.getcwd()), 'hemoglobin_analysis', 'data', 'processed', 'clinvar_training_set.csv')
            ]:
                if os.path.exists(alt_path):
                    print(f"  Found at: {alt_path}")
                    training_file = alt_path
                    break
            else:
                print("  Nessun file trovato. Usando predizioni rule-based.")
                self.ml_model = None
                return None
        
        try:
            ground_truth_df = pd.read_csv(training_file)
            print(f"✅ Caricato file di training con {len(ground_truth_df)} varianti")
        except FileNotFoundError:
            print(f"ERRORE: File di addestramento non trovato: {training_file}")
            self.ml_model = None
            return None

        # Unisci le migliaia di mutazioni che abbiamo generato con le 57 etichette reali
        # Usiamo 'aa_change' (la notazione p.Val6Glu) e il gene come chiave
        training_data = all_generated_mutations.merge(
            ground_truth_df,
            left_on=['gene', 'aa_change'],
            right_on=['gene', 'hgvs_p'],
            how='inner' # 'inner' join per tenere solo le 57 mutazioni che hanno un'etichetta
        )
        
        if len(training_data) < 10: # Se abbiamo meno di 10 esempi, il modello non è affidabile
             print(f"Dati di addestramento insufficienti ({len(training_data)} varianti). Salto l'addestramento.")
             self.ml_model = None
             return None

        print(f"Addestrando il modello su {len(training_data)} varianti con etichette reali da ClinVar.")

        feature_cols = [
            'severity_score', 'conservation_score', 
            'position_criticality', 'chemical_change', 
            'structural_change', 'af_gnomad'
        ]
        
        X = training_data[feature_cols].fillna(0)
        y = training_data['label'] # <-- QUESTA È LA MODIFICA CHE CAMBIA TUTTO!
        
        X_scaled = self.scaler.fit_transform(X)
        
        self.ml_model = RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            random_state=42,
            class_weight='balanced' # Importante per dataset sbilanciati
        )
        self.ml_model.fit(X_scaled, y)
        
        feature_importance = pd.DataFrame({
            'feature': feature_cols,
            'importance': self.ml_model.feature_importances_
        }).sort_values('importance', ascending=False)
        
        print("\nFeature Importance (basata su dati reali):")
        print(feature_importance)
        
        return self.ml_model
    
    def predict_pathogenicity(self, mutation_features: pd.DataFrame) -> np.ndarray:
        """
        Predice patogenicità usando il modello ML
        """
        if self.ml_model is None:
            print("Model not trained. Using rule-based prediction.")
            return (mutation_features['severity_score'] > 5).astype(int).values
        
        feature_cols = [
            'severity_score', 'conservation_score',
            'position_criticality', 'chemical_change',
            'structural_change', 'af_gnomad'
        ]
        
        X = mutation_features[feature_cols].fillna(0)
        X_scaled = self.scaler.transform(X)
        
        # Predici probabilità
        probabilities = self.ml_model.predict_proba(X_scaled)
        
        return probabilities[:, 1]  # Probabilità di essere patogenico


class ComprehensiveMutationAnalyzer:
    """
    Pipeline completa per l'analisi mutazionale sistematica
    
    UTILIZZO PRINCIPALE:
        analyzer = ComprehensiveMutationAnalyzer()
        results = analyzer.run_full_pipeline(genes_dict)
        
    La pipeline esegue:
    1. Generazione di tutte le mutazioni possibili
    2. Addestramento modello ML su dati ClinVar reali
    3. Predizione patogenicità per tutte le mutazioni
    4. Salvataggio risultati e generazione report
    
    IMPORTANTE: Non usare le vecchie funzioni analyze_gene() o batch_analyze()
    che sono state rimosse per evitare logiche errate.
    """
    
    def __init__(self, output_dir="mutation_analysis", vep_data_path=None, reference_fasta=None):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        self.generator = MutationGenerator()
        self.db_integrator = ClinicalDatabaseIntegrator(
            vep_data_path=vep_data_path,
            reference_fasta=reference_fasta
        )
        self.predictor = MutationEffectPredictor()
        
        # GitHub integrations with genomic mapping
        genomic_mapper = getattr(self.db_integrator, 'genomic_mapper', None)
        self.sigprofiler = SigProfilerIntegration(genomic_mapper=genomic_mapper)
        self.hail_integration = None  # Temporaneamente disabilitato
    
    # NOTA: Le vecchie funzioni analyze_gene() e batch_analyze() sono state
    # rimosse perché ridondanti e potenzialmente problematiche.
    # Usare train_and_analyze_all() per l'analisi completa.
    
    def _generate_report(self, gene_name: str, mutations_df: pd.DataFrame):
        """
        Genera report dettagliato dell'analisi
        """
        report_file = f"{self.output_dir}/{gene_name}_report.txt"
        
        with open(report_file, 'w') as f:
            f.write(f"MUTATION ANALYSIS REPORT: {gene_name}\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("="*80 + "\n\n")
            
            # Statistiche generali
            f.write("SUMMARY STATISTICS\n")
            f.write("-"*40 + "\n")
            f.write(f"Total mutations analyzed: {len(mutations_df)}\n")
            f.write(f"Unique positions: {mutations_df['position'].nunique()}\n\n")
            
            # Per tipo di mutazione
            f.write("Mutation types:\n")
            for mut_type, count in mutations_df['type'].value_counts().items():
                percentage = (count / len(mutations_df)) * 100
                f.write(f"  {mut_type}: {count} ({percentage:.1f}%)\n")
            
            # Mutazioni ad alto impatto
            f.write("\n\nHIGH IMPACT MUTATIONS\n")
            f.write("-"*40 + "\n")
            
            high_impact = mutations_df[
                mutations_df['pathogenicity_score'] > 0.8
            ].nlargest(20, 'pathogenicity_score')
            
            for _, row in high_impact.iterrows():
                f.write(f"\n{row['mutation_id']}\n")
                f.write(f"  Type: {row['type']}\n")
                f.write(f"  AA change: {row['aa_change']}\n")
                f.write(f"  Pathogenicity score: {row['pathogenicity_score']:.3f}\n")
                f.write(f"  Functional region: {row.get('functional_region', 'none')}\n")
            
            # Regioni critiche
            f.write("\n\nCRITICAL REGIONS ANALYSIS\n")
            f.write("-"*40 + "\n")
            
            if 'functional_region' in mutations_df.columns:
                region_stats = mutations_df.groupby('functional_region').agg({
                    'pathogenicity_score': 'mean',
                    'mutation_id': 'count'
                }).rename(columns={'mutation_id': 'count'})
                
                region_stats = region_stats.sort_values('pathogenicity_score', ascending=False)
                
                for region, stats in region_stats.iterrows():
                    if region != 'none':
                        f.write(f"\n{region}:\n")
                        f.write(f"  Mutations: {stats['count']}\n")
                        f.write(f"  Avg pathogenicity: {stats['pathogenicity_score']:.3f}\n")
            
            # Hotspots mutazionali
            f.write("\n\nMUTATION HOTSPOTS\n")
            f.write("-"*40 + "\n")
            
            position_stats = mutations_df.groupby('codon_position').agg({
                'pathogenicity_score': 'mean',
                'mutation_id': 'count'
            }).rename(columns={'mutation_id': 'count'})
            
            hotspots = position_stats.nlargest(10, 'pathogenicity_score')
            
            for position, stats in hotspots.iterrows():
                f.write(f"\nCodon {position}:\n")
                f.write(f"  Mutations: {stats['count']}\n")
                f.write(f"  Avg pathogenicity: {stats['pathogenicity_score']:.3f}\n")
        
        print(f"Report saved to: {report_file}")
    
    def run_full_pipeline(self, genes_dict: Dict[str, str]) -> Dict[str, pd.DataFrame]:
        """
        PIPELINE UNIFICATA: Esegue l'intera analisi mutazionale
        
        Passi:
        1. Genera tutte le mutazioni per tutti i geni
        2. Addestra il modello ML UNA SOLA VOLTA usando i dati ClinVar
        3. Applica predizioni a tutti i risultati
        4. Salva file e genera report
        
        Questo approccio evita duplicazioni e garantisce coerenza scientifica
        """
        print("\n" + "="*80)
        print("COMPREHENSIVE ANALYSIS PIPELINE")
        print("="*80)
        
        # STEP 1: Genera tutte le mutazioni per tutti i geni
        print("\n🧬 STEP 1: Generating all mutations...")
        all_gene_results = {}
        all_mutations_for_training = []
        
        for gene_name, cds_sequence in genes_dict.items():
            print(f"  Processing {gene_name}...")
            
            # Genera mutazioni
            gene_mutations = self.generator.generate_all_single_nucleotide_mutations(
                gene_name, cds_sequence
            )
            
            # Aggiungi conservation scores (usa posizione amminoacidica)
            gene_mutations['conservation_score'] = gene_mutations['codon_position'].apply(
                lambda p: self.predictor.calculate_conservation_score(gene_name, p)
            )
            
            # Aggiungi features strutturali
            structural_impacts = []
            for _, row in gene_mutations.iterrows():
                if row['type'] == 'missense':
                    impact = self.predictor.predict_structural_impact(
                        row['wt_aa'], row['mut_aa'], 
                        row['codon_position'], gene_name
                    )
                else:
                    impact = {
                        'position_criticality': 0,
                        'chemical_change': 0,
                        'structural_change': 0,
                        'overall_impact': 0,
                        'functional_region': 'none'
                    }
                structural_impacts.append(impact)
            
            impact_df = pd.DataFrame(structural_impacts)
            gene_mutations = pd.concat([gene_mutations, impact_df], axis=1)
            
            # Aggiungi dati clinici offline 
            gene_mutations['clinvar_found'] = False
            gene_mutations['clinical_significance'] = 'unknown'
            gene_mutations['af_gnomad'] = 0
            
            # Controlla mutazioni note nel database HbVar locale
            for idx, row in gene_mutations.iterrows():
                for variant_name, variant_data in self.db_integrator.hbvar_data.items():
                    if (row['gene'] == variant_data['gene'] and 
                        row['aa_change'] == variant_data['protein']):
                        gene_mutations.at[idx, 'clinvar_found'] = True
                        gene_mutations.at[idx, 'clinical_significance'] = variant_data['severity']
                        freq_map = {'common_in_Africa': 0.1, 'common_in_West_Africa': 0.05, 
                                   'common_in_Southeast_Asia': 0.03}
                        gene_mutations.at[idx, 'af_gnomad'] = freq_map.get(variant_data.get('frequency', ''), 0.001)
            
            # Salva per training e per risultati finali
            all_mutations_for_training.append(gene_mutations)
            all_gene_results[gene_name] = gene_mutations
            
        # STEP 2: Addestra modello ML una volta sola
        print("\n🤖 STEP 2: Training ML model...")
        combined_mutations = pd.concat(all_mutations_for_training, ignore_index=True)
        print(f"  Using {len(combined_mutations)} total mutations for training")
        
        self.predictor.train_ml_model(combined_mutations)
        
        if self.predictor.ml_model is None:
            print("  ⚠️  Using rule-based predictions (no ML model)")
        else:
            print("  ✅ ML model trained successfully!")
        
        # STEP 3: Applica predizioni a tutti i risultati
        print("\n🎯 STEP 3: Applying predictions...")
        for gene_name, gene_mutations in all_gene_results.items():
            print(f"  Predicting for {gene_name}...")
            
            # Predici patogenicità
            gene_mutations['pathogenicity_score'] = self.predictor.predict_pathogenicity(gene_mutations)
            
            # Classifica finale
            gene_mutations['predicted_impact'] = pd.cut(
                gene_mutations['pathogenicity_score'],
                bins=[0, 0.2, 0.5, 0.8, 1.0],
                labels=['benign', 'likely_benign', 'uncertain', 'likely_pathogenic']
            )
            
            # Salva risultati individuali
            output_file = f"{self.output_dir}/{gene_name}_mutation_analysis.csv"
            gene_mutations.to_csv(output_file, index=False)
            
            # Genera report individuale
            self._generate_report(gene_name, gene_mutations)
            
            print(f"    ✅ {len(gene_mutations)} mutations analyzed")
        
        return all_gene_results
    
    
    def _generate_comparative_report(self, results_dict: Dict):
        """
        Genera report comparativo tra geni
        """
        report_file = f"{self.output_dir}/comparative_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("COMPARATIVE MUTATION ANALYSIS REPORT\n")
            f.write("="*80 + "\n\n")
            
            # Confronta statistiche tra geni
            stats_comparison = []
            
            for gene, df in results_dict.items():
                stats_comparison.append({
                    'gene': gene,
                    'total_mutations': len(df),
                    'missense': len(df[df['type'] == 'missense']),
                    'nonsense': len(df[df['type'] == 'nonsense']),
                    'avg_pathogenicity': df['pathogenicity_score'].mean(),
                    'high_impact_count': len(df[df['pathogenicity_score'] > 0.8])
                })
            
            stats_df = pd.DataFrame(stats_comparison)
            
            f.write("GENE COMPARISON\n")
            f.write("-"*40 + "\n")
            f.write(stats_df.to_string(index=False))
            
            # Gene con più mutazioni patogeniche
            f.write("\n\nMOST PATHOGENIC GENES\n")
            f.write("-"*40 + "\n")
            
            for _, row in stats_df.nlargest(5, 'avg_pathogenicity').iterrows():
                f.write(f"\n{row['gene']}:\n")
                f.write(f"  Avg pathogenicity: {row['avg_pathogenicity']:.3f}\n")
                f.write(f"  High impact mutations: {row['high_impact_count']}\n")
        
        print(f"Comparative report saved to: {report_file}")


class SigProfilerIntegration:
    """
    Integrazione con SigProfilerMatrixGenerator per analisi matrici mutazionali
    https://github.com/AlexandrovLab/SigProfilerMatrixGenerator
    """
    
    def __init__(self, output_dir="sigprofiler_output", genomic_mapper=None):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.genomic_mapper = genomic_mapper
        
        if not SIGPROFILER_AVAILABLE:
            print("SigProfilerMatrixGenerator not installed. Skipping matrix analysis.")
            return
        
        # Installa reference genomes se necessario
        try:
            sig_install.install('GRCh38', rsync=False, bash=True)
        except:
            print("Reference genome installation failed. Using existing if available.")
    
    def create_mutation_matrix(self, mutations_df: pd.DataFrame, 
                              genome_build: str = "GRCh38") -> pd.DataFrame:
        """
        Crea matrice mutazionale da dataframe delle mutazioni
        """
        if not SIGPROFILER_AVAILABLE:
            print("SigProfiler not available. Returning empty matrix.")
            return pd.DataFrame()
        
        print("Creating mutational matrix with SigProfiler...")
        
        # Prepara dati nel formato SigProfiler
        # Converte mutazioni in formato VCF-like
        vcf_data = []
        
        for _, row in mutations_df.iterrows():
            if row['type'] in ['missense', 'nonsense']:
                # Use real genomic coordinates if genomic mapper available
                if self.genomic_mapper and ('aa_change' in row or 'hgvs_p' in row):
                    try:
                        # Extract HGVSp notation - prefer aa_change field
                        hgvsp = None
                        if 'aa_change' in row and row['aa_change']:
                            aa_change = str(row['aa_change'])
                            # Clean aa_change: remove "p." prefix if present
                            if aa_change.startswith('p.'):
                                hgvsp = aa_change
                            else:
                                hgvsp = f"p.{aa_change}"
                        elif 'hgvs_p' in row and row['hgvs_p']:
                            hgvsp = str(row['hgvs_p'])
                        
                        if hgvsp:
                            variant_info = self.genomic_mapper.get_genomic_coords(
                                row['gene'], 
                                hgvsp
                            )
                            if variant_info:
                                vcf_entry = {
                                    'CHROM': variant_info.chrom,
                                    'POS': variant_info.pos,
                                    'REF': variant_info.ref,
                                    'ALT': variant_info.alt,
                                    'SAMPLE': row['gene'],
                                    'ID': f"{row['gene']}_{hgvsp.replace('p.', '').replace('(', '').replace(')', '')}"
                                }
                                vcf_data.append(vcf_entry)
                                continue
                    except Exception as e:
                        aa_change = row.get('aa_change', row.get('hgvs_p', ''))
                        print(f"Error getting genomic coords for {row['gene']}:{aa_change}: {e}")
                
                # Fallback to legacy approximation (deprecated) 
                print(f"Warning: Using legacy coordinates for {row['gene']} variant")
                chrom = "11" if "HBB" in row['gene'] else "16"  
                pos = row['position'] * 3 if 'position' in row else 1000  # Approximation
                
                vcf_entry = {
                    'CHROM': chrom,
                    'POS': pos,
                    'REF': row['ref_nt'] if 'ref_nt' in row else 'A',
                    'ALT': row['alt_nt'] if 'alt_nt' in row else 'T', 
                    'SAMPLE': row['gene'],
                    'ID': f"{row['gene']}_pos{pos}"
                }
                vcf_data.append(vcf_entry)
        
        if not vcf_data:
            print("No suitable mutations for matrix generation.")
            return pd.DataFrame()
        
        # Crea file VCF temporaneo con header completo
        vcf_df = pd.DataFrame(vcf_data)
        vcf_file = f"{self.output_dir}/mutations.vcf"
        
        # Scrivi VCF completo con header
        with open(vcf_file, 'w') as f:
            # Header VCF
            f.write("##fileformat=VCFv4.2\n")
            f.write("##reference=GRCh38\n")
            f.write("##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene symbol\">\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            
            # Dati VCF
            for _, row in vcf_df.iterrows():
                line = f"{row['CHROM']}\t{row['POS']}\t.\t{row['REF']}\t{row['ALT']}\t.\tPASS\tGENE={row['SAMPLE']}\tGT\t0/1\n"
                f.write(line)
        
        try:
            # Genera matrice con SigProfiler
            matrices = sig.SigProfilerMatrixGeneratorFunc(
                "test_project", genome_build, self.output_dir,
                plot=True, exome=False
            )
            
            # Carica matrice SBS96
            sbs96_file = f"{self.output_dir}/output/SBS/test_project.SBS96.all"
            if os.path.exists(sbs96_file):
                matrix_df = pd.read_csv(sbs96_file, sep='\t', index_col=0)
                print(f"Mutational matrix created: {matrix_df.shape}")
                return matrix_df
            else:
                print("SBS96 matrix file not found.")
                return pd.DataFrame()
                
        except Exception as e:
            print(f"Error creating mutational matrix: {str(e)}")
            return pd.DataFrame()
    
    def analyze_mutation_signatures(self, matrix_df: pd.DataFrame) -> Dict:
        """
        Analizza signatures mutazionali dalla matrice
        """
        if matrix_df.empty:
            return {}
        
        print("Analyzing mutation signatures...")
        
        # Calcola statistiche di base
        signatures = {
            'total_mutations': matrix_df.sum().sum(),
            'signature_count': len(matrix_df.columns),
            'dominant_signatures': matrix_df.sum().nlargest(5).to_dict(),
            'mutation_types': matrix_df.index.tolist()
        }
        
        # Identifica signatures più comuni nelle emoglobinopatie
        hemoglobin_signatures = {
            'SBS1': 'clock-like aging',
            'SBS2': 'APOBEC activity', 
            'SBS5': 'unknown etiology',
            'SBS13': 'APOBEC activity',
            'SBS18': 'ROS damage'
        }
        
        signatures['hemoglobin_relevant'] = {}
        for sig, desc in hemoglobin_signatures.items():
            if sig in matrix_df.columns:
                contribution = matrix_df[sig].sum()
                signatures['hemoglobin_relevant'][sig] = {
                    'contribution': contribution,
                    'description': desc
                }
        
        return signatures


class HailIntegration:
    """
    Integrazione con Hail per analisi genomiche su larga scala
    https://github.com/hail-is/hail
    """
    
    def __init__(self):
        if not HAIL_AVAILABLE:
            print("Hail not available. Install with: pip install hail")
            return
        
        # Inizializza Hail
        hl.init()
        print("Hail initialized for large-scale genomic analysis")
    
    def create_hail_table(self, mutations_df: pd.DataFrame):
        """
        Converte DataFrame pandas in Hail Table per analisi scalabile
        """
        if not HAIL_AVAILABLE:
            return None
        
        print("Converting to Hail Table...")
        
        # Prepara dati per Hail
        hail_data = []
        for _, row in mutations_df.iterrows():
            hail_data.append({
                'gene': row['gene'],
                'position': row['position'],
                'ref_nt': row.get('ref_nt', 'N'),
                'alt_nt': row.get('alt_nt', 'N'),
                'mutation_type': row['type'],
                'pathogenicity_score': row.get('pathogenicity_score', 0.0),
                'severity_score': row.get('severity_score', 0.0)
            })
        
        # Crea Hail Table
        ht = hl.Table.from_pandas(pd.DataFrame(hail_data))
        
        # Aggiungi annotazioni
        ht = ht.annotate(
            locus = hl.locus(hl.if_else(ht.gene.contains('HBB'), '11', '16'), ht.position),
            alleles = [ht.ref_nt, ht.alt_nt]
        )
        
        print(f"Hail Table created with {ht.count()} variants")
        return ht
    
    def perform_population_analysis(self, ht) -> Dict:
        """
        Analisi popolazionale con Hail
        """
        if not HAIL_AVAILABLE or ht is None:
            return {}
        
        print("Performing population analysis with Hail...")
        
        # Calcola statistiche per gene
        gene_stats = ht.group_by(ht.gene).aggregate(
            total_variants = hl.agg.count(),
            avg_pathogenicity = hl.agg.mean(ht.pathogenicity_score),
            high_impact_count = hl.agg.count_where(ht.pathogenicity_score > 0.8)
        )
        
        # Converti in dict
        stats_dict = {}
        for row in gene_stats.collect():
            stats_dict[row.gene] = {
                'total_variants': row.total_variants,
                'avg_pathogenicity': row.avg_pathogenicity,
                'high_impact_count': row.high_impact_count
            }
        
        return stats_dict
    
    def export_for_genome_browser(self, ht, output_file: str):
        """
        Esporta dati per visualizzazione in genome browser
        """
        if not HAIL_AVAILABLE or ht is None:
            return
        
        print(f"Exporting to genome browser format: {output_file}")
        
        # Esporta in formato BED
        ht.export(f"{output_file}.bed", header=False)


def load_real_gene_sequences(fasta_dir=None):
    """
    Carica automaticamente tutte le sequenze FASTA reali dalla cartella
    """
    gene_sequences = {}
    
    # Se non specificato, usa percorso relativo al progetto
    if fasta_dir is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)  # Up from scripts/ to hemoglobin_analysis/
        fasta_dir = os.path.join(project_root, 'data', 'hemoglobin_genes', 'fasta')
    
    if not os.path.exists(fasta_dir):
        print(f"Directory {fasta_dir} not found. Using hardcoded sequences.")
        return get_hardcoded_sequences()
    
    # Leggi tutti i file FASTA nella directory
    for filename in os.listdir(fasta_dir):
        if filename.endswith('_cds.fasta'):
            gene_name = filename.replace('_cds.fasta', '')
            filepath = os.path.join(fasta_dir, filename)
            
            try:
                # Leggi il file FASTA
                with open(filepath, 'r') as f:
                    lines = f.readlines()
                
                # Estrai la sequenza (salta l'header)
                sequence = ''.join(line.strip() for line in lines[1:] if not line.startswith('>'))
                sequence = sequence.replace('\n', '').replace(' ', '')
                
                if sequence:  # Solo se la sequenza non è vuota
                    gene_sequences[gene_name] = sequence
                    print(f"Loaded {gene_name}: {len(sequence)} bp")
                    
            except Exception as e:
                print(f"Error loading {filename}: {str(e)}")
    
    if not gene_sequences:
        print("No sequences loaded. Using hardcoded sequences.")
        return get_hardcoded_sequences()
    
    print(f"Successfully loaded {len(gene_sequences)} gene sequences")
    return gene_sequences

def get_hardcoded_sequences():
    """
    Fallback: sequenze hardcoded per i geni principali
    """
    return {
        'HBB': 'ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA',
        'HBA1': 'ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGGATGTTCCTGTCCTTCCCCACCACCAAGACCTACTTCCCGCACTTCGACCTGAGCCACGGCTCTGCCCAGGTTAAGGGCCACGGCAAGAAGGTGGCCGACGCGCTGACCAACGCCGTGGCGCACGTGGACGACATGCCCAACGCGCTGTCCGCCCTGAGCGACCTGCACGCGCACAAGCTTCGGGTGGACCCGGTCAACTTCAAGCTCCTAAGCCACTGCCTGCTGGTGACCCTGGCCGCCCACCTCCCCGCCGAGTTCACCCCTGCGGTGCACGCCTCCCTGGACAAGTTCCTGGCTTCTGTGAGCACCGTGCTGACCTCCAAATACCGTTAA',
        'HBA2': 'ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGGATGTTCCTGTCCTTCCCCACCACCAAGACCTACTTCCCGCACTTCGACCTGAGCCACGGCTCTGCCCAGGTTAAGGGCCACGGCAAGAAGGTGGCCGACGCGCTGACCAACGCCGTGGCGCACGTGGACGACATGCCCAACGCGCTGTCCGCCCTGAGCGACCTGCACGCGCACAAGCTTCGGGTGGACCCGGTCAACTTCAAGCTCCTAAGCCACTGCCTGCTGGTGACCCTGGCCGCCCACCTCCCCGCCGAGTTCACCCCTGCGGTGCACGCCTCCCTGGACAAGTTCCTGGCTTCTGTGAGCACCGTGCTGACCTCCAAATACCGTTAA',
        'GATA1': 'ATGGAGTTCCCTGGCCTGGGGTCCCTGGGGACCTCAGAGCCCCTCCCCCAGTTTGTGGATCCTGCTCTGGTGTCCTCCACACCAGAATCAGGGGTTTTCTTCCCCTCTGGGCCTGAGGGCTTGGATGCAGCAGCTTCCTCCACTGCCCCGAGCACAGCCACCGCTGCAGCTGCGGCACTGGCCTACTACAGGGACGCTGAGGCCTACAGACACTCCCCAGTCTTTCAGGTGTACCCATTGCTCAACTGTATGGAGGGGATCCCAGGGGGCTCACCATATGCCGGCTGGGCCTACGGCAAGACGGGGCTCTACCCTGCCTCAACTGTGTGTCCCACCCGCGAGGACTCTCCTCCCCAGGCCGTGGAAGATCTGGATGGAAAAGGCAGCACCAGCTTCCTGGAGACTTTGAAGACAGAGCGGCTGAGCCCAGACCTCCTGACCCTGGGACCTGCACTGCCTTCATCACTCCCTGTCCCCAATAGTGCTTATGGGGGCCCTGACTTTTCCAGTACCTTCTTTTCTCCCACCGGGAGCCCCCTC',
        'KLF1': 'ATGGCCACAGCCGAGACCGCCTTGCCCTCCATCAGCACACTGACCGCCCTGGGCCCCTTCCCGGACACACAGGATGACTTCCTCAAGTGGTGGCGCTCCGAAGAGGCGCAGGACATGGGCCCGGGTCCTCCTGACCCCACGGAGCCGCCCCTCCACGTGAAGTCTGAGGACCAGCCCGGGGAGGAAGAGGACGATGAGAGGGGCGCGGACGCCACCTGGGACCTGGATCTCCTCCTCACCAACTTCTCGGGCCCGGAGCCCGGTGGCGCGCCCCAGACCTGCGCTCTGGCGCCCAGCGAGGCCTCCGGGGCGCAATATCCGCCGCCGCCCGAGACTCTGGGCGCATATGCTGGCGGCCCGGGGCTGGTGGCTGGGCTTTTGGGTTCGGAGGATCACTCGGGTTGGGTGCGCCCTGCCCTGCGAGCCCGGGCTCCCGACGCCTTCGTGGGCCCAGCCCTGGCTCCAGCCCCGGCCCCCGAGCCCAAGGCGCTGGCGCTGCAACCGGTGTACCCGGGGCCCGGCGCCGGCTCCTCGGGTGGC'
    }

def setup_reference_fasta(fasta_path: str) -> Optional[str]:
    """
    Setup reference FASTA: decompress if needed and create index
    Returns path to usable FASTA or None if failed
    """
    import gzip
    import subprocess
    
    if not fasta_path:
        return None
    
    # Check for compressed version first
    gz_path = fasta_path + '.gz' if not fasta_path.endswith('.gz') else fasta_path
    base_path = fasta_path[:-3] if fasta_path.endswith('.gz') else fasta_path
    
    # Try to find the FASTA file (compressed or uncompressed)
    if os.path.exists(gz_path) and not os.path.exists(base_path):
        print(f"📦 Found compressed FASTA: {gz_path}")
        print(f"🔓 Decompressing to: {base_path}")
        
        try:
            # Decompress the FASTA file
            with gzip.open(gz_path, 'rb') as f_in:
                with open(base_path, 'wb') as f_out:
                    f_out.write(f_in.read())
            print(f"✅ Decompression complete")
        except Exception as e:
            print(f"❌ Error decompressing FASTA: {e}")
            return None
    
    elif os.path.exists(base_path):
        print(f"📁 Found uncompressed FASTA: {base_path}")
    else:
        print(f"❌ FASTA not found at either {base_path} or {gz_path}")
        return None
    
    # Check if index exists, create if not
    fai_path = base_path + '.fai'
    if not os.path.exists(fai_path):
        print(f"🔨 Creating FASTA index: {fai_path}")
        try:
            # Use samtools faidx or pysam to create index
            try:
                import pysam
                pysam.faidx(base_path)
                print(f"✅ Index created with pysam")
            except:
                # Fallback to samtools command line
                result = subprocess.run(['samtools', 'faidx', base_path], 
                                      capture_output=True, text=True)
                if result.returncode == 0:
                    print(f"✅ Index created with samtools")
                else:
                    print(f"❌ Error creating index: {result.stderr}")
                    return None
        except Exception as e:
            print(f"❌ Error creating FASTA index: {e}")
            print(f"💡 Please run: samtools faidx {base_path}")
            return None
    else:
        print(f"✅ FASTA index already exists: {fai_path}")
    
    # Validate the FASTA can be opened with pysam
    try:
        import pysam
        with pysam.FastaFile(base_path) as fasta:
            refs = fasta.references
            print(f"✅ FASTA validated: {len(refs)} chromosomes/contigs")
            # Show sample of chromosome names for verification
            sample_refs = list(refs)[:5]
            print(f"📋 Sample references: {sample_refs}")
            
            # Check for common chromosome naming issues
            has_chr_prefix = any(ref.startswith('chr') for ref in refs)
            print(f"🏷️ Chromosome naming: {'chr-prefixed' if has_chr_prefix else 'no-prefix'}")
            
        return base_path
    except Exception as e:
        print(f"❌ Cannot open FASTA with pysam: {e}")
        return None


def main():
    """
    Pipeline principale
    """
    print("╔" + "═" * 78 + "╗")
    print("║" + " COMPREHENSIVE HEMOGLOBIN MUTATION ANALYZER ".center(78) + "║")
    print("║" + " Systematic Analysis of All Possible Mutations ".center(78) + "║")
    print("╚" + "═" * 78 + "╝")
    
    # Carica sequenze reali dalla cartella
    real_genes = load_real_gene_sequences()
    
    print(f"Caricate {len(real_genes)} sequenze geniche reali")
    
    # Usa tutti i geni disponibili per l'analisi completa
    selected_genes = real_genes
    
    if not selected_genes:
        print("ERRORE: Nessun gene trovato!")
        return
    
    print(f"Analisi di {len(selected_genes)} geni disponibili: {list(selected_genes.keys())}")
    
    # Configura percorsi per mapping genomico robusto
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    vep_data_path = os.path.join(project_root, 'data', 'processed', 'vep_extracted.tsv')
    
    # Setup reference FASTA with automatic decompression and indexing
    reference_fasta_path = os.path.join(project_root, 'data', 'reference', 'GRCh38.fa')
    print(f"\n🧬 Reference Genome Setup:")
    reference_fasta = setup_reference_fasta(reference_fasta_path)
    
    if reference_fasta:
        print(f"  ✅ Reference FASTA ready: {reference_fasta}")
    else:
        print(f"  ⚠️ No reference FASTA available")
        print(f"  💡 Download GRCh38 and place at: {reference_fasta_path}")
        print(f"  📥 wget -O {reference_fasta_path}.gz https://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
    
    # Verifica disponibilità dati VEP
    has_vep_data = os.path.exists(vep_data_path)
    print(f"\n🗺️ Genomic Mapping Configuration:")
    print(f"  - VEP data: {'✅ Available' if has_vep_data else '❌ Missing'} ({vep_data_path})")
    print(f"  - Reference FASTA: {'✅ Available' if reference_fasta and os.path.exists(reference_fasta) else '❌ Not provided'}")
    
    if has_vep_data:
        print(f"  ✅ Using robust genomic mapping with MANE transcript priority")
        print(f"  ✅ Real gnomAD variant IDs will be generated")
        print(f"  ✅ VEP REST API fallback available for new variants")
    else:
        print(f"  ⚠️ Falling back to legacy coordinate approximation")
        print(f"  💡 Place VEP data at {vep_data_path} for full functionality")
    
    # Inizializza analyzer con mapping genomico robusto
    analyzer = ComprehensiveMutationAnalyzer(
        vep_data_path=vep_data_path if has_vep_data else None,
        reference_fasta=reference_fasta
    )
    
    # PIPELINE UNIFICATA: genera, addestra, analizza tutto in una volta
    all_results = analyzer.run_full_pipeline(selected_genes)
    
    # Combina tutti i risultati
    combined_results = pd.concat(all_results.values(), ignore_index=True)
    
    # Salva risultati combinati
    combined_file = f"{analyzer.output_dir}/all_genes_mutation_analysis.csv"
    combined_results.to_csv(combined_file, index=False)
    print(f"\n💾 Combined results saved to: {combined_file}")
    
    # Genera report comparativo
    analyzer._generate_comparative_report(all_results)
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print(f"\nTotal mutations analyzed: {len(combined_results)}")
    print(f"Genes analyzed: {list(all_results.keys())}")
    print(f"Output directory: {analyzer.output_dir}/")
    print("\nGenerated files:")
    print("  - Individual gene analyses: [gene]_mutation_analysis.csv")
    print("  - Combined analysis: all_genes_mutation_analysis.csv")
    print("  - Reports: [gene]_report.txt, comparative_report.txt")
    print("\nNext steps:")
    print("  1. Review high-impact mutations in reports")
    print("  2. Validate predictions with experimental data")
    print("  3. Cross-reference with clinical databases")
    print("  4. Design targeted experiments for validation")


if __name__ == "__main__":
    main()