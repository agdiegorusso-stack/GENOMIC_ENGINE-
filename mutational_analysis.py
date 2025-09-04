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
from Bio.SubsMat import MatrixInfo
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler

# Configurazione
Entrez.email = "your.email@example.com"  # MODIFICA CON LA TUA EMAIL
warnings.filterwarnings('ignore')

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
    
    def __init__(self, cache_dir="cache"):
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        self.clinvar_cache = {}
        self.hbvar_data = self._load_hbvar_database()
    
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
            
            time.sleep(0.5)  # Rate limiting
            
            return result
            
        except Exception as e:
            print(f"Error querying ClinVar for {gene_name} {mutation}: {str(e)}")
            return {
                'found': False,
                'error': str(e)
            }
    
    def _extract_clinical_significance(self, parsed_xml):
        """Estrae significato clinico da XML ClinVar"""
        try:
            # Naviga nella struttura XML (può variare)
            return parsed_xml.get('ClinicalSignificance', {}).get('Description', 'unknown')
        except:
            return 'unknown'
    
    def _extract_condition(self, parsed_xml):
        """Estrae condizione clinica da XML ClinVar"""
        try:
            return parsed_xml.get('TraitSet', {}).get('Trait', {}).get('Name', 'unknown')
        except:
            return 'unknown'
    
    def _extract_review_status(self, parsed_xml):
        """Estrae stato di review da XML ClinVar"""
        try:
            return parsed_xml.get('ReviewStatus', 'unknown')
        except:
            return 'unknown'
    
    def query_gnomad(self, gene_name: str, position: int, ref: str, alt: str) -> Dict:
        """
        Query gnomAD per frequenze alleliche popolazione
        """
        # gnomAD API endpoint
        url = f"https://gnomad.broadinstitute.org/api"
        
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
        
        # Costruisci variant ID nel formato gnomAD
        # Questo è semplificato - in produzione dovresti mappare alle coordinate genomiche
        variant_id = f"{gene_name}-{position}-{ref}-{alt}"
        
        try:
            response = requests.post(
                url,
                json={'query': query, 'variables': {'variantId': variant_id}},
                timeout=10
            )
            
            if response.status_code == 200:
                data = response.json()
                if data.get('data', {}).get('variant'):
                    variant = data['data']['variant']
                    return {
                        'found': True,
                        'af_exome': variant.get('exome', {}).get('af', 0),
                        'af_genome': variant.get('genome', {}).get('af', 0),
                        'populations': variant.get('populations', [])
                    }
            
            return {'found': False, 'af_exome': 0, 'af_genome': 0}
            
        except Exception as e:
            print(f"Error querying gnomAD: {str(e)}")
            return {'found': False, 'error': str(e)}
    
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
                return conservation[position]
        
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
    
    def train_ml_model(self, training_data: pd.DataFrame):
        """
        Addestra modello ML per predizione patogenicità
        """
        print("\nTraining machine learning model...")
        
        # Prepara features
        feature_cols = [
            'severity_score', 'conservation_score', 
            'position_criticality', 'chemical_change', 
            'structural_change', 'af_gnomad'
        ]
        
        # Simula labels (in produzione useresti dati reali da ClinVar)
        # 0 = benigno, 1 = patogenico
        training_data['pathogenic'] = (
            training_data['severity_score'] > 5
        ).astype(int)
        
        X = training_data[feature_cols].fillna(0)
        y = training_data['pathogenic']
        
        # Normalizza features
        X_scaled = self.scaler.fit_transform(X)
        
        # Addestra Random Forest
        self.ml_model = RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            random_state=42
        )
        self.ml_model.fit(X_scaled, y)
        
        # Feature importance
        feature_importance = pd.DataFrame({
            'feature': feature_cols,
            'importance': self.ml_model.feature_importances_
        }).sort_values('importance', ascending=False)
        
        print("Feature importance:")
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
    Pipeline completa per l'analisi mutazionale
    """
    
    def __init__(self, output_dir="mutation_analysis"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        self.generator = MutationGenerator()
        self.db_integrator = ClinicalDatabaseIntegrator()
        self.predictor = MutationEffectPredictor()
    
    def analyze_gene(self, gene_name: str, cds_sequence: str, 
                    include_clinical: bool = True,
                    include_double: bool = False) -> pd.DataFrame:
        """
        Analisi completa di tutte le mutazioni per un gene
        """
        print(f"\n{'='*80}")
        print(f"COMPREHENSIVE MUTATION ANALYSIS: {gene_name}")
        print(f"{'='*80}")
        
        # 1. Genera tutte le mutazioni singole
        mutations_df = self.generator.generate_all_single_nucleotide_mutations(
            gene_name, cds_sequence
        )
        
        # 2. Aggiungi score di conservazione (simulato per questo esempio)
        print("\nCalculating conservation scores...")
        mutations_df['conservation_score'] = mutations_df['position'].apply(
            lambda p: self.predictor.calculate_conservation_score(gene_name, p)
        )
        
        # 3. Predici impatto strutturale
        print("Predicting structural impact...")
        structural_impacts = []
        for _, row in mutations_df.iterrows():
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
        
        # Aggiungi colonne impatto strutturale
        impact_df = pd.DataFrame(structural_impacts)
        mutations_df = pd.concat([mutations_df, impact_df], axis=1)
        
        # 4. Integra database clinici (opzionale)
        if include_clinical:
            print("\nQuerying clinical databases...")
            clinical_data = []
            
            # Campiona per evitare troppe query
            sample_size = min(100, len(mutations_df))
            sample = mutations_df.sample(sample_size)
            
            for _, row in sample.iterrows():
                # Query ClinVar
                clinvar_result = self.db_integrator.query_clinvar(
                    gene_name, row['aa_change']
                )
                
                # Query gnomAD (simulato)
                gnomad_result = {'af_exome': np.random.exponential(0.0001)}
                
                clinical_data.append({
                    'mutation_id': row['mutation_id'],
                    'clinvar_found': clinvar_result.get('found', False),
                    'clinical_significance': clinvar_result.get('clinical_significance', 'unknown'),
                    'af_gnomad': gnomad_result.get('af_exome', 0)
                })
            
            clinical_df = pd.DataFrame(clinical_data)
            mutations_df = mutations_df.merge(
                clinical_df, on='mutation_id', how='left'
            )
            mutations_df['af_gnomad'].fillna(0, inplace=True)
        else:
            mutations_df['af_gnomad'] = 0
        
        # 5. Predici patogenicità con ML
        print("\nPredicting pathogenicity...")
        
        # Prima addestra il modello su un subset
        training_sample = mutations_df.sample(min(1000, len(mutations_df)))
        self.predictor.train_ml_model(training_sample)
        
        # Poi predici su tutto il dataset
        mutations_df['pathogenicity_score'] = self.predictor.predict_pathogenicity(mutations_df)
        
        # 6. Classifica finale
        mutations_df['predicted_impact'] = pd.cut(
            mutations_df['pathogenicity_score'],
            bins=[0, 0.2, 0.5, 0.8, 1.0],
            labels=['benign', 'likely_benign', 'uncertain', 'likely_pathogenic']
        )
        
        # 7. Analizza mutazioni doppie (opzionale)
        if include_double:
            double_muts = self.generator.generate_double_mutations(
                gene_name, cds_sequence, max_combinations=100
            )
            double_file = f"{self.output_dir}/{gene_name}_double_mutations.csv"
            double_muts.to_csv(double_file, index=False)
            print(f"\nDouble mutations saved to: {double_file}")
        
        # 8. Salva risultati
        output_file = f"{self.output_dir}/{gene_name}_mutation_analysis.csv"
        mutations_df.to_csv(output_file, index=False)
        print(f"\nResults saved to: {output_file}")
        
        # 9. Genera report
        self._generate_report(gene_name, mutations_df)
        
        return mutations_df
    
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
    
    def batch_analyze(self, genes_dict: Dict[str, str], parallel: bool = True):
        """
        Analizza multipli geni in batch
        """
        print(f"\nBATCH ANALYSIS: {len(genes_dict)} genes")
        print("="*80)
        
        all_results = {}
        
        if parallel and len(genes_dict) > 1:
            # Analisi parallela
            with ThreadPoolExecutor(max_workers=4) as executor:
                futures = {
                    executor.submit(self.analyze_gene, gene, seq): gene
                    for gene, seq in genes_dict.items()
                }
                
                for future in as_completed(futures):
                    gene = futures[future]
                    try:
                        result = future.result()
                        all_results[gene] = result
                    except Exception as e:
                        print(f"Error analyzing {gene}: {str(e)}")
        else:
            # Analisi sequenziale
            for gene, seq in genes_dict.items():
                try:
                    result = self.analyze_gene(gene, seq)
                    all_results[gene] = result
                except Exception as e:
                    print(f"Error analyzing {gene}: {str(e)}")
        
        # Combina risultati
        combined_df = pd.concat(all_results.values(), ignore_index=True)
        
        # Salva risultati combinati
        combined_file = f"{self.output_dir}/all_genes_mutation_analysis.csv"
        combined_df.to_csv(combined_file, index=False)
        
        # Genera report comparativo
        self._generate_comparative_report(all_results)
        
        return combined_df
    
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


def main():
    """
    Pipeline principale
    """
    print("╔" + "═" * 78 + "╗")
    print("║" + " COMPREHENSIVE HEMOGLOBIN MUTATION ANALYZER ".center(78) + "║")
    print("║" + " Systematic Analysis of All Possible Mutations ".center(78) + "║")
    print("╚" + "═" * 78 + "╝")
    
    # Esempio di utilizzo
    # Assumendo che hai già scaricato le sequenze con lo script precedente
    
    # Sequenze di esempio (sostituisci con quelle reali)
    example_genes = {
        'HBB': 'ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGG',
        'HBA1': 'ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGG'
    }
    
    # Inizializza analyzer
    analyzer = ComprehensiveMutationAnalyzer()
    
    # Analizza tutti i geni
    results = analyzer.batch_analyze(example_genes, parallel=True)
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print(f"\nTotal mutations analyzed: {len(results)}")
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