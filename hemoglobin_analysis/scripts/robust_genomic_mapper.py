#!/usr/bin/env python3
"""
Robust Genomic Mapper for HGVS annotations
Handles VEP CSQ parsing, MANE transcript priority, variant normalization
"""

import re
import pandas as pd
import json
import os
import time
from typing import Dict, List, Tuple, Optional, Set
import requests
import pysam
from dataclasses import dataclass


def normalize_hgvsp(hgvsp: str) -> str:
    """
    Normalize HGVSp notation for consistent lookup
    Handles parentheses variations: p.Glu6Val vs p.(Glu6Val)
    """
    if not hgvsp:
        return hgvsp
    
    # Remove leading/trailing whitespace
    hgvsp = hgvsp.strip()
    
    # Skip if not a protein notation
    if not hgvsp.startswith('p.'):
        return hgvsp
    
    # Normalize parentheses - prefer without parentheses for compatibility
    # p.(Glu6Val) -> p.Glu6Val
    if hgvsp.startswith('p.(') and hgvsp.endswith(')'):
        # Extract content within parentheses
        content = hgvsp[3:-1]  # Remove 'p.(' and ')'
        hgvsp = f"p.{content}"
    
    # Handle special cases and common variations
    normalizations = {
        # Common amino acid name variations
        'Ter': '*',  # Termination codon
        'Stop': '*',
        'X': '*',
        
        # Handle insertion/deletion notation
        'ins': 'ins',
        'del': 'del',
        'dup': 'dup',
        'delins': 'delins',
        
        # Frameshifts
        'fs': 'fs',
        'frameshift': 'fs',
    }
    
    # Apply normalizations
    for old, new in normalizations.items():
        hgvsp = hgvsp.replace(old, new)
    
    # Additional cleanup: remove extra spaces around operators
    hgvsp = re.sub(r'\s*([>_])\s*', r'\1', hgvsp)
    
    return hgvsp


def hgvsp_variants(hgvsp: str) -> List[str]:
    """
    Generate common HGVSp notation variants for lookup
    Returns list of possible representations
    """
    variants = [hgvsp]
    normalized = normalize_hgvsp(hgvsp)
    
    if normalized != hgvsp:
        variants.append(normalized)
    
    # Generate parentheses variants
    if normalized.startswith('p.') and not normalized.startswith('p.('):
        # Add parentheses version: p.Glu6Val -> p.(Glu6Val)
        content = normalized[2:]  # Remove 'p.'
        variants.append(f"p.({content})")
    
    # Remove duplicates while preserving order
    seen = set()
    unique_variants = []
    for v in variants:
        if v not in seen:
            seen.add(v)
            unique_variants.append(v)
    
    return unique_variants


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
        self.session.close()


@dataclass
class VariantInfo:
    """Structured variant information"""
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    transcript: str
    hgvsc: str = ""
    hgvsp: str = ""
    allele_num: int = 0
    is_mane: bool = False
    mane_type: str = ""  # "MANE Select", "MANE Plus Clinical"


class CSQParser:
    """Parse VEP CSQ header and annotations"""
    
    def __init__(self, vcf_path: str = None, tsv_path: str = None):
        self.field_map = {}
        self.csq_format = ""
        
        if vcf_path:
            self._parse_vcf_header(vcf_path)
        elif tsv_path:
            # Try to parse TSV header for CSQ format
            self._parse_tsv_header(tsv_path)
    
    def _parse_vcf_header(self, vcf_path: str):
        """Parse CSQ format from VCF header with robust field detection"""
        try:
            with pysam.VariantFile(vcf_path) as vcf:
                if 'CSQ' in vcf.header.info:
                    desc = vcf.header.info['CSQ'].description
                    # Try multiple patterns to extract Format string
                    format_patterns = [
                        r'Format:\s*([^"]+?)(?:\s*"|$)',  # Standard VEP format
                        r'Format:\s*([^"]+)',             # Without closing quote
                        r'Consequence\s+annotations.*?Format:\s*([^"]+)',  # Extended description
                        r'Format=([^"]+)'                 # Alternative format
                    ]
                    
                    fields = None
                    for pattern in format_patterns:
                        match = re.search(pattern, desc, re.IGNORECASE | re.DOTALL)
                        if match:
                            self.csq_format = match.group(1).strip()
                            fields = [f.strip() for f in self.csq_format.split('|')]
                            break
                    
                    if fields:
                        self.field_map = {field: idx for idx, field in enumerate(fields)}
                        print(f"✅ Parsed CSQ format: {len(fields)} fields")
                        
                        # Check for critical fields and report availability
                        critical_fields = ['Allele', 'ALLELE_NUM', 'HGVSp', 'HGVSc', 'SYMBOL', 'MANE_SELECT']
                        available = [f for f in critical_fields if f in self.field_map]
                        missing = [f for f in critical_fields if f not in self.field_map]
                        
                        print(f"🔍 Critical fields available: {available}")
                        if missing:
                            print(f"⚠️  Missing fields: {missing}")
                        
                        # Special handling for ALLELE_NUM
                        if 'ALLELE_NUM' not in self.field_map:
                            print(f"💡 ALLELE_NUM not available - will use Allele field for filtering")
                        else:
                            print(f"✅ ALLELE_NUM available at position {self.field_map['ALLELE_NUM']}")
                            
                    else:
                        print("❌ Could not parse CSQ Format string")
                        self._set_default_format()
                else:
                    print("❌ No CSQ field found in VCF header")
                    self._set_default_format()
        except Exception as e:
            print(f"❌ Error parsing VCF header: {e}")
            self._set_default_format()
    
    def _parse_tsv_header(self, tsv_path: str):
        """Parse CSQ format from TSV file header or first data line with enhanced detection"""
        try:
            with open(tsv_path, 'r') as f:
                # Read first few lines to look for CSQ format info
                lines = [f.readline().strip() for _ in range(20)]  # Increased to 20 lines
                
            # Look for lines starting with # that might contain format info
            format_info = None
            for line in lines:
                if line.startswith('#') and ('CSQ' in line or 'Format' in line):
                    # Extract format information with multiple patterns
                    patterns = [
                        r'Format[:\s]*([^"]+?)(?:\s*"|$)',  # Standard
                        r'Format[:\s]*([^"]+)',             # Without quote
                        r'CSQ[:\s]*([^"]+)',                # Direct CSQ
                        r'Format=([^"]+)'                   # Alternative
                    ]
                    
                    for pattern in patterns:
                        match = re.search(pattern, line, re.IGNORECASE)
                        if match:
                            format_info = match.group(1).strip()
                            break
                    
                    if format_info:
                        break
            
            if format_info:
                fields = [f.strip() for f in format_info.split('|')]
                self.field_map = {field: idx for idx, field in enumerate(fields)}
                self.csq_format = format_info
                print(f"✅ Parsed TSV CSQ format: {len(fields)} fields")
                
                # Enhanced field verification
                required_fields = ['Allele', 'Consequence', 'SYMBOL']
                optional_fields = ['ALLELE_NUM', 'HGVSp', 'HGVSc', 'MANE_SELECT', 'CANONICAL']
                
                available_required = [f for f in required_fields if f in self.field_map]
                available_optional = [f for f in optional_fields if f in self.field_map]
                missing_required = [f for f in required_fields if f not in self.field_map]
                
                print(f"🔍 Required fields: {available_required}")
                if available_optional:
                    print(f"🔍 Optional fields: {available_optional}")
                if missing_required:
                    print(f"⚠️  Missing required fields: {missing_required}")
                    print("💡 Consider regenerating VEP output with --fields to include missing fields")
                    
            else:
                print("⚠️  No CSQ format found in TSV - using default")
                self._set_default_format()
                
        except Exception as e:
            print(f"❌ Error parsing TSV header: {e}")
            self._set_default_format()
    
    def _set_default_format(self):
        """Set standard VEP format as fallback"""
        # Standard VEP output format
        fields = [
            'Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type',
            'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position',
            'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation',
            'DISTANCE', 'STRAND', 'FLAGS', 'VARIANT_CLASS', 'SYMBOL_SOURCE', 'HGNC_ID',
            'CANONICAL', 'MANE_SELECT', 'MANE_PLUS_CLINICAL', 'TSL', 'APPRIS', 'CCDS',
            'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'UNIPROT_ISOFORM', 'GENE_PHENO',
            'SIFT', 'PolyPhen', 'DOMAINS', 'miRNA', 'HGVS_OFFSET', 'AF', 'AFR_AF',
            'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'AA_AF', 'EA_AF', 'gnomAD_AF',
            'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF',
            'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF',
            'MAX_AF', 'MAX_AF_POPS', 'CLIN_SIG', 'SOMATIC', 'PHENO', 'PUBMED', 'MOTIF_NAME',
            'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'TRANSCRIPTION_FACTORS'
        ]
        self.field_map = {field: idx for idx, field in enumerate(fields)}
        self.csq_format = '|'.join(fields)
        print(f"Using default CSQ format: {len(fields)} fields")
    
    def parse_csq_annotation(self, csq_string: str, target_allele: str = None, target_allele_index: int = None) -> List[Dict]:
        """Parse single CSQ annotation string into structured data with ALT-specific filtering"""
        annotations = []
        
        for transcript in csq_string.split(','):
            cols = transcript.split('|')
            if len(cols) < max(self.field_map.values()) + 1:
                continue
                
            # Extract key fields with safe indexing
            annotation = {
                'allele': cols[self.field_map.get('Allele', 0)] if 'Allele' in self.field_map and len(cols) > self.field_map.get('Allele', 0) else '',
                'consequence': cols[self.field_map.get('Consequence', 1)] if 'Consequence' in self.field_map and len(cols) > self.field_map.get('Consequence', 1) else '',
                'symbol': cols[self.field_map.get('SYMBOL', 3)] if 'SYMBOL' in self.field_map and len(cols) > self.field_map.get('SYMBOL', 3) else '',
                'gene': cols[self.field_map.get('Gene', 4)] if 'Gene' in self.field_map and len(cols) > self.field_map.get('Gene', 4) else '',
                'feature': cols[self.field_map.get('Feature', 6)] if 'Feature' in self.field_map and len(cols) > self.field_map.get('Feature', 6) else '',
                'biotype': cols[self.field_map.get('BIOTYPE', 7)] if 'BIOTYPE' in self.field_map and len(cols) > self.field_map.get('BIOTYPE', 7) else '',
                'hgvsc': cols[self.field_map.get('HGVSc', 10)] if 'HGVSc' in self.field_map and len(cols) > self.field_map.get('HGVSc', 10) else '',
                'hgvsp': cols[self.field_map.get('HGVSp', 11)] if 'HGVSp' in self.field_map and len(cols) > self.field_map.get('HGVSp', 11) else '',
                'canonical': cols[self.field_map.get('CANONICAL', 24)] if 'CANONICAL' in self.field_map and len(cols) > self.field_map.get('CANONICAL', 24) else '',
                'mane_select': cols[self.field_map.get('MANE_SELECT', 25)] if 'MANE_SELECT' in self.field_map and len(cols) > self.field_map.get('MANE_SELECT', 25) else '',
                'mane_plus': cols[self.field_map.get('MANE_PLUS_CLINICAL', 26)] if 'MANE_PLUS_CLINICAL' in self.field_map and len(cols) > self.field_map.get('MANE_PLUS_CLINICAL', 26) else '',
                
                # Additional fields for ALT-specific filtering
                'allele_num': cols[self.field_map.get('ALLELE_NUM', -1)] if 'ALLELE_NUM' in self.field_map and len(cols) > self.field_map.get('ALLELE_NUM', -1) else '',
                'variant_class': cols[self.field_map.get('VARIANT_CLASS', -1)] if 'VARIANT_CLASS' in self.field_map and len(cols) > self.field_map.get('VARIANT_CLASS', -1) else ''
            }
            
            # Critical ALT-specific filtering for multi-allelic records
            if target_allele:
                ann_allele = annotation['allele'].strip()
                allele_num = annotation['allele_num'].strip()
                
                # Method 1: Direct allele match (most reliable)
                if ann_allele and ann_allele == target_allele:
                    # Perfect match - keep this annotation
                    pass
                elif ann_allele and ann_allele != target_allele:
                    # Different allele - skip unless it's a placeholder
                    if ann_allele not in ['-', '.', '*']:
                        continue
                
                # Method 2: ALLELE_NUM match for multiallelic variants  
                # ALLELE_NUM is 1-based (1=first ALT, 2=second ALT, etc.)
                if allele_num.isdigit() and target_allele_index is not None:
                    expected_allele_num = target_allele_index + 1  # Convert 0-based to 1-based
                    actual_allele_num = int(allele_num)
                    if actual_allele_num != expected_allele_num:
                        # This annotation is for a different ALT allele
                        continue
                
                # Method 3: Additional validation for complex cases
                # Skip annotations where HGVSp/HGVSc don't make sense for this allele
                hgvsp = annotation['hgvsp'].strip()
                hgvsc = annotation['hgvsc'].strip()
                
                # If we have HGVS annotations but no allele match, be more stringent
                if (hgvsp or hgvsc) and ann_allele and ann_allele != target_allele and ann_allele not in ['-', '.', '*']:
                    continue
                
            annotations.append(annotation)
        
        return annotations


class MANETranscriptPrioritizer:
    """Prioritize transcripts based on MANE status and other criteria"""
    
    @staticmethod
    def prioritize_annotations(annotations: List[Dict]) -> List[Dict]:
        """Sort annotations by priority: MANE Select > MANE Plus > Canonical > Others"""
        def priority_score(ann):
            score = 0
            
            # MANE Select has highest priority
            if ann.get('mane_select'):
                score += 1000
            
            # MANE Plus Clinical second
            elif ann.get('mane_plus'):
                score += 900
            
            # Canonical transcripts third
            elif ann.get('canonical') == 'YES':
                score += 800
            
            # Protein coding preferred over non-coding
            if ann.get('biotype') == 'protein_coding':
                score += 100
            
            # Prefer transcripts with HGVSp annotations
            if ann.get('hgvsp') and 'p.' in ann.get('hgvsp', ''):
                score += 50
            
            return score
        
        return sorted(annotations, key=priority_score, reverse=True)


class VariantNormalizer:
    """Normalize variants for consistent representation"""
    
    def __init__(self, reference_fasta: str = None):
        self.fasta = None
        self.chrom_aliases = {}  # For handling alternative contigs
        
        if reference_fasta:
            try:
                self.fasta = pysam.FastaFile(reference_fasta)
                print(f"✅ Loaded reference FASTA: {reference_fasta}")
                
                # Check for primary_assembly vs toplevel contigs
                if 'primary_assembly' in reference_fasta.lower():
                    print("✅ Using primary_assembly - chromosome names should be standard (1-22,X,Y,MT)")
                elif 'toplevel' in reference_fasta.lower():
                    print("⚠️  Using toplevel FASTA - may include alternative contigs")
                    print("💡 Consider using primary_assembly for standard chromosome names")
                    # Could add contig alias mapping here if needed
                    
            except Exception as e:
                print(f"❌ Could not load reference FASTA: {e}")
                print("💡 Ensure FASTA is indexed (.fai file present)")
                print("💡 For best results, use GRCh38.primary_assembly.fa.gz")
    
    def normalize_variant(self, chrom: str, pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
        """
        Normalize variant using left-alignment and minimal representation
        Returns (normalized_pos, normalized_ref, normalized_alt)
        """
        if not self.fasta:
            # Return as-is if no reference available
            return pos, ref, alt
        
        try:
            # Ensure chromosome format consistency
            if not chrom.startswith('chr') and chrom.isdigit():
                chrom_lookup = f"chr{chrom}"
            else:
                chrom_lookup = chrom
            
            if chrom_lookup not in self.fasta.references:
                chrom_lookup = chrom.replace('chr', '') if 'chr' in chrom else f"chr{chrom}"
            
            if chrom_lookup not in self.fasta.references:
                print(f"Chromosome {chrom} not found in reference")
                return pos, ref, alt
            
            # Verify reference allele matches
            if len(ref) == 1 and len(alt) == 1:
                # SNV - check reference base
                ref_base = self.fasta.fetch(chrom_lookup, pos-1, pos).upper()
                if ref_base != ref.upper():
                    print(f"Reference mismatch at {chrom}:{pos} - expected {ref}, found {ref_base}")
                    # Try to correct position
                    for offset in [-1, 1, -2, 2]:
                        test_pos = pos + offset
                        if test_pos > 0:
                            test_ref = self.fasta.fetch(chrom_lookup, test_pos-1, test_pos).upper()
                            if test_ref == ref.upper():
                                print(f"Corrected position to {test_pos}")
                                return test_pos, ref, alt
                return pos, ref, alt
            
            else:
                # Indel - normalize using left-alignment
                return self._left_align(chrom_lookup, pos, ref, alt)
                
        except Exception as e:
            print(f"Error normalizing variant {chrom}:{pos}:{ref}>{alt}: {e}")
            return pos, ref, alt
    
    def _left_align(self, chrom: str, pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
        """Left-align indel variant"""
        try:
            # Basic left-alignment implementation
            while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
                ref = ref[1:]
                alt = alt[1:]
                pos += 1
            
            # Ensure at least one base
            if not ref or not alt:
                pos -= 1
                context_base = self.fasta.fetch(chrom, pos-1, pos).upper()
                ref = context_base + ref
                alt = context_base + alt
            
            return pos, ref, alt
            
        except Exception as e:
            print(f"Error in left-alignment: {e}")
            return pos, ref, alt


class RobustGenomicMapper:
    """
    Robust genomic mapping with proper CSQ parsing, MANE prioritization,
    and variant normalization
    """
    
    def __init__(self, vep_data_path: str, reference_fasta: str = None, vcf_path: str = None, cache_dir: str = "transcript_cache"):
        self.vep_df = None
        self.variant_index = {}
        self.hgvs_index = {}
        
        # Initialize components
        self.csq_parser = CSQParser(vcf_path=vcf_path, tsv_path=vep_data_path)
        self.prioritizer = MANETranscriptPrioritizer()
        self.normalizer = VariantNormalizer(reference_fasta)
        self.transcript_mapper = TranscriptProteinMapper(cache_dir=cache_dir)
        
        # Initialize gnomAD cache
        self.gnomad_cache = {}
        self.load_gnomad_cache()
        
        # Load VEP data and pre-populate transcript mappings
        self.load_vep_data(vep_data_path)
        self._populate_transcript_cache()
    
    def __del__(self):
        """Save gnomAD cache when object is destroyed"""
        try:
            self.save_gnomad_cache()
        except:
            pass  # Ignore errors during cleanup
    
    def load_vep_data(self, vep_data_path: str):
        """Load and index VEP data"""
        try:
            self.vep_df = pd.read_csv(
                vep_data_path,
                sep='\t',
                names=['CHROM', 'POS', 'REF', 'ALT', 'CLNSIG', 'GENEINFO', 'CSQ'],
                dtype={'CHROM': str, 'POS': int}
            )
            print(f"Loaded {len(self.vep_df)} VEP annotations")
            self._build_indices()
            
        except Exception as e:
            print(f"Error loading VEP data: {e}")
            self.vep_df = pd.DataFrame()
    
    def _build_indices(self):
        """Build comprehensive indices for fast lookup"""
        print("Building genomic variant indices...")
        
        for idx, row in self.vep_df.iterrows():
            chrom = str(row['CHROM'])
            pos = int(row['POS'])
            ref = str(row['REF'])
            alt = str(row['ALT'])
            csq = str(row['CSQ'])
            
            # Parse CSQ annotations
            annotations = self.csq_parser.parse_csq_annotation(csq)
            prioritized = self.prioritizer.prioritize_annotations(annotations)
            
            if not prioritized:
                continue
            
            # Create variant info for each ALT (handle multiallelic)
            alt_alleles = [a.strip() for a in alt.split(',')]
            for alt_idx, alt_allele in enumerate(alt_alleles):
                # Parse CSQ with ALT-specific filtering using both allele string and index
                alt_specific_annotations = self.csq_parser.parse_csq_annotation(
                    csq, 
                    target_allele=alt_allele,
                    target_allele_index=alt_idx
                )
                
                if not alt_specific_annotations:
                    # Fallback: use all annotations if no ALT-specific match
                    alt_specific_annotations = annotations
                
                # Prioritize annotations for this specific ALT
                prioritized_alt = self.prioritizer.prioritize_annotations(alt_specific_annotations)
                
                if not prioritized_alt:
                    continue
                
                # Use the best annotation for this ALT
                matching_ann = prioritized_alt[0]
                
                # Create variant info
                variant_info = VariantInfo(
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt_allele,
                    gene=matching_ann.get('symbol', ''),
                    transcript=matching_ann.get('feature', ''),
                    hgvsc=matching_ann.get('hgvsc', ''),
                    hgvsp=matching_ann.get('hgvsp', ''),
                    is_mane=bool(matching_ann.get('mane_select') or matching_ann.get('mane_plus')),
                    mane_type="MANE Select" if matching_ann.get('mane_select') else 
                              "MANE Plus Clinical" if matching_ann.get('mane_plus') else ""
                )
                
                # Index by genomic coordinates
                var_key = f"{chrom}:{pos}:{ref}:{alt_allele}"
                if var_key not in self.variant_index:
                    self.variant_index[var_key] = []
                self.variant_index[var_key].append(variant_info)
                
                # Index by HGVS annotations with normalization
                if variant_info.gene and variant_info.hgvsp:
                    # Clean HGVSp
                    hgvsp_clean = variant_info.hgvsp
                    if ':' in hgvsp_clean:
                        hgvsp_clean = hgvsp_clean.split(':')[-1]
                    
                    # Index all HGVSp variants for robust lookup
                    for hgvsp_variant in hgvsp_variants(hgvsp_clean):
                        hgvs_key = f"{variant_info.gene}:{hgvsp_variant}"
                        if hgvs_key not in self.hgvs_index:
                            self.hgvs_index[hgvs_key] = []
                        self.hgvs_index[hgvs_key].append(variant_info)
        
        print(f"Built indices: {len(self.variant_index)} genomic variants, {len(self.hgvs_index)} HGVS annotations")
    
    def _populate_transcript_cache(self):
        """Pre-populate transcript->protein mappings for transcripts found in VEP data"""
        if not hasattr(self, 'vep_df') or self.vep_df is None:
            return
        
        print("🔄 Extracting transcripts from VEP data for cache population...")
        
        # Extract all unique transcripts from VEP data
        all_transcripts = set()
        
        try:
            for _, row in self.vep_df.iterrows():
                csq = str(row.get('CSQ', ''))
                if csq and csq != 'nan':
                    annotations = self.csq_parser.parse_csq_annotation(csq)
                    
                    for ann in annotations:
                        transcript = ann.get('feature', '')
                        if transcript and transcript.startswith('ENST'):
                            all_transcripts.add(transcript)
        
            print(f"📊 Found {len(all_transcripts)} unique transcripts in VEP data")
            
            # Batch process transcripts for efficiency
            if all_transcripts:
                transcript_list = list(all_transcripts)
                # Process in smaller batches to avoid overwhelming API
                batch_size = 20
                total_mapped = 0
                
                for i in range(0, len(transcript_list), batch_size):
                    batch = transcript_list[i:i+batch_size]
                    print(f"🔍 Processing transcript batch {i//batch_size + 1}/{(len(transcript_list) + batch_size - 1)//batch_size}")
                    
                    mappings = self.transcript_mapper.bulk_map_transcripts(batch)
                    mapped_count = sum(1 for v in mappings.values() if v)
                    total_mapped += mapped_count
                    
                    print(f"  ✅ Mapped {mapped_count}/{len(batch)} transcripts in this batch")
                    
                    # Small delay between batches
                    time.sleep(0.5)
                
                print(f"🎯 Total transcript mappings cached: {total_mapped}/{len(all_transcripts)}")
                
        except Exception as e:
            print(f"⚠️ Error populating transcript cache: {e}")
    
    def get_genomic_coords(self, gene: str, hgvsp: str, prefer_mane: bool = True) -> Optional[VariantInfo]:
        """
        Get genomic coordinates for gene and HGVSp annotation with HGVSp normalization
        Returns best matching VariantInfo based on MANE priority
        Falls back to VEP REST API if not found in local data
        """
        # Clean and normalize HGVSp notation
        hgvsp_clean = hgvsp
        if ':' in hgvsp_clean:
            hgvsp_clean = hgvsp_clean.split(':')[-1]
        
        # Try all HGVSp variants for lookup
        variants = []
        for hgvsp_variant in hgvsp_variants(hgvsp_clean):
            key = f"{gene}:{hgvsp_variant}"
            found_variants = self.hgvs_index.get(key, [])
            variants.extend(found_variants)
        
        # Remove duplicates
        unique_variants = []
        seen_ids = set()
        for v in variants:
            variant_id = f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}"
            if variant_id not in seen_ids:
                seen_ids.add(variant_id)
                unique_variants.append(v)
        
        variants = unique_variants
        
        if variants:
            if prefer_mane:
                # Prioritize MANE transcripts
                mane_variants = [v for v in variants if v.is_mane]
                if mane_variants:
                    # Prefer MANE Select over MANE Plus Clinical
                    mane_select = [v for v in mane_variants if "Select" in v.mane_type]
                    if mane_select:
                        return mane_select[0]
                    return mane_variants[0]
            
            return variants[0]  # Return first available
        
        # Fallback to VEP REST API for unmapped variants
        normalized_hgvsp = normalize_hgvsp(hgvsp_clean)
        print(f"Variant {gene}:{normalized_hgvsp} not found in local VEP data, trying REST API...")
        return self._query_vep_rest_api(gene, normalized_hgvsp)
    
    def _query_vep_rest_api(self, gene: str, hgvsp: str) -> Optional[VariantInfo]:
        """
        Query VEP REST API to resolve HGVS protein notation to genomic coordinates
        Uses comprehensive transcript->protein mapping
        """
        try:
            # Get canonical transcript for the gene
            canonical_mappings = self.transcript_mapper.get_canonical_transcript_mapping(gene)
            
            if not canonical_mappings:
                print(f"❌ No canonical transcript found for gene {gene}")
                return None
            
            # Try each transcript->protein mapping
            for transcript, protein in canonical_mappings.items():
                if not protein:
                    continue
                    
                hgvs_notation = f"{protein}:{hgvsp}"
                print(f"🔍 Trying {transcript} -> {protein}: {hgvs_notation}")
                
                result = self._query_vep_api_with_hgvs(hgvs_notation, gene, transcript)
                if result:
                    return result
            
            # Fallback: try to get any transcript for this gene from VEP data
            print(f"🔄 Trying fallback transcript discovery for {gene}")
            return self._try_fallback_transcripts(gene, hgvsp)
                
        except Exception as e:
            print(f"❌ Error in VEP REST API query: {e}")
            return None
    
    def _query_vep_api_with_hgvs(self, hgvs_notation: str, gene: str, transcript: str) -> Optional[VariantInfo]:
        """Query VEP API with specific HGVS notation"""
        try:
            # VEP REST API endpoint
            url = f"https://rest.ensembl.org/vep/human/hgvs/{hgvs_notation}"
            
            headers = {
                'Content-Type': 'application/json',
                'Accept': 'application/json'
            }
            
            response = self.transcript_mapper.session.get(url, headers=headers, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                
                if data and len(data) > 0:
                    result = data[0]  # Take first result
                    
                    # Extract genomic coordinates
                    if 'seq_region_name' in result and 'start' in result:
                        chrom = str(result['seq_region_name'])
                        pos = int(result['start'])
                        
                        # Extract alleles
                        allele_string = result.get('allele_string', '')
                        if '/' in allele_string:
                            ref, alt = allele_string.split('/', 1)
                        else:
                            # Default for SNVs
                            ref, alt = 'A', 'T'
                        
                        variant_info = VariantInfo(
                            chrom=chrom,
                            pos=pos,
                            ref=ref,
                            alt=alt,
                            gene=gene,
                            transcript=transcript,
                            hgvsc=result.get('hgvsc', ''),
                            hgvsp=hgvs_notation.split(':')[-1],
                            is_mane=False,  # Unknown from REST API
                            mane_type=""
                        )
                        
                        print(f"✅ Resolved via VEP REST API: {chrom}:{pos}:{ref}>{alt}")
                        return variant_info
            
            elif response.status_code == 429:
                print(f"⏱️ Rate limited, waiting 2 seconds...")
                time.sleep(2)
                # Single retry
                return self._query_vep_api_with_hgvs(hgvs_notation, gene, transcript)
            
            else:
                print(f"❌ VEP REST API returned status {response.status_code}")
                
        except Exception as e:
            print(f"❌ Error querying VEP REST API: {e}")
        
        return None
    
    def _try_fallback_transcripts(self, gene: str, hgvsp: str) -> Optional[VariantInfo]:
        """Try to find transcripts from existing VEP data for fallback mapping"""
        try:
            # Extract unique transcripts for this gene from VEP data
            gene_transcripts = set()
            if hasattr(self, 'vep_df') and self.vep_df is not None:
                for _, row in self.vep_df.iterrows():
                    csq = str(row.get('CSQ', ''))
                    annotations = self.csq_parser.parse_csq_annotation(csq)
                    
                    for ann in annotations:
                        if ann.get('symbol', '').upper() == gene.upper():
                            transcript = ann.get('feature', '')
                            if transcript.startswith('ENST'):
                                gene_transcripts.add(transcript)
            
            if gene_transcripts:
                print(f"🔍 Found {len(gene_transcripts)} transcripts for {gene} in VEP data")
                
                # Try each transcript
                for transcript in list(gene_transcripts)[:5]:  # Limit to 5 attempts
                    protein = self.transcript_mapper.get_protein_id(transcript)
                    if protein:
                        hgvs_notation = f"{protein}:{hgvsp}"
                        result = self._query_vep_api_with_hgvs(hgvs_notation, gene, transcript)
                        if result:
                            return result
            
        except Exception as e:
            print(f"❌ Error in fallback transcript discovery: {e}")
        
        return None
    
    def close(self):
        """Clean up resources"""
        if hasattr(self, 'transcript_mapper'):
            self.transcript_mapper.close()
    
    def get_gnomad_variant_id(self, gene: str, hgvsp: str, normalize: bool = True) -> Optional[str]:
        """
        Get gnomAD-compatible variant ID (chr-pos-ref-alt)
        Ensures proper chromosome format (no 'chr' prefix) and 1-based coordinates
        """
        variant_info = self.get_genomic_coords(gene, hgvsp)
        if not variant_info:
            return None
        
        chrom = variant_info.chrom
        pos = variant_info.pos
        ref = variant_info.ref
        alt = variant_info.alt
        
        # Remove 'chr' prefix for gnomAD compatibility
        if chrom.startswith('chr'):
            chrom = chrom[3:]
        
        # Normalize variant coordinates if requested
        if normalize:
            pos, ref, alt = self.normalizer.normalize_variant(chrom, pos, ref, alt)
        
        # Ensure 1-based coordinates (gnomAD standard)
        return f"{chrom}-{pos}-{ref}-{alt}"
    
    def query_gnomad_with_mapping(self, gene: str, hgvsp: str) -> Dict:
        """
        Query gnomAD using robust genomic mapping with caching and backoff
        """
        variant_id = self.get_gnomad_variant_id(gene, hgvsp)
        if not variant_id:
            return {
                'found': False, 
                'error': f'No genomic mapping found for {gene}:{hgvsp}',
                'variant_id': None
            }
        
        # Check cache first (if we add a cache dict)
        cache_key = f"gnomad_{variant_id}"
        if hasattr(self, 'gnomad_cache') and cache_key in self.gnomad_cache:
            cached_result = self.gnomad_cache[cache_key]
            print(f"📋 Using cached gnomAD result for {variant_id}")
            return cached_result
        
        # gnomAD GraphQL query
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
        
        max_retries = 3
        backoff_delay = 1
        
        for attempt in range(max_retries):
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
                        result = {
                            'found': True,
                            'variant_id': variant_id,
                            'af_exome': variant.get('exome', {}).get('af', 0),
                            'af_genome': variant.get('genome', {}).get('af', 0),
                            'populations': variant.get('populations', [])
                        }
                        
                        # Cache successful results
                        if not hasattr(self, 'gnomad_cache'):
                            self.gnomad_cache = {}
                        self.gnomad_cache[cache_key] = result
                        
                        return result
                
                elif response.status_code == 429:
                    # Rate limited - exponential backoff
                    delay = backoff_delay * (2 ** attempt)
                    print(f"⏱️ Rate limited, waiting {delay}s before retry {attempt + 1}/{max_retries}")
                    time.sleep(delay)
                    continue
                    
                elif response.status_code == 404:
                    # Variant not found - cache negative result
                    result = {
                        'found': False, 
                        'variant_id': variant_id,
                        'af_exome': 0, 
                        'af_genome': 0,
                        'error': 'Variant not found in gnomAD'
                    }
                    
                    if not hasattr(self, 'gnomad_cache'):
                        self.gnomad_cache = {}
                    self.gnomad_cache[cache_key] = result
                    
                    return result
                    
                else:
                    print(f"❌ gnomAD API error {response.status_code} for {variant_id}")
                    if attempt < max_retries - 1:
                        time.sleep(backoff_delay)
                        continue
            
            except Exception as e:
                print(f"❌ Error querying gnomAD for {variant_id}: {e}")
                if attempt < max_retries - 1:
                    time.sleep(backoff_delay)
                    continue
        
        # All retries failed
        return {
            'found': False, 
            'variant_id': variant_id,
            'af_exome': 0, 
            'af_genome': 0,
            'error': 'Failed to query gnomAD after retries'
        }
    
    def get_variant_summary(self, gene: str, hgvsp: str) -> Dict:
        """Get comprehensive variant information"""
        variant_info = self.get_genomic_coords(gene, hgvsp)
        if not variant_info:
            return {'error': f'No mapping found for {gene}:{hgvsp}'}
        
        return {
            'gene': variant_info.gene,
            'hgvsp': variant_info.hgvsp,
            'hgvsc': variant_info.hgvsc,
            'genomic_coords': f"{variant_info.chrom}:{variant_info.pos}:{variant_info.ref}>{variant_info.alt}",
            'transcript': variant_info.transcript,
            'is_mane': variant_info.is_mane,
            'mane_type': variant_info.mane_type,
            'gnomad_id': self.get_gnomad_variant_id(gene, hgvsp)
        }


    def save_gnomad_cache(self, cache_file: str = None):
        """Save gnomAD cache to disk for persistence"""
        if not hasattr(self, 'gnomad_cache') or not self.gnomad_cache:
            print("📋 No gnomAD cache to save")
            return
        
        if cache_file is None:
            cache_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'cache', 'gnomad_cache.json')
        
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        
        try:
            with open(cache_file, 'w') as f:
                json.dump(self.gnomad_cache, f, indent=2)
            print(f"💾 Saved {len(self.gnomad_cache)} gnomAD cache entries to {cache_file}")
        except Exception as e:
            print(f"❌ Error saving gnomAD cache: {e}")
    
    def load_gnomad_cache(self, cache_file: str = None):
        """Load gnomAD cache from disk"""
        if cache_file is None:
            cache_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'cache', 'gnomad_cache.json')
        
        if not os.path.exists(cache_file):
            print("📋 No existing gnomAD cache file found")
            return
        
        try:
            with open(cache_file, 'r') as f:
                self.gnomad_cache = json.load(f)
            print(f"📋 Loaded {len(self.gnomad_cache)} gnomAD cache entries from {cache_file}")
        except Exception as e:
            print(f"❌ Error loading gnomAD cache: {e}")
            self.gnomad_cache = {}


# Example usage and testing
if __name__ == "__main__":
    # Initialize with paths
    vep_path = "../data/processed/vep_extracted.tsv"
    # reference_fasta = "/path/to/GRCh38.fasta"  # Optional
    
    mapper = RobustGenomicMapper(
        vep_data_path=vep_path,
        # reference_fasta=reference_fasta
    )
    
    # Test with hemoglobin variants
    test_variants = [
        ("HBB", "p.(Glu6Val)"),  # Sickle cell
        ("HBB", "p.(Glu6Lys)"),  # HbC
        ("HBA1", "p.(Phe8Leu)") # Alpha variant
    ]
    
    for gene, hgvsp in test_variants:
        print(f"\nTesting {gene}:{hgvsp}")
        summary = mapper.get_variant_summary(gene, hgvsp)
        print(f"Summary: {summary}")
        
        gnomad_result = mapper.query_gnomad_with_mapping(gene, hgvsp)
        print(f"gnomAD: {gnomad_result}")