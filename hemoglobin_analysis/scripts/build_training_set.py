#!/usr/bin/env python3
# scripts/build_training_set.py (VERSIONE 6.2 - Correzioni finali)

import os
import re
import pandas as pd
from cyvcf2 import VCF

# Proprietà AA di base
def _load_aa_properties():
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
        '*': {'hydrophobic': 0, 'size': 0, 'charge': 0, 'polar': 0},
    }

AA_PROPERTIES = _load_aa_properties()

# Mappa 3-lettere -> 1-lettera
AA3_TO1 = {
    'Ala': 'A','Arg': 'R','Asn': 'N','Asp': 'D','Cys': 'C','Glu': 'E','Gln': 'Q','Gly': 'G',
    'His': 'H','Ile': 'I','Leu': 'L','Lys': 'K','Met': 'M','Phe': 'F','Pro': 'P','Ser': 'S',
    'Thr': 'T','Trp': 'W','Tyr': 'Y','Val': 'V','Ter': '*','Stop': '*'
}

# Regex per HGVSp (3-lettere e 1-lettera); asterisco escapato
RE_HGVSP_3 = re.compile(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)')
RE_HGVSP_1 = re.compile(r'p\.([A-Z\*])(\d+)([A-Z\*])')

def _calculate_missense_severity(wt_aa, mut_aa):
    if wt_aa not in AA_PROPERTIES or mut_aa not in AA_PROPERTIES:
        return 5.0
    wt_props, mut_props = AA_PROPERTIES[wt_aa], AA_PROPERTIES[mut_aa]
    hydro_diff = abs(wt_props['hydrophobic'] - mut_props['hydrophobic'])
    size_diff = abs(wt_props['size'] - mut_props['size']) / 100
    charge_diff = abs(wt_props['charge'] - mut_props['charge']) * 2
    polar_diff = abs(wt_props['polar'] - mut_props['polar']) * 2
    return round(min(10, hydro_diff + size_diff + charge_diff + polar_diff), 2)

def _normalize_info(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]

def _get_csq_fields(vcf_reader):
    h = vcf_reader.get_header_type('CSQ')
    if not h:
        return []
    desc = h.get('Description', '')
    m = re.search(r'Format:\s*([^"]+)', desc)
    return m.group(1).split('|') if m else []

def _iter_csq_items(info_value):
    if not info_value:
        return []
    if isinstance(info_value, (list, tuple)):
        items = []
        for s in info_value:
            items.extend(str(s).split(','))
        return items
    return str(info_value).split(',')

def _extract_hgvsp_from_texts(texts):
    for s in texts:
        m = RE_HGVSP_3.search(s) or RE_HGVSP_1.search(s)
        if not m:
            continue
        if m.re is RE_HGVSP_3:
            wt3, pos, mut3 = m.groups()
            wt = AA3_TO1.get(wt3)
            mut = '*' if mut3 == '*' else AA3_TO1.get(mut3)
            if wt and mut:
                return m.group(0), wt, mut
        else:
            wt, pos, mut = m.groups()
            if wt in AA_PROPERTIES and mut in AA_PROPERTIES:
                return m.group(0), wt, mut
    return None, None, None

def build_feature_rich_training_set(vcf_path, output_csv):
    print(f"Building feature-rich training set from {vcf_path}...")
    vcf_reader = VCF(vcf_path)

    # Prepara parsing CSQ
    csq_fields = _get_csq_fields(vcf_reader)
    csq_index = {name: idx for idx, name in enumerate(csq_fields)} if csq_fields else {}
    idx_consequence = csq_index.get('Consequence')
    idx_hgvsp = csq_index.get('HGVSp')
    idx_symbol = csq_index.get('SYMBOL')

    training_variants = []

    PATHOGENIC_TERMS = ['pathogenic', 'likely_pathogenic', 'pathogenic/likely_pathogenic']
    BENIGN_TERMS = ['benign', 'likely_benign', 'benign/likely_benign']

    for variant in vcf_reader:
        info = variant.INFO

        # Missense via MC (SO:0001583) se disponibile
        mc_values = _normalize_info(info.get('MC'))
        is_missense_mc = any(('SO:0001583' in v) or ('missense_variant' in v) for v in mc_values)

        # CLNSIG/CLNREVSTAT
        clnsig_s = '|'.join(_normalize_info(info.get('CLNSIG', ''))).lower()
        revstat_s = '|'.join(_normalize_info(info.get('CLNREVSTAT', ''))).lower()

        label = -1
        if any(term in clnsig_s for term in PATHOGENIC_TERMS) and 'conflict' not in clnsig_s:
            label = 1
        elif any(term in clnsig_s for term in BENIGN_TERMS) and 'conflict' not in clnsig_s:
            label = 0

        has_good_review = (
            ('reviewed_by_expert_panel' in revstat_s) or
            ('practice_guideline' in revstat_s) or
            ('multiple_submitters' in revstat_s) or
            ('single_submitter' in revstat_s)
        )
        if not (has_good_review and label != -1):
            continue

        # Estrazione HGVSp e gene
        hgvsp = wt_aa = mut_aa = None
        gene_symbol = '.'

        # CSQ strutturato (VEP)
        csq_raw = info.get('CSQ')
        csq_items = _iter_csq_items(csq_raw)
        if csq_items and idx_consequence is not None:
            for item in csq_items:
                cols = item.split('|')
                consequence = cols[idx_consequence] if idx_consequence < len(cols) else ''
                if 'missense_variant' not in consequence:
                    continue
                hsp = cols[idx_hgvsp] if (idx_hgvsp is not None and idx_hgvsp < len(cols)) else ''
                if hsp:
                    p_part = hsp[hsp.find('p.'):] if 'p.' in hsp else ''
                    if p_part:
                        m = RE_HGVSP_3.search(p_part) or RE_HGVSP_1.search(p_part)
                        if m:
                            if m.re is RE_HGVSP_3:
                                wt3, pos, mut3 = m.groups()
                                wt_aa = AA3_TO1.get(wt3)
                                mut_aa = '*' if mut3 == '*' else AA3_TO1.get(mut3)
                            else:
                                wt_aa, pos, mut_aa = m.groups()
                            if wt_aa and mut_aa:
                                hgvsp = m.group(0)
                                if idx_symbol is not None and idx_symbol < len(cols):
                                    gene_symbol = cols[idx_symbol] or '.'
                                break

        # Fallback ANN / CLNHGVS
        if not hgvsp:
            hgvsp, wt_aa, mut_aa = _extract_hgvsp_from_texts(_normalize_info(info.get('ANN')))
        if not hgvsp:
            hgvsp, wt_aa, mut_aa = _extract_hgvsp_from_texts(_normalize_info(info.get('CLNHGVS')))

        if not (hgvsp and wt_aa and mut_aa):
            continue

        # Gene fallback da GENEINFO se non trovato in CSQ
        if not gene_symbol or gene_symbol == '.':
            geneinfo = info.get('GENEINFO', '')
            if geneinfo:
                first_pair = str(geneinfo).split('|')[0]  # Prima coppia
                gene_symbol = first_pair.split(':')[0] if ':' in first_pair else (first_pair or '.')
            else:
                gene_symbol = '.'

        # Richiedi missense in MC o p. presente
        if not (is_missense_mc or 'p.' in hgvsp):
            continue

        severity = _calculate_missense_severity(wt_aa, mut_aa)

        # AF_EXAC normalizzato
        af_exac_val = info.get('AF_EXAC', 0.0) or 0.0
        try:
            if isinstance(af_exac_val, (list, tuple)):
                af_exac = float(af_exac_val[0]) if af_exac_val else 0.0
            else:
                af_exac = float(af_exac_val)
        except Exception:
            af_exac = 0.0

        training_variants.append({
            'gene': gene_symbol,
            'hgvs_p': hgvsp,
            'severity_score': severity,
            'conservation_score': 0.5,      # placeholder
            'position_criticality': 0.5,    # placeholder
            'chemical_change': abs(AA_PROPERTIES[wt_aa]['charge'] - AA_PROPERTIES[mut_aa]['charge']),
            'structural_change': 0.5,       # placeholder
            'af_exac': af_exac,
            'label': label
        })

    df = pd.DataFrame(training_variants).drop_duplicates(subset=['gene', 'hgvs_p'])
    df.to_csv(output_csv, index=False)

    print(f"\nTraining set salvato in {output_csv}")
    print(f"Totale varianti INCLUSE nel dataset: {len(df)}")
    if not df.empty:
        print("Distribuzione delle etichette:\n", df['label'].value_counts())

if __name__ == "__main__":
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    vcf_path = os.path.join(project_root, 'data', 'processed', 'clinvar_annotated_final.vcf')
    output_csv = os.path.join(project_root, 'data', 'processed', 'clinvar_training_set_features.csv')
    build_feature_rich_training_set(vcf_path, output_csv)
