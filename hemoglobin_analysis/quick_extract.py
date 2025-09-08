#!/usr/bin/env python3
import pandas as pd
import re

# Estrazione rapida dai dati TSV  
df = pd.read_csv('data/processed/vep_extracted.tsv', sep='\t', names=['CHROM', 'POS', 'REF', 'ALT', 'CLNSIG', 'GENEINFO', 'CSQ'])

results = []
for _, row in df.iterrows():
    clnsig = str(row['CLNSIG']).lower()
    
    # Label semplificato
    label = -1
    if 'pathogenic' in clnsig and 'conflict' not in clnsig:
        label = 1
    elif 'benign' in clnsig and 'conflict' not in clnsig:
        label = 0
    
    if label == -1:
        continue
        
    # Estrai HGVSp dal CSQ (campo 11)
    for transcript in str(row['CSQ']).split(','):
        cols = transcript.split('|')
        if len(cols) > 11 and 'missense_variant' in cols[1]:
            hgvsp = cols[11]
            symbol = cols[3] if len(cols) > 3 else '.'
            
            if 'p.(' in hgvsp and symbol and symbol != '.':
                match = re.search(r'p\.\(([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)', hgvsp)
                if match:
                    results.append({
                        'gene': symbol,
                        'hgvs_p': hgvsp.split(':')[-1] if ':' in hgvsp else hgvsp,
                        'label': label
                    })
                break

results_df = pd.DataFrame(results).drop_duplicates()
results_df.to_csv('data/processed/clinvar_training_set_features.csv', index=False)
print(f'Varianti trovate: {len(results_df)}')
print('Distribuzione etichette:')
print(results_df['label'].value_counts())
print('Prime varianti:')
print(results_df.head())