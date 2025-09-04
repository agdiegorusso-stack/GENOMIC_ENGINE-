#!/usr/bin/env python3
"""
Script per analisi mutazionale batch di tutti i geni dell'emoglobina
Processa un gene alla volta per evitare timeout e gestire errori
"""

import os
import sys
import time
import glob
from datetime import datetime
from mutational_analysis import ComprehensiveMutationAnalyzer

def log_message(message):
    """Log con timestamp"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def get_gene_list():
    """Ottieni lista di tutti i geni da analizzare"""
    fasta_files = sorted(glob.glob('hemoglobin_genes/fasta/*_cds.fasta'))
    return [os.path.basename(f).replace('_cds.fasta', '') for f in fasta_files]

def load_sequence(gene_name):
    """Carica sequenza di un gene"""
    fasta_file = f'hemoglobin_genes/fasta/{gene_name}_cds.fasta'
    
    if not os.path.exists(fasta_file):
        return None
    
    with open(fasta_file, 'r') as f:
        content = f.read().strip()
        lines = content.split('\n')
        sequence = lines[1] if len(lines) > 1 else ''
    
    return sequence if len(sequence) > 0 else None

def analyze_single_gene(analyzer, gene_name, progress_info=""):
    """Analizza singolo gene"""
    log_message(f"{progress_info}Iniziando analisi {gene_name}")
    
    # Carica sequenza
    sequence = load_sequence(gene_name)
    if not sequence:
        log_message(f"[ERROR] Sequenza non valida per {gene_name}")
        return False
    
    log_message(f"  Sequenza {gene_name}: {len(sequence)} bp")
    
    try:
        # Analizza gene
        start_time = time.time()
        results = analyzer.analyze_gene(gene_name, sequence, 
                                      include_clinical=True,
                                      include_double=False)
        
        elapsed = time.time() - start_time
        
        # Conteggio per tipo
        mutation_types = {}
        for mut in results:
            mut_type = mut.get('mutation_type', 'unknown')
            mutation_types[mut_type] = mutation_types.get(mut_type, 0) + 1
        
        log_message(f"[OK] {gene_name} COMPLETATO: {len(results)} mutazioni in {elapsed:.1f}s")
        
        # Log breakdown
        for mut_type, count in sorted(mutation_types.items()):
            log_message(f"    {mut_type}: {count}")
        
        return True
        
    except Exception as e:
        log_message(f"[ERROR] in {gene_name}: {str(e)}")
        return False

def main():
    """Main batch analysis"""
    log_message("="*80)
    log_message("BATCH MUTATION ANALYSIS - TUTTI I GENI EMOGLOBINA")
    log_message("="*80)
    
    # Ottieni lista geni
    gene_list = get_gene_list()
    total_genes = len(gene_list)
    
    log_message(f"Geni da analizzare: {total_genes}")
    log_message(f"Primi 10: {gene_list[:10]}")
    
    # Inizializza analyzer
    log_message("Inizializzazione ComprehensiveMutationAnalyzer...")
    try:
        analyzer = ComprehensiveMutationAnalyzer()
        log_message("[OK] Analyzer inizializzato")
    except Exception as e:
        log_message(f"[ERROR] inizializzazione: {str(e)}")
        return
    
    # Contatori
    completed = 0
    failed = 0
    start_time = time.time()
    
    # Analizza ogni gene
    for i, gene_name in enumerate(gene_list, 1):
        progress_info = f"[{i}/{total_genes}] "
        
        if analyze_single_gene(analyzer, gene_name, progress_info):
            completed += 1
        else:
            failed += 1
        
        # Progress report ogni 10 geni
        if i % 10 == 0 or i == total_genes:
            elapsed = time.time() - start_time
            rate = i / elapsed * 60  # geni per minuto
            eta_minutes = (total_genes - i) / rate if rate > 0 else 0
            
            log_message(f"PROGRESS: {i}/{total_genes} geni | Completati: {completed} | Falliti: {failed}")
            log_message(f"         Rate: {rate:.1f} geni/min | ETA: {eta_minutes:.1f} min")
        
        # Pausa per evitare sovraccarico
        time.sleep(0.5)
    
    # Report finale
    total_time = time.time() - start_time
    
    log_message("="*80)
    log_message("ANALISI BATCH COMPLETATA!")
    log_message("="*80)
    log_message(f"Totale geni processati: {total_genes}")
    log_message(f"Completati con successo: {completed}")
    log_message(f"Falliti: {failed}")
    log_message(f"Tempo totale: {total_time/60:.1f} minuti")
    log_message(f"Rate medio: {total_genes/(total_time/60):.1f} geni/min")
    
    # Lista files generati
    log_message("\nFile generati:")
    output_files = sorted(glob.glob('mutation_analysis/*.csv'))
    for f in output_files:
        size = os.path.getsize(f)
        log_message(f"  {os.path.basename(f)} ({size} bytes)")
    
    log_message(f"\nDirectory output: mutation_analysis/")
    log_message("Prossimo passo: Integrazione database clinici")

if __name__ == "__main__":
    main()