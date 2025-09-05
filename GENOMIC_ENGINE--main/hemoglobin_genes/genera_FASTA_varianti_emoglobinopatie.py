#!/usr/bin/env python3
"""
Algoritmo per generare un file FASTA di varianti di emoglobinopatie
a partire dalle sequenze di riferimento e da un file JSON di mutazioni.

Autore: Assistant
Data: 2025

Prerequisiti:
1. Eseguire prima `hemoglobin_downloader.py` per scaricare le sequenze di riferimento.
2. Assicurarsi che la cartella `hemoglobin_genes/fasta/` sia presente.
3. Posizionare `IthaGenes_variations_export (1).json` nella stessa directory.
4. Installare Biopython: pip install biopython
"""

import os
import json
import re
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord

# --- CONFIGURAZIONE ---
FASTA_REF_DIR = "fasta"
MUTATION_JSON_FILE = "IthaGenes_variations_export.json"
OUTPUT_FASTA_FILE = "hemoglobinopathy_variants.fasta"

# Mappatura per gestire nomi di geni non standard nella colonna 'Genes' del JSON
GENE_MAP = {
    "β": ["HBB"],
    "α2": ["HBA2"],
    "α1": ["HBA1"],
    "α1 or α2": ["HBA1", "HBA2"],
    "δ": ["HBD"],
    "Gγ": ["HBG2"],
    "Aγ": ["HBG1"],
    "ε": ["HBE1"],
    "ζ": ["HBZ"],
}

def load_reference_sequences(fasta_dir):
    """Carica tutte le sequenze di riferimento FASTA in un dizionario."""
    sequences = {}
    if not os.path.isdir(fasta_dir):
        print(f"ERRORE: La directory '{fasta_dir}' non esiste. Esegui prima lo script 'hemoglobin_downloader.py'.")
        exit()
        
    print(f"Caricamento delle sequenze di riferimento da '{fasta_dir}'...")
    for filename in os.listdir(fasta_dir):
        if filename.endswith(".fasta"):
            gene_name = filename.replace("_cds.fasta", "")
            record = SeqIO.read(os.path.join(fasta_dir, filename), "fasta")
            sequences[gene_name] = record.seq
    print(f"Caricate {len(sequences)} sequenze di riferimento.\n")
    return sequences

def parse_hgvs(hgvs_string):
    """Estrae informazioni utili dalla notazione HGVS."""
    if not isinstance(hgvs_string, str) or ":" not in hgvs_string:
        return None, None

    gene_part, mutation_part = hgvs_string.split(":", 1)
    
    # Ignora le mutazioni complesse/non-CDS per ora
    if not mutation_part.startswith("c."):
        return gene_part, None
        
    return gene_part, mutation_part

def apply_mutation(ref_seq, mutation_str):
    """
    Applica una singola mutazione (formato HGVS 'c.') a una sequenza.
    Restituisce la sequenza mutata o None se la mutazione non è applicabile.
    """
    try:
        mut_seq = MutableSeq(str(ref_seq))
        
        # Sostituzione: c.79G>A
        match = re.match(r"c\.(\d+)(\w)>(\w)", mutation_str)
        if match:
            pos, ref, alt = match.groups()
            pos = int(pos) - 1  # Converte in indice 0-based
            if pos < len(mut_seq) and mut_seq[pos].upper() == ref.upper():
                mut_seq[pos] = alt
                return Seq(mut_seq)
            return None # Ref non corrisponde

        # Delezione singola/multipla: c.126_129delCTTT o c.20delA
        match = re.match(r"c\.(\d+)_?(\d*)del(\w*)", mutation_str)
        if match:
            start, end, deleted_seq = match.groups()
            start = int(start) - 1
            end = int(end) if end else start
            
            # Verifica che la sequenza da eliminare corrisponda
            ref_segment = str(mut_seq[start:end+1]).upper()
            if not deleted_seq or ref_segment == deleted_seq.upper():
                 del mut_seq[start:end+1]
                 return Seq(mut_seq)
            return None

        # Duplicazione: c.27dupG o c.143_146dupATCT
        match = re.match(r"c\.(\d+)_?(\d*)dup(\w*)", mutation_str)
        if match:
            start, end, duplicated_seq = match.groups()
            start = int(start) - 1
            end = int(end) if end else start
            
            ref_segment = str(mut_seq[start:end+1])
            mut_seq[end+1:end+1] = ref_segment
            return Seq(mut_seq)

        # Inserzione: c.27_28insAGAA
        match = re.match(r"c\.(\d+)_(\d+)ins(\w+)", mutation_str)
        if match:
            start_pos, _, ins_seq = match.groups()
            pos = int(start_pos) # Inserisce dopo la base specificata
            mut_seq[pos:pos] = ins_seq.upper()
            return Seq(mut_seq)
            
    except (ValueError, IndexError):
        return None # Errore nel parsing o posizione fuori dai limiti

    return None # Nessun pattern corrispondente

def main():
    """Funzione principale per generare il file FASTA delle varianti."""
    print("Avvio dell'algoritmo per la generazione del file FASTA delle varianti.")
    
    # Carica dati di base
    ref_sequences = load_reference_sequences(FASTA_REF_DIR)
    
    try:
        with open(MUTATION_JSON_FILE, 'r', encoding='utf-8') as f:
            mutation_data = json.load(f)
    except FileNotFoundError:
        print(f"ERRORE: File delle mutazioni '{MUTATION_JSON_FILE}' non trovato.")
        return
        
    variants_to_write = []
    processed_count = 0
    skipped_count = 0

    # Itera su ogni variante nel file JSON
    for entry in mutation_data.get("body", []):
        try:
            # Estrai informazioni rilevanti
            itha_id = entry[0]
            common_name = entry[1].replace(" ", "_").replace("/", "-") if entry[1] != "N/A" else ""
            hgvs_name = entry[3]
            gene_symbols_str = entry[4]
            
            # Determina il gene target
            target_genes = []
            hgvs_gene, hgvs_mutation = parse_hgvs(hgvs_name)
            
            if hgvs_gene and hgvs_gene != "N/A":
                target_genes.append(hgvs_gene)
            elif gene_symbols_str in GENE_MAP:
                 target_genes.extend(GENE_MAP[gene_symbols_str])
            
            if not hgvs_mutation or not target_genes:
                skipped_count += 1
                continue

            # Applica la mutazione per ogni gene target
            for gene in set(target_genes): # Usa set per evitare duplicati
                if gene in ref_sequences:
                    ref_seq = ref_sequences[gene]
                    mutated_seq = apply_mutation(ref_seq, hgvs_mutation)
                    
                    if mutated_seq:
                        # Crea un'intestazione FASTA descrittiva
                        header_parts = [itha_id, gene, common_name, hgvs_name]
                        header = "|".join(filter(None, header_parts)) # Filtra parti vuote
                        
                        record = SeqRecord(mutated_seq, id=header, description="")
                        variants_to_write.append(record)
                        processed_count += 1
                    else:
                        skipped_count += 1
                else:
                    skipped_count += 1

        except Exception as e:
            #print(f"Errore nell'elaborazione della riga {entry}: {e}")
            skipped_count += 1

    # Scrivi il file FASTA di output
    if variants_to_write:
        SeqIO.write(variants_to_write, OUTPUT_FASTA_FILE, "fasta")
        print("\n" + "="*50)
        print("PROCESSO COMPLETATO")
        print(f"✓ Generato il file '{OUTPUT_FASTA_FILE}' con successo.")
        print(f"  - Varianti create: {processed_count}")
        print(f"  - Varianti saltate (formato non supportato o dati incompleti): {skipped_count}")
        print("="*50)
    else:
        print("\nNessuna variante è stata generata. Controlla i file di input e i percorsi.")

if __name__ == "__main__":
    main()