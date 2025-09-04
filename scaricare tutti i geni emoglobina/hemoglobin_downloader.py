#!/usr/bin/env python3
"""
Script completo per scaricare TUTTE le sequenze geniche del sistema emoglobinico umano
Autore: Assistant
Data: 2025
Requisiti: pip install biopython pandas requests
"""

import os
import time
import json
from datetime import datetime
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import requests

# Configurazione Entrez (IMPORTANTE: inserisci la tua email)
Entrez.email = "agdiegorusso@gmail.com"  # MODIFICA CON LA TUA EMAIL
Entrez.api_key = None  # Opzionale: inserisci API key NCBI per aumentare rate limit

class HemoglobinGeneDownloader:
    def __init__(self, output_dir="hemoglobin_genes"):
        """
        Inizializza il downloader con tutti i geni del sistema emoglobinico
        """
        self.output_dir = output_dir
        self.create_output_dirs()
        
        # Database completo di tutti i geni coinvolti
        self.gene_database = {
            "globin_structural": {
                "alpha_cluster": {
                    "HBZ": {"refseq": "NM_005332.3", "chr": "16p13.3", "desc": "Zeta globin (embryonic)"},
                    "HBM": {"refseq": "NM_001003938.3", "chr": "16p13.3", "desc": "Mu globin (embryonic)"},
                    "HBA2": {"refseq": "NM_000517.6", "chr": "16p13.3", "desc": "Alpha-2 globin (adult)"},
                    "HBA1": {"refseq": "NM_000558.5", "chr": "16p13.3", "desc": "Alpha-1 globin (adult)"},
                    "HBQ1": {"refseq": "NM_005331.4", "chr": "16p13.3", "desc": "Theta globin (embryonic)"}
                },
                "beta_cluster": {
                    "HBE1": {"refseq": "NM_005330.4", "chr": "11p15.4", "desc": "Epsilon globin (embryonic)"},
                    "HBG2": {"refseq": "NM_000184.3", "chr": "11p15.4", "desc": "Gamma-G globin (fetal)"},
                    "HBG1": {"refseq": "NM_000559.3", "chr": "11p15.4", "desc": "Gamma-A globin (fetal)"},
                    "HBD": {"refseq": "NM_000519.4", "chr": "11p15.4", "desc": "Delta globin (minor adult)"},
                    "HBB": {"refseq": "NM_000518.5", "chr": "11p15.4", "desc": "Beta globin (adult)"}
                }
            },
            "transcription_factors": {
                "primary": {
                    "GATA1": {"refseq": "NM_002049.4", "chr": "Xp11.23", "desc": "GATA binding protein 1"},
                    "KLF1": {"refseq": "NM_006563.5", "chr": "19p13.13", "desc": "Kruppel-like factor 1 (EKLF)"},
                    "TAL1": {"refseq": "NM_003189.3", "chr": "1p33", "desc": "TAL1/SCL transcription factor"},
                    "NFE2": {"refseq": "NM_001136023.3", "chr": "12q13.13", "desc": "Nuclear factor erythroid 2"},
                    "NFE2L2": {"refseq": "NM_006164.5", "chr": "2q31.2", "desc": "NFE2-like 2 (NRF2)"},
                    "ZFPM1": {"refseq": "NM_002625.3", "chr": "16q24.2", "desc": "FOG1 - GATA1 cofactor"}
                },
                "switch_regulators": {
                    "BCL11A": {"refseq": "NM_022718.4", "chr": "2p16.1", "desc": "BCL11 transcription factor A"},
                    "ZBTB7A": {"refseq": "NM_015898.3", "chr": "19p13.3", "desc": "LRF - fetal globin repressor"},
                    "MYB": {"refseq": "NM_001130172.2", "chr": "6q23.3", "desc": "c-Myb proto-oncogene"},
                    "SOX6": {"refseq": "NM_033326.4", "chr": "11p15.2", "desc": "SRY-box 6"},
                    "HBS1L": {"refseq": "NM_006620.3", "chr": "6q23.3", "desc": "HBS1-like translational GTPase"}
                },
                "additional": {
                    "LDB1": {"refseq": "NM_001113407.2", "chr": "10q24.32", "desc": "LIM domain binding 1"},
                    "GATA2": {"refseq": "NM_001145661.2", "chr": "3q21.3", "desc": "GATA binding protein 2"},
                    "TCF3": {"refseq": "NM_001136139.3", "chr": "19p13.3", "desc": "E2A transcription factor"},
                    "RUNX1": {"refseq": "NM_001754.5", "chr": "21q22.12", "desc": "Runt-related TF1"}
                }
            },
            "heme_biosynthesis": {
                "enzymes": {
                    "ALAS2": {"refseq": "NM_000032.5", "chr": "Xp11.21", "desc": "5-aminolevulinate synthase 2"},
                    "ALAD": {"refseq": "NM_000031.6", "chr": "9q32", "desc": "Delta-aminolevulinic acid dehydratase"},
                    "HMBS": {"refseq": "NM_000190.4", "chr": "11q23.3", "desc": "Hydroxymethylbilane synthase"},
                    "UROS": {"refseq": "NM_000375.3", "chr": "10q26.13", "desc": "Uroporphyrinogen III synthase"},
                    "UROD": {"refseq": "NM_000374.5", "chr": "1p34.1", "desc": "Uroporphyrinogen decarboxylase"},
                    "CPOX": {"refseq": "NM_000097.6", "chr": "3q11.2", "desc": "Coproporphyrinogen oxidase"},
                    "PPOX": {"refseq": "NM_000309.5", "chr": "1q23.3", "desc": "Protoporphyrinogen oxidase"},
                    "FECH": {"refseq": "NM_000140.5", "chr": "18q21.31", "desc": "Ferrochelatase"}
                },
                "regulation": {
                    "ALAS1": {"refseq": "NM_000688.6", "chr": "3p21.2", "desc": "5-aminolevulinate synthase 1"},
                    "ABCB6": {"refseq": "NM_005689.4", "chr": "2q35", "desc": "ATP-binding cassette B6"},
                    "SLC25A38": {"refseq": "NM_017875.3", "chr": "3p22.1", "desc": "Glycine transporter"}
                }
            },
            "iron_metabolism": {
                "transport": {
                    "TFR1": {"refseq": "NM_003234.4", "chr": "3q29", "desc": "Transferrin receptor 1"},
                    "TFR2": {"refseq": "NM_003227.3", "chr": "7q22.1", "desc": "Transferrin receptor 2"},
                    "SLC11A2": {"refseq": "NM_001174125.2", "chr": "12q13.12", "desc": "DMT1 - iron transporter"},
                    "SLC40A1": {"refseq": "NM_014585.6", "chr": "2q32.2", "desc": "Ferroportin"},
                    "SLC25A37": {"refseq": "NM_016612.4", "chr": "8p21.2", "desc": "Mitoferrin-1"},
                    "SLC25A28": {"refseq": "NM_001319043.2", "chr": "10q24.2", "desc": "Mitoferrin-2"}
                },
                "regulation": {
                    "HAMP": {"refseq": "NM_021175.4", "chr": "19q13.12", "desc": "Hepcidin"},
                    "HFE": {"refseq": "NM_000410.4", "chr": "6p22.2", "desc": "Hemochromatosis protein"},
                    "HJV": {"refseq": "NM_213653.4", "chr": "1q21.1", "desc": "Hemojuvelin"},
                    "BMP6": {"refseq": "NM_001718.6", "chr": "6p24.3", "desc": "Bone morphogenetic protein 6"},
                    "TMPRSS6": {"refseq": "NM_153609.4", "chr": "22q12.3", "desc": "Matriptase-2"},
                    "ACO1": {"refseq": "NM_002197.3", "chr": "9p21.1", "desc": "IRP1 - iron regulatory protein"}
                },
                "storage": {
                    "FTL": {"refseq": "NM_000146.4", "chr": "19q13.33", "desc": "Ferritin light chain"},
                    "FTH1": {"refseq": "NM_002032.3", "chr": "11q12.3", "desc": "Ferritin heavy chain"},
                    "FTMT": {"refseq": "NM_177478.3", "chr": "5q23.1", "desc": "Mitochondrial ferritin"}
                }
            },
            "chaperones_cofactors": {
                "hemoglobin_specific": {
                    "AHSP": {"refseq": "NM_016633.3", "chr": "16p11.2", "desc": "Alpha hemoglobin stabilizing protein"},
                    "GCLC": {"refseq": "NM_001498.4", "chr": "6p12.1", "desc": "Glutamate-cysteine ligase catalytic"},
                    "GCLM": {"refseq": "NM_002061.3", "chr": "1p22.1", "desc": "Glutamate-cysteine ligase modifier"},
                    "GSR": {"refseq": "NM_000637.5", "chr": "8p12", "desc": "Glutathione reductase"},
                    "G6PD": {"refseq": "NM_001360016.2", "chr": "Xq28", "desc": "Glucose-6-phosphate dehydrogenase"}
                },
                "protein_quality": {
                    "HSPA8": {"refseq": "NM_006597.6", "chr": "11q24.1", "desc": "Heat shock 70kDa protein 8"},
                    "HSP90AA1": {"refseq": "NM_005348.4", "chr": "14q32.31", "desc": "Heat shock protein 90 alpha"},
                    "DNAJB1": {"refseq": "NM_006145.3", "chr": "19p13.2", "desc": "DnaJ homolog B1"}
                }
            },
            "epigenetic_modifiers": {
                "chromatin": {
                    "KDM1A": {"refseq": "NM_001009999.3", "chr": "1p36.12", "desc": "Lysine demethylase 1A"},
                    "HDAC1": {"refseq": "NM_004964.3", "chr": "1p35.2", "desc": "Histone deacetylase 1"},
                    "HDAC2": {"refseq": "NM_001527.4", "chr": "6q21", "desc": "Histone deacetylase 2"},
                    "DNMT1": {"refseq": "NM_001130823.3", "chr": "19p13.2", "desc": "DNA methyltransferase 1"},
                    "DNMT3A": {"refseq": "NM_022552.5", "chr": "2p23.3", "desc": "DNA methyltransferase 3A"},
                    "DNMT3B": {"refseq": "NM_006892.4", "chr": "20q11.21", "desc": "DNA methyltransferase 3B"}
                }
            }
        }
        
    def create_output_dirs(self):
        """Crea le directory di output"""
        dirs = [
            self.output_dir,
            f"{self.output_dir}/fasta",
            f"{self.output_dir}/genbank",
            f"{self.output_dir}/analysis",
            f"{self.output_dir}/mutations"
        ]
        for dir_path in dirs:
            os.makedirs(dir_path, exist_ok=True)
    
    def download_sequence(self, gene_name, accession, retries=3):
        """
        Scarica una singola sequenza da NCBI
        """
        for attempt in range(retries):
            try:
                print(f"  Downloading {gene_name} ({accession})... ", end="")
                
                # Scarica il record GenBank completo
                handle = Entrez.efetch(
                    db="nucleotide", 
                    id=accession, 
                    rettype="gb", 
                    retmode="text"
                )
                record = SeqIO.read(handle, "genbank")
                handle.close()
                
                # Estrai la sequenza CDS
                cds_seq = None
                cds_location = None
                
                for feature in record.features:
                    if feature.type == "CDS":
                        cds_seq = feature.extract(record.seq)
                        cds_location = feature.location
                        break
                
                if cds_seq:
                    print(f"✓ ({len(cds_seq)} bp)")
                    return {
                        "gene": gene_name,
                        "accession": accession,
                        "full_sequence": str(record.seq),
                        "cds_sequence": str(cds_seq),
                        "cds_location": str(cds_location),
                        "record": record
                    }
                else:
                    print("✗ (No CDS found)")
                    return None
                    
            except Exception as e:
                print(f"✗ (Attempt {attempt+1}/{retries}: {str(e)})")
                if attempt < retries - 1:
                    time.sleep(2)  # Pausa tra i tentativi
                else:
                    return None
        
        return None
    
    def download_all_genes(self):
        """
        Scarica tutte le sequenze dal database
        """
        all_sequences = {}
        stats = {"total": 0, "success": 0, "failed": []}
        
        print("=" * 80)
        print("DOWNLOADING HEMOGLOBIN SYSTEM GENES")
        print("=" * 80)
        
        for category, subcategories in self.gene_database.items():
            print(f"\n[{category.upper()}]")
            all_sequences[category] = {}
            
            for subcat, genes in subcategories.items():
                print(f"\n  {subcat}:")
                all_sequences[category][subcat] = {}
                
                for gene_name, gene_info in genes.items():
                    stats["total"] += 1
                    result = self.download_sequence(gene_name, gene_info["refseq"])
                    
                    if result:
                        all_sequences[category][subcat][gene_name] = result
                        stats["success"] += 1
                        
                        # Salva la sequenza
                        self.save_sequence(result, category, subcat)
                    else:
                        stats["failed"].append(f"{gene_name} ({gene_info['refseq']})")
                    
                    time.sleep(0.5)  # Rispetta i rate limits di NCBI
        
        # Report finale
        print("\n" + "=" * 80)
        print(f"DOWNLOAD COMPLETED: {stats['success']}/{stats['total']} successful")
        if stats["failed"]:
            print(f"Failed downloads: {', '.join(stats['failed'])}")
        print("=" * 80)
        
        return all_sequences, stats
    
    def save_sequence(self, seq_data, category, subcategory):
        """
        Salva la sequenza in formato FASTA
        """
        # Crea il record FASTA per la sequenza CDS
        fasta_record = SeqRecord(
            Seq(seq_data["cds_sequence"]),
            id=f"{seq_data['gene']}|{seq_data['accession']}",
            description=f"{category}|{subcategory}|CDS"
        )
        
        # Salva in file FASTA individuale
        filename = f"{self.output_dir}/fasta/{seq_data['gene']}_cds.fasta"
        SeqIO.write([fasta_record], filename, "fasta")
        
        # Salva anche il record GenBank completo
        gb_filename = f"{self.output_dir}/genbank/{seq_data['gene']}.gb"
        SeqIO.write([seq_data["record"]], gb_filename, "genbank")
    
    def create_combined_fasta(self):
        """
        Crea un file FASTA combinato con tutte le sequenze
        """
        print("\nCreating combined FASTA file...")
        
        all_records = []
        for fasta_file in os.listdir(f"{self.output_dir}/fasta"):
            if fasta_file.endswith("_cds.fasta"):
                record = SeqIO.read(f"{self.output_dir}/fasta/{fasta_file}", "fasta")
                all_records.append(record)
        
        # Salva il file combinato
        combined_file = f"{self.output_dir}/all_hemoglobin_genes_cds.fasta"
        SeqIO.write(all_records, combined_file, "fasta")
        print(f"Combined FASTA saved: {combined_file} ({len(all_records)} sequences)")
    
    def generate_mutation_analysis_template(self):
        """
        Genera un template per l'analisi delle mutazioni
        """
        print("\nGenerating mutation analysis template...")
        
        template = {
            "metadata": {
                "created": datetime.now().isoformat(),
                "total_genes": sum(
                    len(genes) 
                    for cat in self.gene_database.values() 
                    for genes in cat.values()
                ),
                "categories": list(self.gene_database.keys())
            },
            "genes": {}
        }
        
        # Aggiungi informazioni per ogni gene
        for category, subcats in self.gene_database.items():
            for subcat, genes in subcats.items():
                for gene_name, info in genes.items():
                    template["genes"][gene_name] = {
                        "category": category,
                        "subcategory": subcat,
                        "refseq": info["refseq"],
                        "chromosome": info["chr"],
                        "description": info["desc"],
                        "known_mutations": [],
                        "analysis_pending": True
                    }
        
        # Salva il template
        template_file = f"{self.output_dir}/analysis/mutation_analysis_template.json"
        with open(template_file, 'w') as f:
            json.dump(template, f, indent=2)
        print(f"Template saved: {template_file}")
    
    def analyze_sequence_properties(self):
        """
        Analizza le proprietà basiche delle sequenze scaricate
        """
        print("\nAnalyzing sequence properties...")
        
        analysis_results = []
        
        for fasta_file in os.listdir(f"{self.output_dir}/fasta"):
            if fasta_file.endswith("_cds.fasta"):
                gene_name = fasta_file.replace("_cds.fasta", "")
                record = SeqIO.read(f"{self.output_dir}/fasta/{fasta_file}", "fasta")
                
                seq = record.seq
                
                # Calcola proprietà
                analysis_results.append({
                    "gene": gene_name,
                    "length_bp": len(seq),
                    "length_aa": len(seq) // 3,
                    "gc_content": round((seq.count("G") + seq.count("C")) / len(seq) * 100, 2),
                    "start_codon": str(seq[:3]),
                    "stop_codon": str(seq[-3:]),
                    "has_valid_start": str(seq[:3]) == "ATG",
                    "has_valid_stop": str(seq[-3:]) in ["TAA", "TAG", "TGA"]
                })
        
        # Salva i risultati
        df = pd.DataFrame(analysis_results)
        df = df.sort_values("gene")
        
        # Salva come CSV
        csv_file = f"{self.output_dir}/analysis/sequence_properties.csv"
        df.to_csv(csv_file, index=False)
        print(f"Analysis saved: {csv_file}")
        
        # Stampa statistiche
        print(f"\nSequence Statistics:")
        print(f"  Total genes analyzed: {len(df)}")
        print(f"  Average CDS length: {df['length_bp'].mean():.0f} bp")
        print(f"  Average protein length: {df['length_aa'].mean():.0f} aa")
        print(f"  Average GC content: {df['gc_content'].mean():.1f}%")
        print(f"  Valid start codons: {df['has_valid_start'].sum()}/{len(df)}")
        print(f"  Valid stop codons: {df['has_valid_stop'].sum()}/{len(df)}")
        
        return df

def main():
    """
    Funzione principale
    """
    print("╔" + "═" * 78 + "╗")
    print("║" + " HEMOGLOBIN GENE SYSTEM DOWNLOADER ".center(78) + "║")
    print("║" + " Complete Database of All Related Genes ".center(78) + "║")
    print("╚" + "═" * 78 + "╝")
    
    # Verifica email Entrez
    if Entrez.email == "your.email@example.com":
        print("\n⚠️  WARNING: Please set your email in the script (line 16)")
        email = input("Enter your email for NCBI Entrez: ")
        if "@" in email:
            Entrez.email = email
        else:
            print("Invalid email. Exiting.")
            return
    
    # Inizializza il downloader
    downloader = HemoglobinGeneDownloader()
    
    # Download delle sequenze
    print(f"\nStarting download at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    sequences, stats = downloader.download_all_genes()
    
    # Crea file combinato
    downloader.create_combined_fasta()
    
    # Genera template per analisi mutazionali
    downloader.generate_mutation_analysis_template()
    
    # Analizza proprietà delle sequenze
    properties_df = downloader.analyze_sequence_properties()
    
    # Report finale
    print("\n" + "=" * 80)
    print("PROCESS COMPLETED SUCCESSFULLY!")
    print("=" * 80)
    print(f"\nOutput directory: {downloader.output_dir}/")
    print(f"  - Individual FASTA files: {downloader.output_dir}/fasta/")
    print(f"  - GenBank files: {downloader.output_dir}/genbank/")
    print(f"  - Combined FASTA: {downloader.output_dir}/all_hemoglobin_genes_cds.fasta")
    print(f"  - Analysis results: {downloader.output_dir}/analysis/")
    
    print("\nNext steps:")
    print("  1. Use the combined FASTA for multiple alignment")
    print("  2. Analyze mutations using the template")
    print("  3. Import into your mutation prediction pipeline")
    print("\nFor mutation analysis, you can use:")
    print("  - SIFT: https://sift.bii.a-star.edu.sg/")
    print("  - PolyPhen-2: http://genetics.bwh.harvard.edu/pph2/")
    print("  - MutPred: http://mutpred.mutdb.org/")
    print("  - PROVEAN: http://provean.jcvi.org/")

if __name__ == "__main__":
    main()