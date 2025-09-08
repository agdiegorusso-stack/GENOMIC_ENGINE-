# in scripts/deep_inspect_vcf.py
from cyvcf2 import VCF

def deep_inspect_variants(vcf_path, num_variants_to_inspect=20):
    """
    Esegue un'ispezione approfondita delle prime N varianti missense
    e stampa i loro campi INFO critici per il debugging.
    """
    print(f"Eseguendo un'ispezione approfondita di: {vcf_path}\n")
    vcf_reader = VCF(vcf_path)
    
    inspected_count = 0
    
    for variant in vcf_reader:
        info = variant.INFO
        mc_info = info.get('MC', '')
        
        # Controlliamo se è una variante missense
        is_missense = 'missense_variant' in mc_info or 'SO:0001583' in mc_info
        
        if is_missense:
            inspected_count += 1
            
            print(f"--- Variante Missense Trovata #{inspected_count} ({variant.CHROM}:{variant.POS}) ---")
            
            # Estraiamo e stampiamo i campi di interesse
            clnsig_value = info.get('CLNSIG')
            clnrevstat_value = info.get('CLNREVSTAT')
            
            # Stampa il valore e il TIPO di dato, che è cruciale
            print(f"CLNSIG       -> Valore: {clnsig_value} | Tipo: {type(clnsig_value)}")
            print(f"CLNREVSTAT   -> Valore: {clnrevstat_value} | Tipo: {type(clnrevstat_value)}")
            print(f"MC (completo)-> Valore: {mc_info}")
            print("-" * 50)
            
            if inspected_count >= num_variants_to_inspect:
                break # Fermiamoci dopo averne ispezionate abbastanza

    if inspected_count == 0:
        print("Non è stata trovata nessuna variante missense per l'ispezione.")
    else:
        print(f"\nIspezione completata. Mostrate le prime {inspected_count} varianti missense.")

if __name__ == "__main__":
    deep_inspect_variants('data/raw/clinvar_hemoglobin.vcf')