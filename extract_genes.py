import json

# Carica il file JSON
with open('GENOMIC_ENGINE--main/hemoglobin_genes/IthaGenes_variations_export.json', 'r') as f:
    data = json.load(f)

# Conta le variazioni totali
total_variants = len(data["body"])
print(f"Numero totale di variazioni: {total_variants}")

# Conta variazioni con "N/A" nel campo Genes
na_count = 0
for variant in data["body"]:
    if len(variant) > 4 and variant[4] == "N/A":
        na_count += 1

print(f"Variazioni con 'N/A' nel campo Genes: {na_count}")

# Estrai i nomi dei geni unici (escludendo N/A)
genes = set()
for variant in data["body"]:
    if len(variant) > 4:  # Assicurati che ci siano abbastanza elementi
        gene = variant[4]  # Campo "Genes" è il quinto elemento (indice 4)
        if gene and gene != "N/A":
            genes.add(gene)

print(f"Numero di geni unici (escludendo N/A): {len(genes)}")
print(f"Nomi dei geni unici: {sorted(list(genes))}")