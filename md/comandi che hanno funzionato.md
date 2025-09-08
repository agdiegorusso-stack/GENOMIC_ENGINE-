 ~/ensembl-vep/vep \
      -i data/raw/clinvar_hemoglobin_exact.vcf \
      -o data/processed/clinvar_annotated_final.vcf \
      --vcf --offline --dir_cache ~/.vep --force_overwrite \
      --symbol --hgvsp --hgvs --variant_class \
      --fasta ~/.vep/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz