mindmap
  root((Risorse del Progetto Hemoglobin Analysis))
    Librerie Python
      pandas
        ::icon(fa:database)
        Manipolazione dati e gestione CSV
      cyvcf2
        ::icon(fa:dna)
        Parsing di file VCF per dati genomici
      scikit-learn
        ::icon(fa:brain)
        Machine learning per analisi predittive
      re
        ::icon(fa:search)
        Espressioni regolari per pattern matching
    Librerie di Sistema Linux
      libhts-dev
        ::icon(fa:file-code)
        Libreria HTS per cyvcf2
      libcurl4-gnutls-dev
        ::icon(fa:globe)
        Libreria per connessioni HTTP
      libdeflate-dev
        ::icon(fa:compress)
        Libreria per compressione dati
      libc6-dev
        ::icon(fa:cogs)
        Libreria C standard
      zlib1g-dev
        ::icon(fa:archive)
        Libreria per compressione zlib
    Dipendenze Python
      coloredlogs
        ::icon(fa:palette)
        Logging colorato
      click
        ::icon(fa:terminal)
        Interfaccia a linea di comando
      numpy
        ::icon(fa:calculator)
        Calcoli numerici
      scipy
        ::icon(fa:chart-line)
        Scienza e ingegneria
      joblib
        ::icon(fa:tasks)
        Parallelizzazione
      threadpoolctl
        ::icon(fa:users)
        Controllo thread pool
    File di Dati Esterni
      clinvar.vcf.gz
        ::icon(fa:file-archive)
        Database ClinVar compresso (fonte principale)
      clinvar.vcf.gz.tbi
        ::icon(fa:index)
        Indice per ClinVar VCF
      clinvar_hemoglobin.vcf
        ::icon(fa:dna)
        VCF filtrato per geni dell'emoglobina
      hemoglobin_loci.bed
        ::icon(fa:map-marker)
        Loci genomici per geni dell'emoglobina
    Script e Codice
      build_training_set.py
        ::icon(fa:python)
        Script per costruire il dataset di training
    Output del Modello
      clinvar_training_set.csv
        ::icon(fa:table)
        Dataset di training generato
