flowchart TD
    A[Inizio: Avvio dell'App di Analisi Emoglobina] --> B[Caricamento Dati ClinVar VCF]
    B --> C[Filtro Varianti Missense nei Geni dell'Emoglobina]
    C --> D[Etichettatura: Patogeniche/Benigne]
    D --> E[Estrazione Informazioni HGVS e Gene]
    E --> F[Salvataggio Dataset di Training in CSV]
    F --> G[Fine: Dataset Pronto per ML]

    H[Obiettivo: Predire Patogenicità Varianti Genetiche]
    I[Importanza: Diagnosi Malattie Genetiche Emoglobina]
    J[Applicazioni: Sickle Cell, Talassemia, Anemie]

    A --> H
    H --> I
    I --> J

    style A fill:#e1f5fe
    style G fill:#c8e6c9
    style H fill:#fff3e0
    style I fill:#ffebee
    style J fill:#f3e5f5
