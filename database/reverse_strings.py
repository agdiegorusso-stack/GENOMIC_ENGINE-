import re

# Script di reverse engineering: estrazione avanzata di stringhe e pattern dal file binario
FILE_PATH = 'database/zephyr.db'
OUT_PATH = 'zephyr_strings.txt'

# Pattern per stringhe ASCII stampabili di almeno 6 caratteri
ascii_pattern = re.compile(rb'[\x20-\x7E]{6,}')

# Pattern per delimitatori binari noti
delimiters = [b'fe fe fe fe', b'eXc\r', b'Bio-Rad', b'V2_BThal', b'Arial', b'Courier New']

with open(FILE_PATH, 'rb') as f, open(OUT_PATH, 'w') as out:
    data = f.read()
    # Estrazione stringhe ASCII
    for match in ascii_pattern.finditer(data):
        out.write(match.group().decode(errors='ignore') + '\n')
    out.write('\n--- Delimitatori trovati ---\n')
    # Ricerca pattern/delimitatori
    for delim in delimiters:
        idx = 0
        while True:
            idx = data.find(delim, idx)
            if idx == -1:
                break
            out.write(f"Delimitatore '{delim.decode(errors='ignore')}' trovato a offset {idx}\n")
            idx += len(delim)
print(f"Estrazione completata. Risultati in {OUT_PATH}")
