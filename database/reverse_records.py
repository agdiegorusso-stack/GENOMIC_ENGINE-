import struct
import re

FILE_PATH = 'database/zephyr.db'
OUT_PATH = 'zephyr_records.txt'

# Pattern per stringhe ASCII stampabili di almeno 6 caratteri
ascii_pattern = re.compile(rb'[\x20-\x7E]{6,}')

# Funzione per estrarre possibili record binari
# Cerca blocchi che iniziano con 'eXc\r' e contengono stringhe ASCII

def extract_complex_records():
    with open(FILE_PATH, 'rb') as f, open(OUT_PATH, 'w') as out:
        data = f.read()
        offset = 0
        while offset < len(data):
            # Cerca header 'eXc\r' (come inizio record)
            idx = data.find(b'eXc\r', offset)
            if idx == -1:
                break
            # Prova a leggere 512 byte dal record
            record = data[idx:idx+512]
            # Estrai tutte le stringhe ASCII dal record
            strings = ascii_pattern.findall(record)
            out.write(f'Record trovato a offset {idx} (lunghezza {len(record)}):\n')
            for s in strings:
                out.write('  ' + s.decode(errors='ignore') + '\n')
            out.write('\n')
            offset = idx + 512
    print(f"Estrazione record avanzata completata. Risultati in {OUT_PATH}")

if __name__ == '__main__':
    extract_complex_records()
