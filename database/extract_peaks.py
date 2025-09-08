
import re
import struct
import csv

FILE_PATH = 'database/zephyr.db'
OUT_PATH = 'zephyr_peaks.csv'

# Nomi dei picchi e campi da cercare (come nello screenshot)
peak_names = [b'Unknown', b'P2', b'P3', b'A0', b'A2', b'S-window']
sample_fields = [b'Sample ID', b'Lot Number', b'Injection Time']

# Pattern per stringhe ASCII
ascii_pattern = re.compile(rb'[\x20-\x7E]{3,}')

def extract_peaks():
    with open(FILE_PATH, 'rb') as f, open(OUT_PATH, 'w', newline='') as out:
        data = f.read()
        writer = csv.writer(out)
        writer.writerow(['Peak Name', 'RT', 'Area', 'Area %', 'Concentration', 'Offset', 'Sample ID', 'Lot Number', 'Injection Time'])
        for peak in peak_names:
            idx = 0
            while True:
                idx = data.find(peak, idx)
                if idx == -1:
                    break
                chunk = data[idx:idx+128]
                floats = []
                doubles = []
                # Prova a leggere float e double
                for i in range(0, len(chunk)-4, 4):
                    try:
                        val = struct.unpack('<f', chunk[i:i+4])[0]
                        if 0 < val < 1e7:
                            floats.append(val)
                    except Exception:
                        continue
                for i in range(0, len(chunk)-8, 8):
                    try:
                        val = struct.unpack('<d', chunk[i:i+8])[0]
                        if 0 < val < 1e7:
                            doubles.append(val)
                    except Exception:
                        continue
                # Cerca Sample ID, Lot Number, Injection Time
                sample_id = lot_number = injection_time = ''
                for field in sample_fields:
                    fidx = data.find(field, max(0, idx-256), idx+256)
                    if fidx != -1:
                        # Prendi la stringa dopo il campo
                        end = data.find(b'\x00', fidx)
                        value = data[fidx+len(field):end].decode(errors='ignore').strip()
                        if field == b'Sample ID':
                            sample_id = value
                        elif field == b'Lot Number':
                            lot_number = value
                        elif field == b'Injection Time':
                            injection_time = value
                # Log dettagliato
                print(f"Peak: {peak.decode()} @ {idx} | floats: {floats[:4]} | doubles: {doubles[:4]} | SampleID: {sample_id} | Lot: {lot_number} | Time: {injection_time}")
                # Scrivi solo se trovi almeno 3 valori
                if len(floats) >= 3 or len(doubles) >= 3:
                    writer.writerow([
                        peak.decode(),
                        *(floats[:4] if len(floats)>=4 else ['']*4),
                        idx,
                        sample_id,
                        lot_number,
                        injection_time
                    ])
                idx += len(peak)
    print(f"Estrazione picchi completata. Risultati in {OUT_PATH}")

if __name__ == '__main__':
    extract_peaks()
