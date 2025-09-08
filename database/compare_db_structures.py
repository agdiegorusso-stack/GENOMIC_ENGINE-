import struct
import re

DB_PATH = 'database/zephyr.db'
BAK_PATH = '26092022_10_11_Zephyr.dbbak'
HEADER_PATTERNS = [
    b'ObJeCtStOrEdAtAbAsEbAcKuPiMaGe',
    b'ObJeCtStOrEaRcHiVeImAgE'
]

# Cerca header e blocchi dati

def find_headers(data, patterns):
    found = []
    for pat in patterns:
        idx = 0
        while True:
            idx = data.find(pat, idx)
            if idx == -1:
                break
            found.append(idx)
            idx += len(pat)
    return found

def extract_floats_doubles(chunk):
    floats, doubles = [], []
    for i in range(0, len(chunk)-4, 4):
        try:
            val = struct.unpack('<f', chunk[i:i+4])[0]
            if 0 < abs(val) < 1e8:
                floats.append(val)
        except Exception:
            continue
    for i in range(0, len(chunk)-8, 8):
        try:
            val = struct.unpack('<d', chunk[i:i+8])[0]
            if 0 < abs(val) < 1e8:
                doubles.append(val)
        except Exception:
            continue
    return floats, doubles

def main():
    with open(DB_PATH, 'rb') as f:
        db_data = f.read()
    with open(BAK_PATH, 'rb') as f:
        bak_data = f.read()

    db_headers = find_headers(db_data, HEADER_PATTERNS)
    bak_headers = find_headers(bak_data, HEADER_PATTERNS)

    print(f"zephyr.db: trovati {len(db_headers)} header: {db_headers[:10]}")
    print(f"dbbak: trovati {len(bak_headers)} header: {bak_headers[:10]}")

    # Analizza blocchi dopo header
    for idx in db_headers[:5]:
        chunk = db_data[idx:idx+256]
        floats, doubles = extract_floats_doubles(chunk)
        print(f"zephyr.db @ {idx}: floats={floats[:5]}, doubles={doubles[:5]}")
    for idx in bak_headers[:5]:
        chunk = bak_data[idx:idx+256]
        floats, doubles = extract_floats_doubles(chunk)
        print(f"dbbak @ {idx}: floats={floats[:5]}, doubles={doubles[:5]}")

if __name__ == '__main__':
    main()
