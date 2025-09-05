#!/usr/bin/env python3

import requests
import time
import json
import pandas as pd
from bs4 import BeautifulSoup
from urllib.parse import urljoin, urlencode
import re
import os

class IthaGenesDatabase:
    """
    Interfaccia corretta per IthaGenes con 3.519+ varianti
    """
    
    def __init__(self, cache_dir="ithagenes_cache"):
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        self.base_url = "https://www.ithanet.eu/db/ithagenes"
        self.session = requests.Session()
        # Headers per evitare blocchi
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
            'Accept-Language': 'en-US,en;q=0.5',
            'Accept-Encoding': 'gzip, deflate',
            'Connection': 'keep-alive'
        })
    
    def download_all_variants(self, max_variants=None):
        """
        Scarica tutte le varianti da IthaGenes (3.519+ varianti)
        """
        print("Scaricando tutte le varianti da IthaGenes...")
        
        all_variants = []
        page = 1
        variants_per_page = 50  # IthaGenes mostra ~50 per pagina
        
        while True:
            print(f"  Scaricando pagina {page}...")
            
            # URL per la lista paginata
            params = {
                'action': 'list',
                'page': page
            }
            
            try:
                response = self.session.get(self.base_url, params=params, timeout=30)
                
                if response.status_code == 200:
                    variants = self._parse_variants_page(response.text)
                    
                    if not variants:  # Nessuna variante = fine delle pagine
                        break
                    
                    all_variants.extend(variants)
                    print(f"    Trovate {len(variants)} varianti (totale: {len(all_variants)})")
                    
                    # Controllo limite opzionale
                    if max_variants and len(all_variants) >= max_variants:
                        all_variants = all_variants[:max_variants]
                        break
                    
                    page += 1
                    time.sleep(2)  # Rate limiting più conservativo
                    
                else:
                    print(f"    Errore HTTP {response.status_code}")
                    if response.status_code == 500:
                        print("    Server error - aspetto 10 secondi...")
                        time.sleep(10)
                        continue
                    break
                    
            except requests.exceptions.RequestException as e:
                print(f"    Errore connessione: {e}")
                time.sleep(5)
                continue
            
            except Exception as e:
                print(f"    Errore: {e}")
                break
        
        print(f"\nTotale varianti scaricate da IthaGenes: {len(all_variants)}")
        return all_variants
    
    def _parse_variants_page(self, html_content):
        """
        Parse delle varianti dalla pagina HTML di IthaGenes
        """
        soup = BeautifulSoup(html_content, 'html.parser')
        variants = []
        
        # Cerca la tabella delle varianti
        tables = soup.find_all('table')
        
        for table in tables:
            rows = table.find_all('tr')
            
            # Identifica l'header per capire la struttura
            header_row = None
            for i, row in enumerate(rows):
                if any('gene' in cell.get_text().lower() or 'variant' in cell.get_text().lower() 
                       for cell in row.find_all(['th', 'td'])):
                    header_row = i
                    break
            
            if header_row is not None:
                headers = [th.get_text(strip=True).lower() for th in rows[header_row].find_all(['th', 'td'])]
                
                for row in rows[header_row + 1:]:
                    cells = row.find_all('td')
                    
                    if len(cells) >= 3:  # Minimo 3 colonne
                        variant = {}
                        
                        for i, cell in enumerate(cells):
                            if i < len(headers):
                                header = headers[i]
                                text = cell.get_text(strip=True)
                                
                                # Link details se presente
                                link = cell.find('a')
                                if link and link.get('href'):
                                    variant[f'{header}_link'] = urljoin(self.base_url, link.get('href'))
                                
                                variant[header] = text
                        
                        if variant:  # Solo se ha contenuto
                            variants.append(variant)
        
        return variants
    
    def get_variant_details(self, variant_id):
        """
        Ottiene dettagli completi per una specifica variante
        """
        if not variant_id:
            return None
        
        url = f"{self.base_url}?ithaID={variant_id}"
        
        try:
            response = self.session.get(url, timeout=30)
            
            if response.status_code == 200:
                return self._parse_variant_detail(response.text)
                
        except Exception as e:
            print(f"Errore getting details for {variant_id}: {e}")
        
        return None
    
    def _parse_variant_detail(self, html_content):
        """
        Parse dettagli di una singola variante
        """
        soup = BeautifulSoup(html_content, 'html.parser')
        details = {}
        
        # Cerca informazioni strutturate
        for label in soup.find_all(['dt', 'th', 'strong']):
            label_text = label.get_text(strip=True).lower()
            
            # Trova il valore associato
            value_elem = label.find_next_sibling(['dd', 'td']) or label.find_next()
            if value_elem:
                value = value_elem.get_text(strip=True)
                if value:
                    details[label_text] = value
        
        return details

    def search_by_gene(self, gene_name):
        """
        Cerca varianti per un gene specifico
        """
        print(f"Cercando varianti per gene: {gene_name}")
        
        params = {
            'action': 'search',
            'gene': gene_name
        }
        
        try:
            response = self.session.get(self.base_url, params=params, timeout=30)
            
            if response.status_code == 200:
                variants = self._parse_variants_page(response.text)
                print(f"  Trovate {len(variants)} varianti per {gene_name}")
                return variants
            else:
                print(f"  Errore HTTP {response.status_code}")
                
        except Exception as e:
            print(f"  Errore: {e}")
        
        return []

# Classe HbVar corretta
class HbVarDatabase:
    """
    Interfaccia corretta per HbVar con gestione errori migliorata
    """
    
    def __init__(self, cache_dir="hbvar_cache"):
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        self.base_url = "https://globin.bx.psu.edu/hbvar"
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (compatible; Research Bot 1.0)',
            'Accept': 'text/html,application/xhtml+xml',
        })
    
    def download_all_variants(self):
        """
        Tentativo di download da HbVar con fallback su LOVD
        """
        print("Tentando download da HbVar...")
        
        variants = []
        
        # Prova URL principale
        urls_to_try = [
            f"{self.base_url}/menu.html",
            f"{self.base_url}/",
            "https://lovd.bx.psu.edu/"  # Fallback su LOVD
        ]
        
        for url in urls_to_try:
            try:
                print(f"  Provando: {url}")
                response = self.session.get(url, timeout=30)
                
                if response.status_code == 200:
                    print(f"   Connessione riuscita")
                    # Implementa parsing specifico per questo URL
                    page_variants = self._parse_hbvar_content(response.text, url)
                    variants.extend(page_variants)
                    
                    if variants:
                        break
                        
                else:
                    print(f"   HTTP {response.status_code}")
                    
            except Exception as e:
                print(f"   Errore: {e}")
                continue
        
        print(f"Varianti HbVar trovate: {len(variants)}")
        return variants
    
    def _parse_hbvar_content(self, html_content, source_url):
        """
        Parse adattivo del contenuto HbVar
        """
        soup = BeautifulSoup(html_content, 'html.parser')
        variants = []
        
        # Cerca link a pagine di dati
        data_links = []
        for link in soup.find_all('a', href=True):
            href = link.get('href')
            if any(keyword in href.lower() for keyword in ['variant', 'mutation', 'data', 'search']):
                full_url = urljoin(source_url, href)
                data_links.append(full_url)
        
        # Prova a scaricare da questi link
        for link in data_links[:5]:  # Limita a 5 per test
            try:
                print(f"    Esplorando: {link}")
                response = self.session.get(link, timeout=20)
                
                if response.status_code == 200:
                    link_variants = self._extract_variants_from_page(response.text)
                    variants.extend(link_variants)
                    
                time.sleep(1)
                
            except Exception as e:
                print(f"    Errore su {link}: {e}")
                continue
        
        return variants
    
    def _extract_variants_from_page(self, html_content):
        """
        Estrae varianti da una pagina generica
        """
        soup = BeautifulSoup(html_content, 'html.parser')
        variants = []
        
        # Cerca tabelle con pattern di varianti
        tables = soup.find_all('table')
        
        for table in tables:
            rows = table.find_all('tr')
            
            for row in rows:
                cells = row.find_all(['td', 'th'])
                
                if len(cells) >= 2:
                    text_content = ' '.join(cell.get_text(strip=True) for cell in cells)
                    
                    # Pattern per varianti di emoglobina
                    if re.search(r'[αβγδ]\d+|Hb[A-Z]|[A-Z]\d+[A-Z]', text_content):
                        variant = {
                            'raw_text': text_content,
                            'cells': [cell.get_text(strip=True) for cell in cells]
                        }
                        variants.append(variant)
        
        return variants

# Funzione principale migliorata
def main():
    """
    Esegue download integrato con gestione errori robusta
    """
    print("="*80)
    print("DOWNLOAD DATABASE EMOGLOBINOPATIE - VERSIONE CORRETTA")
    print("="*80)
    
    # Inizializza database
    ithagenes = IthaGenesDatabase()
    hbvar = HbVarDatabase()
    
    # Scarica da IthaGenes (fonte principale)
    print("\n1. DOWNLOAD DA ITHAGENES (3.519+ varianti)")
    ithagenes_variants = ithagenes.download_all_variants(max_variants=4000)
    
    # Tenta HbVar come fonte secondaria
    print("\n2. DOWNLOAD DA HBVAR (fonte secondaria)")
    hbvar_variants = hbvar.download_all_variants()
    
    # Combina risultati
    total_variants = len(ithagenes_variants) + len(hbvar_variants)
    print(f"\n3. RIEPILOGO:")
    print(f"   IthaGenes: {len(ithagenes_variants)} varianti")
    print(f"   HbVar: {len(hbvar_variants)} varianti")
    print(f"   TOTALE: {total_variants} varianti")
    
    if total_variants == 0:
        print("\n NESSUNA VARIANTE SCARICATA")
        print("Possibili soluzioni:")
        print("1. Verifica connessione internet")
        print("2. I siti potrebbero essere temporaneamente down")
        print("3. Prova più tardi")
        return
    
    # Salva risultati
    print(f"\n4. SALVATAGGIO DATI...")
    
    # Salva IthaGenes
    if ithagenes_variants:
        ithagenes_df = pd.DataFrame(ithagenes_variants)
        ithagenes_df.to_csv('ithagenes_variants.csv', index=False)
        ithagenes_df.to_json('ithagenes_variants.json', orient='records', indent=2)
        print(f"    IthaGenes salvato: ithagenes_variants.csv/json")
    
    # Salva HbVar
    if hbvar_variants:
        hbvar_df = pd.DataFrame(hbvar_variants)
        hbvar_df.to_csv('hbvar_variants.csv', index=False)
        print(f"    HbVar salvato: hbvar_variants.csv")
    
    print(f"\n DOWNLOAD COMPLETATO!")
    
    # Esempi di query
    if ithagenes_variants:
        print(f"\n5. ESEMPI DI VARIANTI TROVATE:")
        for i, variant in enumerate(ithagenes_variants[:5]):
            print(f"   {i+1}. {variant}")

if __name__ == "__main__":
    main()
