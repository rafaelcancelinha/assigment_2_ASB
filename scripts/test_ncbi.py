#!/usr/bin/env python3
import sys
import subprocess
import requests
import xml.etree.ElementTree as ET
import time 

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def esearch_sra(accession):
    url = EUTILS + "esearch.fcgi"
    params = {
        "db": "sra",
        "term": accession,
        "retmode": "xml"
    }
    r = requests.get(url, params=params)
    r.raise_for_status()
    root = ET.fromstring(r.text)
    ids = [id_elem.text for id_elem in root.findall(".//Id")]
    return ids

def efetch_runs(uid):
    url = EUTILS + "efetch.fcgi"
    params = {
        "db": "sra",
        "id": uid,
        "rettype": "runinfo",
        "retmode": "text"
    }
    r = requests.get(url, params=params)
    r.raise_for_status()
    lines = r.text.splitlines()
    runs = []
    for line in lines:
        if line.startswith("Run,"):
            continue
        parts = line.split(",")
        if parts and parts[0].strip():
            runs.append(parts[0].strip())
    return runs

def download_fastq(run_id):
    print(f"   Downloading FASTQ for {run_id} ...")
    subprocess.run(["fasterq-dump", run_id], check=True)
    print(f"   Done: {run_id}")

    # Compress all .fastq files produced by this run to .fastq.gz
    for suffix in ["_1.fastq", "_2.fastq", ".fastq"]:
        fastq_file = f"{run_id}{suffix}"
        gz_file    = f"{fastq_file}.gz"
        import os
        if os.path.exists(fastq_file):
            print(f"   Compressing {fastq_file} -> {gz_file} ...")
            subprocess.run(["gzip", fastq_file], check=True)
            print(f"   Compressed: {gz_file}")

def process_accession(accession):
    print(f"\n{'='*50}")
    print(f"Processing accession: {accession}")
    print(f"{'='*50}")

    uids = esearch_sra(accession)
    if not uids:
        print(f"   [!] No SRA entries found for {accession}, skipping.")
        return

    print(f"   Found {len(uids)} SRA entry/entries")

    all_runs = []
    for uid in uids:
        runs = efetch_runs(uid)
        all_runs.extend(runs)
        time.sleep(1) 

    if not all_runs:
        print(f"   [!] No Run IDs found for {accession}, skipping.")
        return

    print(f"   Run IDs: {all_runs}")

    for run in all_runs:
        download_fastq(run)

def main():
    if len(sys.argv) != 2:
        print("Usage: python sra_fetch.py <accessions.txt>")
        sys.exit(1)

    txt_file = sys.argv[1]

    try:
        with open(txt_file, "r") as f:
            accessions = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Error: File '{txt_file}' not found.")
        sys.exit(1)

    if not accessions:
        print("Error: No accessions found in file.")
        sys.exit(1)

    print(f"Loaded {len(accessions)} accession(s) from '{txt_file}'")

    for i, accession in enumerate(accessions, 1):
        print(f"\n[{i}/{len(accessions)}]", end="")
        process_accession(accession)
        
        if i < len(accessions):
            print("   [Info] A aguardar 3 segundos para evitar bloqueio do NCBI...")
            time.sleep(3)

    print(f"\n{'='*50}")
    print("All accessions processed.")

if __name__ == "__main__":
    main()
