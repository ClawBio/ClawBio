#!/usr/bin/env python3
"""
Post cancer gene list to Enrichr API
"""

import requests

# Gene list
genes = [
    "TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "MYC", "APC", "RB1", "PTEN", 
    "PIK3CA", "BRAF", "NRAS", "CDH1", "VHL", "WT1", "NF1", "NF2", "RET", 
    "KIT", "PDGFRA", "ALK", "ERBB2", "FGFR3", "IDH1", "IDH2", "NPM1", "FLT3", 
    "DNMT3A", "TET2", "JAK2", "MPL", "CALR", "SF3B1", "ASXL1", "EZH2", 
    "NOTCH1", "FBXW7", "CTNNB1", "SMAD4", "STK11", "CDKN2A", "CDK4", "MDM2", 
    "MET", "ROS1", "NTRK1", "MAP2K1", "ARID1A", "KMT2A", "CREBBP"
]

# Create newline-separated string
gene_list_str = "\n".join(genes)

# Enrichr API endpoint
url = "https://maayanlab.cloud/Enrichr/addList"

# Prepare data for POST request (multipart/form-data)
files = {
    "list": (None, gene_list_str),
    "description": (None, "Cancer-associated genes")
}

try:
    response = requests.post(url, files=files)
    response.raise_for_status()
    
    result = response.json()
    print(f"Success! List ID: {result.get('userListId')}")
    print(f"Short URL: {result.get('shortUrl')}")
    print(f"Full response: {result}")
    
except requests.exceptions.RequestException as e:
    print(f"Error posting to Enrichr: {e}")
    print(f"Response: {response.text if 'response' in locals() else 'No response'}")
