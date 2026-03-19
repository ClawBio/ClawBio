#!/usr/bin/env python3
"""
Fetch Enrichr enrichment results for cancer gene list
"""

import requests
import json

# Gene list ID from previous upload
user_list_id = 124669542

# Enrichr API endpoint for fetching results
url = f"https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType=KEGG_2021_Human"

try:
    response = requests.get(url)
    response.raise_for_status()
    
    result = response.json()
    
    # Pretty print the results
    print(json.dumps(result, indent=2))
    
    # If there are KEGG results, display them in a more readable format
    if "results" in result and "KEGG_2021_Human" in result["results"]:
        print("\n" + "="*80)
        print("KEGG 2021 Human Pathway Enrichment Results")
        print("="*80)
        
        kegg_results = result["results"]["KEGG_2021_Human"]
        for i, pathway in enumerate(kegg_results[:10], 1):  # Show top 10
            print(f"\n{i}. {pathway[1]}")
            print(f"   p-value: {pathway[2]:.2e}")
            print(f"   Combined Score: {pathway[3]:.2f}")
            print(f"   Genes: {', '.join(pathway[5][:5])}{'...' if len(pathway[5]) > 5 else ''}")
    
except requests.exceptions.RequestException as e:
    print(f"Error fetching results from Enrichr: {e}")
    print(f"Response: {response.text if 'response' in locals() else 'No response'}")
