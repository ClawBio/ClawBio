#!/usr/bin/env python3
"""
Gene Set Enrichment Analysis Skill
Performs enrichment analysis using Enrichr with multiple pathway libraries
"""

import argparse
import json
import os
import sys
from pathlib import Path
import requests
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Tuple

class EnrichrEnrichment:
    """Perform gene set enrichment analysis using Enrichr API"""
    
    # Multiple pathway libraries for comprehensive analysis
    LIBRARIES = [
        "KEGG_2021_Human",
        "GO_Biological_Process_2021",
        "GO_Molecular_Function_2021"
    ]
    
    ENRICHR_URL = "https://maayanlab.cloud/Enrichr"
    
    def __init__(self, output_dir: str = "."):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results_df = None
        
    def load_genes(self, gene_file: str) -> List[str]:
        """Load genes from file (one per line)"""
        with open(gene_file, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        print(f"Loaded {len(genes)} genes from {gene_file}")
        return genes
    
    def post_to_enrichr(self, genes: List[str]) -> int:
        """Upload gene list to Enrichr and return list ID"""
        gene_list_str = "\n".join(genes)
        
        url = f"{self.ENRICHR_URL}/addList"
        files = {
            "list": (None, gene_list_str),
            "description": (None, "Gene set enrichment analysis")
        }
        
        print(f"Uploading {len(genes)} genes to Enrichr...")
        response = requests.post(url, files=files)
        response.raise_for_status()
        
        result = response.json()
        user_list_id = result.get("userListId")
        print(f"Successfully uploaded gene list. ID: {user_list_id}")
        return user_list_id
    
    def fetch_results(self, user_list_id: int, library: str) -> List[List]:
        """Fetch enrichment results for a specific library"""
        url = f"{self.ENRICHR_URL}/enrich?userListId={user_list_id}&backgroundType={library}"
        
        print(f"Fetching results for {library}...")
        response = requests.get(url)
        response.raise_for_status()
        
        result = response.json()
        return result.get(library, [])
    
    def parse_results(self, raw_results: List[List], library: str) -> pd.DataFrame:
        """
        Parse raw Enrichr results into DataFrame
        
        Each result entry format:
        [rank, term_name, p_value, combined_score, overlap, overlap_genes, 
         adjusted_p_value, z_score, ...]
        """
        parsed_data = []
        
        for result in raw_results:
            parsed_data.append({
                'Library': library,
                'Term': result[1],
                'P-value': result[2],
                'Combined Score': result[3],
                'Overlap': result[4],
                'Genes': ','.join(result[5]) if len(result) > 5 else '',
                'Adjusted P-value': result[6] if len(result) > 6 else None,
                'Z-score': result[7] if len(result) > 7 else None,
            })
        
        return pd.DataFrame(parsed_data)
    
    def fetch_all_libraries(self, user_list_id: int) -> pd.DataFrame:
        """Fetch and combine results from all libraries"""
        all_results = []
        
        for library in self.LIBRARIES:
            try:
                raw_results = self.fetch_results(user_list_id, library)
                df = self.parse_results(raw_results, library)
                all_results.append(df)
                print(f"Retrieved {len(df)} terms from {library}")
            except Exception as e:
                print(f"Warning: Failed to fetch {library}: {e}")
        
        # Combine all results
        if all_results:
            combined_df = pd.concat(all_results, ignore_index=True)
            
            # Sort by adjusted p-value (or p-value if adjusted not available)
            sort_col = 'Adjusted P-value' if 'Adjusted P-value' in combined_df.columns else 'P-value'
            combined_df = combined_df.dropna(subset=[sort_col]).sort_values(by=sort_col)
            
            print(f"\nCombined {len(combined_df)} results from {len(self.LIBRARIES)} libraries")
            return combined_df
        else:
            raise ValueError("No results retrieved from any library")
    
    def generate_chart(self, df: pd.DataFrame, top_n: int = 15):
        """Generate visualization of top enriched terms"""
        # Get top N results by adjusted p-value
        top_df = df.head(top_n).copy()
        top_df['Log P-value'] = -np.log10(top_df['P-value'])
        top_df = top_df.sort_values('Log P-value', ascending=True)
        
        # Create figure with subplots
        fig, axes = plt.subplots(1, 2, figsize=(16, 8))
        
        # Plot 1: P-value by term
        colors = {'KEGG_2021_Human': '#1f77b4', 
                  'GO_Biological_Process_2021': '#ff7f0e',
                  'GO_Molecular_Function_2021': '#2ca02c'}
        bar_colors = [colors.get(lib, '#999999') for lib in top_df['Library']]
        
        axes[0].barh(range(len(top_df)), top_df['Log P-value'], color=bar_colors)
        axes[0].set_yticks(range(len(top_df)))
        axes[0].set_yticklabels(top_df['Term'], fontsize=9)
        axes[0].set_xlabel('-log10(P-value)', fontsize=11)
        axes[0].set_title(f'Top {top_n} Enriched Terms by P-value', fontsize=12, fontweight='bold')
        axes[0].grid(axis='x', alpha=0.3)
        
        # Plot 2: Combined Score by term
        top_df_score = df.head(top_n).sort_values('Combined Score', ascending=True)
        bar_colors_score = [colors.get(lib, '#999999') for lib in top_df_score['Library']]
        
        axes[1].barh(range(len(top_df_score)), top_df_score['Combined Score'], color=bar_colors_score)
        axes[1].set_yticks(range(len(top_df_score)))
        axes[1].set_yticklabels(top_df_score['Term'], fontsize=9)
        axes[1].set_xlabel('Combined Score', fontsize=11)
        axes[1].set_title(f'Top {top_n} Enriched Terms by Combined Score', fontsize=12, fontweight='bold')
        axes[1].grid(axis='x', alpha=0.3)
        
        plt.tight_layout()
        chart_file = self.output_dir / "enrichment_chart.png"
        plt.savefig(chart_file, dpi=300, bbox_inches='tight')
        print(f"Chart saved to {chart_file}")
        plt.close()
        
        return chart_file
    
    def save_results(self, df: pd.DataFrame, filename: str = "enrichment_results.csv"):
        """Save results to CSV"""
        output_file = self.output_dir / filename
        df.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")
        return output_file
    
    def generate_report(self, genes: List[str], df: pd.DataFrame, chart_file: Path):
        """Generate markdown report"""
        report_file = self.output_dir / "report.md"
        
        with open(report_file, 'w') as f:
            f.write("# Gene Set Enrichment Analysis Report\n\n")
            
            f.write("## Summary\n")
            f.write(f"- **Input genes**: {len(genes)}\n")
            f.write(f"- **Unique genes analyzed**: {len(set(genes))}\n")
            f.write(f"- **Pathways identified**: {len(df)}\n")
            f.write(f"- **Libraries queried**: {', '.join(self.LIBRARIES)}\n\n")
            
            f.write("## Top 10 Enriched Pathways\n\n")
            f.write("| Library | Term | P-value | Combined Score | Overlap | Key Genes |\n")
            f.write("|---------|------|---------|-----------------|---------|----------|\n")
            
            for idx, row in df.head(10).iterrows():
                genes_preview = ', '.join(row['Genes'].split(',')[:3])
                if len(row['Genes'].split(',')) > 3:
                    genes_preview += ", ..."
                f.write(f"| {row['Library']} | {row['Term']} | {row['P-value']:.2e} | "
                       f"{row['Combined Score']:.2f} | {row['Overlap']} | {genes_preview} |\n")
            
            f.write(f"\n## Full Results\n")
            f.write(f"See `enrichment_results.csv` for complete results ({len(df)} pathways)\n\n")
            
            f.write(f"## Visualization\n")
            f.write(f"![Enrichment Chart](enrichment_chart.png)\n\n")
            
            f.write("## Methods\n")
            f.write("Gene set enrichment analysis performed using Enrichr (Kuleshov et al., 2016) "
                   "against the following pathway libraries:\n")
            for lib in self.LIBRARIES:
                f.write(f"- {lib}\n")
            f.write("\nResults combined and sorted by adjusted p-value.\n")
        
        print(f"Report saved to {report_file}")
        return report_file
    
    def run(self, gene_file: str):
        """Execute full enrichment analysis workflow"""
        try:
            # Load genes
            genes = self.load_genes(gene_file)
            
            # Upload to Enrichr
            user_list_id = self.post_to_enrichr(genes)
            
            # Fetch and combine results from all libraries
            df = self.fetch_all_libraries(user_list_id)
            
            # Save results
            self.save_results(df)
            
            # Generate chart
            chart_file = self.generate_chart(df)
            
            # Generate report
            self.generate_report(genes, df, chart_file)
            
            print(f"\n✓ Analysis complete. Results saved to {self.output_dir}")
            return True
            
        except Exception as e:
            print(f"✗ Error during analysis: {e}", file=sys.stderr)
            raise

def main():
    parser = argparse.ArgumentParser(
        description="Gene Set Enrichment Analysis using Enrichr"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input gene list file (one gene per line)"
    )
    parser.add_argument(
        "--output", "-o",
        default="enrichment_results",
        help="Output directory for results (default: enrichment_results)"
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run with demo cancer gene set"
    )
    
    args = parser.parse_args()
    
    # Handle demo mode
    if args.demo:
        demo_genes = [
            "TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "MYC", "APC", "RB1",
            "PTEN", "PIK3CA", "BRAF", "NRAS", "CDH1", "VHL", "WT1"
        ]
        demo_file = Path(args.output) / "demo_genes.txt"
        demo_file.parent.mkdir(parents=True, exist_ok=True)
        with open(demo_file, 'w') as f:
            f.write("\n".join(demo_genes))
        args.input = str(demo_file)
        print("Running demo mode with cancer genes...")
    
    # Run analysis
    enricher = EnrichrEnrichment(output_dir=args.output)
    enricher.run(args.input)

if __name__ == "__main__":
    import numpy as np
    main()
