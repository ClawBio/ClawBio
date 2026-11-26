# Drug Photo Demo Images

Drag these images into an agent chat with a vision model for live demos.
Results are checked against Manuel's real 23andMe genotype data.

## Priority Demo Drugs

| File | Drug (Generic) | Brand | Gene | Expected Result |
|------|----------------|-------|------|-----------------|
| clopidogrel.jpg | Clopidogrel | Plavix | CYP2C19 | Depends on genotype |
| codeine.jpg | Codeine | Tylenol #3 | CYP2D6 | Depends on genotype |
| warfarin.jpg | Warfarin | Coumadin | CYP2C9+VKORC1 | Multi-gene |
| sertraline.jpg | Sertraline | Zoloft | CYP2C19 | Depends on genotype |
| omeprazole.jpg | Omeprazole | Prilosec | CYP2C19 | Depends on genotype |
| simvastatin.jpg | Simvastatin | Zocor | SLCO1B1 | Depends on genotype |
| tamoxifen.jpg | Tamoxifen | Nolvadex | CYP2D6 | Depends on genotype |
| fluoxetine.jpg | Fluoxetine | Prozac | CYP2D6 | Depends on genotype |
| metoprolol.jpg | Metoprolol | Lopressor | CYP2D6 | Depends on genotype |
| tacrolimus.jpg | Tacrolimus | Prograf | CYP3A5 | Depends on genotype |

## Image Sources

Download from DailyMed (dailymed.nlm.nih.gov) — the official FDA drug label database.
Search by drug name, download the highest-quality pill/packaging image.

## Hackathon Demo Script (2 min)

1. Show Sertraline 50mg box photo → card shows genotype-specific guidance
2. Show Plavix 75mg photo → card contextualises dose against CYP2C19
3. Show codeine photo → dramatic result if CYP2D6 is poor metabolizer
4. Pitch: "51 drugs, 12 genes, 31 SNPs, checked against my actual genome — all local, FDA pharmacogenomic biomarkers, running in 3 seconds"
