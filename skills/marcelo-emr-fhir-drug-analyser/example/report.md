# 🏥 FHIR EMR Pharmacogenomic Medication Review

**Generated**: 2026-03-19 11:40  
**FHIR Server**: https://hapi.fhir.org/baseR4  
**Patient ID**: 131299660  

---

## Patient Context

| Field | Value |
|---|---|
| Name | Hungdung Vo |
| Date of Birth | 1992-03-02 |
| Gender | male |

## 🟢 Standard PGx Profile

| Medication | Gene | Phenotype | Note |
|---|---|---|---|
| simvastatin 40mg | SLCO1B1 | — | No actionable PGx concern |

## ⚪ Medications Without PGx Coverage

| Medication | Status |
|---|---|
| acetaminophen 325 mg oral tablet | current |
| ibuprofen 200 mg oral tablet [advil] | current |
| amoxicillin 200 mg oral tablet [amoxi-tabs] | current |
| 20 ml morphine sulfate 10 mg/ml injection [infumorph] | current |
| salbutamol 100micrograms/inhaler | current |
| aspirin 81mg e/c tablet | current |
| beclomethasone dipropionate 40mcg aerosol | current |
| carvedilol 25 mg | current |
| chlortalidone 15 mg | current |
| gabapentin 600mg | current |
| insulin glargine 100units/ml injection | current |
| losartan 100mg + hydrochlorothiazide 12.5mg | current |
| metformin 1000mg | current |
| naproxen 500mg e/c tablet | current |
| prednisone 20mg | current |
| terbinafine 250mg | current |
| zolpidem 5mg sublingual tablet | current |
| lisinopril 20mg | current |
| isosorbide mononitrate 20mg | current |
| buspirone 5mg | current |

## ⚠️ Allergy Conflicts

The following medications appear in both the prescription list and the patient's AllergyIntolerance record:

- **ibuprofen 200 mg oral tablet [advil]** — review allergy record before any PGx guidance

## Notes for Clinician

- Medications marked 🔴 have a pharmacogenomic signal suggesting potential harm or reduced efficacy.
- Medications marked 🟡 may benefit from dose review or closer monitoring.
- Candidate medication alternatives, if any, should be confirmed against the full clinical picture including allergies, renal/hepatic function, and co-morbidities before any prescribing decision.

## Limitations

- PGx analysis covers only the genes and SNPs included in the pharmgx-reporter database (12 genes, ~31 SNPs).
- Medication names are normalized from FHIR coded entries; brand-to-generic mapping may be incomplete.
- ClinPGx evidence is retrieved only for flagged gene-drug pairs; other medications may have unpublished PGx associations.
- FHIR data completeness depends on what the source EMR has recorded.

---

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*