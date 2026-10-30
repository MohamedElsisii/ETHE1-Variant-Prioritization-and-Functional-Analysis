# ETHE1 Variant Prioritization and Functional Analysis

<img width="1536" height="1024" alt="Untitled-1" src="https://github.com/user-attachments/assets/dc63d64b-1863-47ae-b500-6173163306c3" />

An integrated computational framework for the annotation, prioritization, and functional assessment of **ETHE1** variants associated with ethylmalonic encephalopathy (EE). The project combines transcript-level variant processing, protein sequence mutagenesis, functional annotation, CADD-based scoring, and structural proximity analysis to support systematic investigation of ETHE1 variation.

## Research Team

**Horizon Science Communication LLC — Bioinformatics Team**

- **Mohamed Elsisi**
- **Dr. Hesham Abdullah**

## Project Background

Ethylmalonic encephalopathy is a rare metabolic disorder associated with pathogenic variants in **ETHE1**, which encodes mitochondrial persulfide dioxygenase. Characterizing the potential functional consequences of ETHE1 variants requires the integration of sequence-level, functional, computational, and structural evidence.

This project establishes a computational workflow to organize ETHE1 variant information, generate mutant protein sequences, and prioritize variants using a combination of molecular consequence, functional annotation, CADD scores, and protein structural context.

**Reference gene:** [ETHE1 — NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/23474)  
**Reference protein:** [UniProt O95571](https://www.uniprot.org/uniprotkb/O95571/entry)  
**Reference transcript:** `NM_014297.5`

## Project Objective

The objective of this project is to develop an integrated computational framework for ETHE1 variant annotation and prioritization.

The workflow is designed to:

- Organize ETHE1 variants using transcript-level HGVS cDNA descriptions.
- Use the reference ETHE1 transcript to generate and annotate mutant protein sequences.
- Integrate variant consequences and functional annotations.
- Incorporate CADD PHRED scores into a rule-based prioritization model.
- Assess structural proximity to active-site and iron-binding residues.
- Support the interpretation of prioritized variants using supplied structural analysis results.

## Project Overview

| Attribute | Description |
|---|---|
| Target gene | ETHE1 |
| Associated disorder | Ethylmalonic encephalopathy |
| Reference transcript | NM_014297.5 |
| Reference protein | UniProt O95571 |
| Variant records in prioritization dataset | 364 |
| Variant representation | HGVS cDNA |
| Protein sequence analysis | R / Biostrings |
| Variant prioritization | Rule-based scoring model |
| Structural proximity analysis | PyMOL / Python |
| CADD annotation | CADD PHRED scores |
| Structure-related results | AlphaFold prediction and truncated-variant analyses |

The repository contains reference data, analysis scripts, and selected result tables. It does not contain every external resource or intermediate file used during the wider research process.

## Why Integrate Sequence, Functional, and Structural Evidence?

### Variant annotation and protein mutagenesis

Transcript-level HGVS descriptions provide a standardized representation of nucleotide variants. Applying these variants to the reference transcript allows the corresponding mutant protein sequences and protein-level consequences to be examined.

The R mutagenesis script uses the full reference transcript, including untranslated regions, and the annotated coding sequence coordinates. This is relevant for variants that may affect translation or transcript regions beyond the coding sequence.

### Rule-based variant prioritization

A single annotation or prediction score does not capture every aspect of a variant's potential functional relevance. The prioritization model therefore combines multiple evidence categories, including molecular consequence, active-site and iron-binding annotations, epitope annotation, CADD PHRED score, initiator methionine position, and second-order structural positions.

### Structural proximity analysis

The PyMOL script evaluates side-chain distances between selected active-site or iron-binding residues and other residues in the protein structure. This provides an additional structural-context measure for interpreting residue positions.

## Workflow Overview

The project follows the analysis sequence below:

```text
ETHE1 reference protein and transcript
                 |
                 v
Protein localization, topology
and signal peptide prediction
                 |
                 v
Functional annotation
                 |
                 v
STRING functional interaction network
                 |
                 v
ETHE1 variant collection
                 |
                 v
HGVS standardization
                 |
                 v
CADD annotation
                 |
                 v
Mutant protein sequence generation
                 |
                 v
Rule-based variant prioritization
                 |
                 v
Clinical significance-based validation
                 |
                 v
Selection of High and Very High variants
                 |
                 v
Truncated variant analysis
                 |
                 v
AlphaFold structural prediction
                 |
                 v
GSSH substrate preparation
```

The workflow diagram is intended to summarize the study's analysis sequence. Not every stage above is implemented as an executable script in this repository; some stages rely on external databases, prediction servers, or separately performed analyses.

## Repository Structure

```text
ETHE1-Variant-Prioritization-and-Functional-Analysis/
├── README.md
├── LICENSE
├── .gitignore
│
├── data/
│   └── reference/
│       ├── ETHE1_transcript_NM_014297.5.fasta
│       └── ETHE1_hgvs_cdna.txt
│
├── scripts/
│   ├── R/
│   │   ├── ETHE1 Variant Mutagenesis Script.R
│   │   └── Risk Model Scores Script.R
│   │
│   └── python/
│       └── PyMol Proximity Residues Script.py
│
└── results/
    ├── cadd/
    │   └── CADD Results.tsv
    │
    ├── prioritization/
    │   ├── ETHE1 Variants Risk Model Results.csv
    │   └── ETHE1_NEW_MODEL.csv
    │
    └── structural/
        ├── ETHE1 Predicted AlphaFold Structures pTM & Validation Results.csv
        └── ETHE1 Truncated Variants Analysis.csv
```

The repository includes selected result tables and reference inputs. Raw structural files, complete external-tool workspaces, and all intermediate files are not included.

## Highlights

- **Transcript-aware mutagenesis:** Generates mutant protein sequences from transcript-level HGVS cDNA variants using the reference transcript.
- **Variant-level quality checks:** The mutagenesis functions include reference-sequence checks and a QC summary.
- **Integrated prioritization:** Combines consequence, functional annotation, CADD, and residue-position evidence in a revised scoring model.
- **Structural context:** Evaluates proximity to active-site and iron-binding residues using PyMOL.
- **Retained analysis outputs:** Includes CADD, prioritization, AlphaFold validation, and truncated-variant result tables.

## Results Showcase

### 1. Variant prioritization

The revised prioritization model was applied to 364 variant records.

| Priority category | Variant count |
|---|---:|
| Very High | 19 |
| High | 116 |
| Moderate | 210 |
| Low | 19 |
| **Total** | **364** |

The revised model classified 135 variants as High or Very High. These are computational prioritization categories and should not be interpreted as clinical classifications.

The original and revised prioritization tables are both retained in the repository to distinguish the earlier results from the revised scoring model.

### 2. CADD annotation

The CADD results are provided in:

`results/cadd/CADD Results.tsv`

The revised scoring model uses CADD PHRED score bands as one component of the overall prioritization score.

### 3. Structural analysis

The repository includes result tables documenting predicted AlphaFold structure pTM and validation metrics, as well as truncated-variant analysis:

- `results/structural/ETHE1 Predicted AlphaFold Structures pTM & Validation Results.csv`
- `results/structural/ETHE1 Truncated Variants Analysis.csv`

These are supplied analysis outputs. The repository does not include all structure files or an automated workflow to regenerate every structural prediction.

## Analysis Scripts

### 1. ETHE1 variant mutagenesis

**Script:** `scripts/R/ETHE1 Variant Mutagenesis Script.R`

This R file defines functions for processing transcript-level HGVS cDNA variants, applying variants to the reference transcript, generating mutant protein sequences, and producing variant annotations and QC outputs.

The script uses:

- Reference transcript: `NM_014297.5`
- CDS start: nucleotide 25
- CDS end: nucleotide 789

**R dependencies:**

- Biostrings
- dplyr
- readr
- stringr
- tibble

The function library can be sourced from the repository root:

```r
source("scripts/R/ETHE1 Variant Mutagenesis Script.R")
```

The main function is `run_ethe1_pipeline()`. It accepts a transcript FASTA file, a variants file, an output directory, and optional CDS and HGVS-column settings.

Example invocation:

```r
run_ethe1_pipeline(
  transcript_fasta = "data/reference/ETHE1_transcript_NM_014297.5.fasta",
  variants_file = "data/reference/ETHE1_hgvs_cdna.txt",
  outdir = "results/mutagenesis"
)
```

Expected outputs include:

- `ETHE1_variant_sequences.csv`
- `ETHE1_QC_summary.csv`
- `ETHE1_AlphaFold.fasta`
- `ETHE1_run_info.txt`

### 2. Revised ETHE1 variant prioritization model

**Script:** `scripts/R/Risk Model Scores Script.R`

This standalone R script reads the existing variant annotation and prioritization table, calculates revised scores, assigns priority categories, and exports the revised results.

**Input:**

```text
results/prioritization/ETHE1 Variants Risk Model Results.csv
```

**Output:**

```text
results/prioritization/ETHE1_NEW_MODEL.csv
```

**R dependencies:**

- dplyr
- readr

Run from the repository root:

```bash
Rscript "scripts/R/Risk Model Scores Script.R"
```

#### Scoring scheme

| Scoring component | Rule | Score |
|---|---|---:|
| Frameshift variant | `frameshift_variant` | 5 |
| Stop gained | `stop_gained` | 5 |
| Start lost | `start_lost` | 4 |
| Stop lost | `stop_lost` | 4 |
| In-frame insertion/deletion | `inframe_insertion`, `inframe_deletion` | 4 |
| Missense variant | `missense_variant` | 2 |
| Synonymous variant | `synonymous_variant` | 0 |
| Active-site annotation | `Active_Site = TRUE` | 4 |
| Iron-binding annotation | `Iron_Binding = TRUE` | 4 |
| Epitope annotation | `Epitope = TRUE` | 1 |
| CADD PHRED | ≥ 29 | 5 |
| CADD PHRED | ≥ 20 and < 29 | 4 |
| CADD PHRED | ≥ 10 and < 20 | 3 |
| CADD PHRED | < 10 | 0 |
| Initiator methionine | Position 1 | 3 |
| Second-order structural positions | Position-specific lists in the script | 0–2 |

Priority categories are assigned using the total revised score:

| Total score | Priority |
|---:|---|
| ≥ 12 | Very High |
| 8–11 | High |
| 5–7 | Moderate |
| 0–4 | Low |

### 3. PyMOL structural proximity analysis

**Script:** `scripts/python/PyMol Proximity Residues Script.py`

This script evaluates side-chain contacts between selected ETHE1 active-site or iron-binding residues and other protein residues.

The script uses:

- Chain A.
- Nine active-site reference residues: 79, 83, 135, 154, 157, 161, 162, 163, and 197.
- Three iron-binding reference residues: 79, 135, and 154.
- Maximum contact distance: 6 Å.
- Strong-contact threshold: 4 Å.

Contacts at or below 4 Å receive a structural score of 2. Contacts above 4 Å and up to 6 Å receive a score of 1.

**Execution requirements:**

1. Open the appropriate ETHE1 structure in PyMOL.
2. Confirm that the structure uses chain A.
3. Confirm that its residue numbering matches the reference positions used in the script.
4. Run the script in the PyMOL command line:

```python
run scripts/python/PyMol Proximity Residues Script.py
```

The script prints the active-site network, iron-binding network, combined structural network, and final score summary in the PyMOL console. It does not currently export a separate results file.

## Software

| Tool / package | Role |
|---|---|
| R | Variant sequence processing and prioritization |
| Biostrings | DNA and protein sequence handling |
| dplyr | Data manipulation and scoring |
| readr | Tabular data import and export |
| Python | Structural analysis scripting |
| PyMOL | Protein structure inspection and atom-coordinate access |
| CADD | Computational deleteriousness annotation |
| AlphaFold | Protein structure prediction |

External resources and software may require separate installation, access, or configuration. Versions and execution settings should be recorded when reproducing analyses.

## Reproducibility and Scope

The scripts in this repository cover selected parts of the overall research workflow. In particular:

- The R mutagenesis file is a function library and requires the supplied reference transcript and variant input file.
- The revised risk model is executable from the repository root and reads the supplied original prioritization table.
- The PyMOL script requires a compatible ETHE1 structure to be loaded before execution.
- CADD annotation, external functional predictions, AlphaFold structure prediction, and truncated-variant analyses are represented by supplied result tables; the repository does not currently automate every upstream external analysis.
- Exact reproduction of external analyses may require the original database versions, server settings, software versions, and input structures.

## Interpretation and Limitations

The prioritization score is a research-oriented evidence-integration framework. It is not a validated clinical classification system and should not be used alone to establish pathogenicity, diagnosis, or treatment decisions.

Structural proximity is based on the loaded structure, chain selection, residue numbering, and distance thresholds defined in the PyMOL script. Results may differ if a different structure or residue mapping is used.

## Contributors and Affiliation

This project is associated with **Horizon Science Communication LLC**.

### Contributors

**Mohamed Elsisi**
- Email: `M.Elsisi@horizon-sci-comm.us`
- GitHub: [MohamedElsisii](https://github.com/MohamedElsisii)
- LinkedIn: [Mohamed Elsisi](https://www.linkedin.com/in/mohamed-elsisii)

**Dr. Hesham Abdullah**
- Email: `ceo@horizon-sci-comm.us`
- GitHub: [HABIOTECH](https://github.com/HABIOTECH)
- LinkedIn: [Dr. Hesham Abdullah](https://www.linkedin.com/in/hesham-abdullah-40bab1190/)

### Organization

**Horizon Science Communication LLC**

<img width="1448" height="1086" alt="Screenshot 2026-06-21 211033" src="https://github.com/user-attachments/assets/beb83a84-6581-4c7e-a142-bc490b616b6f" />

[Website](https://www.horizon-sci-comm.us/)

This repository documents the computational methods, selected analysis outputs, and research workflow developed for ETHE1 variant prioritization and functional analysis.

---

For questions about the workflow or collaboration, please contact the contributors using the information above.
