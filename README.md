# NLRcraft

**NLRcraft** is a structure-based pipeline for NLR identification and classification.


## ✨ Features

- Structure-based domain detection
- Quality control of predicted structures
- Machine-learning-assisted false positive removal
- Structural classification of NLRs based on N-terminal domain clustering


## 🧭 Pipeline Overview

NLRcraft consists of **four major steps**:

1. **Domain database building**  
   Build a reference database for domain detection.

2. **pLDDT-based structure filtering**  
   Remove low-confidence regions from predicted protein structures.

3. **NLR identification by structural alignment**  
   - Perform domain-level structural alignments
   - Infer protein-level NLR status (N / T / TN / Na)
   - Remove false positives using a Random Forest model

4. **NLR classification by structural clustering**  
   - Split N-terminal domains
   - Cluster N-terminal domains for classification


## 📂 Directory Structure

```text
NLRcraft/
├── NLRcraft.py
├── domains_pdb/               # Reference NLR domain structures
├── scripts/
│   ├── run_plddt_filter.py
│   ├── extract_align_results.py
│   ├── rm_FPs.py
│   ├── extract_NBS_pos.py
│   ├── split_pdb_by_NBS.py
│   └── R/
│       ├── rf_predict.R
│       └── final_rf_model_undersampling.rds
```
## 🚀 Usage

```text
python NLRcraft.py \
  --structs predicted_structures \
  --ids protein_ids.txt \
  --plddt 60
```

## 📄 Key Output Files

| File                     | Description                                        |
| ------------------------ | -------------------------------------------------- |
| `aln_all.txt`            | Raw structural alignment results        |
| `aln_filtered.txt`       | Filtered structural alignment results                   |
| `results_all.txt`        | NLR inference (N / T / TN / Na)          |
| `results_all_rm_FPs.txt` | NLR inference results after false positive removal |
| `NLR_annotation.tsv`     | Classification of NLRs                   |

#### Notes on NLR Inference

Foldseek was used to perform structure-based alignments between target protein structures and a curated reference database of NLR-related domains. The alignment results capture structural similarity at the domain level.

NLR inference was subsequently performed by integrating structural hits to the TIR and NB-ARC domains. Each protein was classified into one of the following categories:

- **N**: NLR proteins containing a structurally supported NB-ARC domain.
- **T**: TX proteins lacking an NB-ARC domain but carrying a TIR domain.
- **TN**: TIR-NLR proteins containing both TIR and NB-ARC domains.
- **Na**: Non-NLR proteins lacking structural support for NLR-related domains.

## 🔧 Dependencies

This is a pipeline to be run on unix based machines. The following software must be available in your path. 

-   [Foldseek](https://github.com/steineggerlab/foldseek)
-   Python ≥ 3.8
-   R version ≥ 4.4.1
-   Python package
    -   [biopython](https://github.com/biopython/biopython)
    -   [pandas](https://github.com/pandas-dev/pandas)
-   R package
    -   [tidyverse](https://www.tidyverse.org/)
    -   [randomForest](https://cran.r-project.org/web/packages/randomForest/randomForest.pdf)


## 📖 Citation

If you use **NLRcraft** in your research, please cite:

> Jiang, Q. *et al.*  
> **Structural conservation reveals cryptic homology and reshapes the evolutionary history of plant NLR immune receptors**.  
> *Manuscript in preparation.*

## 📬 Contact
Qian Jiang
Email: jiangqian@hnu.edu.cn