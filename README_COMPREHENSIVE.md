# ACF Immune Gene Expression Analysis

**Complete Reconstruction & Publication Package**

---

## 📋 Project Overview

This repository contains the complete analysis pipeline for investigating immune-related gene expression changes in aberrant crypt foci (ACF) - the earliest detectable precancerous lesions in colorectal cancer.

### Key Information

- **Project Title:** The Epithelial-Stromal Microenvironment in Early Colonic Neoplasia
- **Lead Investigators:** Dr. Daniel W. Rosenberg, Dr. Charles Giardina, Boyang Li, Takayasu Ideta
- **Institution:** University of Connecticut School of Medicine
- **Status:** Manuscript in preparation (85% complete)
- **Target Journals:** Clinical Cancer Research, Cancer Research

### Study Design

- **Patients:** 10 with dysplastic proximal colon ACF
- **Samples:** 40 total (10 patients × 2 conditions × 2 compartments)
  - ACF vs Matched Normal
  - Epithelial vs Stromal compartments
- **Technology:** Ion Torrent sequencing with Oncomine Immune Response Assay
- **Genes Profiled:** ~395 immune-related genes

### Main Findings

1. **PD1 is the most upregulated gene** in ACF epithelium (FC=-4.24, p<0.001)
2. **Chemokine signaling activated** in both compartments with distinct patterns
3. **Epithelial compartment:** Neutrophil chemotaxis genes upregulated (CXCL8, CXCL1, S100A8/A9)
4. **Stromal compartment:** CXCL13, CXCR4, and immune regulatory genes upregulated
5. **Early immune evasion:** PD1+ intraepithelial lymphocytes suggest immunosuppression

---

## 📁 Repository Structure

```
ACF-Immune-Analysis/
├── README.md                          # This file
├── todo.md                            # Project task list
├── LICENSE                            # MIT License
│
├── docs/                              # Documentation
│   ├── COMPREHENSIVE_PROJECT_ANALYSIS.md  # Complete project analysis
│   ├── DATA_DICTIONARY.md             # Data variable descriptions
│   └── ANALYSIS_GUIDE.md              # Step-by-step analysis guide
│
├── data/                              # Data files
│   ├── raw/                           # Raw sequencing data (not included)
│   ├── processed/                     # Processed data files
│   │   ├── ACF normalized epi.csv     # Normalized expression (epithelial)
│   │   ├── ACF normalized stro.csv    # Normalized expression (stromal)
│   │   ├── DGE_epithelial.xlsx        # Differential expression (epithelial)
│   │   ├── DGE_stromal.xlsx           # Differential expression (stromal)
│   │   ├── PCA_cell.csv               # PCA sample coordinates
│   │   └── PCA_loading.csv            # PCA gene loadings
│   └── metadata/                      # Sample and gene metadata
│       └── Gene panel.xlsx            # Gene annotations
│
├── scripts/                           # Analysis scripts
│   ├── ACF_master_analysis.R          # Complete analysis pipeline
│   ├── figure_generation.R            # Publication figure scripts
│   └── helper_functions.R             # Utility functions
│
├── analysis/                          # Analysis outputs
│   └── outputs/
│       ├── figures/                   # Generated figures
│       ├── tables/                    # Generated tables
│       ├── data/                      # Intermediate data
│       └── jmp_export/                # JMP-compatible exports
│
├── manuscript/                        # Manuscript files
│   ├── drafts/                        # Manuscript versions
│   └── figures/                       # Publication figures
│       ├── publication/               # Main figures
│       └── supplementary/             # Supplementary figures
│
└── references/                        # Reference materials
    ├── Gene-List_Oncomine-Immune-Response-Assay.pdf
    └── [Other reference PDFs]
```

---

## 🚀 Quick Start

### Prerequisites

**R Environment:**
- R version 4.0+ (tested on 3.6.3 and 4.3.0)
- RStudio (recommended)

**Required R Packages:**
```r
# CRAN packages
install.packages(c("tidyverse", "readxl", "writexl", "data.table",
                   "ggplot2", "pheatmap", "ggrepel", "cowplot",
                   "RColorBrewer", "viridis", "ggsci",
                   "scales", "reshape2", "gridExtra"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "edgeR", "limma",
                       "clusterProfiler", "enrichplot", "DOSE",
                       "org.Hs.eg.db", "ComplexHeatmap"))
```

### Running the Analysis

1. **Clone the repository:**
```bash
git clone https://github.com/RustemKurkcu/ACF-Immune-Analysis.git
cd ACF-Immune-Analysis
```

2. **Open R/RStudio and set working directory:**
```r
setwd("/path/to/ACF-Immune-Analysis")
```

3. **Run the master analysis script:**
```r
source("scripts/ACF_master_analysis.R")
```

The script will:
- Load all data files
- Perform quality control
- Generate differential expression results
- Create publication-quality figures
- Export JMP-compatible data
- Generate summary report

**Expected runtime:** 30-60 minutes (depending on system)

---

## 📊 Data Description

### Expression Matrices

**Files:** `data/processed/ACF normalized epi.csv`, `ACF normalized stro.csv`

- **Format:** CSV with genes as rows, samples as columns
- **Normalization:** DESeq2 size factor normalization
- **Values:** Normalized expression counts
- **Dimensions:** 395 genes × 20 samples (per compartment)

**Sample Naming Convention:**
- Format: `PatientID-Replicate-Compartment`
- Example: `110-1E` = Patient 110, Replicate 1, Epithelial
- Compartments: E (Epithelial), S (Stromal)
- Conditions: 1 (ACF), 2 (Normal)

### Differential Expression Results

**Files:** `data/processed/DGE_epithelial.xlsx`, `DGE_stromal.xlsx`

- **Method:** DESeq2 (version 1.26.0)
- **Design:** ~Condition (ACF vs Normal)
- **Test:** Wald test
- **Multiple testing:** Benjamini-Hochberg FDR
- **Significance:** FC >1.5 or <0.67, padj <0.05

**Columns:**
- `Gene`: Gene symbol
- `baseMean`: Mean normalized expression
- `log2FoldChange`: Log2 fold change (ACF vs Normal)
- `FC`: Fold change (linear scale)
- `lfcSE`: Standard error of log2FC
- `stat`: Wald statistic
- `pvalue`: Raw p-value
- `padj`: Adjusted p-value (FDR)

### Gene Panel

**File:** `data/metadata/Gene panel.xlsx`

- **Total genes:** 398
- **Categories:** 41 functional categories
- **Top categories:**
  - Lymphocyte_infiltrate (46 genes)
  - Checkpoint_pathway (30 genes)
  - Tumor_marker (23 genes)
  - Chemokine_signaling (10 genes)

---

## 📈 Key Results

### Epithelial Compartment (ACF vs Normal)

**Significant genes:** 113 (64 upregulated, 49 downregulated)

**Top Upregulated:**
1. ARG1 (FC=2,176, p=1.5e-04)
2. CCR7 (FC=260, p=8.7e-11)
3. CRTAM (FC=98, p=2.7e-11)
4. TDO2 (FC=33, p=8.6e-06)
5. MADCAM1 (FC=26, p=9.4e-05)

**Top Downregulated:**
1. **PDCD1 (PD1)** (FC=0.053, p=6.1e-07) ⭐ **MOST SIGNIFICANT**
2. IL2 (FC=0.17, p=4.9e-02)
3. S100A9 (FC=0.18, p=1.9e-04)
4. FCGR3B (FC=0.18, p=2.6e-02)
5. CXCL8 (FC=0.18, p=3.1e-02)

**Enriched Pathways:**
- Neutrophil chemotaxis
- Chemokine production
- Leukocyte chemotaxis
- NF-κB signaling

### Stromal Compartment (ACF vs Normal)

**Significant genes:** 93 (54 upregulated, 39 downregulated)

**Top Upregulated:**
1. TWIST1 (FC=105, p=4.1e-06)
2. CEACAM8 (FC=104, p=2.3e-05)
3. KLRF1 (FC=73, p=2.1e-09)
4. FCRLA (FC=51, p=2.2e-02)
5. CMKLR1 (FC=10.5, p=1.6e-05)

**Top Downregulated:**
1. CXCL13 (FC=0.054, p=9.7e-05)
2. CDKN2A (FC=0.059, p=1.3e-02)
3. EGR3 (FC=0.073, p=1.1e-04)
4. CXCL8 (FC=0.085, p=1.2e-04)
5. CCR7 (FC=0.147, p=5.9e-03)

**Enriched Pathways:**
- Chemokine-mediated signaling
- Cell chemotaxis
- Immune cell recruitment
- JAK/STAT3 signaling

---

## 📊 Figures

### Main Figures

1. **Figure 1:** Study design and ACF histology
2. **Figure 2:** PCA analysis (all samples, epithelial, stromal)
3. **Figure 3:** Volcano plots and expression heatmaps
4. **Figure 4:** GSEA pathway enrichment
5. **Figure 5:** PD1/CD8 immunofluorescence validation

### Supplementary Figures

- **Figure S1:** Additional PCA views
- **Figure S2:** Gene expression patterns
- **Figure S3:** Pathway networks
- **Figure S4:** Correlation heatmaps
- **Figure S5:** Quality control metrics

All figures are generated by `scripts/ACF_master_analysis.R` and saved to `analysis/outputs/figures/`

---

## 📝 Manuscript Status

**Current Version:** May 11, 2025 draft  
**Completion:** 85%  
**Word Count:** 8,413 words

### Section Status

| Section | Status | Completion |
|---------|--------|------------|
| Abstract | ✅ Complete | 100% |
| Introduction | ✅ Complete | 95% |
| Methods | ⚠️ Needs revision | 90% |
| Results | ⚠️ Needs revision | 85% |
| Discussion | ⚠️ Major revision needed | 70% |
| References | ✅ Complete | 95% |

### Critical Tasks

- [ ] Rewrite Discussion section
- [ ] Clarify Methods (replicate strategy, DESeq2 parameters)
- [ ] Fix terminology issues (PD1 up vs down)
- [ ] Remake Figure 3 (volcano plots, heatmaps)
- [ ] Remake Figure 4 (pathway enrichment)
- [ ] Create supplementary figures
- [ ] Format tables properly

---

## 🔬 Methods Summary

### Sample Processing

1. **Colonoscopy:** High-definition chromoendoscopy with indigo carmine dye
2. **Biopsy collection:** ACF and matched normal from same colon region
3. **Tissue processing:** OCT freezing, storage at -80°C
4. **Laser capture microdissection (LCM):**
   - 9-μm frozen sections on PEN membrane slides
   - Rapid dehydration with RNase inhibitor
   - Separate capture of epithelial and stromal compartments

### RNA Sequencing

1. **RNA extraction:** Arcturus PicoPure kit (1-10ng input)
2. **cDNA synthesis:** SuperScript VILO
3. **Amplification:** Oncomine Immune Response Assay (395 genes)
4. **Sequencing:** Ion S5XL, Ion 540 chip, 2-3M reads/sample
5. **Replicates:** Two technical replicates per sample

### Data Analysis

1. **Alignment & Quantification:** Torrent Suite "ImmuneResponseRNA" plugin
2. **Normalization:** DESeq2 size factors
3. **Differential Expression:** DESeq2 (Wald test, FDR correction)
4. **PCA:** Variance stabilizing transformation, top 50 variable genes
5. **Pathway Enrichment:** GeneAnalytics (GO:BP terms, p<0.0001)

### Statistical Thresholds

- **Fold change:** >1.5 or <0.67
- **Adjusted p-value:** <0.05
- **Detection rate:** >50% in upregulated group

---

## 💾 JMP Data Export

For Dr. Nelson's independent analysis, all data is exported in JMP-compatible formats:

**Location:** `analysis/outputs/jmp_export/`

**Files:**
1. `expression_epithelial.csv` - Normalized expression (epithelial)
2. `expression_stromal.csv` - Normalized expression (stromal)
3. `metadata_epithelial.csv` - Sample information (epithelial)
4. `metadata_stromal.csv` - Sample information (stromal)
5. `DGE_epithelial.csv` - Differential expression results (epithelial)
6. `DGE_stromal.csv` - Differential expression results (stromal)
7. `PCA_coordinates.csv` - PCA sample coordinates
8. `PCA_loadings.csv` - PCA gene loadings
9. `gene_annotations.csv` - Gene panel information
10. `ACF_Complete_Dataset.xlsx` - Comprehensive Excel workbook with all data

---

## 📚 Documentation

### Comprehensive Guides

1. **[COMPREHENSIVE_PROJECT_ANALYSIS.md](docs/COMPREHENSIVE_PROJECT_ANALYSIS.md)**
   - Complete project overview
   - Detailed data analysis
   - Biological interpretation
   - Manuscript assessment
   - Publication recommendations

2. **[DATA_DICTIONARY.md](docs/DATA_DICTIONARY.md)** (to be created)
   - All variables explained
   - File format specifications
   - Naming conventions

3. **[ANALYSIS_GUIDE.md](docs/ANALYSIS_GUIDE.md)** (to be created)
   - Step-by-step instructions
   - Troubleshooting tips
   - Expected outputs

---

## 🤝 Contributing

This is a research project repository. For questions or collaborations:

- **Primary Contact:** Dr. Daniel W. Rosenberg (daniel.rosenberg@uconn.edu)
- **Data Analysis:** Boyang Li (boyang.li@uconn.edu)
- **Repository Maintainer:** Rustem Kurkcu

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- **Funding:** [Add funding sources]
- **Tissue Collection:** Dr. Rosenberg's lab, UConn Health
- **Sequencing:** ThermoFisher Scientific (Geoffrey Lowman, Tim Looney)
- **Data Analysis:** Boyang Li, NinjaTech AI team
- **Clinical Support:** Dr. Thomas J Devers, John Dempsey Hospital

---

## 📊 Citation

**Manuscript in preparation:**

Ideta T*, Li B*, Lemos BS, Lowman G, Looney T, Devers TJ, Giardina C, Rosenberg DW. The Epithelial-Stromal Microenvironment in Early Colonic Neoplasia. *In preparation for Clinical Cancer Research*, 2025.

*These authors contributed equally to this work.

---

## 📞 Contact

For questions about this analysis or collaboration opportunities:

- **Dr. Daniel W. Rosenberg:** daniel.rosenberg@uconn.edu
- **Dr. Charles Giardina:** charles.giardina@uconn.edu
- **Boyang Li:** boyang.li@uconn.edu

**Institution:**  
Center for Molecular Oncology  
University of Connecticut School of Medicine  
Farmington, CT 06030

---

**Last Updated:** October 19, 2025  
**Repository:** https://github.com/RustemKurkcu/ACF-Immune-Analysis  
**Status:** Active development - manuscript preparation phase