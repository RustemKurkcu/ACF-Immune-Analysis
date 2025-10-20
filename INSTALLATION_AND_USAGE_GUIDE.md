# ACF IMMUNE ANALYSIS - COMPLETE INSTALLATION & USAGE GUIDE

**For Windows with WSL2 Ubuntu**

---

## 📋 TABLE OF CONTENTS

1. [System Requirements](#system-requirements)
2. [R Installation & Setup](#r-installation--setup)
3. [Package Installation](#package-installation)
4. [Running the Analysis](#running-the-analysis)
5. [Generating Figures](#generating-figures)
6. [Generating Tables](#generating-tables)
7. [Troubleshooting](#troubleshooting)
8. [Expected Outputs](#expected-outputs)

---

## 1. SYSTEM REQUIREMENTS

### Your System (✅ All Requirements Met!)
- **OS:** Windows with WSL2 Ubuntu
- **CPU:** i9 processor ✅
- **RAM:** 128 GB ✅ (only need 8-16 GB)
- **Storage:** ~20 TB ✅ (only need ~10 GB for project)

### Software Requirements
- **R:** Version 4.0 or higher (recommended: 4.3+)
- **RStudio:** Latest version (optional but recommended)
- **Git:** For repository management

---

## 2. R INSTALLATION & SETUP

### Option A: Install R on Windows (Recommended)

1. **Download R:**
   - Go to: https://cran.r-project.org/bin/windows/base/
   - Download the latest version (e.g., R-4.3.2-win.exe)
   - Run the installer with default settings

2. **Download RStudio:**
   - Go to: https://posit.co/download/rstudio-desktop/
   - Download RStudio Desktop (free version)
   - Install with default settings

3. **Verify Installation:**
   - Open RStudio
   - In the Console, type: `R.version`
   - Should show R version 4.x.x

### Option B: Install R on WSL2 Ubuntu

```bash
# Update package list
sudo apt update

# Install R
sudo apt install -y r-base r-base-dev

# Install system dependencies for R packages
sudo apt install -y \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libfontconfig1-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype6-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev

# Verify installation
R --version
```

---

## 3. PACKAGE INSTALLATION

### Step 1: Open R or RStudio

**On Windows:** Open RStudio  
**On WSL2:** Type `R` in terminal

### Step 2: Install Required Packages

Copy and paste this entire code block into R console:

```r
# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager first
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# CRAN packages
cran_packages <- c(
  # Data manipulation
  "tidyverse", "readxl", "writexl", "data.table",
  
  # Visualization
  "ggplot2", "pheatmap", "ggrepel", "cowplot",
  "RColorBrewer", "viridis", "ggsci",
  
  # Utilities
  "scales", "reshape2", "gridExtra", "knitr", "kableExtra"
)

# Install CRAN packages
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg)
  }
}

# Bioconductor packages
bioc_packages <- c(
  "DESeq2", "edgeR", "limma",
  "clusterProfiler", "enrichplot", "DOSE",
  "org.Hs.eg.db", "ComplexHeatmap"
)

# Install Bioconductor packages
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

cat("\n✅ All packages installed successfully!\n")
```

**Expected time:** 10-30 minutes (depending on internet speed)

### Step 3: Verify Installation

```r
# Test that all packages load
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)

cat("✅ All packages loaded successfully!\n")
```

---

## 4. RUNNING THE ANALYSIS

### Step 1: Clone the Repository

**On Windows (Git Bash or Command Prompt):**
```bash
cd C:\Users\YourUsername\Documents
git clone https://github.com/RustemKurkcu/ACF-Immune-Analysis.git
cd ACF-Immune-Analysis
```

**On WSL2 Ubuntu:**
```bash
cd ~
git clone https://github.com/RustemKurkcu/ACF-Immune-Analysis.git
cd ACF-Immune-Analysis
```

### Step 2: Open RStudio and Set Working Directory

**In RStudio:**
1. File → Open Project
2. Navigate to `ACF-Immune-Analysis` folder
3. Open `ACF-Immune-Analysis.Rproj` (if exists) or just set working directory

**Or in R Console:**
```r
# Windows
setwd("C:/Users/YourUsername/Documents/ACF-Immune-Analysis")

# WSL2
setwd("~/ACF-Immune-Analysis")

# Verify
getwd()  # Should show the ACF-Immune-Analysis path
```

### Step 3: Run the Master Analysis Script

```r
# Run the complete analysis
source("scripts/ACF_master_analysis.R")
```

**What this does:**
- Loads all data files
- Performs quality control
- Analyzes differential expression
- Creates PCA plots
- Performs pathway enrichment
- Generates figures
- Exports JMP data
- Creates summary report

**Expected runtime:** 30-60 minutes

**Progress:** The script will print progress messages as it runs:
```
================================================================================
ACF IMMUNE ANALYSIS - MASTER PIPELINE
================================================================================

Loading required libraries...
All packages loaded successfully!

================================================================================
SECTION 2: DATA LOADING
================================================================================

Loading normalized expression matrices...
  Epithelial: 395 genes × 20 samples
  Stromal: 395 genes × 20 samples
...
```

---

## 5. GENERATING FIGURES

### Option 1: Run Figure Generation Script

```r
# Generate all figures
source("scripts/generate_all_figures.R")
```

**Generates:**
- Figure 1: Study Design Schematic
- Figure 2: PCA Analysis (A-C panels)
- Figure 3: Volcano Plots (A-B) & Heatmaps (C-D)
- Figure 4: Pathway Enrichment (requires GSEA data)
- Supplementary Figure S1: QC Metrics

**Output location:** `manuscript/figures/`

### Option 2: Generate Individual Figures

```r
# Load required libraries
library(tidyverse)
library(ggplot2)
library(pheatmap)

# Load data
expr_epi <- read.csv("data/processed/ACF normalized epi.csv", row.names = 1)
dge_epi <- readxl::read_excel("data/processed/DGE_epithelial.xlsx")

# Create a volcano plot
volcano_data <- dge_epi %>%
  mutate(
    neg_log10_padj = -log10(padj),
    significant = padj < 0.05 & (FC > 1.5 | FC < 0.67)
  )

ggplot(volcano_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("gray70", "red")) +
  theme_bw() +
  labs(title = "Epithelial: ACF vs Normal",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted P-value)")

# Save
ggsave("my_volcano_plot.pdf", width = 10, height = 8)
```

---

## 6. GENERATING TABLES

### Run Table Generation Script

```r
# Generate all tables
source("scripts/generate_all_tables.R")
```

**Generates:**
- Table 1: Patient Demographics
- Table 2: Chemokine-Related Genes
- Table S1: Complete DGE Results (Epithelial)
- Table S2: Complete DGE Results (Stromal)
- Table S3: Gene Panel Categories
- Table S4: Sample Sequencing Metrics
- Table S5: Top Differentially Expressed Genes
- Master Workbook: All tables in one Excel file

**Output location:** `manuscript/tables/`

---

## 7. TROUBLESHOOTING

### Problem: Package Installation Fails

**Solution 1: Install system dependencies (WSL2 only)**
```bash
sudo apt install -y \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libfontconfig1-dev
```

**Solution 2: Install packages one at a time**
```r
# Try installing problematic package individually
install.packages("tidyverse")

# If that fails, try from source
install.packages("tidyverse", type = "source")
```

**Solution 3: Use different CRAN mirror**
```r
# Try a different mirror
options(repos = "https://cran.rstudio.com/")
install.packages("tidyverse")
```

### Problem: "Cannot find file" Error

**Solution:**
```r
# Check working directory
getwd()

# List files in current directory
list.files()

# Set correct working directory
setwd("path/to/ACF-Immune-Analysis")

# Verify data files exist
file.exists("data/processed/ACF normalized epi.csv")
```

### Problem: Out of Memory Error

**Solution:**
```r
# Clear workspace
rm(list = ls())
gc()

# Increase memory limit (Windows only)
memory.limit(size = 16000)  # 16 GB
```

### Problem: Plots Don't Display

**Solution:**
```r
# For RStudio: Check Plots pane
# For R console: Save directly to file

# Save plot to file
pdf("my_plot.pdf", width = 10, height = 8)
plot(1:10)
dev.off()
```

### Problem: Script Runs Slowly

**This is normal!** The complete analysis takes 30-60 minutes.

**To speed up:**
```r
# Run only specific sections
# Comment out sections you don't need with #

# For example, skip pathway enrichment:
# source("scripts/ACF_master_analysis.R")
# Then manually comment out Section 9 in the script
```

---

## 8. EXPECTED OUTPUTS

### After Running Master Analysis Script

**Directory structure:**
```
ACF-Immune-Analysis/
├── analysis/outputs/
│   ├── figures/
│   │   ├── QC_Epithelial_expression_boxplot.pdf
│   │   ├── QC_Stromal_expression_boxplot.pdf
│   │   ├── Volcano_Epithelial.pdf
│   │   ├── Volcano_Stromal.pdf
│   │   ├── Heatmap_Epithelial_top50.pdf
│   │   ├── Heatmap_Stromal_top50.pdf
│   │   ├── PCA_all_samples.pdf
│   │   └── GO_BP_*.pdf
│   │
│   ├── tables/
│   │   ├── DGE_results_processed.xlsx
│   │   ├── QC_Epithelial_sample_summary.csv
│   │   ├── QC_Stromal_sample_summary.csv
│   │   └── GO_BP_*.csv
│   │
│   ├── data/
│   │   ├── metadata_epithelial.csv
│   │   ├── metadata_stromal.csv
│   │   └── analysis_summary.rds
│   │
│   └── jmp_export/
│       ├── expression_epithelial.csv
│       ├── expression_stromal.csv
│       ├── DGE_epithelial.csv
│       ├── DGE_stromal.csv
│       ├── PCA_coordinates.csv
│       └── ACF_Complete_Dataset.xlsx
```

### After Running Figure Generation Script

**Additional outputs:**
```
manuscript/figures/
├── publication/
│   ├── Figure1_Study_Design.pdf
│   ├── Figure2_PCA_Analysis.pdf
│   ├── Figure3AB_Volcano_Plots.pdf
│   ├── Figure3C_Epithelial_Heatmap.pdf
│   └── Figure3D_Stromal_Heatmap.pdf
│
└── supplementary/
    └── FigureS1_QC_Metrics.pdf
```

### After Running Table Generation Script

**Additional outputs:**
```
manuscript/tables/
├── Table1_Demographics.xlsx
├── Table2_Chemokine_Genes.xlsx
├── TableS1_DGE_Epithelial_Complete.xlsx
├── TableS2_DGE_Stromal_Complete.xlsx
├── TableS3_Gene_Panel_Categories.xlsx
├── TableS4_Sample_Metrics.xlsx
├── TableS5_Top_DE_Genes.xlsx
└── ACF_All_Tables_Master.xlsx
```

---

## 9. QUICK START CHECKLIST

- [ ] Install R (version 4.0+)
- [ ] Install RStudio (optional but recommended)
- [ ] Clone GitHub repository
- [ ] Open RStudio and set working directory
- [ ] Install required packages (run package installation code)
- [ ] Run master analysis script: `source("scripts/ACF_master_analysis.R")`
- [ ] Run figure generation: `source("scripts/generate_all_figures.R")`
- [ ] Run table generation: `source("scripts/generate_all_tables.R")`
- [ ] Check outputs in `analysis/outputs/` and `manuscript/`
- [ ] Share JMP data from `analysis/outputs/jmp_export/` with Dr. Nelson

---

## 10. GETTING HELP

### If You Encounter Issues:

1. **Check the error message carefully**
   - Most errors are self-explanatory
   - Google the error message for solutions

2. **Verify your setup**
   - R version: `R.version`
   - Package versions: `packageVersion("DESeq2")`
   - Working directory: `getwd()`

3. **Check file paths**
   - Use forward slashes `/` not backslashes `\`
   - Use relative paths, not absolute paths

4. **Restart R session**
   - In RStudio: Session → Restart R
   - Clears memory and resets environment

5. **Contact for help**
   - Email: boyang.li@uconn.edu
   - Include: Error message, R version, what you were trying to do

---

## 11. NEXT STEPS AFTER RUNNING ANALYSIS

1. **Review all outputs**
   - Check figures in `manuscript/figures/`
   - Review tables in `manuscript/tables/`
   - Examine QC metrics

2. **Share JMP data with Dr. Nelson**
   - Send entire `analysis/outputs/jmp_export/` folder
   - Or just send `ACF_Complete_Dataset.xlsx`

3. **Begin manuscript revisions**
   - Use generated figures
   - Update tables
   - Revise Discussion section

4. **Prepare for publication**
   - Format figures for journal
   - Create figure legends
   - Finalize supplementary materials

---

## 12. ESTIMATED TIME REQUIREMENTS

| Task | Time Required |
|------|---------------|
| R Installation | 10-15 minutes |
| Package Installation | 10-30 minutes |
| Clone Repository | 1-2 minutes |
| Run Master Analysis | 30-60 minutes |
| Generate Figures | 5-10 minutes |
| Generate Tables | 2-5 minutes |
| Review Outputs | 30-60 minutes |
| **TOTAL** | **~2-3 hours** |

---

## 13. TIPS FOR SUCCESS

1. **Be patient** - Package installation and analysis take time
2. **Read error messages** - They usually tell you what's wrong
3. **Save your work** - RStudio auto-saves, but save scripts manually
4. **Use RStudio** - Much easier than command-line R
5. **Check outputs** - Verify figures and tables look correct
6. **Ask for help** - Don't struggle alone if stuck

---

**Good luck! You've got this! 🚀**

---

**Last Updated:** October 19, 2025  
**Maintained by:** SuperNinja AI Agent  
**For:** Rustem Kurkcu, UConn Health