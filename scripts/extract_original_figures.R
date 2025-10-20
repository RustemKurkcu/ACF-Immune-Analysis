################################################################################
# EXTRACT ORIGINAL FIGURES FROM BO'S ANALYSIS
################################################################################
#
# This script finds and copies all original figures from Bo's analysis
# to the manuscript figures folder for reference
#
# Author: SuperNinja AI Agent
# Date: October 19, 2025
#
################################################################################

cat("\n")
cat("================================================================================\n")
cat("EXTRACTING ORIGINAL FIGURES FROM BO'S ANALYSIS\n")
cat("================================================================================\n\n")

# Create output directory
dir.create("manuscript/figures/original_from_bo", recursive = TRUE, showWarnings = FALSE)

# Define source directory (adjust path as needed)
source_dir <- "../ACF paper revised"

if (!dir.exists(source_dir)) {
  cat("ERROR: Source directory not found!\n")
  cat("Please adjust the source_dir path in this script.\n")
  cat("Current path:", source_dir, "\n")
  stop("Source directory not found")
}

cat("Source directory:", source_dir, "\n\n")

# Find all image files
cat("Searching for image files...\n")

image_extensions <- c("pdf", "png", "jpg", "jpeg", "tiff", "tif")
all_images <- c()

for (ext in image_extensions) {
  pattern <- paste0("\\.", ext, "$")
  images <- list.files(source_dir, pattern = pattern, 
                      recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  all_images <- c(all_images, images)
}

cat("Found", length(all_images), "image files\n\n")

# Categorize images
cat("Categorizing images...\n")

# PCA figures
pca_images <- all_images[grepl("PCA|pca", all_images, ignore.case = TRUE)]
cat("  PCA figures:", length(pca_images), "\n")

# Heatmap figures
heatmap_images <- all_images[grepl("heatmap|Heatmap", all_images, ignore.case = TRUE)]
cat("  Heatmap figures:", length(heatmap_images), "\n")

# Barchart/GSEA figures
gsea_images <- all_images[grepl("barchart|GSEA|gsea", all_images, ignore.case = TRUE)]
cat("  GSEA/Barchart figures:", length(gsea_images), "\n")

# Expression plot figures
expr_images <- all_images[grepl("expression|Expression", all_images, ignore.case = TRUE)]
cat("  Expression plot figures:", length(expr_images), "\n")

# Figure S (supplementary)
supp_images <- all_images[grepl("Fig S|FigS|Figure S", all_images, ignore.case = TRUE)]
cat("  Supplementary figures:", length(supp_images), "\n")

cat("\n")

# Function to copy files with organized naming
copy_organized <- function(files, category) {
  if (length(files) == 0) return()
  
  cat("Copying", category, "figures...\n")
  
  dest_dir <- file.path("manuscript/figures/original_from_bo", category)
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (i in seq_along(files)) {
    src <- files[i]
    
    # Get filename
    filename <- basename(src)
    
    # Get parent directory name for context
    parent_dir <- basename(dirname(src))
    
    # Create descriptive name
    new_name <- paste0(sprintf("%03d", i), "_", parent_dir, "_", filename)
    
    dest <- file.path(dest_dir, new_name)
    
    # Copy file
    file.copy(src, dest, overwrite = TRUE)
    
    cat("  ", i, "/", length(files), ":", filename, "\n")
  }
  
  cat("  Copied", length(files), category, "figures\n\n")
}

# Copy all categorized images
copy_organized(pca_images, "PCA")
copy_organized(heatmap_images, "Heatmaps")
copy_organized(gsea_images, "GSEA_Pathways")
copy_organized(expr_images, "Expression_Plots")
copy_organized(supp_images, "Supplementary")

# Copy all other images
other_images <- setdiff(all_images, c(pca_images, heatmap_images, gsea_images, 
                                      expr_images, supp_images))
if (length(other_images) > 0) {
  copy_organized(other_images, "Other")
}

# Create index file
cat("Creating index file...\n")

index_content <- c(
  "# ORIGINAL FIGURES FROM BO'S ANALYSIS",
  "",
  paste("Extracted on:", Sys.time()),
  paste("Total figures:", length(all_images)),
  "",
  "## Categories:",
  "",
  paste("- PCA figures:", length(pca_images)),
  paste("- Heatmap figures:", length(heatmap_images)),
  paste("- GSEA/Pathway figures:", length(gsea_images)),
  paste("- Expression plots:", length(expr_images)),
  paste("- Supplementary figures:", length(supp_images)),
  paste("- Other figures:", length(other_images)),
  "",
  "## Directory Structure:",
  "",
  "```",
  "original_from_bo/",
  "├── PCA/",
  "├── Heatmaps/",
  "├── GSEA_Pathways/",
  "├── Expression_Plots/",
  "├── Supplementary/",
  "└── Other/",
  "```",
  "",
  "## Notes:",
  "",
  "- Files are numbered sequentially within each category",
  "- Original parent directory name is preserved in filename",
  "- Use these as reference for recreating figures",
  "- Compare with newly generated figures for consistency",
  "",
  "## Key Figures to Reference:",
  "",
  "### For Figure 2 (PCA):",
  "- Look in PCA/ folder",
  "- Files from 'Figur2 PCA' or 'Figure S1' directories",
  "",
  "### For Figure 3 (Heatmaps):",
  "- Look in Heatmaps/ folder",
  "- Files from 'heatmap' directories",
  "",
  "### For Figure 4 (GSEA):",
  "- Look in GSEA_Pathways/ folder",
  "- Files with 'barchart' in name",
  "",
  "### For Supplementary Figures:",
  "- Look in Supplementary/ folder",
  "- Files with 'Fig S' or 'FigS' in name"
)

writeLines(index_content, "manuscript/figures/original_from_bo/INDEX.md")

cat("  Index file created\n\n")

# Create summary report
summary_df <- data.frame(
  Category = c("PCA", "Heatmaps", "GSEA/Pathways", "Expression Plots", 
               "Supplementary", "Other", "TOTAL"),
  Count = c(length(pca_images), length(heatmap_images), length(gsea_images),
            length(expr_images), length(supp_images), length(other_images),
            length(all_images))
)

write.csv(summary_df, "manuscript/figures/original_from_bo/summary.csv", 
          row.names = FALSE)

cat("================================================================================\n")
cat("EXTRACTION COMPLETE!\n")
cat("================================================================================\n\n")

cat("Summary:\n")
print(summary_df)

cat("\nAll original figures copied to:\n")
cat("  manuscript/figures/original_from_bo/\n\n")

cat("Next steps:\n")
cat("  1. Review original figures in each category\n")
cat("  2. Compare with newly generated figures\n")
cat("  3. Use as reference for figure recreation\n")
cat("  4. Ensure consistency in style and content\n\n")

cat("================================================================================\n")