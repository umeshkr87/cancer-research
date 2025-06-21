#!/usr/bin/env Rscript
library(TCGAmutations)
cat("📥 Downloading PRAD mutations...\n")
prad_data <- tcga_load(study = "PRAD", source = "MC3")
mutations <- prad_data@data
output_dir <- "/u/aa107/uiuc-cancer-research/data/processed/tcga_prad_prostate"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(mutations, file.path(output_dir, "tcga_prad_mutations.csv"), row.names = FALSE)
cat("✅ Saved", nrow(mutations), "mutations to", file.path(output_dir, "tcga_prad_mutations.csv"), "\n")
