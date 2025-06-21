#!/usr/bin/env Rscript
input_file <- "/u/aa107/uiuc-cancer-research/data/processed/tcga_prad_prostate/tcga_prad_mutations.csv"
output_file <- "/u/aa107/uiuc-cancer-research/data/processed/tcga_prad_prostate/tcga_prad_summary.txt"

data <- read.csv(input_file)
cat("🔍 Generating TCGA-PRAD summary report...\n")

# Generate summary report
sink(output_file)
cat("TCGA-PRAD Mutation Analysis Summary\n")
cat("===================================\n\n")
cat("Total mutations:", nrow(data), "\n")
cat("Unique genes:", length(unique(data$Hugo_Symbol)), "\n")
cat("Unique samples:", length(unique(data$Tumor_Sample_Barcode)), "\n")
cat("Chromosomes:", paste(sort(unique(data$Chromosome)), collapse = ", "), "\n\n")

cat("Mutation Type Distribution:\n")
variant_counts <- sort(table(data$Variant_Classification), decreasing = TRUE)
for (i in 1:length(variant_counts)) { cat("  ", names(variant_counts)[i], ":", variant_counts[i], "\n") }

cat("\nTop 20 Genes:\n")
gene_counts <- sort(table(data$Hugo_Symbol), decreasing = TRUE)
for (i in 1:min(20, length(gene_counts))) { cat("  ", names(gene_counts)[i], ":", gene_counts[i], "\n") }

cat("\nChromosome Distribution:\n")
chr_counts <- sort(table(data$Chromosome))
for (i in 1:length(chr_counts)) { cat("  Chr", names(chr_counts)[i], ":", chr_counts[i], "\n") }

prostate_genes <- c("AR", "PTEN", "TP53", "BRCA1", "BRCA2", "ATM", "PIK3CA", "MTOR", "TSC1", "TSC2")
prostate_muts <- sum(data$Hugo_Symbol %in% prostate_genes)
cat("\nKey Prostate Cancer Genes:\n")
cat("  Total mutations in key genes:", prostate_muts, sprintf("(%.1f%%)\n", 100*prostate_muts/nrow(data)))
sink()

cat("✅ Summary saved to:", output_file, "\n")
