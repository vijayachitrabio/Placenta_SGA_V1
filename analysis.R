
# ============================================================================
#
#   PLACENTA DNA METHYLATION ANALYSIS — SGA vs CONTROL
#   Teaching Script for Students
#
#   Project : S-2025-138
#   Array   : Illumina EPIC 850K (~850,000 CpG sites)
#   Tissue  : Placenta
#   Groups  : SGA (Small for Gestational Age) vs AGA (Control)
#   Samples : n = 72 (35 SGA + 37 Control)
#
#   This script walks through a complete DNA methylation analysis:
#     1.  Load and check data
#     2.  Quality control
#     3.  Normalisation
#     4.  Filter unreliable probes
#     5.  Exploratory plots (PCA, heatmap)
#     6.  Find differentially methylated positions (DMPs)
#     7.  Find differentially methylated regions (DMRs)
#     8.  Gene ontology enrichment
#     9.  Epigenetic clock
#     10. Save all results
#
#   HOW TO USE THIS SCRIPT
#   ──────────────────────
#   - Read every comment before running each section
#   - Run one section at a time (Ctrl+Enter in RStudio)
#   - Check the output in the Console before moving on
#   - If something fails, read the error carefully — it usually tells you
#     exactly what went wrong
#
#   WHAT YOU NEED BEFORE STARTING
#   ──────────────────────────────
#   - R (version >= 4.2)  : https://cran.r-project.org
#   - RStudio             : https://posit.co/download/rstudio-desktop/
#   - The iDAT files folder (provided by supervisor)
#   - The sample sheet Excel file (provided by supervisor)
#
# ============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 0 : INSTALL PACKAGES
# ─────────────────────────────────────────────────────────────────────────────
#
# You only need to run this section ONCE on a new computer.
# After packages are installed, you can skip straight to Section 1.
#
# This will take 15–30 minutes on first run — that is normal.
# If asked "Update all/some/none?" type: n (for none) and press Enter.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

# Install packages from CRAN (the standard R repository)
cran_needed <- c("ggplot2", "pheatmap", "RColorBrewer",
                 "openxlsx", "dplyr", "tidyr", "ggrepel")
install.packages(
  setdiff(cran_needed, rownames(installed.packages())),
  repos = "https://cloud.r-project.org"
)

# Install packages from Bioconductor (specialised biology packages)
bioc_needed <- c(
  "minfi",                                          # core methylation analysis
  "limma",                                          # linear models (differential analysis)
  "DMRcate",                                        # find differentially methylated regions
  "missMethyl",                                     # gene ontology for methylation
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",  # EPIC array probe annotation
  "IlluminaHumanMethylationEPICmanifest",           # EPIC array manifest
  "ENmix",                                          # quality control tools
  "sva",                                            # surrogate variable analysis
  "EpiDISH"                                         # cell-type deconvolution
)
BiocManager::install(
  setdiff(bioc_needed, rownames(installed.packages())),
  update = FALSE,
  ask    = FALSE
)

# methylclock is only available from GitHub (not Bioconductor)
if (!requireNamespace("methylclock", quietly = TRUE))
  remotes::install_github("isglobal-brge/methylclock", quiet = TRUE)

cat("All packages installed successfully!\n")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1 : LOAD PACKAGES AND SET PATHS
# ─────────────────────────────────────────────────────────────────────────────
#
# Every time you start a new R session you need to load the packages.
# Installation (Section 0) only needs to happen once.

library(minfi)
library(limma)
library(DMRcate)
library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(ENmix)
library(sva)
library(EpiDISH)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)
library(dplyr)

# ── Tell R where your files are ──────────────────────────────────────────────
#
# IMPORTANT: Change these three paths to match your own computer.
# Use forward slashes / even on Windows.
# Example Windows path: "C:/Users/YourName/Documents/SGA/idats"

IDAT_DIR     <- path.expand(
  "~/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/SGA/S-2025-138/idats"
)
SAMPLE_SHEET <- path.expand(
  "~/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/SGA/S-2025-138/Placenta_SampleSheet_minfi_n72.xlsx"
)
OUT_DIR <- path.expand(
  "~/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/SGA/S-2025-138/results"
)

# Create the results folder if it does not exist yet
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Check that the paths exist — this will stop with a clear error if not
if (!dir.exists(IDAT_DIR))
  stop("Cannot find iDAT folder: ", IDAT_DIR,
       "\nPlease check the path above.")
if (!file.exists(SAMPLE_SHEET))
  stop("Cannot find sample sheet: ", SAMPLE_SHEET,
       "\nPlease check the path above.")

cat("Paths OK\n")
cat("Results will be saved to:", OUT_DIR, "\n")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2 : LOAD THE SAMPLE SHEET
# ─────────────────────────────────────────────────────────────────────────────
#
# The sample sheet is an Excel file that tells R which iDAT files belong to
# which sample, and what group each sample is in (SGA or Control).
#
# Key columns:
#   Sample_Name      : unique identifier e.g. KR001-P1
#   Sample_Group     : SGA or Control
#   Foetal_Sex       : F or M  (derived from sex chromosome check — NOT maternal sex)
#   Gestational_Age  : weeks+days format e.g. "39+2"
#   Birth_Weight_g   : birth weight in grams
#   Sentrix_ID       : 10-digit chip barcode
#   Sentrix_Position : position on chip e.g. R04C01
#   Basename         : full path to iDAT files (built below)

targets <- read.xlsx(SAMPLE_SHEET, sheet = "Placenta_Analysis_n72")

# Convert columns to the correct data type
targets$Sample_Group     <- factor(targets$Sample_Group, levels = c("Control", "SGA"))
targets$Foetal_Sex       <- factor(targets$Foetal_Sex,   levels = c("F", "M"))
targets$Sentrix_ID       <- as.character(targets$Sentrix_ID)
targets$Sentrix_Position <- as.character(targets$Sentrix_Position)

# Convert gestational age from "39+2" string to decimal weeks (39.29)
targets$GA_weeks <- sapply(strsplit(targets$Gestational_Age, "\\+"), function(x) {
  as.numeric(x[1]) + as.numeric(x[2]) / 7
})

# Build the full file path to each iDAT pair
# minfi looks for files named: <Sentrix_ID>_<Position>_Red.idat
#                          and: <Sentrix_ID>_<Position>_Grn.idat
targets$Basename <- file.path(
  IDAT_DIR,
  targets$Sentrix_ID,
  paste0(targets$Sentrix_ID, "_", targets$Sentrix_Position)
)

# Verify every iDAT file actually exists before trying to load
missing <- !file.exists(paste0(targets$Basename, "_Red.idat")) |
           !file.exists(paste0(targets$Basename, "_Grn.idat"))

if (any(missing)) {
  cat("Missing iDAT files for these samples:\n")
  print(targets$Sample_Name[missing])
  stop("Fix missing files before continuing.")
}

# Print a summary so you can confirm everything looks right
cat("\n=== Sample Sheet Summary ===\n")
cat("Total samples:", nrow(targets), "\n")
cat("\nGroup breakdown:\n")
print(table(targets$Sample_Group))
cat("\nFoetal sex by group:\n")
print(table(targets$Foetal_Sex, targets$Sample_Group))
cat("\nBirth weight (grams) — mean by group:\n")
print(tapply(targets$Birth_Weight_g, targets$Sample_Group, mean, na.rm = TRUE))
cat("\nAll iDAT files found ✓\n")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3 : READ iDAT FILES INTO R
# ─────────────────────────────────────────────────────────────────────────────
#
# iDAT files contain the raw fluorescence intensities from the array scanner.
# read.metharray.exp() reads all iDAT files and creates an RGChannelSet object
# (RG = Red/Green, the two fluorescent dyes used on the array).
#
# This step will take a few minutes — there are 144 iDAT files to read.

cat("Reading iDAT files — this may take a few minutes...\n")



# Attach our sample sheet as phenotype data so we can access it later
pData(RGSet) <- DataFrame(targets)
sampleNames(RGSet) <- targets$Sample_Name

cat("\niDAT files loaded successfully!\n")
cat("Array type detected:", annotation(RGSet)["array"], "\n")
cat("Number of samples:", ncol(RGSet), "\n")
cat("Number of probes:", nrow(RGSet), "\n")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4 : QUALITY CONTROL
# ─────────────────────────────────────────────────────────────────────────────
#
# Quality control (QC) checks whether each sample and each probe behaved
# reliably on the array.
#
# Detection p-value: for each probe in each sample, minfi tests whether the
# signal is significantly above background noise. A p-value > 0.01 means
# the probe did not produce a reliable signal for that sample.
#
# A good sample should have <1% of probes failing this test.

# Calculate detection p-values for all probes in all samples
detP <- detectionP(RGSet)

cat("=== Detection P-value QC ===\n")
cat("Probes per sample failing detP > 0.01:\n")
failed_pct <- colMeans(detP > 0.01) * 100
print(round(sort(failed_pct, decreasing = TRUE), 3))

# Flag samples where more than 1% of probes failed
bad_samples <- names(failed_pct[failed_pct > 1])
if (length(bad_samples) > 0) {
  cat("\nWARNING — these samples have >1% failed probes and will be removed:\n")
  print(bad_samples)
} else {
  cat("\nAll samples passed QC ✓\n")
}

# Generate a full PDF QC report (opens in your results folder)
qcReport(
  rgSet      = RGSet,
  sampNames  = targets$Sample_Name,
  sampGroups = targets$Sample_Group,
  pdf        = file.path(OUT_DIR, "QC_report.pdf")
)
cat("QC report saved to: QC_report.pdf\n")

# Plot methylation vs unmethylation signal (should cluster together)
mSet_raw <- preprocessRaw(RGSet)
qc_metrics <- getQC(mSet_raw)

pdf(file.path(OUT_DIR, "QC_signal_plot.pdf"))
plotQC(qc_metrics)
title("Sample QC — Methylated vs Unmethylated Signal")
dev.off()
cat("QC signal plot saved\n")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5 : NORMALISATION
# ─────────────────────────────────────────────────────────────────────────────
#
# Normalisation corrects for technical variation between arrays so that
# biological differences are what drive the results, not batch effects.
#
# We use preprocessFunnorm() — functional normalisation.
# It is the best choice when:
#   - Comparing different tissue types or groups with large methylation differences
#   - You have many samples across multiple chips (we have 27 chips)
#
# After normalisation we have a GenomicRatioSet (GRSet) containing
# beta values (0-1 scale, easy to interpret) and
# M values (log-scale, better for statistics).

cat("Normalising data with preprocessFunnorm...\n")

GRSet <- preprocessFunnorm(
   RGSet,
  nPCs     = 2,       # number of principal components of control probes to use
  bgCorr   = TRUE,    # correct for background fluorescence
  dyeCorr  = TRUE,    # correct for dye intensity differences between arrays
  verbose  = TRUE
)

# Remove any bad samples identified in QC
if (length(bad_samples) > 0) {
  GRSet  <- GRSet[, !sampleNames(GRSet) %in% bad_samples]
  detP   <- detP[, !colnames(detP) %in% bad_samples]
  targets <- targets[!targets$Sample_Name %in% bad_samples, ]
  cat("Removed", length(bad_samples), "failed sample(s)\n")
}

cat("Normalisation complete\n")
cat("Samples remaining:", ncol(GRSet), "\n")

# Plot beta value density — should look smooth after normalisation
pdf(file.path(OUT_DIR, "Beta_density_normalised.pdf"))
densityPlot(
  getBeta(GRSet),
  sampGroups = pData(GRSet)$Sample_Group,
  main       = "Beta Value Density After Normalisation",
  legend     = FALSE,
  pal        = c("#27AE60", "#E67E22")
)
legend("topright", legend = c("Control", "SGA"),
       col = c("#27AE60", "#E67E22"), lty = 1)
dev.off()
cat("Beta density plot saved\n")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6 : PROBE FILTERING
# ─────────────────────────────────────────────────────────────────────────────
#
# Not all ~850,000 probes should be included in the analysis.
# We remove probes that are:
#
#   1. Failed (detP > 0.01 in any sample) — unreliable signal
#   2. Cross-reactive — bind to multiple genomic locations, so we cannot
#      be sure which site is being measured
#   3. SNP-overlapping — a genetic variant at the CpG site means methylation
#      and genotype are confounded
#   4. On sex chromosomes — we have mixed foetal sexes, so X/Y probes would
#      create a sex effect that is not related to SGA biology

# Start with the annotation (information about each probe's genomic location)
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# --- Filter 1: Failed probes ---
keep_detP <- rowSums(detP[rownames(detP) %in% featureNames(GRSet), ] > 0.01) == 0
cat("Probes removed (failed detection p-value):", sum(!keep_detP), "\n")
GRSet <- GRSet[keep_detP, ]

# --- Filter 2: Cross-reactive probes ---
# Download this file from: https://github.com/sirselim/illumina_array_cross-reactive_probes
# Save it to your project folder and update the path below
cross_react_file <- file.path(dirname(SAMPLE_SHEET), "cross_reactive_EPIC.csv")
if (file.exists(cross_react_file)) {
  cross_react <- read.csv(cross_react_file, header = FALSE)$V1
  keep_cr     <- !featureNames(GRSet) %in% cross_react
  cat("Probes removed (cross-reactive):", sum(!keep_cr), "\n")
  GRSet <- GRSet[keep_cr, ]
} else {
  cat("Cross-reactive probe list not found — skipping this filter.\n")
  cat("Download from: https://github.com/sirselim/illumina_array_cross-reactive_probes\n")
}

# --- Filter 3: SNP-overlapping probes ---
# Remove probes where a common SNP (minor allele freq > 5%) overlaps
# the CpG site or the single-base extension site
GRSet <- dropLociWithSnps(GRSet, snps = c("SBE", "CpG"), maf = 0.05)
cat("Probes after SNP filter:", nrow(GRSet), "\n")

# --- Filter 4: Sex chromosome probes ---
keep_auto <- !featureNames(GRSet) %in% ann$Name[ann$chr %in% c("chrX", "chrY")]
cat("Probes removed (sex chromosomes):", sum(!keep_auto), "\n")
GRSet <- GRSet[keep_auto, ]

cat("\nProbes remaining after all filters:", nrow(GRSet), "\n")
cat("(Started with ~850,000; expect ~700,000–750,000 after filtering)\n")

# Extract beta values (0–1) and M values (log-scale) for downstream analysis
beta  <- getBeta(GRSet)    # use for: visualisation, heatmaps, clocks
M_val <- getM(GRSet)       # use for: statistical testing (limma)

cat("\nBeta matrix dimensions:", nrow(beta), "probes x", ncol(beta), "samples\n")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7 : EXPLORATORY ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────
#
# Before running any statistical tests, we look at the data visually.
# This helps us:
#   - Spot outlier samples that should be investigated
#   - Check whether SGA and Control samples separate as expected
#   - Identify any batch effects (chip, plate, sex)
#
# PCA = Principal Component Analysis
# It reduces ~750,000 dimensions down to 2 so we can plot the samples.
# Samples that are methylation-similar will cluster together.

# --- PCA ---
pca     <- prcomp(t(beta), center = TRUE, scale. = FALSE)
pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_df <- data.frame(
  PC1    = pca$x[, 1],
  PC2    = pca$x[, 2],
  Group  = pData(GRSet)$Sample_Group,
  Sex    = pData(GRSet)$Foetal_Sex,
  Sample = sampleNames(GRSet)
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                             colour = Group, shape = Sex, label = Sample)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_colour_manual(values = c("Control" = "#27AE60", "SGA" = "#E67E22")) +
  labs(
    title    = "PCA — Placenta EPIC Methylation (n=72)",
    subtitle = "Each point = one sample. Colour = group, Shape = foetal sex",
    x        = paste0("PC1 (", pct_var[1], "% variance)"),
    y        = paste0("PC2 (", pct_var[2], "% variance)")
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

ggsave(file.path(OUT_DIR, "Fig1_PCA.pdf"), p_pca, width = 7, height = 5)
ggsave(file.path(OUT_DIR, "Fig1_PCA.png"), p_pca, width = 7, height = 5, dpi = 300)
cat("PCA plot saved\n")

# --- MDS plot (alternative to PCA, built into minfi) ---
pdf(file.path(OUT_DIR, "MDS_plot.pdf"))
mdsPlot(
  dat        = beta,
  numPositions = 1000,
  sampGroups = pData(GRSet)$Sample_Group,
  sampNames  = sampleNames(GRSet),
  pal        = c("#27AE60", "#E67E22"),
  main       = "MDS — Top 1000 Variable Probes"
)
dev.off()
cat("MDS plot saved\n")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8 : CELL-TYPE DECONVOLUTION
# ─────────────────────────────────────────────────────────────────────────────
#
# The placenta contains multiple cell types: trophoblasts, stromal cells,
# Hofbauer cells (macrophages), and endothelial cells.
#
# If SGA and control placentas have different cell-type compositions,
# the methylation differences we find could reflect cell-type differences
# rather than true epigenetic changes.
#
# EpiDISH estimates the proportion of each cell type in each sample.
# We then include these estimates as covariates in our statistical model
# to "adjust" for cell-type composition.

if (requireNamespace("EpiDISH", quietly = TRUE)) {
  library(EpiDISH)

  # Note: Use a placenta-specific reference when available (PlaNET, Yuan 2021)
  # For now we use the built-in blood reference as a demonstration
  # Replace centDHSbloodDMC.m with a placenta reference matrix if you have one
  data(centDHSbloodDMC.m)

  epi_out    <- epidish(beta.m = beta, ref.m = centDHSbloodDMC.m, method = "RPC")
  cell_props <- as.data.frame(epi_out$estF)

  cat("Cell-type proportions estimated for", nrow(cell_props), "samples\n")
  cat("Cell types:", colnames(cell_props), "\n")
  cat("Mean proportions:\n")
  print(round(colMeans(cell_props), 3))

  write.csv(cell_props, file.path(OUT_DIR, "CellType_proportions.csv"))
  cat("Cell-type proportions saved\n")

} else {
  cat("EpiDISH not available — proceeding without cell-type deconvolution\n")
  cell_props <- NULL
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9 : DIFFERENTIAL METHYLATION — DMPs
# ─────────────────────────────────────────────────────────────────────────────
#
# We use limma (Linear Models for Microarray/Methylation data) to find
# CpG sites (probes) that are significantly differently methylated between
# SGA and Control placentas.
#
# We work with M values (not beta values) for statistics because M values
# are approximately normally distributed, which is an assumption of limma.
# Beta values are better for interpretation and visualisation.
#
# The design matrix tells limma what to model:
#   ~ Sample_Group   : the effect we care about (SGA vs Control)
#     + Foetal_Sex   : adjust for sex differences
#     + GA_weeks     : adjust for gestational age differences
#
# After fitting, we apply the Benjamini-Hochberg (BH) correction to control
# the False Discovery Rate (FDR) — this reduces false positives when testing
# hundreds of thousands of probes simultaneously.

pd <- pData(GRSet)  # phenotype data shortcut

# Build the design matrix
if (!is.null(cell_props)) {
  # Include cell-type proportions as covariates (drop one column to avoid
  # perfect collinearity — the dropped column is implicit in the intercept)
  cell_cov <- cell_props[rownames(cell_props) %in% pd$Sample_Name,
                         -ncol(cell_props), drop = FALSE]
  design <- model.matrix(
    ~ Sample_Group + Foetal_Sex + GA_weeks,
    data = cbind(pd, cell_cov)
  )
} else {
  design <- model.matrix(~ Sample_Group + Foetal_Sex + GA_weeks, data = pd)
}

cat("Design matrix (first 3 rows):\n")
print(head(design, 3))
cat("Columns:", colnames(design), "\n")
cat("The coefficient we will test: 'Sample_GroupSGA'\n\n")

# Fit the linear model to every probe
fit  <- lmFit(M_val, design)
fit2 <- eBayes(fit)     # empirical Bayes moderation — improves power

# Extract results for the SGA vs Control comparison
DMPs <- topTable(
  fit2,
  coef      = "Sample_GroupSGA",  # the comparison of interest
  number    = Inf,                 # return all probes, not just top 10
  adjust    = "BH",               # Benjamini-Hochberg FDR correction
  sort.by   = "p"                 # sort by p-value
)

# Add genomic annotation to the results
ann_cols <- c("chr", "pos", "strand", "UCSC_RefGene_Name",
              "UCSC_RefGene_Group", "Relation_to_Island")
DMPs <- cbind(DMPs, ann[match(rownames(DMPs), ann$Name), ann_cols])

# How many probes are significant?
sig <- DMPs[DMPs$adj.P.Val < 0.05, ]
cat("=== DMP Results ===\n")
cat("Total probes tested:", nrow(DMPs), "\n")
cat("Significant DMPs (FDR < 0.05):", nrow(sig), "\n")
cat("  Hypermethylated in SGA (logFC > 0):", sum(sig$logFC > 0), "\n")
cat("  Hypomethylated in SGA  (logFC < 0):", sum(sig$logFC < 0), "\n")

# Save results
write.xlsx(DMPs, file.path(OUT_DIR, "DMPs_all.xlsx"),         rowNames = TRUE)
write.xlsx(sig,  file.path(OUT_DIR, "DMPs_FDR05_sig.xlsx"),   rowNames = TRUE)
cat("DMP results saved\n")

# --- Volcano plot ---
# Shows every probe: x-axis = effect size, y-axis = significance
# Red = hypermethylated in SGA, Blue = hypomethylated in SGA
DMPs$Direction <- "Not significant"
DMPs$Direction[DMPs$adj.P.Val < 0.05 & DMPs$logFC > 0] <- "Hyper in SGA"
DMPs$Direction[DMPs$adj.P.Val < 0.05 & DMPs$logFC < 0] <- "Hypo in SGA"

# FDR threshold line
fdr_line <- if (any(DMPs$adj.P.Val < 0.05)) {
  -log10(max(DMPs$P.Value[DMPs$adj.P.Val < 0.05]))
} else { Inf }

p_volc <- ggplot(DMPs, aes(x = logFC, y = -log10(P.Value), colour = Direction)) +
  geom_point(size = 0.6, alpha = 0.5) +
  scale_colour_manual(values = c(
    "Not significant" = "grey75",
    "Hyper in SGA"    = "#C0392B",
    "Hypo in SGA"     = "#2471A3"
  )) +
  geom_hline(yintercept = fdr_line,
             linetype = "dashed", colour = "black", linewidth = 0.4) +
  annotate("text", x = max(DMPs$logFC) * 0.8, y = fdr_line + 0.2,
           label = "FDR = 0.05", size = 3) +
  labs(
    title    = "Volcano Plot — Placenta DMPs (SGA vs Control)",
    subtitle = paste0("Significant: ", nrow(sig), " probes"),
    x        = "log2 Fold Change (M-value)",
    y        = expression(-log[10](p-value)),
    colour   = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(OUT_DIR, "Fig2_Volcano.pdf"), p_volc, width = 7, height = 5)
ggsave(file.path(OUT_DIR, "Fig2_Volcano.png"), p_volc, width = 7, height = 5, dpi = 300)
cat("Volcano plot saved\n")

# --- Heatmap of top significant DMPs ---
if (nrow(sig) >= 2) {
  n_heatmap <- min(100, nrow(sig))
  top_probes <- rownames(sig)[1:n_heatmap]

  # Annotation bars above the heatmap columns
  col_ann <- data.frame(
    Group      = pData(GRSet)$Sample_Group,
    Foetal_Sex = pData(GRSet)$Foetal_Sex,
    row.names  = sampleNames(GRSet)
  )
  ann_colours <- list(
    Group      = c(Control = "#27AE60", SGA = "#E67E22"),
    Foetal_Sex = c(F = "#E91E8C", M = "#1E90FF")
  )

  pdf(file.path(OUT_DIR, "Fig3_Heatmap_topDMPs.pdf"), width = 10, height = 10)
  pheatmap(
    beta[top_probes, ],
    annotation_col    = col_ann,
    annotation_colors = ann_colours,
    show_rownames     = FALSE,
    show_colnames     = FALSE,
    scale             = "row",              # scale each probe to mean=0, sd=1
    clustering_method = "ward.D2",
    color = colorRampPalette(c("#2471A3", "white", "#C0392B"))(100),
    main  = paste0("Top ", n_heatmap, " DMPs — Placenta SGA vs Control"),
    fontsize = 9
  )
  dev.off()
  cat("Heatmap saved\n")
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 10 : DIFFERENTIALLY METHYLATED REGIONS (DMRs)
# ─────────────────────────────────────────────────────────────────────────────
#
# A DMR is a stretch of genomic DNA where multiple neighbouring CpGs are
# all differentially methylated in the same direction.
# DMRs are more biologically meaningful than individual DMPs because:
#   - They are less likely to be noise
#   - They more likely reflect genuine gene regulation changes
#
# DMRcate works by:
#   1. Identifying probes with differential methylation signal
#   2. Smoothing the signal across neighbouring probes
#   3. Calling regions where the smoothed signal exceeds a threshold
#
# Parameters:
#   lambda   = bandwidth (bp) for smoothing — 1000bp is standard
#   C        = scaling factor — 2 is recommended for EPIC arrays
#   min.cpgs = minimum CpGs in a DMR — we require at least 3

cat("Annotating CpGs for DMRcate...\n")

cpg_ann <- cpg.annotate(
  datatype      = "array",
  object        = M_val,
  what          = "M",
  arraytype     = "EPIC",
  analysis.type = "differential",
  design        = design,
  coef          = "Sample_GroupSGA",
  fdr           = 0.05
)

cat("Finding DMRs...\n")

DMR_obj <- dmrcate(
  object   = cpg_ann,
  lambda   = 1000,
  C        = 2,
  min.cpgs = 3
)

# Extract as a data frame
DMR_df <- as.data.frame(extractRanges(DMR_obj, genome = "hg19"))
DMR_df <- DMR_df[order(DMR_df$Stouffer), ]  # sort by significance

cat("=== DMR Results ===\n")
cat("Total DMRs found:", nrow(DMR_df), "\n")
cat("DMRs with Stouffer p < 0.05:", sum(DMR_df$Stouffer < 0.05, na.rm = TRUE), "\n")

if (nrow(DMR_df) > 0) {
  cat("\nTop 5 DMRs:\n")
  print(head(DMR_df[, c("seqnames","start","end","no.cpgs","meandiff","Stouffer")], 5))
}

write.xlsx(DMR_df, file.path(OUT_DIR, "DMRs_all.xlsx"), rowNames = FALSE)
cat("DMR results saved\n")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 11 : GENE ONTOLOGY ENRICHMENT
# ─────────────────────────────────────────────────────────────────────────────
#
# Once we have a list of significant DMPs, we ask: what biological processes
# are these genes involved in?
#
# We use missMethyl::gometh() which is specifically designed for methylation
# data. It accounts for the fact that some genes have many more probes than
# others on the EPIC array (otherwise large genes would be artificially
# enriched in GO results).
#
# GO = Gene Ontology (a standardised vocabulary for biological processes,
#      molecular functions, and cellular components)
#
# We test three GO categories:
#   BP = Biological Process  (e.g. "regulation of transcription")
#   MF = Molecular Function  (e.g. "DNA binding")
#   CC = Cellular Component  (e.g. "nucleus")

if (nrow(sig) >= 10) {

  sig_probes <- rownames(sig)
  all_probes <- rownames(DMPs)

  cat("Running GO enrichment on", length(sig_probes), "significant probes...\n")

  # GO enrichment
  GO_res <- gometh(
    sig.cpg    = sig_probes,
    all.cpg    = all_probes,
    collection = "GO",
    array.type = "EPIC",
    prior.prob = TRUE   # correct for number of probes per gene
  )
  GO_res  <- GO_res[order(GO_res$FDR), ]
  GO_sig  <- GO_res[GO_res$FDR < 0.05, ]

  cat("Significant GO terms (FDR < 0.05):", nrow(GO_sig), "\n")
  if (nrow(GO_sig) > 0) {
    cat("\nTop 10 GO terms:\n")
    print(head(GO_sig[, c("ONTOLOGY","TERM","N","DE","FDR")], 10))
  }
  write.xlsx(GO_res,  file.path(OUT_DIR, "GO_enrichment_all.xlsx"),  rowNames = TRUE)
  write.xlsx(GO_sig,  file.path(OUT_DIR, "GO_enrichment_FDR05.xlsx"), rowNames = TRUE)

  # KEGG pathway enrichment
  KEGG_res <- gometh(
    sig.cpg    = sig_probes,
    all.cpg    = all_probes,
    collection = "KEGG",
    array.type = "EPIC",
    prior.prob = TRUE
  )
  KEGG_res <- KEGG_res[order(KEGG_res$FDR), ]
  write.xlsx(KEGG_res, file.path(OUT_DIR, "KEGG_enrichment.xlsx"), rowNames = TRUE)
  cat("KEGG enrichment saved\n")

  # Dotplot of top GO Biological Process terms
  top_bp <- head(GO_sig[GO_sig$ONTOLOGY == "BP", ], 20)

  if (nrow(top_bp) > 0) {
    top_bp$TERM <- factor(top_bp$TERM, levels = rev(top_bp$TERM))

    p_go <- ggplot(top_bp, aes(x = -log10(FDR), y = TERM,
                                size = DE, colour = -log10(FDR))) +
      geom_point() +
      scale_colour_gradient(low = "#AED6F1", high = "#1A5276") +
      scale_size_continuous(range = c(2, 8)) +
      labs(
        title    = "GO Biological Process Enrichment",
        subtitle = "Significant DMPs (FDR < 0.05)",
        x        = expression(-log[10](FDR)),
        y        = NULL,
        size     = "DE probes",
        colour   = expression(-log[10](FDR))
      ) +
      theme_classic(base_size = 11) +
      theme(axis.text.y = element_text(size = 9))

    ggsave(file.path(OUT_DIR, "Fig4_GO_dotplot.pdf"), p_go, width = 9, height = 7)
    ggsave(file.path(OUT_DIR, "Fig4_GO_dotplot.png"), p_go, width = 9, height = 7, dpi = 300)
    cat("GO dotplot saved\n")
  }

} else {
  cat("Fewer than 10 significant DMPs — skipping GO enrichment.\n")
  cat("This can happen with a small sample size or stringent FDR threshold.\n")
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 12 : EPIGENETIC CLOCK
# ─────────────────────────────────────────────────────────────────────────────
#
# Epigenetic clocks estimate "biological age" from DNA methylation patterns.
# For placenta we use the Lee 2019 clock, which was trained specifically on
# placental tissue and estimates gestational age.
#
# "Epigenetic age acceleration" = difference between the clock's predicted age
# and the actual gestational age. Positive = biologically older than expected.
#
# If SGA placentas show significant age acceleration it suggests premature
# or accelerated epigenetic ageing, which has implications for fetal development.

if (requireNamespace("methylclock", quietly = TRUE)) {
  library(methylclock)

  cat("Running epigenetic clock (Lee 2019 placenta clock)...\n")

  # DNAmAge requires CpG names as rows, samples as columns
  clock_out <- DNAmAge(
    x      = beta,
    clocks = "Lee2019"
  )

  clock_df <- data.frame(
    Sample      = sampleNames(GRSet),
    Group       = pData(GRSet)$Sample_Group,
    Foetal_Sex  = pData(GRSet)$Foetal_Sex,
    GA_actual   = pData(GRSet)$GA_weeks,
    EpiAge      = clock_out$Lee2019,
    AgeAccel    = clock_out$Lee2019 - pData(GRSet)$GA_weeks
  )

  write.xlsx(clock_df, file.path(OUT_DIR, "Epigenetic_clock.xlsx"), rowNames = FALSE)

  # Compare age acceleration between groups
  t_res <- t.test(AgeAccel ~ Group, data = clock_df)
  cat("\nEpigenetic age acceleration — t-test:\n")
  cat("  SGA mean:     ", round(mean(clock_df$AgeAccel[clock_df$Group == "SGA"],     na.rm=TRUE), 2), "weeks\n")
  cat("  Control mean: ", round(mean(clock_df$AgeAccel[clock_df$Group == "Control"], na.rm=TRUE), 2), "weeks\n")
  cat("  p-value:      ", round(t_res$p.value, 4), "\n")

  # Boxplot
  p_clock <- ggplot(clock_df, aes(x = Group, y = AgeAccel, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.12, alpha = 0.6, size = 1.8, aes(colour = Foetal_Sex)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values  = c("Control" = "#27AE60", "SGA" = "#E67E22")) +
    scale_colour_manual(values = c("F" = "#E91E8C", "M" = "#1E90FF"),
                        name  = "Foetal sex") +
    labs(
      title    = "Placental Epigenetic Age Acceleration (Lee 2019)",
      subtitle = paste0("p = ", round(t_res$p.value, 4)),
      x        = NULL,
      y        = "Age Acceleration (weeks)",
      fill     = "Group"
    ) +
    theme_classic(base_size = 12)

  ggsave(file.path(OUT_DIR, "Fig5_EpiClock.pdf"), p_clock, width = 5, height = 5)
  ggsave(file.path(OUT_DIR, "Fig5_EpiClock.png"), p_clock, width = 5, height = 5, dpi = 300)
  cat("Epigenetic clock plot saved\n")

} else {
  cat("methylclock not installed. Install with:\n")
  cat("  remotes::install_github('isglobal-brge/methylclock')\n")
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 13 : SENSITIVITY ANALYSES
# ─────────────────────────────────────────────────────────────────────────────
#
# Sensitivity analyses test whether our findings hold in subgroups.
# Because foetal sex strongly influences methylation (X-inactivation etc.),
# we repeat the DMP analysis separately in female and male foetuses.
#
# If the same DMPs appear in both sexes, the finding is more robust.
# If DMPs only appear in one sex, the biology may be sex-specific.

run_dmp <- function(grset, label) {
  pd_sub  <- pData(grset)
  des_sub <- model.matrix(~ Sample_Group + GA_weeks, data = pd_sub)
  fit_sub <- eBayes(lmFit(getM(grset), des_sub))
  res_sub <- topTable(fit_sub, coef = "Sample_GroupSGA",
                      number = Inf, adjust = "BH", sort.by = "p")
  cat(label, "— significant DMPs (FDR<0.05):", sum(res_sub$adj.P.Val < 0.05), "\n")
  write.xlsx(res_sub,
             file.path(OUT_DIR, paste0("DMPs_", label, ".xlsx")),
             rowNames = TRUE)
  return(res_sub)
}

# Female foetuses
GRSet_F <- GRSet[, pData(GRSet)$Foetal_Sex == "F"]
cat("Female foetus analysis (n=", ncol(GRSet_F), "):\n")
DMPs_F  <- run_dmp(GRSet_F, "female_foetus")

# Male foetuses
GRSet_M <- GRSet[, pData(GRSet)$Foetal_Sex == "M"]
cat("Male foetus analysis (n=", ncol(GRSet_M), "):\n")
DMPs_M  <- run_dmp(GRSet_M, "male_foetus")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 14 : SAVE EVERYTHING AND WRAP UP
# ─────────────────────────────────────────────────────────────────────────────

# Save the main R objects so you can reload them without re-running the whole
# pipeline (useful for making additional plots or checking specific results)
save(
  GRSet, beta, M_val,
  DMPs, DMR_df,
  file = file.path(OUT_DIR, "Placenta_analysis_objects.RData")
)
cat("R objects saved (reload with: load('Placenta_analysis_objects.RData'))\n")

# Save session info — important for reproducibility and methods section
sink(file.path(OUT_DIR, "session_info.txt"))
cat("Analysis completed:", format(Sys.time()), "\n\n")
sessionInfo()
sink()

# ── Final summary ─────────────────────────────────────────────────────────────
cat("\n")
cat(strrep("═", 60), "\n")
cat("  ANALYSIS COMPLETE\n")
cat(strrep("═", 60), "\n")
cat("  Results folder:", OUT_DIR, "\n\n")
cat("  Key output files:\n")
cat("  ├── QC_report.pdf\n")
cat("  ├── Fig1_PCA.png\n")
cat("  ├── Fig2_Volcano.png\n")
cat("  ├── Fig3_Heatmap_topDMPs.pdf\n")
cat("  ├── Fig4_GO_dotplot.png\n")
cat("  ├── Fig5_EpiClock.png\n")
cat("  ├── DMPs_all.xlsx           (all probes with stats)\n")
cat("  ├── DMPs_FDR05_sig.xlsx     (significant DMPs only)\n")
cat("  ├── DMRs_all.xlsx\n")
cat("  ├── GO_enrichment_FDR05.xlsx\n")
cat("  ├── KEGG_enrichment.xlsx\n")
cat("  ├── Epigenetic_clock.xlsx\n")
cat("  ├── DMPs_female_foetus.xlsx\n")
cat("  ├── DMPs_male_foetus.xlsx\n")
cat("  └── Placenta_analysis_objects.RData\n")
cat(strrep("═", 60), "\n")
