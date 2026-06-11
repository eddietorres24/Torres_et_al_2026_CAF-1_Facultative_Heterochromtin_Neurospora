#!/usr/bin/env Rscript

## --- set working directory to this script's folder ---
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa   <- grep("^--file=", args, value = TRUE)
  if (length(fa) == 1) return(dirname(normalizePath(sub("^--file=", "", fa))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
}
setwd(get_script_dir()); cat("Working directory set to:", getwd(), "\n")

## ================== user settings (edit) ==================
counts_file      <- "readcounts_All_CAF1paper.merged.txt"   # featureCounts output
coldata_file     <- "coldata.csv"                            # must have columns: sample_id, condition
names_file       <- "names.txt"                              # mapping columns -> desired names
gene_map_file    <- "gene_annotation.csv"                    # optional: gene_id,gene_symbol
k27_bed_file     <- "H3K27me3_methylated_genes_FINAL.bed"    # optional: IDs in col 10 (not required)
outdir           <- "RNAseq_QC_out"
reference_level  <- "WT"                                     # reference in 'condition'
min_count_filter <- 10                                       # keep genes with rowSums >= this
ma_ylim          <- c(-8, 8)                                 # MA plot Y limits
simplify_names   <- TRUE                                     # turn BAM paths into SRR… sample names
## ==========================================================

## Packages must already be installed locally
suppressPackageStartupMessages({
  library(DESeq2); library(pheatmap); library(data.table); library(dplyr)
  library(tibble); library(ggplot2); library(tidyr); library(readr); library(stringr)
})
have_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

## --------- helpers ----------
read_featureCounts <- function(path, simplify_names = TRUE) {
  dt <- data.table::fread(path, header = TRUE, sep = "\t", data.table = FALSE, check.names = FALSE)
  if (!"Geneid" %in% names(dt)) stop("Counts file must have a 'Geneid' column (featureCounts style).")
  meta <- c("Geneid","Chr","Start","End","Strand","Length","Gene.Symbol","Symbol","GeneSymbol")
  mat  <- dt[, setdiff(names(dt), meta), drop = FALSE]
  rownames(mat) <- dt$Geneid
  
  if (simplify_names) {
    newn <- basename(colnames(mat))
    newn <- sub("_Aligned.sortedByCoord.out.bam$", "", newn)
    newn <- sub("\\.bam$", "", newn)
    colnames(mat) <- newn
  }
  
  to_num <- function(x) as.integer(round(as.numeric(gsub(",", "", trimws(x)))))
  num <- as.data.frame(lapply(mat, to_num), check.names = FALSE, stringsAsFactors = FALSE)
  rownames(num) <- rownames(mat)
  if (anyNA(num)) {
    bad <- names(num)[colSums(is.na(num)) > 0]
    stop("Non-numeric values detected in: ", paste(bad, collapse=", "), ".")
  }
  num
}

parse_names_file <- function(file, n_expected) {
  if (!file.exists(file)) stop("names file not found: ", file)
  lines <- trimws(readLines(file, warn = FALSE))
  lines <- lines[nzchar(lines)]
  if (length(lines) > 0 && grepl('^"?[Xx]"?$', lines[1])) lines <- lines[-1]
  m <- regmatches(lines, gregexpr('"[^"]+"|\\S+', lines))
  if (all(lengths(m) == 2)) {
    nm <- vapply(m, function(tokens) gsub('^"|"$', "", tokens[2]), character(1))
  } else {
    nm <- gsub('^"|"$', "", lines)
  }
  if (length(nm) != n_expected) {
    stop("names.txt has ", length(nm), " entries, but counts has ", n_expected, " samples.\n",
         "One target name per sample column, same order as counts.")
  }
  nm
}
## -----------------------------

## I/O setup
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
figdir <- file.path(outdir, "figures"); dir.create(figdir, FALSE, TRUE)
tbldir <- file.path(outdir, "tables");  dir.create(tbldir, FALSE, TRUE)
qcdir  <- file.path(outdir, "qc");      dir.create(qcdir,  FALSE, TRUE)

cat("\n=== Loading coldata ===\n")
coldata <- read.csv(coldata_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!"sample_id" %in% names(coldata)) stop("coldata.csv must have a 'sample_id' column.")
if ("condition" %in% names(coldata)) {
  coldata$condition <- factor(coldata$condition)
  if (reference_level %in% levels(coldata$condition)) {
    coldata$condition <- stats::relevel(coldata$condition, ref = reference_level)
  }
} else {
  warning("'condition' missing; will use ~1 design (no contrasts).")
}
rownames(coldata) <- coldata$sample_id

cat("=== Loading counts (featureCounts) ===\n")
counts_raw <- read_featureCounts(counts_file, simplify_names = simplify_names)

cat("=== Renaming counts columns from names.txt ===\n")
desired_names <- parse_names_file(names_file, n_expected = ncol(counts_raw))
colnames(counts_raw) <- desired_names

missing_in_counts <- setdiff(coldata$sample_id, colnames(counts_raw))
if (length(missing_in_counts)) {
  stop("These coldata sample_id are not present in counts after renaming: ",
       paste(missing_in_counts, collapse=", "),
       "\n• Check names.txt order/entries or update coldata$sample_id.")
}
counts_df <- counts_raw[, coldata$sample_id, drop = FALSE]

## ----- normalize/limit samples & drop old cac-1 -----
# make condition strings uniform
coldata$condition <- as.character(coldata$condition)
coldata$condition <- gsub("-", "", coldata$condition)   # set-7 -> set7, naf-1 -> naf1, etc.
coldata$condition <- tolower(coldata$condition)
coldata$condition[coldata$condition == "wt"] <- "WT"     # keep WT uppercase for labels

# groups to keep in this figure
target_groups <- c("WT","set7","cac1_new","cac2","cac3","naf1","naf2","asf1","atrx","ash1")

# drop legacy cac1 (old) entirely
drop_idx <- which(coldata$condition == "cac1")
if (length(drop_idx)) {
  message("Dropping legacy cac-1 samples: ", paste(rownames(coldata)[drop_idx], collapse=", "))
  coldata   <- coldata[-drop_idx, , drop = FALSE]
}

# keep only the strains of interest
keep_idx <- which(coldata$condition %in% target_groups)
coldata   <- coldata[keep_idx, , drop = FALSE]
counts_df <- counts_df[, rownames(coldata), drop = FALSE]

# tidy factor order
coldata$condition <- factor(coldata$condition,
                            levels = target_groups[target_groups %in% unique(coldata$condition)])

## Optional: gene symbol map
gene_map <- NULL
if (file.exists(gene_map_file)) {
  gm <- read.csv(gene_map_file, stringsAsFactors = FALSE)
  if (all(c("gene_id","gene_symbol") %in% names(gm))) {
    gene_map <- gm[!duplicated(gm$gene_id), c("gene_id","gene_symbol")]
  }
}

## Optional: K27 set (not used for these all-genes QC plots)
k27_ids <- NULL
if (file.exists(k27_bed_file)) {
  bed <- tryCatch(read.delim(k27_bed_file, header = FALSE, stringsAsFactors = FALSE), error=function(e) NULL)
  if (!is.null(bed)) {
    k27_ids <- unique(if (ncol(bed) >= 10) bed[[10]] else bed[[ncol(bed)]])
    k27_ids <- k27_ids[!is.na(k27_ids)]
  }
}

## =================== Build DESeq2 & QC ===================
cat("=== Build DESeq2 & QC ===\n")
keep <- rowSums(counts_df) >= min_count_filter
counts_f <- counts_df[keep, , drop=FALSE]
design_formula <- if ("condition" %in% names(coldata)) ~ condition else ~ 1
dds <- DESeqDataSetFromMatrix(countData = counts_f, colData = coldata, design = design_formula)
dds <- estimateSizeFactors(dds)

# sample QC metrics
lib_sizes <- colSums(counts_df)
detected  <- colSums(counts_df > 0)
pct_zero  <- colSums(counts_df == 0) / nrow(counts_df) * 100
qc_tbl <- tibble(
  sample_id     = names(lib_sizes),
  library_size  = as.numeric(lib_sizes),
  detected_genes= as.numeric(detected),
  pct_zero      = round(pct_zero, 2),
  size_factor   = as.numeric(sizeFactors(dds))
)
write.csv(qc_tbl, file.path(qcdir, "sample_qc_metrics.csv"), row.names = FALSE)

# library size / detected genes
p_lib <- ggplot(qc_tbl, aes(reorder(sample_id, library_size), library_size, fill = coldata$condition)) +
  geom_col() + coord_flip() + labs(x="Sample", y="Library size", fill="condition") + theme_minimal(11)
ggsave(file.path(figdir,"library_sizes.pdf"), p_lib, width=7.5, height=6, useDingbats=FALSE)

p_det <- ggplot(qc_tbl, aes(reorder(sample_id, detected_genes), detected_genes, fill = coldata$condition)) +
  geom_col() + coord_flip() + labs(x="Sample", y="Detected genes (>0)", fill="condition") + theme_minimal(11)
ggsave(file.path(figdir,"detected_genes.pdf"), p_det, width=7.5, height=6, useDingbats=FALSE)

# export matrices
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file.path(tbldir, "normalized_counts.csv"))
write.csv(as.data.frame(counts(dds)), file.path(tbldir, "raw_counts_filtered.csv"))
write.csv(tibble(sample_id = names(sizeFactors(dds)), size_factor = sizeFactors(dds)),
          file.path(tbldir, "size_factors.csv"), row.names = FALSE)

# normalized-count violin
log_norm <- log10(norm_counts + 1)
ln_df <- as.data.frame(log_norm) |>
  tibble::rownames_to_column("gene_id") |>
  tidyr::pivot_longer(-gene_id, names_to="sample_id", values_to="log10_norm")
ln_df$condition <- coldata[ln_df$sample_id, "condition", drop=TRUE]
p_nd <- ggplot(ln_df, aes(sample_id, log10_norm, fill = condition)) +
  geom_violin(scale="width", trim=TRUE) + coord_flip() +
  labs(x="Sample", y="log10(normalized counts + 1)", fill="condition") + theme_minimal(11)
ggsave(file.path(figdir,"normalized_count_violin.pdf"), p_nd, width=7.5, height=6.5, useDingbats=FALSE)

## =================== VST, PCA, correlation ===================
cat("=== VST, PCA, correlation ===\n")
vsd <- vst(dds, blind = FALSE)
write.csv(as.data.frame(assay(vsd)), file.path(tbldir, "vst_matrix.csv"))

## ----- PCA (square aspect, symmetric limits, no point labels) -----
pca_df <- DESeq2::plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))

# Make axes symmetric around 0 and equal-scaled
xrange <- range(pca_df$PC1, na.rm = TRUE)
yrange <- range(pca_df$PC2, na.rm = TRUE)
lim <- max(abs(c(xrange, yrange)))             # common half-range
xlims <- c(-lim, lim)
ylims <- c(-lim, lim)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  scale_x_continuous(limits = xlims, expand = expansion(mult = 0.05),
                     breaks = scales::pretty_breaks()) +
  scale_y_continuous(limits = ylims, expand = expansion(mult = 0.05),
                     breaks = scales::pretty_breaks()) +
  coord_equal() +  # ensures a true square aspect (1 unit on x == 1 unit on y)
  labs(
    x = paste0("PC1: ", percentVar[1], "%"),
    y = paste0("PC2: ", percentVar[2], "%"),
    color = ""
  ) +
  theme_minimal(12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# save as a square page so nothing looks squished
ggsave(file.path(figdir, "PCA_vst_chaperones.pdf"),
       p_pca, width = 6.5, height = 6.5, useDingbats = FALSE)


# ---------- Correlation heatmap (original colors; legend ON; auto-size) ----------
cm  <- cor(assay(vsd), method = "pearson")
ann <- data.frame(condition = coldata$condition); rownames(ann) <- rownames(coldata)

# Auto-size so labels + legend fit; DO NOT set colors/breaks (keeps original look)
n <- ncol(cm)
pdf_w <- max(7.5, 0.35 * n + 3.5)
pdf_h <- max(6.0, 0.35 * n + 2.5)

pheatmap(cm,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         legend            = TRUE,          # scale bar back on
         annotation_legend = TRUE,
         annotation_col    = ann,
         annotation_row    = ann,
         display_numbers   = FALSE,
         show_rownames     = TRUE,
         show_colnames     = TRUE,
         fontsize          = 10,
         fontsize_col      = 8,
         angle_col         = 45,
         main              = "Sample–sample correlation (Pearson, VST)",
         filename          = file.path(figdir, "correlation_heatmap_vst_chaperones.pdf"),
         width             = pdf_w,
         height            = pdf_h
)

## ======= DESeq2 fit, dispersion, MA plots & p-value hist (no shrinkage) =======
cat("=== DESeq2 fit (dispersion, MA, p-hist; no shrinkage) ===\n")
dds <- DESeq(dds)

# Dispersion plot (Illustrator-safe export)
pdf(file.path(figdir, "dispersion_plot.pdf"),
    width = 6, height = 5, useDingbats = FALSE)
plotDispEsts(dds, cex = 0.6)   # no pch=; just slightly larger points
dev.off()

if ("condition" %in% names(coldata)) {
  ref <- levels(coldata$condition)[1]  # after relevel should equal reference_level
  for (lvl in setdiff(levels(coldata$condition), ref)) {
    res <- results(dds, contrast = c("condition", lvl, ref), alpha = 0.05)
    res_df <- as.data.frame(res) |> tibble::rownames_to_column("gene_id")
    
    if (!is.null(gene_map)) res_df <- dplyr::left_join(res_df, gene_map, by="gene_id")
    write.csv(res_df, file.path(tbldir, paste0("DE_results_", lvl, "_vs_", ref, ".csv")), row.names = FALSE)
    
    pdf(file.path(figdir, paste0("MA_", lvl, "_vs_", ref, ".pdf")),
        width = 6.2, height = 5.2, useDingbats = FALSE)
    plotMA(res, ylim = ma_ylim, cex = 0.8, main = paste(lvl, "vs", ref))
    dev.off()
    
    pdf(file.path(figdir, paste0("pvalue_hist_", lvl, "_vs_", ref, ".pdf")), width=5.5, height=4.5, useDingbats=TRUE)
    hist(res$pvalue, breaks=50, col="grey", main=paste("p-value distribution:", lvl, "vs", ref), xlab="p-value"); dev.off()
  }
} else {
  message("No 'condition' column; skipping contrasts/MA/p-hist.")
}

writeLines(c(capture.output(sessionInfo())), con = file.path(outdir, "sessionInfo.txt"))
cat("\nAll done. Outputs in: ", normalizePath(outdir), "\n")

