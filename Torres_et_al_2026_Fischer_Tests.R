install.packages('readxl')
install.packages('dplyr')
library(readxl)
library(dplyr)

#Fill in with path location of Dataset 1 and 2
#dataset1 = 
#dataset2 = 
sheetnames1 = excel_sheets(dataset1)
sheetnames = excel_sheets(dataset2)

#read in Dataset 1
rna <- lapply(sheetnames1, function(sheet) {
  read_excel(dataset1, sheet = sheet, skip = 1)
})

names(rna) <- sheetnames1

#read in Dataset 2
dfs <- lapply(sheetnames, function(sheet) {
  read_excel(dataset2, sheet = sheet, skip = 1)
})

names(dfs) <- sheetnames

###Get Genes of Interest###
#K27 genes is our "Background" for tests (n = 571)
K27genes = dfs$`Dataset S2K`
#using correct K36 data tables (REPLACE w/ DATASET 2 K36 Data to replicate)
K36data = read.csv(file = "bed_csv_txt_files/H3K36me3_Results_allgenespromTSS_vWT.csv")
K36datainput = read.csv(file = "bed_csv_txt_files/H3K36me3_Results_allgenespromTSS_vinput.csv")
#K27 data from Dataset 2
K27data = dfs$`Dataset S2C`
###########################

###Genes that lose H3K36me3 & H3K27me3 in CAF-1 mutants###
#First, subset to H3K27me3-marked genes also marked by K36 (vs input)
K36enriched = K36datainput %>%
  filter(if_any(all_of(c(6, 9, 12, 15, 18, 24)), ~ coalesce(.x > 0, FALSE)))
#genes that lose K36 in promoters
K36datasub = subset(K36data, K36data$gene_id %in% K36enriched$gene_id)
K36lost = subset(K36datasub, (K36datasub$logFC_cac1 < -0.75 | K36datasub$logFC_cac2 < -0.75) & K36datasub$gene_id %in% K27genes$Name)
#Lose K27
K27lost = subset(K27data, (K27data$logFC_cac1 < -0.75 | K27data$logFC_cac2 < -0.75) & K27data$Name %in% K27genes$Name)
##########################################################

###Genes that lose H3K27me3 in CAF mutants & ASH-1 mutant###
#lose H3K27me3 in ash-1
K27lost_ash1 = subset(K27data, K27data$logFC_ash1 < -0.75 & K27data$Name %in% K27genes$Name)
#lose H3K27me3 in ash-1 and caf-1
############################################################

###K27-marked Genes upregulated in CAF-1 and ASH-1###
#First, get H3K27me3-marked genes that are actually protein coding for "background" (n = 545)
expressed_genes = rna$`Datset S1B`
#get relevant dataframes
cac1_data = subset(expressed_genes, expressed_genes$Name %in% K27genes$Name)
cac2_data = rna$`Dataset S1C`
ash1_data = rna$`Dataset S1H`
cac1_data$padj <- as.numeric(as.character(cac1_data$padj))
cac2_data$padj <- as.numeric(as.character(cac2_data$padj))
ash1_data$padj <- as.numeric(as.character(ash1_data$padj))
cac2_data_sub = subset(cac2_data, cac2_data$Name %in% cac1_data$Name)
ash1_data_sub = subset(ash1_data, ash1_data$Name %in% cac1_data$Name)
#Genes Upregulated in CAF-1
CAF1_UP = subset(cac1_data, (cac1_data$log2FoldChange > 1 & cac1_data$padj < 0.05) | (cac2_data_sub$log2FoldChange > 1 & cac2_data_sub$padj < 0.05))
#Genes upregulated in ash-1
ASH1_UP = subset(ash1_data_sub, ash1_data_sub$log2FoldChange > 1 & ash1_data_sub$padj < 0.05)

#####################################################

###Fishcer Test Function###

run_fisher_overlap <- function(background_df, setA_df, setB_df,
                               background_col = "Name",
                               setA_col = "Name",
                               setB_col = "Name",
                               test_name = "Fisher test") {
  
  # Get unique gene names
  background <- unique(background_df[[background_col]])
  setA <- unique(setA_df[[setA_col]])
  setB <- unique(setB_df[[setB_col]])
  
  # Restrict sets to background
  setA <- intersect(setA, background)
  setB <- intersect(setB, background)
  
  # Counts
  overlap <- length(intersect(setA, setB))
  setA_only <- length(setdiff(setA, setB))
  setB_only <- length(setdiff(setB, setA))
  neither <- length(setdiff(background, union(setA, setB)))
  
  # 2x2 contingency table
  mat <- matrix(
    c(overlap, setA_only,
      setB_only, neither),
    nrow = 2,
    byrow = TRUE
  )
  
  rownames(mat) <- c("SetA_yes", "SetA_no")
  colnames(mat) <- c("SetB_yes", "SetB_no")
  
  # Fisher's exact test for enrichment
  fisher_result <- fisher.test(mat, alternative = "greater")
  
  # Clean summary
  summary_df <- data.frame(
    test = test_name,
    background_n = length(background),
    setA_n = length(setA),
    setB_n = length(setB),
    overlap_n = overlap,
    overlap_percent_of_setA = round(100 * overlap / length(setA), 2),
    odds_ratio = unname(fisher_result$estimate),
    p_value = fisher_result$p.value
  )
  
  return(list(
    table = mat,
    fisher = fisher_result,
    summary = summary_df
  ))
}
###########################

###Run Tests###
#K27 and K36
test_K36_K27_loss <- run_fisher_overlap(
  background_df = K27genes,
  setA_df = K36lost,
  setB_df = K27lost,
  background_col = "Name",
  setA_col = "gene_id",
  setB_col = "Name",
  test_name = "CAF-1 K36 loss vs K27 loss"
)

#K27 CAF-1 vs ASH-1
test_CAF1_ASH1_K27_loss <- run_fisher_overlap(
  background_df = K27genes,
  setA_df = K27lost,
  setB_df = K27lost_ash1,
  background_col = "Name",
  setA_col = "Name",
  setB_col = "Name",
  test_name = "CAF-1 K27 loss vs ash-1 K27 loss"
)

#Upregulated CAF-1 vs ASH-1
test_CAF1_ASH1_UP <- run_fisher_overlap(
  background_df = cac1_data,
  setA_df = CAF1_UP,
  setB_df = ASH1_UP,
  background_col = "Name",
  setA_col = "Name",
  setB_col = "Name",
  test_name = "CAF-1 upregulated vs ash-1 upregulated"
)

#Combine Test data and write to excel
all_fisher_summaries <- bind_rows(
  test_K36_K27_loss$summary,
  test_CAF1_ASH1_K27_loss$summary,
  test_CAF1_ASH1_UP$summary
)

write.csv(all_fisher_summaries, file = "fisher_exact_test_summaries.csv", row.names = FALSE)

###############
