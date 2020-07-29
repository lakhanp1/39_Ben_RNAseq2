suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(org.HSapiens.gencodev30.eg.db))

## plot combined log2FoldChange heatmap of all SM OE polII strains

rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################

analysisName <- "DEG_summary"
outDir <- here::here("analysis", "03_DEG_summary")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

if(!dir.exists(outDir)){
  dir.create(outDir, recursive = T)
}

diffDataPath <- here::here("analysis", "02_DESeq2_diff")
file_sampleInfo <- here::here("data", "reference_data", "sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "DESeq2_DEG_info.txt")

cutoff_fdr <- 0.05
cutoff_lfc <- 0.585
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1

orgDb <- org.HSapiens.gencodev30.eg.db

###########################################################################
rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath)

rowId <- 1
degData <- NULL

for (rowId in 1:nrow(rnaseqInfo)) {
  
  tmpDf <- suppressMessages(readr::read_tsv(file = rnaseqInfo$deg[rowId])) %>% 
    dplyr::select(geneId, contrast, log2FoldChange, padj, GENE_NAME) %>% 
    dplyr::mutate(comparison = rnaseqInfo$comparison[rowId])
  
  degData <- dplyr::bind_rows(degData, tmpDf)
}


degData <- dplyr::mutate(
  degData,
  # contrast = stringr::str_replace(string = contrast, pattern = "(condition_|SM_TF_)", replacement = ""),
  diff = dplyr::case_when(
    padj <= cutoff_fdr & log2FoldChange <= cutoff_down ~ "down",
    padj <= cutoff_fdr & log2FoldChange >= cutoff_up ~ "up",
    TRUE ~ "noDEG"
  )
) %>% 
  dplyr::filter(diff != "noDEG")

degStats <- dplyr::group_by(degData, comparison, diff) %>% 
  dplyr::count() %>% 
  tidyr::pivot_wider(
    names_from = diff, values_from = n, values_fill = list(n = 0)
  )

readr::write_tsv(x = degStats, path = paste(outPrefix, ".stats.tab", sep = ""))

plotData <- dplyr::left_join(x = degStats, y = rnaseqInfo, by = "comparison") %>% 
  dplyr::mutate_at(.vars = c("drug", "concentration", "time"), .funs = forcats::as_factor) %>%
  dplyr::mutate(
    drug = forcats::fct_relevel(.f = drug, "Tet", "ADTet"),
    concentration = forcats::fct_relevel(.f = concentration, "5uM", "15uM"),
    time = forcats::fct_relevel(.f = time, "12hr", "24hr")
  ) %>% 
  tidyr::unite(col = "treatment", drug, concentration, sep = ":", remove = FALSE)


pt_stats <- ggplot(
  data = plotData,
  mapping = aes(y = forcats::fct_rev(treatment), fill = time)
) +
  geom_bar(
    mapping = aes(x = up), color = "black", stat = "identity",
    position = position_dodge(), width = 0.7, size = 1
  ) +
  geom_bar(
    mapping = aes(x = -down), color = "black", stat = "identity",
    position = position_dodge(), width = 0.7, size = 1
  ) +
  geom_vline(xintercept = 0, size = 1) +
  labs(
    title = "Significant DEG stats",
    x = "Significant differentially expressed genes"
  ) +
  scale_fill_manual(
    values = c("12hr" = "black", "24hr" = "white")
  ) +
  scale_x_continuous(
    breaks = seq(-1500, 1500, by = 500),
    labels = abs(seq(-1500, 1500, by = 500))
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(2, "mm"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    legend.key.size = unit(2, "cm"),
    legend.text = element_text(size = 16, face = "bold", color = "black"),
    legend.title = element_text(size = 18, face = "bold")
  )

ggsave(filename = paste(outPrefix, ".stats.pdf", sep = ""),
       plot = pt_stats, width = 10, height = 6)


# 
# ## a dataframe with DEG status column for each comparison
# diffDf <- tidyr::pivot_wider(
#   data = degData,
#   id_cols = c(geneId),
#   names_from = comparison,
#   values_from = diff,
#   values_fill = list(diff = "noDEG")
# )
# 
# readr::write_tsv(x = diffDf, path = paste(outPrefix, ".diff.tab", sep = ""))
# 
# 
# ###########################################################################
# ## heatmap of log2FoldChanges
# lfcDf <- tidyr::pivot_wider(
#   data = degData,
#   id_cols = c(geneId),
#   names_from = comparison,
#   values_from = log2FoldChange,
#   values_fill = list(log2FoldChange = 0)
# )
# 
# readr::write_tsv(x = lfcDf, path = paste(outPrefix, ".log2FoldChange.tab", sep = ""))
# 
# lfcMat <- as.data.frame(lfcDf) %>% 
#   tibble::column_to_rownames(var = "geneId") %>% 
#   as.matrix()
# 
# colAnDf <- tibble::tibble(comparison = colnames(lfcMat)) %>% 
#   dplyr::left_join(y = degStats, by = "comparison") %>% 
#   dplyr::mutate(
#     total = up + down
#   ) %>% 
#   as.data.frame() %>% 
#   tibble::column_to_rownames(var = "comparison")
# 
# 
# htAn <- ComplexHeatmap::HeatmapAnnotation(
#   copy = colAnDf$copy,
#   up = anno_text(
#     x = colAnDf$up, rot = 90, just = "center", location = 0.5,
#     gp = gpar(fill = "#ffb3c4", col = "black", border = "black"),
#     height = max_text_width(text = colAnDf$up)*1.5
#   ),
#   down = anno_text(
#     x = colAnDf$down, rot = 90, just = "center", location = 0.5,
#     gp = gpar(fill = "#d8daf3", col = "black", border = "black"),
#     height = max_text_width(text = colAnDf$down)*1.5
#   ),
#   fraction = anno_text(
#     x = colAnDf$fraction, rot = 90, just = "center", location = 0.5,
#     gp = gpar(col = "black", border = "black"),
#     height = max_text_width(text = colAnDf$fraction)*1.2
#   ),
#   which = "column",
#   col = list(
#     copy = c("sCopy" = "#ff7f00", "msCopy" = "#8dd3c7")
#   ),
#   annotation_label = c(
#     "copy number", "DEG Up", "DEG down", "fraction"
#   ),
#   show_annotation_name = TRUE,
#   annotation_name_side = "right",
#   annotation_legend_param = list(
#     copy = list(
#       title = "Copy number", ncol = 1, title_position = "leftcenter-rot",
#       title_gp = gpar(fontsize = 16, fontface = "bold"),
#       labels_gp = gpar(fontsize = 14)
#     )
#   )
# )
# 
# 
# ht_lfc <- ComplexHeatmap::Heatmap(
#   matrix = lfcMat,
#   name = "log2FoldChange",
#   col = colorRamp2(breaks = c(-3, 0, 3), colors = c("#313695", "white", "#a50026")),
#   top_annotation = htAn,
#   row_title = "genes",
#   column_title = "log2FoldChange heatmap of all significant DEGs in all OE/WT comparisons", 
#   column_title_gp = gpar(fontsize = 18, fontface = "bold"),
#   show_row_names = FALSE,
#   row_dend_reorder = TRUE, column_dend_reorder = TRUE,
#   column_labels = stringr::str_replace(
#     string = colnames(lfcMat), pattern = "_vs_MH11036", replacement = ""
#   ),
#   heatmap_legend_param = list(
#     title = "log2(fold-change)",  title_position = "leftcenter-rot",
#     legend_height = unit(4, "cm"), title_gp = gpar(fontsize = 16, fontface = "bold"),
#     labels_gp = gpar(fontsize = 14)
#   ),
#   use_raster = FALSE
# )
# 
# htList <- ht_lfc
# 
# png(filename = paste(outPrefix, ".lfc_heatmap.png", sep = ""), width = 6000, height = 4000, res = 400)
# draw(
#   htList,
#   merge_legends = TRUE,
#   legend_gap = unit(3, "cm"),
#   padding = unit(c(0.5, 0.5, 0.5, 1), "cm")
# )
# dev.off()
# 
# ###########################################################################




