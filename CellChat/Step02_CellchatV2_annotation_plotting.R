#================================== NOTE ==================================
# NOTE ## CellChat --- Plotting ReinstatementvsNR (Dataset1: NR; Dataset2: Reinstatement)
#================================== NOTE ==================================

set.seed(123456)

library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(UCell)
library(scCustomize)
library(viridis)
theme_set(theme_cowplot())
library(ggplot2)
library(ggrastr)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

# Add this line to print warnings as they occur
options(warn = 1)

umap_theme <- theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

#
# Set working directory
setwd('/dfs7/swaruplab/shared_lab/Collaborations/Wood/Reinstatement_Project/reinstatement_2022/CellChat/annotation/')

fig_dir <- "/dfs7/swaruplab/shared_lab/Collaborations/Wood/Reinstatement_Project/reinstatement_2022/CellChat/annotation/figures/"
data_dir <- '/dfs7/swaruplab/shared_lab/Collaborations/Wood/Reinstatement_Project/reinstatement_2022/CellChat/annotation/data/'

compare <- "NRvsReinstatement"


##========================================
# NOTE loading the data
# Define a vector with the paths to the .rds files
rds_files <- c(
  '/dfs7/swaruplab/shared_lab/Collaborations/Wood/Reinstatement_Project/reinstatement_2022/CellChat/annotation/data/NR_annotation_cellchat_filtered_n_processed_Reinstatement_0930_2024.rds',
  # '/dfs7/swaruplab/zechuas/Projects/AD_GM_WM_2022/snRNA_seq/Analysis/Cellchat/PFC/ADvsControl/data/AD_PFC_cellchat_filtered_n_processed_AD_GMWM_0924_2024.rds',
  '/dfs7/swaruplab/shared_lab/Collaborations/Wood/Reinstatement_Project/reinstatement_2022/CellChat/annotation/data/Reinstatement_annotation_cellchat_filtered_n_processed_Reinstatement_0930_2024.rds'
)


# # Assign condition names corresponding to the files
conditions <- c("NR", "Reinstatement")

# Create an empty list to hold the data
cellchat_list <- list()

# Use a loop or lapply to read and assign the files to the list with the condition names
for (i in seq_along(rds_files)) {
  cellchat_list[[conditions[i]]] <- readRDS(rds_files[i])
}

cellchat_list

# Assuming the previous steps have been completed, and cellchat_list is filled with the data
# Now merge the datasets into one cellchat object:
cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
cellchat
# > cellchat

# # #=================
# NOTE change for later plotting -- make it easier
object.list <- cellchat_list

# NOTE Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), color.use = c("darkorchid3", "seagreen"))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), measure = "weight", color.use = c("darkorchid3", "seagreen"))

pdf(paste0(fig_dir, compare, '_', 'cellchat_compareInteractions_1003_2024.pdf'), width=10, height=6)
gg1 + gg2
dev.off()



# #=================
# NOTE Compare the number of interactions and interaction strength among different cell populations
pdf(paste0(fig_dir, compare, '_', 'cellchat_netVisual_diffInteraction_1003_2024.pdf'), width=10, height=8)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat, comparison = c(1, 2), weight.scale = T, vertex.weight = 5)
netVisual_diffInteraction(cellchat, comparison = c(1, 2), weight.scale = T, vertex.weight = 5, measure = "weight")
dev.off()

pdf(paste0(fig_dir, compare, '_', 'cellchat_netVisual_diffInteraction_1003_2024_colorUpdates.pdf'), width=10, height=8)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat, comparison = c(2, 1), weight.scale = T, vertex.weight = 5, color.edge = c("darkorchid3", "seagreen"))
netVisual_diffInteraction(cellchat, comparison = c(2, 1), weight.scale = T, vertex.weight = 5, measure = "weight", color.edge = c("darkorchid3", "seagreen"))
dev.off()


# (B) Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets
gg1 <- netVisual_heatmap(cellchat, color.heatmap = c("darkorchid3", "seagreen"))
gg2 <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = c("darkorchid3", "seagreen"))

pdf(paste0(fig_dir, compare, '_', 'cellchat_netVisual_heatmap_1003_2024.pdf'), width=14, height=8)
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

pdf(paste0(fig_dir, compare, '_', 'cellchat_netVisual_heatmap_1003_2024_colorUpdates.pdf'), width=14, height=8)
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

gg1 <- netVisual_heatmap(cellchat, sources.use = c(1:11)) # , targets.use = c(1:11, 15:17, 19:22)
gg2 <- netVisual_heatmap(cellchat, sources.use = c(1:11), measure = "weight")

pdf(paste0(fig_dir, compare, '_', 'cellchat_netVisual_heatmap_wsourcesMHbLHb_1003_2024.pdf'), width=14, height=8)
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

gg1 <- netVisual_heatmap(cellchat, sources.use = c(15:17, 19:22))
gg2 <- netVisual_heatmap(cellchat, sources.use = c(15:17, 19:22), measure = "weight")

pdf(paste0(fig_dir, compare, '_', 'cellchat_netVisual_heatmap_wsourcesGlia_1003_2024.pdf'), width=14, height=8)
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()


# (C) Circle plot showing the number of interactions or interaction strength among different cell populations across multiple datasets
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))

pdf(paste0(fig_dir, compare, '_', 'cellchat_netVisual_circle_1003_2024.pdf'), width=30, height=26)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()


# #=================
# (A) Identify cell populations with significant changes in sending or receiving signals

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf(paste0(fig_dir, compare, '_', 'cellchat_netAnalysis_signalingRole_scatter_1003_2024.pdf'), width=8, height=6)
patchwork::wrap_plots(plots = gg)
dev.off()

# #=================
# (B) Identify the signaling changes of specific cell populations

# NOTE  ODC1
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "ODC1", comparison = c(2, 1),
  color.use = c("grey10", "#00BFC4", "#F8766D")) # , signaling.exclude = "MIF"
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "ODC1", comparison = c(2, 1),
  color.use = c("grey10", "#00BFC4", "#F8766D")) # , signaling.exclude = c("MIF")

pdf(paste0(fig_dir, compare, '_ODC1_', 'cellchat_netAnalysis_signalingChanges_scatter_1003_2024.pdf'), width=18, height=6)
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()



# Compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", color.use = c("darkorchid3", "seagreen"), sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", color.use = c("darkorchid3", "seagreen"), sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

pdf(paste0(fig_dir, compare, '_', 'cellchat_rankNet_1003_2024.pdf'), width=8, height=8)
gg1 + gg2
dev.off()



library(ComplexHeatmap)
# (B) Compare outgoing (or incoming) signaling patterns associated with each cell population

i = 1
# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 16)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 16)


pdf(paste0(fig_dir,  compare, '_', 'cellchat_netAnalysis_signalingRole_heatmap_outgoing_1003_2024.pdf'), width=12, height=16)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 16, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 16, color.heatmap = "GnBu")

pdf(paste0(fig_dir,  compare, '_', 'cellchat_netAnalysis_signalingRole_heatmap_incoming_1003_2024.pdf'), width=12, height=16)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



# ####
# NOTE Saving the data
saveRDS(cellchat_list, file=paste0(data_dir, compare, '_', "merged_cellchat_list.rda"))
saveRDS(cellchat, file=paste0(data_dir, compare, '_', "merged_cellchat.rda"))
