

# conda activate SeuratV4
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(cowplot)
library(ggrepel)
library(ggpubr)
library(ggrastr)
library(patchwork)
library(MetBrewer)
theme_set(theme_cowplot())

# source('/dfs7/swaruplab/smorabit/analysis/scWGCNA/bin/spatial_functions.R')

# set random seed for reproducibility
set.seed(12345)

setwd('/dfs7/swaruplab/shared_lab/Collaborations/Wood/Reinstatement_Project/reinstatement_2022/')
fig_dir <- 'MELD/'
data_dir <- 'data/'

# load seurat object
seurat_obj <- readRDS(paste0(data_dir, 'Reinstatement_2022_seurat.rds'))

head(seurat_obj@meta.data)

#
# # subset seurat obj to just Nurr2c & GFP
# seurat_full <- seurat_obj
# seurat_obj <- subset(seurat_full, Group %in% c('Nurr2c', 'GFP'))

# add MELD score to the seurat obj
# meld_scores <- read.csv(paste0(fig_dir, 'MELD_Reinstatement_likelihood_behavior.csv'))
Update_Path <- '/dfs7/swaruplab/zechuas/Collaborations/Wood/reinstatement_2022/Analysis/MELD/'
meld_scores <- read.csv(paste0(Update_Path, 'MELD_Reinstatement_likelihood_behavior_Update.csv'))

head(meld_scores)

#########
# check stuff
#########
identical(seurat_obj@meta.data$barcode.1, rownames(seurat_obj@meta.data))

identical(meld_scores$barcode, rownames(seurat_obj@meta.data))

identical(meld_scores$barcode, meld_scores$barcode.1)


identical(seurat_obj@meta.data$barcode.1, meld_scores$barcode.1)
# [1] TRUE

# Find the first mismatch
mismatch_index <- which(seurat_obj@meta.data$barcode.1 != rownames(seurat_obj@meta.data))[1]
# Print the mismatch
cat("Mismatch at index:", mismatch_index, "\n")
cat("barcode.1:", seurat_obj@meta.data$barcode.1[mismatch_index], "\n")
cat("rownames:", rownames(seurat_obj@meta.data)[mismatch_index], "\n")

# > cat("Mismatch at index:", mismatch_index, "\n")
# Mismatch at index: 31133
# > cat("barcode.1:", seurat_obj@meta.data$barcode.1[mismatch_index], "\n")
# barcode.1: TCATCATCAATGAAAC-1-5
# > cat("rownames:", rownames(seurat_obj@meta.data)[mismatch_index], "\n")
# rownames: TCATCATCAATGAAAC-1-5-1


head(meld_scores[meld_scores$barcode.1 == "TCATCATCAATGAAAC-1-5", ])


head(seurat_obj@meta.data[14272, ])
head(seurat_obj@meta.data[31133, ])



# Check if all values in barcode.1 are unique
is_unique <- length(unique(seurat_obj@meta.data$barcode.1)) == length(seurat_obj@meta.data$barcode.1)

# Print the result
if (is_unique) {
  cat("The column seurat_obj@meta.data$barcode.1 has unique values.\n")
} else {
  cat("The column seurat_obj@meta.data$barcode.1 does not have unique values.\n")
}

################################################################################
# add MELD likelihood
################################################################################


seurat_obj@meta.data['Reinstatement_likelihood'] <- meld_scores$Reinstatement_likelihood

quantile(seurat_obj@meta.data$Reinstatement_likelihood)
quantile(seurat_obj@meta.data$Reinstatement_likelihood[seurat_obj@meta.data$annotation == "MHb3"])

################################################################################
# MELD Score vlnplot sorted by median value
################################################################################

plot_df <- seurat_obj@meta.data %>% select(c(annotation, Reinstatement_likelihood, Group, cell_type))

table(plot_df$cell_type)


color2 <- 'seagreen'; color1 <- 'darkorchid3'
# color2 <- 'darkgoldenrod3'; color1 <- 'hotpink3'


# order groups by Neuronal vs non-neuronal, and by desceinding likelihood:
neuron_order <- plot_df %>% subset(cell_type %in% c('MHb', 'LHb', 'PHb')) %>% # 'MHb-Neuron', 'LHb-Neuron', 'PHb-Neuron'
  mutate(annotation = droplevels(annotation)) %>%
  mutate(annotation = forcats::fct_reorder(annotation, Reinstatement_likelihood, .desc=TRUE)) %>%
  .$annotation %>% levels

nonneuron_order <- plot_df %>% subset(!(cell_type %in% c('MHb', 'LHb', 'PHb'))) %>% # 'MHb-Neuron', 'LHb-Neuron', 'PHb-Neuron'
  mutate(annotation = droplevels(annotation)) %>%
  mutate(annotation = forcats::fct_reorder(annotation, Reinstatement_likelihood, .desc=TRUE)) %>%
  .$annotation %>% levels

plot_df$annotation <- factor(as.character(plot_df$annotation), levels=c(neuron_order, nonneuron_order))

p <- ggplot(plot_df, aes(x = annotation, y = Reinstatement_likelihood)) +
  ggrastr::rasterise(
    geom_jitter(
      shape=16, position=position_jitter(0.2), aes(color=Reinstatement_likelihood), alpha=0.8), dpi=400, scale=0.3) +
  scale_color_gradient2(low=color1, mid='lightgrey', high=color2, midpoint=0.5) +
  geom_boxplot(outlier.colour="black", outlier.shape=NA, notch=TRUE, fill='grey', alpha=0.25) +
  geom_hline(yintercept = 0.5, linetype='dashed', color='darkgrey') +
  RotatedAxis() + xlab('') + ylab('MELD Reinstatement Likelihood') +
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(size=1, fill=NA, color='black'),
) +scale_y_continuous(expand = c(0, 0), limits = c(0, NA))


# pdf(paste0(fig_dir, 'Reinstatement_likelihood_boxplot.pdf'), width=8, height=2.5)
# pdf(paste0(fig_dir, 'Reinstatement_likelihood_boxplot.pdf'), width=8, height=5)
pdf(paste0(fig_dir, 'Reinstatement_likelihood_boxplot_Update.pdf'), width=8, height=5)
p
dev.off()

# p <- ggplot(plot_df, aes(y = fct_rev(annotation), x = 1 - Reinstatement_likelihood)) +
p <- ggplot(plot_df, aes(y = fct_rev(annotation), x = Reinstatement_likelihood)) +
  ggrastr::rasterise(
    geom_jitter(
      shape=16, aes(color=Reinstatement_likelihood), alpha=0.8, width=0, height=0.2), dpi=400, scale=0.3) +
  scale_color_gradient2(low=color1, mid='lightgrey', high=color2, midpoint=0.5) +
  geom_boxplot(outlier.colour="black", outlier.shape=NA, notch=TRUE, fill='grey', alpha=0.25) +
  geom_vline(xintercept = 0.5, linetype='dashed', color='darkgrey') +
  ylab('') + xlab('MELD Reinstatement Likelihood') +
  RotatedAxis() +
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(size=1, fill=NA, color='black'),
) +scale_x_continuous(trans = "reverse", breaks = c(0, 0.25, 0.5, 0.75, 1)) # c(0,0.25,0.5,0.75,1) # reverse the order


#
# p <- ggplot(plot_df, aes(y = fct_rev(annotation), x = 1 - Reinstatement_likelihood)) +
#   ggrastr::rasterise(
#     geom_jitter(
#       shape=16, aes(color=Reinstatement_likelihood), alpha=0.8, width=0, height=0.2), dpi=400, scale=0.3) +
#   scale_color_gradient2(low=color1, mid='lightgrey', high=color2, midpoint=0.5) +
#   geom_boxplot(outlier.colour="black", outlier.shape=NA, notch=TRUE, fill='grey', alpha=0.25) +
#   geom_vline(xintercept = 0.5, linetype='dashed', color='darkgrey') +
#   ylab('') + xlab('MELD Reinstatement Likelihood') +
#   RotatedAxis() +
#   theme(
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank(),
#     panel.border = element_rect(size=1, fill=NA, color='black'),
# ) +scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1))



# pdf(paste0(fig_dir, 'Reinstatement_likelihood_boxplot_tall.pdf'), width=4, height=6)
# pdf(paste0(fig_dir, 'Reinstatement_likelihood_boxplot_tall.pdf'), width=5, height=6)
pdf(paste0(fig_dir, 'Reinstatement_likelihood_boxplot_tall_Update.pdf'), width=5, height=6)
p
dev.off()
