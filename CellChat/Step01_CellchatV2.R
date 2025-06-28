#================================== NOTE ==================================
# NOTE # CellChat
#================================== NOTE ==================================

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


# Set working directory
setwd('/dfs7/swaruplab/shared_lab/Collaborations/Wood/Reinstatement_Project/reinstatement_2022/CellChat/')

# List of cell types and annotations to loop over
identifiers <- c("cell_type", "annotation")  # Loop through these identifiers in Seurat object

for (id in identifiers) {
  cat("Processing identifier:", id, "\n")
  start_id <- Sys.time()

  # Set paths for output data and figures
  data_dir <- file.path(id, 'data')
  fig_dir <- file.path(id, 'figures')
  data_in <- file.path("/dfs7/swaruplab/shared_lab/Collaborations/Wood/Reinstatement_Project/reinstatement_2022/data/")

  # Create directories if they don't exist
  dir.create(data_dir, recursive = TRUE)
  dir.create(fig_dir, recursive = TRUE)

  # Load Seurat object (same for all identifiers)
  cat("Loading Seurat object\n")
  start_load <- Sys.time()
  NucSeq <- readRDS(file.path(data_in, "Reinstatement_2022_seurat.rds"))
  seurat_obj <- NucSeq

  cat("Time to load Seurat object:", as.numeric(Sys.time() - start_load, units = "mins"), "minutes\n")

  # Set Idents dynamically based on the identifier (cell_type or annotation)
  Idents(seurat_obj) <- id
  data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  labels <- Idents(seurat_obj)
  meta <- data.frame(group = labels, row.names = names(labels))

  # Check if rownames in seurat_obj@meta.data match those in meta
  cat("Checking if rownames in seurat_obj@meta.data match those in meta:\n")
  rownames_match <- all(rownames(seurat_obj@meta.data) == rownames(meta))
  cat("Rownames match result:", rownames_match, "\n")

  if (!rownames_match) {
    stop("Rownames are NOT aligned. Please check the data!")
  }

  seurat_obj$samples <- seurat_obj@meta.data$Sample
  conditions <- unique(seurat_obj$Group)

  # Load the mouse version of CellChatDB
  CellChatDB <- CellChatDB.mouse
  cellchat_list <- list()

  # Loop through conditions (AD and Control)
  for (cond in conditions) {
    cat("Creating CellChat object for condition:", cond, "and identifier:", id, "\n")
    start_cond <- Sys.time()

    # Create CellChat object
    cellchat_list[[cond]] <- createCellChat(
      object = data.input[, seurat_obj$Group == cond],
      meta = seurat_obj@meta.data %>% subset(Group == cond),
      group.by = id
    )
    cellchat_list[[cond]]@DB <- CellChatDB

    cat("Time to create CellChat object for condition", cond, ":", as.numeric(Sys.time() - start_cond, units = "mins"), "minutes\n")
  }

  # Process each condition for prefiltering, filtering, and further analysis
  for (cond in conditions) {
    cat("Processing condition:", cond, "for identifier:", id, "\n")

    # Prefiltering steps
    start_prefilter <- Sys.time()
    cat("Prefiltering for condition:", cond, "\n")

    cellchat_list[[cond]] <- subsetData(cellchat_list[[cond]])
    cellchat_list[[cond]] <- identifyOverExpressedGenes(cellchat_list[[cond]])
    cellchat_list[[cond]] <- identifyOverExpressedInteractions(cellchat_list[[cond]])
    cellchat_list[[cond]] <- computeCommunProb(cellchat_list[[cond]], type = "triMean")

    saveRDS(cellchat_list[[cond]], file = file.path(data_dir, paste0(gsub(' ', '_', cond), '_', id, '_cellchat_prefiltering_Reinstatement_0930_2024.rds')))

    cat("Time for prefiltering condition", cond, ":", as.numeric(Sys.time() - start_prefilter, units = "mins"), "minutes\n")

    # Filtering and further processing
    start_filter <- Sys.time()
    cat("Filtering and processing condition:", cond, "\n")

    cellchat_list[[cond]] <- filterCommunication(cellchat_list[[cond]], min.cells = 50)
    df.net <- subsetCommunication(cellchat_list[[cond]])
    cellchat_list[[cond]] <- computeCommunProbPathway(cellchat_list[[cond]])
    cellchat_list[[cond]] <- aggregateNet(cellchat_list[[cond]])
    cellchat_list[[cond]] <- netAnalysis_computeCentrality(cellchat_list[[cond]], slot.name = "netP")

    saveRDS(cellchat_list[[cond]], file = file.path(data_dir, paste0(gsub(' ', '_', cond), '_', id, '_cellchat_filtered_n_processed_Reinstatement_0930_2024.rds')))

    cat("Time for filtering and processing condition", cond, ":", as.numeric(Sys.time() - start_filter, units = "mins"), "minutes\n")
  }

  cat("Finished processing identifier:", id, "\n")
  cat("Total time for identifier", id, ":", as.numeric(Sys.time() - start_id, units = "mins"), "minutes\n\n")

  # Clean up and free memory
  rm(NucSeq, seurat_obj, cellchat_list)
  gc()  # Garbage collection
}
