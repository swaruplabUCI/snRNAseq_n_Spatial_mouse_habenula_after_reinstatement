#Installation of high dimensional WGCNA
# create new conda environment for R
conda create -n hdWGCNA -c conda-forge r-base r-essentials

# activate conda environment
conda activate hdWGCNA


devtools::install_github('smorabit/hdWGCNA', ref='dev')

BiocManager::install("MAST")
BiocManager::install("conflicted")
BiocManager::install("qlcMatrix")
BiocManager::install("UCell")
BiocManager::install("gganimate")
BiocManager::install("enrichR")
BiocManager::install("GeneOverlap")
BiocManager::install("dittoSeq")
BiocManager::install("scCustomize")
BiocManager::install("plotly")
BiocManager::install("Matrix")
BiocManager::install("WGCNA")
BiocManager::install("Seurat")
BiocManager::install("UCell")
BiocManager::install("MetBrewer")
BiocManager::install("plotly")


# single-cell analysis package
library(Seurat)
library(SeuratObject)
# plotting and data science packages
library(tidyverse)
library(Matrix)
library(MAST)
library(devtools)
library(cowplot)
library(patchwork)
library(igraph)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
library(qlcMatrix)
library(UCell)
library(gganimate)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())
# gene enrichment packages
library(enrichR)
library(GeneOverlap)
library(dittoSeq)
library(scCustomize)
library(qs)
library(RColorBrewer)
library(MetBrewer)
library(viridis)
library(grid)
library(gridExtra)
library(EnhancedVolcano)




# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# load the our snRNA-seq dataset
seurat_obj <- readRDS('Reinstatement_2022_seurat.rds')
head(seurat_obj)

p <- DimPlot(seurat_obj, group.by='cell_type', label=TRUE) +
   umap_theme() + ggtitle('UMAP') + NoLegend()

png('hdWGCNA_figures/DimPlot.png', width=4, height=4, res=200, units='in')
print(p)
dev.off()

seurat_obj <- seurat_obj[,!seurat_obj$cell_type %in% c('PER', 'END', 'EPD')]
head(seurat_obj)

#ODC celltype for hdWGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  group.by = 'cell_type',
  wgcna_name = "ODC"
)
length(GetWGCNAGenes(seurat_obj))

test <- SelectNetworkGenes(
    subset(seurat_obj, cell_type == 'ODC'),
    gene_select = 'fraction',
    fraction = 0.05
)
length(GetWGCNAGenes(test))
cur_genes <- GetWGCNAGenes(test)
seurat_obj <- SetWGCNAGenes(seurat_obj, cur_genes)

test <- MetacellsByGroups(
  seurat_obj = test,
  group.by = c("cell_type","Sample", "Group"),
  k = 25, # k=25 # this was used on the first try
  max_shared=15,
  ident.group = 'cell_type',
  reduction='harmony',
  target_metacells=100,
  min_cells = 25
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type","Sample", "Group"),
  k = 25, # k=25 # this was used on the first try
  max_shared=15,
  ident.group = 'cell_type',
  reduction='harmony',
  target_metacells=100,
  min_cells = 50
)


seurat_obj <- NormalizeMetacells(seurat_obj)


seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "ODC",
  group.by = "cell_type"
)

# test the soft power threshold
seurat_obj <- TestSoftPowers(seurat_obj)

plot_list <- PlotSoftPowers(seurat_obj)

#assemble with patchwork
pdf(paste0('hdWGCNA_figures/softpower_odc.pdf'), width=12, height=8)
wrap_plots(plot_list, ncol=2)
dev.off()


seurat_obj <- ConstructNetwork(
    seurat_obj, 
    soft_power=6,
    tom_name='ODC', 
    overwrite_tom=TRUE
)

# plot the dendrogram
pdf(paste0("hdWGCNA_figures/dendro_odc.pdf"),height=3, width=6)
PlotDendrogram(seurat_obj, main='ODC hdWGCNA Dendrogram')
dev.off()

modules <- GetModules(seurat_obj)
table(modules$module)
head(modules)


seurat_obj <- ScaleData(seurat_obj, features=rownames(seurat_obj)[1:100])
seurat_obj$Assignment <- droplevels(seurat_obj$Assignment)
seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars="Assignment" # snRNAseq batch
)

# compute module connectivity:
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type', group_name = 'ODC'
)

########################################
# Feature Plot
########################################


plot_list <- ModuleFeaturePlot(
  seurat_obj, 
  order=TRUE, 
  raster=TRUE, 
  restrict_range=FALSE, 
  raster_dpi=400, 
  raster_scale=0.5, 
  point_size=1
  )

plot_list <- lapply(plot_list, function(x){
  x + NoLegend() + theme(plot.margin=margin(0,0,0,0))
})

pdf(paste0("hdWGCNA_figures/ODC_featureplot_MEs.pdf"),height=6, width=12)
wrap_plots(plot_list, ncol=5)
dev.off()

########################################
# Dot Plot
########################################


MEs <- GetMEs(seurat_obj)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module)

write.csv(modules , file='ODC_modules.csv')
mods <- mods[mods!='grey']
mods

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)


# make dotplot
p <- DotPlot(
  seurat_obj,
  group.by='annotation',
  features = rev(mods)
#  scale=FALSE, col.max=5
) + coord_flip() + RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') + xlab('') + ylab('') +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

pdf(paste0('hdWGCNA_figures/ODC_MEs_dotplot.pdf'), width=9, height=4)
p
dev.off()

# remove from seurat obj
seurat_obj@meta.data <- seurat_obj@meta.data[,!(colnames(seurat_obj@meta.data) %in% colnames(MEs))]


################################################################################
# Hubgene circle plots:
################################################################################

library(igraph)

# individual module networks
ModuleNetworkPlot(
  seurat_obj,
  mods = "all",
  #label_center=TRUE,
  outdir = paste0('hdWGCNA_figures/ODC_hubNetworks/')
)

################################################################################
# UMAP:
################################################################################


seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 3,
  n_neighbors=10,
  min_dist=0.3,
  spread=3,
  supervised=TRUE,
  target_weight=0.5
)


# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

write.csv(umap_df , file='ODC_modules_UMAP.csv')
# plot with ggplot
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
   color=umap_df$color,
   size=umap_df$kME*2
  ) +
  umap_theme()

pdf(paste0('hdWGCNA_figures/ODC_hubgene_umap_ggplot2.pdf'), width=5, height=5)
p
dev.off()


library(reshape2)
library(igraph)
pdf(paste0('hdWGCNA_figures/ODC_hubgene_umap_igraph.pdf'), width=10, height=10)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.5,
  sample_edges=TRUE,
  keep_grey_edges=FALSE,
  edge_prop=0.075, # taking the top 20% strongest edges in each module
  label_hubs=4, # how many hub genes to plot per module?
)
dev.off()

saveRDS(seurat_obj, file = "ODC_hdWGCNA.rds")

######## Enrichment analysis ##################

dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')
head(dbs)
# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)
head(enrich_df)
# make GO term plots:
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(9,12), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

# enrichr dotplot
pdf(paste0('hdWGCNA_figures/ODC_BP.pdf'), width=8, height=6)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

pdf(paste0('hdWGCNA_figures/ODC_CC.pdf'), width=8, height=6)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Cellular_Component_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

pdf(paste0('hdWGCNA_figures/ODC_MF.pdf'), width=8, height=6)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Molecular_Function_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

#Marker gene overlap analysis
# compute cell-type marker genes with Seurat:
Idents(seurat_obj) <- seurat_obj$cell_type
markers <- Seurat::FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  logfc.threshold=1
)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)

# overlap barplot, produces a plot for each cell type
plot_list <- OverlapBarPlot(overlap_df)

# stitch plots with patchwork
pdf(paste0('hdWGCNA_figures/ODC_Overlap_barplot.pdf'), width=5, height=6)
wrap_plots(plot_list, ncol=3)
dev.off()

# plot odds ratio of the overlap as a dot plot
pdf(paste0('hdWGCNA_figures/ODC_Overlap_dotplot.pdf'), width=5, height=4)
OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cell-type markers')
dev.off()



############## DME analysis comparing two groups ###############
library(ggforestplot)
table(seurat_obj@meta.data$Group)

group1 <- seurat_obj@meta.data %>% subset(cell_type == 'ODC' & Group == "Reinstatement") %>% rownames
group2 <- seurat_obj@meta.data %>% subset(cell_type == 'ODC' & Group == "NR") %>% rownames
head(group1)
head(group2)

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='ODC'
)
head(DMEs)

pdf(paste0('hdWGCNA_figures/ODC_PlotDMEsLollipop.pdf'), width=4, height=5)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name='ODC', 
  pvalue = "p_val_adj"
)
dev.off()

pdf(paste0('hdWGCNA_figures/ODC_PlotDMEsVolcanop.pdf'), width=4, height=4)
PlotDMEsVolcano(
  seurat_obj,
  DMEs,
  wgcna_name = 'ODC'
)
dev.off()

#One-versus-all DME analysis
group.by = 'cell_type'

DMEs_all <- FindAllDMEs(
  seurat_obj,
  group.by = 'cell_type',
  wgcna_name = 'ODC'
)

head(DMEs_all)

p <- PlotDMEsVolcano(
  seurat_obj,
  DMEs_all,
  wgcna_name = 'ODC',
  plot_labels=FALSE,
  show_cutoff=FALSE
)

# facet wrap by each cell type
pdf(paste0('hdWGCNA_figures/ODC_PlotDMEsVolcano_One-versus-all.pdf'), width=6, height=7)
p + facet_wrap(~group, ncol=3)
dev.off()


##############################################################################################
#Astrocytes (ASC)
# load the our snRNA-seq dataset
seurat_obj <- readRDS('Reinstatement_2022_seurat.rds')
head(seurat_obj)

p <- DimPlot(seurat_obj, group.by='cell_type', label=TRUE) +
   umap_theme() + ggtitle('UMAP') + NoLegend()

png('hdWGCNA_figures/DimPlot.png', width=4, height=4, res=200, units='in')
print(p)
dev.off()

seurat_obj <- seurat_obj[,!seurat_obj$cell_type %in% c('PER', 'END', 'EPD')]
head(seurat_obj)

#ASC celltype for hdWGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  group.by = 'cell_type',
  wgcna_name = "ASC"
)
length(GetWGCNAGenes(seurat_obj))

test <- SelectNetworkGenes(
    subset(seurat_obj, cell_type == 'ASC'),
    gene_select = 'fraction',
    fraction = 0.05
)
length(GetWGCNAGenes(test))
cur_genes <- GetWGCNAGenes(test)
seurat_obj <- SetWGCNAGenes(seurat_obj, cur_genes)

test <- MetacellsByGroups(
  seurat_obj = test,
  group.by = c("cell_type","Sample", "Group"),
  k = 25, # k=25 # this was used on the first try
  max_shared=15,
  ident.group = 'cell_type',
  reduction='harmony',
  target_metacells=100,
  min_cells = 25
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type","Sample", "Group"),
  k = 25, # k=25 # this was used on the first try
  max_shared=15,
  ident.group = 'cell_type',
  reduction='harmony',
  target_metacells=100,
  min_cells = 50
)


seurat_obj <- NormalizeMetacells(seurat_obj)


seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "ASC",
  group.by = "cell_type"
)

# test the soft power threshold
seurat_obj <- TestSoftPowers(seurat_obj)

plot_list <- PlotSoftPowers(seurat_obj)

#assemble with patchwork
pdf(paste0('hdWGCNA_figures/softpower_asc.pdf'), width=12, height=8)
wrap_plots(plot_list, ncol=2)
dev.off()


seurat_obj <- ConstructNetwork(
    seurat_obj, 
    soft_power=6,
    tom_name='ASC', 
    overwrite_tom=TRUE
)

# plot the dendrogram
pdf(paste0("hdWGCNA_figures/dendro_asc.pdf"),height=3, width=6)
PlotDendrogram(seurat_obj, main='ASC hdWGCNA Dendrogram')
dev.off()

modules <- GetModules(seurat_obj)
table(modules$module)
head(modules)


seurat_obj <- ScaleData(seurat_obj, features=rownames(seurat_obj)[1:100])
seurat_obj$Assignment <- droplevels(seurat_obj$Assignment)
seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars="Assignment" # snRNAseq batch
)

# compute module connectivity:
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type', group_name = 'ASC'
)

########################################
# Feature Plot
########################################


plot_list <- ModuleFeaturePlot(
  seurat_obj, 
  order=TRUE, 
  raster=TRUE, 
  restrict_range=FALSE, 
  raster_dpi=400, 
  raster_scale=0.5, 
  point_size=1
  )

plot_list <- lapply(plot_list, function(x){
  x + NoLegend() + theme(plot.margin=margin(0,0,0,0))
})

pdf(paste0("hdWGCNA_figures/ASC_featureplot_MEs.pdf"),height=6, width=12)
wrap_plots(plot_list, ncol=5)
dev.off()

########################################
# Dot Plot
########################################


MEs <- GetMEs(seurat_obj)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module)

write.csv(modules , file='ASC_modules.csv')
mods <- mods[mods!='grey']
mods

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)


# make dotplot
p <- DotPlot(
  seurat_obj,
  group.by='annotation',
  features = rev(mods)
#  scale=FALSE, col.max=5
) + coord_flip() + RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') + xlab('') + ylab('') +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

pdf(paste0('hdWGCNA_figures/ASC_MEs_dotplot.pdf'), width=9, height=4)
p
dev.off()

# remove from seurat obj
seurat_obj@meta.data <- seurat_obj@meta.data[,!(colnames(seurat_obj@meta.data) %in% colnames(MEs))]


################################################################################
# Hubgene circle plots:
################################################################################

library(igraph)

# individual module networks
ModuleNetworkPlot(
  seurat_obj,
  mods = "all",
  #label_center=TRUE,
  outdir = paste0('hdWGCNA_figures/ASC_hubNetworks/')
)

################################################################################
# UMAP:
################################################################################


seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 3,
  n_neighbors=10,
  min_dist=0.3,
  spread=3,
  supervised=TRUE,
  target_weight=0.5
)


# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

write.csv(umap_df , file='ASC_modules_UMAP.csv')
# plot with ggplot
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
   color=umap_df$color,
   size=umap_df$kME*2
  ) +
  umap_theme()

pdf(paste0('hdWGCNA_figures/ASC_hubgene_umap_ggplot2.pdf'), width=5, height=5)
p
dev.off()


library(reshape2)
library(igraph)
pdf(paste0('hdWGCNA_figures/ASC_hubgene_umap_igraph.pdf'), width=10, height=10)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.5,
  sample_edges=TRUE,
  keep_grey_edges=FALSE,
  edge_prop=0.075, # taking the top 20% strongest edges in each module
  label_hubs=4, # how many hub genes to plot per module?
)
dev.off()

saveRDS(seurat_obj, file = "ASC_hdWGCNA.rds")

######## Enrichment analysis ##################

dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')
head(dbs)
# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)
head(enrich_df)
# make GO term plots:
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(9,12), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

# enrichr dotplot
pdf(paste0('hdWGCNA_figures/ASC_BP.pdf'), width=8, height=6)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

pdf(paste0('hdWGCNA_figures/ASC_CC.pdf'), width=8, height=6)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Cellular_Component_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

pdf(paste0('hdWGCNA_figures/ASC_MF.pdf'), width=8, height=6)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Molecular_Function_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

#Marker gene overlap analysis
# compute cell-type marker genes with Seurat:
Idents(seurat_obj) <- seurat_obj$cell_type
markers <- Seurat::FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  logfc.threshold=1
)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)

# overlap barplot, produces a plot for each cell type
plot_list <- OverlapBarPlot(overlap_df)

# stitch plots with patchwork
pdf(paste0('hdWGCNA_figures/ASC_Overlap_barplot.pdf'), width=5, height=6)
wrap_plots(plot_list, ncol=3)
dev.off()

# plot odds ratio of the overlap as a dot plot
pdf(paste0('hdWGCNA_figures/ASC_Overlap_dotplot.pdf'), width=5, height=4)
OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cell-type markers')
dev.off()


############## DME analysis comparing two groups ###############
library(ggforestplot)
table(seurat_obj@meta.data$Group)

group1 <- seurat_obj@meta.data %>% subset(cell_type == 'ASC' & Group == "Reinstatement") %>% rownames
group2 <- seurat_obj@meta.data %>% subset(cell_type == 'ASC' & Group == "NR") %>% rownames
head(group1)
head(group2)

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='ASC'
)
head(DMEs)

pdf(paste0('hdWGCNA_figures/ASC_PlotDMEsLollipop.pdf'), width=4, height=5)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name='ASC', 
  pvalue = "p_val_adj"
)
dev.off()

pdf(paste0('hdWGCNA_figures/ASC_PlotDMEsVolcanop.pdf'), width=4, height=4)
PlotDMEsVolcano(
  seurat_obj,
  DMEs,
  wgcna_name = 'ASC'
)
dev.off()

#One-versus-all DME analysis
group.by = 'cell_type'

DMEs_all <- FindAllDMEs(
  seurat_obj,
  group.by = 'cell_type',
  wgcna_name = 'ASC'
)

head(DMEs_all)

p <- PlotDMEsVolcano(
  seurat_obj,
  DMEs_all,
  wgcna_name = 'ASC',
  plot_labels=FALSE,
  show_cutoff=FALSE
)

# facet wrap by each cell type
pdf(paste0('hdWGCNA_figures/ASC_PlotDMEsVolcano_One-versus-all.pdf'), width=6, height=7)
p + facet_wrap(~group, ncol=3)
dev.off()


##############################################################################################
# MHb neuron(MHb)
# load the our snRNA-seq dataset
seurat_obj <- readRDS('Reinstatement_2022_seurat.rds')
head(seurat_obj)

p <- DimPlot(seurat_obj, group.by='cell_type', label=TRUE) +
   umap_theme() + ggtitle('UMAP') + NoLegend()

png('hdWGCNA_figures/DimPlot.png', width=4, height=4, res=200, units='in')
print(p)
dev.off()

seurat_obj <- seurat_obj[,!seurat_obj$cell_type %in% c('PER', 'END', 'EPD')]
head(seurat_obj)

#MHb celltype for hdWGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  group.by = 'cell_type',
  wgcna_name = "MHb"
)
length(GetWGCNAGenes(seurat_obj))

test <- SelectNetworkGenes(
    subset(seurat_obj, cell_type == 'MHb'),
    gene_select = 'fraction',
    fraction = 0.05
)
length(GetWGCNAGenes(test))
cur_genes <- GetWGCNAGenes(test)
seurat_obj <- SetWGCNAGenes(seurat_obj, cur_genes)

test <- MetacellsByGroups(
  seurat_obj = test,
  group.by = c("cell_type","Sample", "Group"),
  k = 25, # k=25 # this was used on the first try
  max_shared=15,
  ident.group = 'cell_type',
  reduction='harmony',
  target_metacells=100,
  min_cells = 25
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type","Sample", "Group"),
  k = 25, # k=25 # this was used on the first try
  max_shared=15,
  ident.group = 'cell_type',
  reduction='harmony',
  target_metacells=100,
  min_cells = 50
)


seurat_obj <- NormalizeMetacells(seurat_obj)


seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "MHb",
  group.by = "cell_type"
)

# test the soft power threshold
seurat_obj <- TestSoftPowers(seurat_obj)

plot_list <- PlotSoftPowers(seurat_obj)

#assemble with patchwork
pdf(paste0('hdWGCNA_figures/softpower_MHb.pdf'), width=12, height=8)
wrap_plots(plot_list, ncol=2)
dev.off()


seurat_obj <- ConstructNetwork(
    seurat_obj, 
    soft_power=6,
    tom_name='MHb', 
    overwrite_tom=TRUE
)

# plot the dendrogram
pdf(paste0("hdWGCNA_figures/dendro_MHb.pdf"),height=3, width=6)
PlotDendrogram(seurat_obj, main='MHb hdWGCNA Dendrogram')
dev.off()

modules <- GetModules(seurat_obj)
table(modules$module)
head(modules)


seurat_obj <- ScaleData(seurat_obj, features=rownames(seurat_obj)[1:100])
seurat_obj$Assignment <- droplevels(seurat_obj$Assignment)
seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars="Assignment" # snRNAseq batch
)

# compute module connectivity:
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type', group_name = 'MHb'
)

########################################
# Feature Plot
########################################


plot_list <- ModuleFeaturePlot(
  seurat_obj, 
  order=TRUE, 
  raster=TRUE, 
  restrict_range=FALSE, 
  raster_dpi=400, 
  raster_scale=0.5, 
  point_size=1
  )

plot_list <- lapply(plot_list, function(x){
  x + NoLegend() + theme(plot.margin=margin(0,0,0,0))
})

pdf(paste0("hdWGCNA_figures/MHb_featureplot_MEs.pdf"),height=6, width=12)
wrap_plots(plot_list, ncol=5)
dev.off()

########################################
# Dot Plot
########################################


MEs <- GetMEs(seurat_obj)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module)

write.csv(modules , file='MHb_modules.csv')
mods <- mods[mods!='grey']
mods

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)


# make dotplot
p <- DotPlot(
  seurat_obj,
  group.by='annotation',
  features = rev(mods)
#  scale=FALSE, col.max=5
) + coord_flip() + RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') + xlab('') + ylab('') +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

pdf(paste0('hdWGCNA_figures/MHb_MEs_dotplot.pdf'), width=9, height=4)
p
dev.off()

# remove from seurat obj
seurat_obj@meta.data <- seurat_obj@meta.data[,!(colnames(seurat_obj@meta.data) %in% colnames(MEs))]


################################################################################
# Hubgene circle plots:
################################################################################

library(igraph)

# individual module networks
ModuleNetworkPlot(
  seurat_obj,
  mods = "all",
  #label_center=TRUE,
  outdir = paste0('hdWGCNA_figures/MHb_hubNetworks/')
)

################################################################################
# UMAP:
################################################################################


seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 3,
  n_neighbors=10,
  min_dist=0.3,
  spread=3,
  supervised=TRUE,
  target_weight=0.5
)


# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

write.csv(umap_df , file='MHb_modules_UMAP.csv')
# plot with ggplot
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
   color=umap_df$color,
   size=umap_df$kME*2
  ) +
  umap_theme()

pdf(paste0('hdWGCNA_figures/MHb_hubgene_umap_ggplot2.pdf'), width=5, height=5)
p
dev.off()


library(reshape2)
library(igraph)
pdf(paste0('hdWGCNA_figures/MHb_hubgene_umap_igraph.pdf'), width=10, height=10)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.5,
  sample_edges=TRUE,
  keep_grey_edges=FALSE,
  edge_prop=0.075, # taking the top 20% strongest edges in each module
  label_hubs=4, # how many hub genes to plot per module?
)
dev.off()

saveRDS(seurat_obj, file = "MHb_hdWGCNA.rds")

######## Enrichment analysis ##################

dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')
head(dbs)
# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)
head(enrich_df)
# make GO term plots:
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(9,12), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

# enrichr dotplot
pdf(paste0('hdWGCNA_figures/MHb_BP.pdf'), width=8, height=6)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

pdf(paste0('hdWGCNA_figures/MHb_CC.pdf'), width=8, height=6)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Cellular_Component_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

pdf(paste0('hdWGCNA_figures/MHb_MF.pdf'), width=8, height=6)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Molecular_Function_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

#Marker gene overlap analysis
# compute cell-type marker genes with Seurat:
Idents(seurat_obj) <- seurat_obj$cell_type
markers <- Seurat::FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  logfc.threshold=1
)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)

# overlap barplot, produces a plot for each cell type
plot_list <- OverlapBarPlot(overlap_df)

# stitch plots with patchwork
pdf(paste0('hdWGCNA_figures/MHb_Overlap_barplot.pdf'), width=5, height=6)
wrap_plots(plot_list, ncol=3)
dev.off()

# plot odds ratio of the overlap as a dot plot
pdf(paste0('hdWGCNA_figures/MHb_Overlap_dotplot.pdf'), width=5, height=4)
OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cell-type markers')
dev.off()


############## DME analysis comparing two groups ###############
library(ggforestplot)
table(seurat_obj@meta.data$Group)

group1 <- seurat_obj@meta.data %>% subset(cell_type == 'MHb' & Group == "Reinstatement") %>% rownames
group2 <- seurat_obj@meta.data %>% subset(cell_type == 'MHb' & Group == "NR") %>% rownames
head(group1)
head(group2)

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='MHb'
)
head(DMEs)

pdf(paste0('hdWGCNA_figures/MHb_PlotDMEsLollipop.pdf'), width=4, height=5)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name='MHb', 
  pvalue = "p_val_adj"
)
dev.off()

pdf(paste0('hdWGCNA_figures/MHb_PlotDMEsVolcanop.pdf'), width=4, height=4)
PlotDMEsVolcano(
  seurat_obj,
  DMEs,
  wgcna_name = 'MHb'
)
dev.off()

#One-versus-all DME analysis
group.by = 'cell_type'

DMEs_all <- FindAllDMEs(
  seurat_obj,
  group.by = 'cell_type',
  wgcna_name = 'MHb'
)

head(DMEs_all)

p <- PlotDMEsVolcano(
  seurat_obj,
  DMEs_all,
  wgcna_name = 'MHb',
  plot_labels=FALSE,
  show_cutoff=FALSE
)

# facet wrap by each cell type
pdf(paste0('hdWGCNA_figures/MHb_PlotDMEsVolcano_One-versus-all.pdf'), width=6, height=7)
p + facet_wrap(~group, ncol=3)
dev.off()


##############################################################################################
# LHb neuron(LHb)
# load the our snRNA-seq dataset
seurat_obj <- readRDS('Reinstatement_2022_seurat.rds')
head(seurat_obj)

p <- DimPlot(seurat_obj, group.by='cell_type', label=TRUE) +
   umap_theme() + ggtitle('UMAP') + NoLegend()

png('hdWGCNA_figures/DimPlot.png', width=4, height=4, res=200, units='in')
print(p)
dev.off()

seurat_obj <- seurat_obj[,!seurat_obj$cell_type %in% c('PER', 'END', 'EPD')]
head(seurat_obj)

#LHb celltype for hdWGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  group.by = 'cell_type',
  wgcna_name = "LHb"
)
length(GetWGCNAGenes(seurat_obj))

test <- SelectNetworkGenes(
    subset(seurat_obj, cell_type == 'LHb'),
    gene_select = 'fraction',
    fraction = 0.05
)
length(GetWGCNAGenes(test))
cur_genes <- GetWGCNAGenes(test)
seurat_obj <- SetWGCNAGenes(seurat_obj, cur_genes)

test <- MetacellsByGroups(
  seurat_obj = test,
  group.by = c("cell_type","Sample", "Group"),
  k = 25, # k=25 # this was used on the first try
  max_shared=15,
  ident.group = 'cell_type',
  reduction='harmony',
  target_metacells=100,
  min_cells = 25
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type","Sample", "Group"),
  k = 25, # k=25 # this was used on the first try
  max_shared=15,
  ident.group = 'cell_type',
  reduction='harmony',
  target_metacells=100,
  min_cells = 50
)


seurat_obj <- NormalizeMetacells(seurat_obj)


seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "LHb",
  group.by = "cell_type"
)

# test the soft power threshold
seurat_obj <- TestSoftPowers(seurat_obj)

plot_list <- PlotSoftPowers(seurat_obj)

#assemble with patchwork
pdf(paste0('hdWGCNA_figures/softpower_LHb.pdf'), width=12, height=8)
wrap_plots(plot_list, ncol=2)
dev.off()


seurat_obj <- ConstructNetwork(
    seurat_obj, 
    soft_power=6,
    tom_name='LHb', 
    overwrite_tom=TRUE
)

# plot the dendrogram
pdf(paste0("hdWGCNA_figures/dendro_LHb.pdf"),height=3, width=6)
PlotDendrogram(seurat_obj, main='LHb hdWGCNA Dendrogram')
dev.off()

modules <- GetModules(seurat_obj)
table(modules$module)
head(modules)


seurat_obj <- ScaleData(seurat_obj, features=rownames(seurat_obj)[1:100])
seurat_obj$Assignment <- droplevels(seurat_obj$Assignment)
seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars="Assignment" # snRNAseq batch
)

# compute module connectivity:
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type', group_name = 'LHb'
)

########################################
# Feature Plot
########################################


plot_list <- ModuleFeaturePlot(
  seurat_obj, 
  order=TRUE, 
  raster=TRUE, 
  restrict_range=FALSE, 
  raster_dpi=400, 
  raster_scale=0.5, 
  point_size=1
  )

plot_list <- lapply(plot_list, function(x){
  x + NoLegend() + theme(plot.margin=margin(0,0,0,0))
})

pdf(paste0("hdWGCNA_figures/LHb_featureplot_MEs.pdf"),height=6, width=12)
wrap_plots(plot_list, ncol=5)
dev.off()

########################################
# Dot Plot
########################################


MEs <- GetMEs(seurat_obj)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module)

write.csv(modules , file='LHb_modules.csv')
mods <- mods[mods!='grey']
mods

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)


# make dotplot
p <- DotPlot(
  seurat_obj,
  group.by='annotation',
  features = rev(mods)
#  scale=FALSE, col.max=5
) + coord_flip() + RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') + xlab('') + ylab('') +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

pdf(paste0('hdWGCNA_figures/LHb_MEs_dotplot.pdf'), width=9, height=4)
p
dev.off()

# remove from seurat obj
seurat_obj@meta.data <- seurat_obj@meta.data[,!(colnames(seurat_obj@meta.data) %in% colnames(MEs))]


################################################################################
# Hubgene circle plots:
################################################################################

library(igraph)

# individual module networks
ModuleNetworkPlot(
  seurat_obj,
  mods = "all",
  #label_center=TRUE,
  outdir = paste0('hdWGCNA_figures/LHb_hubNetworks/')
)

################################################################################
# UMAP:
################################################################################


seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 3,
  n_neighbors=10,
  min_dist=0.3,
  spread=3,
  supervised=TRUE,
  target_weight=0.5
)


# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

write.csv(umap_df , file='LHb_modules_UMAP.csv')
# plot with ggplot
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
   color=umap_df$color,
   size=umap_df$kME*2
  ) +
  umap_theme()

pdf(paste0('hdWGCNA_figures/LHb_hubgene_umap_ggplot2.pdf'), width=5, height=5)
p
dev.off()


library(reshape2)
library(igraph)
pdf(paste0('hdWGCNA_figures/LHb_hubgene_umap_igraph.pdf'), width=10, height=10)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.5,
  sample_edges=TRUE,
  keep_grey_edges=FALSE,
  edge_prop=0.075, # taking the top 20% strongest edges in each module
  label_hubs=4, # how many hub genes to plot per module?
)
dev.off()

saveRDS(seurat_obj, file = "LHb_hdWGCNA.rds")

######## Enrichment analysis ##################

dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')
head(dbs)
# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)
head(enrich_df)
# make GO term plots:
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(9,12), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

# enrichr dotplot
pdf(paste0('hdWGCNA_figures/LHb_BP.pdf'), width=8, height=6)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

pdf(paste0('hdWGCNA_figures/LHb_CC.pdf'), width=8, height=6)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Cellular_Component_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

pdf(paste0('hdWGCNA_figures/LHb_MF.pdf'), width=8, height=6)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Molecular_Function_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

#Marker gene overlap analysis
# compute cell-type marker genes with Seurat:
Idents(seurat_obj) <- seurat_obj$cell_type
markers <- Seurat::FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  logfc.threshold=1
)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)

# overlap barplot, produces a plot for each cell type
plot_list <- OverlapBarPlot(overlap_df)

# stitch plots with patchwork
pdf(paste0('hdWGCNA_figures/LHb_Overlap_barplot.pdf'), width=5, height=6)
wrap_plots(plot_list, ncol=3)
dev.off()

# plot odds ratio of the overlap as a dot plot
pdf(paste0('hdWGCNA_figures/LHb_Overlap_dotplot.pdf'), width=5, height=4)
OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cell-type markers')
dev.off()


############## DME analysis comparing two groups ###############
library(ggforestplot)
table(seurat_obj@meta.data$Group)

group1 <- seurat_obj@meta.data %>% subset(cell_type == 'LHb' & Group == "Reinstatement") %>% rownames
group2 <- seurat_obj@meta.data %>% subset(cell_type == 'LHb' & Group == "NR") %>% rownames
head(group1)
head(group2)

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='LHb'
)
head(DMEs)

pdf(paste0('hdWGCNA_figures/LHb_PlotDMEsLollipop.pdf'), width=4, height=5)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name='LHb', 
  pvalue = "p_val_adj"
)
dev.off()

pdf(paste0('hdWGCNA_figures/LHb_PlotDMEsVolcanop.pdf'), width=4, height=4)
PlotDMEsVolcano(
  seurat_obj,
  DMEs,
  wgcna_name = 'LHb'
)
dev.off()

#One-versus-all DME analysis
group.by = 'cell_type'

DMEs_all <- FindAllDMEs(
  seurat_obj,
  group.by = 'cell_type',
  wgcna_name = 'LHb'
)

head(DMEs_all)

p <- PlotDMEsVolcano(
  seurat_obj,
  DMEs_all,
  wgcna_name = 'LHb',
  plot_labels=FALSE,
  show_cutoff=FALSE
)

# facet wrap by each cell type
pdf(paste0('hdWGCNA_figures/LHb_PlotDMEsVolcano_One-versus-all.pdf'), width=6, height=7)
p + facet_wrap(~group, ncol=3)
dev.off()



##############################################################################################
# Microglia(MG)
# load the our snRNA-seq dataset
seurat_obj <- readRDS('Reinstatement_2022_seurat.rds')
head(seurat_obj)

p <- DimPlot(seurat_obj, group.by='cell_type', label=TRUE) +
   umap_theme() + ggtitle('UMAP') + NoLegend()

png('hdWGCNA_figures/DimPlot.png', width=4, height=4, res=200, units='in')
print(p)
dev.off()

seurat_obj <- seurat_obj[,!seurat_obj$cell_type %in% c('PER', 'END', 'EPD')]
head(seurat_obj)

#MG celltype for hdWGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  group.by = 'cell_type',
  wgcna_name = "MG"
)
length(GetWGCNAGenes(seurat_obj))

test <- SelectNetworkGenes(
    subset(seurat_obj, cell_type == 'MG'),
    gene_select = 'fraction',
    fraction = 0.05
)
length(GetWGCNAGenes(test))
cur_genes <- GetWGCNAGenes(test)
seurat_obj <- SetWGCNAGenes(seurat_obj, cur_genes)

test <- MetacellsByGroups(
  seurat_obj = test,
  group.by = c("cell_type","Sample", "Group"),
  k = 15, # k=25 # this was used on the first try
  max_shared=15,
  ident.group = 'cell_type',
  reduction='harmony',
  target_metacells=100,
  min_cells = 15
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type","Sample", "Group"),
  k = 25, # k=25 # this was used on the first try
  max_shared=15,
  ident.group = 'cell_type',
  reduction='harmony',
  target_metacells=100,
  min_cells = 50
)


seurat_obj <- NormalizeMetacells(seurat_obj)


seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "MG",
  group.by = "cell_type"
)

# test the soft power threshold
seurat_obj <- TestSoftPowers(seurat_obj)

plot_list <- PlotSoftPowers(seurat_obj)

#assemble with patchwork
pdf(paste0('hdWGCNA_figures/softpower_MG.pdf'), width=12, height=8)
wrap_plots(plot_list, ncol=2)
dev.off()


seurat_obj <- ConstructNetwork(
    seurat_obj, 
    soft_power=6,
    tom_name='MG', 
    overwrite_tom=TRUE
)

# plot the dendrogram
pdf(paste0("hdWGCNA_figures/dendro_MG.pdf"),height=3, width=6)
PlotDendrogram(seurat_obj, main='MG hdWGCNA Dendrogram')
dev.off()

modules <- GetModules(seurat_obj)
table(modules$module)
head(modules)


seurat_obj <- ScaleData(seurat_obj, features=rownames(seurat_obj)[1:100])
seurat_obj$Assignment <- droplevels(seurat_obj$Assignment)
seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars="Assignment" # snRNAseq batch
)

# compute module connectivity:
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type', group_name = 'MG'
)

########################################
# Feature Plot
########################################


plot_list <- ModuleFeaturePlot(
  seurat_obj, 
  order=TRUE, 
  raster=TRUE, 
  restrict_range=FALSE, 
  raster_dpi=400, 
  raster_scale=0.5, 
  point_size=1
  )

plot_list <- lapply(plot_list, function(x){
  x + NoLegend() + theme(plot.margin=margin(0,0,0,0))
})

pdf(paste0("hdWGCNA_figures/MG_featureplot_MEs.pdf"),height=6, width=12)
wrap_plots(plot_list, ncol=5)
dev.off()

########################################
# Dot Plot
########################################


MEs <- GetMEs(seurat_obj)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module)

write.csv(modules , file='MG_modules.csv')
mods <- mods[mods!='grey']
mods

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)


# make dotplot
p <- DotPlot(
  seurat_obj,
  group.by='annotation',
  features = rev(mods)
#  scale=FALSE, col.max=5
) + coord_flip() + RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') + xlab('') + ylab('') +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

pdf(paste0('hdWGCNA_figures/MG_MEs_dotplot.pdf'), width=9, height=4)
p
dev.off()

# remove from seurat obj
seurat_obj@meta.data <- seurat_obj@meta.data[,!(colnames(seurat_obj@meta.data) %in% colnames(MEs))]


################################################################################
# Hubgene circle plots:
################################################################################



# individual module networks
ModuleNetworkPlot(
  seurat_obj,
  mods = "all",
  #label_center=TRUE,
  outdir = paste0('hdWGCNA_figures/MG_hubNetworks/')
)

################################################################################
# UMAP:
################################################################################


seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 3,
  n_neighbors=10,
  min_dist=0.3,
  spread=3,
  supervised=TRUE,
  target_weight=0.5
)


# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

write.csv(umap_df , file='MG_modules_UMAP.csv')
# plot with ggplot
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
   color=umap_df$color,
   size=umap_df$kME*2
  ) +
  umap_theme()

pdf(paste0('hdWGCNA_figures/MG_hubgene_umap_ggplot2.pdf'), width=5, height=5)
p
dev.off()


library(reshape2)
library(igraph)
pdf(paste0('hdWGCNA_figures/MG_hubgene_umap_igraph.pdf'), width=10, height=10)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.5,
  sample_edges=TRUE,
  keep_grey_edges=FALSE,
  edge_prop=0.075, # taking the top 20% strongest edges in each module
  label_hubs=4, # how many hub genes to plot per module?
)
dev.off()

saveRDS(seurat_obj, file = "MG_hdWGCNA.rds")

######## Enrichment analysis ##################

dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')
head(dbs)
# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)
head(enrich_df)
# make GO term plots:
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(9,12), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

# enrichr dotplot
pdf(paste0('hdWGCNA_figures/MG_BP.pdf'), width=8, height=6)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

pdf(paste0('hdWGCNA_figures/MG_CC.pdf'), width=8, height=6)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Cellular_Component_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

pdf(paste0('hdWGCNA_figures/MG_MF.pdf'), width=8, height=6)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Molecular_Function_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)
dev.off()

#Marker gene overlap analysis
# compute cell-type marker genes with Seurat:
Idents(seurat_obj) <- seurat_obj$cell_type
markers <- Seurat::FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  logfc.threshold=1
)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)

# overlap barplot, produces a plot for each cell type
plot_list <- OverlapBarPlot(overlap_df)

# stitch plots with patchwork
pdf(paste0('hdWGCNA_figures/MG_Overlap_barplot.pdf'), width=5, height=6)
wrap_plots(plot_list, ncol=3)
dev.off()

# plot odds ratio of the overlap as a dot plot
pdf(paste0('hdWGCNA_figures/MG_Overlap_dotplot.pdf'), width=5, height=4)
OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cell-type markers')
dev.off()


############## DME analysis comparing two groups ###############
library(ggforestplot)
table(seurat_obj@meta.data$Group)

group1 <- seurat_obj@meta.data %>% subset(cell_type == 'MG' & Group == "Reinstatement") %>% rownames
group2 <- seurat_obj@meta.data %>% subset(cell_type == 'MG' & Group == "NR") %>% rownames
head(group1)
head(group2)

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='MG'
)
head(DMEs)

pdf(paste0('hdWGCNA_figures/MG_PlotDMEsLollipop.pdf'), width=4, height=5)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name='MG', 
  pvalue = "p_val_adj"
)
dev.off()

pdf(paste0('hdWGCNA_figures/MG_PlotDMEsVolcanop.pdf'), width=4, height=4)
PlotDMEsVolcano(
  seurat_obj,
  DMEs,
  wgcna_name = 'MG'
)
dev.off()

#One-versus-all DME analysis
group.by = 'cell_type'

DMEs_all <- FindAllDMEs(
  seurat_obj,
  group.by = 'cell_type',
  wgcna_name = 'MG'
)

head(DMEs_all)

p <- PlotDMEsVolcano(
  seurat_obj,
  DMEs_all,
  wgcna_name = 'MG',
  plot_labels=FALSE,
  show_cutoff=FALSE
)

# facet wrap by each cell type
pdf(paste0('hdWGCNA_figures/MG_PlotDMEsVolcano_One-versus-all.pdf'), width=6, height=7)
p + facet_wrap(~group, ncol=3)
dev.off()
















