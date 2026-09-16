# Single nucleus RNA sequencing of mouse habenula after reinstatement of cocaine self-administration

This study used single-nucleus RNA sequencing to characterize transcriptional states in the mouse habenula associated with cocaine-seeking behavior after withdrawal and reinstatement. Male wild-type C57BL/6J mice underwent intravenous cocaine self-administration followed by 30 days of withdrawal. Bilateral habenula tissue was collected for single-nucleus RNA sequencing using the 10x Genomics Chromium Next GEM Single Cell 3' platform. The dataset was used to characterize habenular cell populations and transcriptional states associated with reinstatement of cocaine-seeking behavior.

The analyses compare two experimental groups:

- **Reinstatement**: mice subjected to reinstatement testing
- **NR**: no-reinstatement controls

> **Note:** This repository contains analysis code rather than a fully packaged software application. Several scripts retain project-specific HPC paths and therefore require path updates before reuse. Raw and processed data are not stored in this code repository, but in the GEO repository.

## Data availability

Raw and processed data are not included in this repository. Add the public accession and a description of the deposited files here when available:

- **GEO accession:** 
  - Reinstatement Project: [GSE344264](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE344264) (Unpublished, Reviewer access only)
  - NR4A2/NURR1 Relapse Project: [GSE208081](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208081) (Published)
- **Website:** [Habenula Atlas](https://swaruplab.bio.uci.edu/habenula_ivsa_reinstatement/overview)

## Citation

If you use this code, please cite the associated manuscript:

> Childs J.*, Shi Z.*, et al. Single-nucleus transcriptomic analysis of the mouse habenula after reinstatement of cocaine seeking. 

## Repository structure

```text
.
├── CellChat/
│   ├── Step01_CellchatV2.R
│   └── Step02_CellchatV2_annotation_plotting.R
├── Data_processing/
│   ├── notebooks/
│   │   └── reinstatement_scanpy_processing.ipynb
│   ├── scanpy_to_seurat.Rmd
│   ├── plotting_for_paper.Rmd
│   └── vln.R
├── DEG/
│   ├── DEG_analysis.Rmd
│   └── DGE_scripts/
│       ├── celltype_markers.sub
│       ├── cluster_markers.sub
│       ├── celltype_Reinstatement_vs_WT.sub
│       ├── cluster_Reinstatement_vs_WT.sub
│       └── plotting_cluster_Reinstatement_vs_WT.sub
├── MELD/
│   ├── Step01ab_run_MELD_reinstatement.ipynb
│   └── Step02_plotting.R
└── hdWGCNA/
    └── hdWGCNA_reinstatement.r
```

## Analysis workflow

### 1. snRNA-seq preprocessing and annotation

`Data_processing/notebooks/reinstatement_scanpy_processing.ipynb`

The Scanpy workflow:

1. Loads CellBender-processed count matrices from multiple sequencing batches.
2. Adds sample and experimental metadata.
3. Performs quality control and removes low-quality nuclei and outlier clusters.
4. Normalizes and log-transforms counts.
5. Selects highly variable genes and regresses out library-size and mitochondrial-count effects.
6. Performs PCA and Harmony batch correction using `Assignment` as the batch variable.
7. Constructs the neighbor graph, calculates UMAP embeddings, and performs Leiden clustering.
8. Assigns broad cell classes and finer habenular and glial annotations.
9. Exports the count matrix, metadata, Harmony embeddings, and processed AnnData object.

The annotated cell populations include medial habenula (MHb), lateral habenula (LHb), parahabenular (PHb), astrocyte (ASC), oligodendrocyte precursor (OPC), oligodendrocyte (ODC), microglia (MG), pericyte (PER), endothelial (END), and ependymal (EPD) populations.

### 2. Conversion to Seurat and descriptive visualization

`Data_processing/scanpy_to_seurat.Rmd` imports the Scanpy outputs into a Seurat object, attaches the Harmony and UMAP reductions, sets annotation levels, normalizes the RNA assay, and saves `Reinstatement_2022_seurat.rds`.

`Data_processing/plotting_for_paper.Rmd` generates dataset-level quality-control and summary figures, including:

- UMAPs grouped by cell type, annotation, sample, sequencing batch, and condition
- cell-composition plots by sample and batch
- quality-control violin plots
- selected gene and signature visualizations

### 3. Differential expression

The SLURM scripts in `DEG/DGE_scripts/` run parallel differential-expression analyses at two resolutions:

- broad cell type (`cell_type`)
- fine-grained cell annotation (`annotation`)

They calculate both cluster markers and Reinstatement-versus-NR differential expression using MAST, with `total_counts` and sequencing assignment included as latent variables. `DEG/DEG_analysis.Rmd` combines the per-cluster output files and generates heatmaps, dot plots, volcano plots, and selected-gene violin plots.

The `.sub` files call lab-specific helper scripts (`parallel_DEGs.R` and `parallel_DEG_plotting.R`) that are not included here. To reproduce these analyses independently, replace those calls with an equivalent Seurat `FindMarkers()`/MAST workflow or add the required helper scripts.

### 4. MELD analysis

`MELD/Step01ab_run_MELD_reinstatement.ipynb` estimates a continuous, graph-based Reinstatement likelihood for each nucleus. The notebook benchmarks MELD parameters using the Harmony embedding, fits the selected model, and exports per-nucleus likelihood scores.

`MELD/Step02_plotting.R` adds the likelihood scores to the Seurat metadata and compares their distributions across annotated cell populations.

### 5. Cell–cell communication

`CellChat/Step01_CellchatV2.R` uses the mouse CellChat database to infer signaling networks separately in the Reinstatement and NR groups at both broad cell-type and fine-annotation resolutions. Interactions represented by fewer than 50 cells are filtered before pathway aggregation and centrality analysis.

`CellChat/Step02_CellchatV2_annotation_plotting.R` merges the condition-specific CellChat objects and compares interaction numbers, interaction strengths, signaling pathways, and outgoing and incoming signaling roles.

### 6. Co-expression network analysis

`hdWGCNA/hdWGCNA_reinstatement.r` performs hdWGCNA separately for major cell populations, including ODC, ASC, MHb, LHb, and MG. The workflow includes:

- metacell construction
- soft-power selection
- co-expression module construction
- module eigengene harmonization
- hub-gene and module visualization
- functional enrichment and marker-overlap analysis
- differential module eigengene analysis between Reinstatement and NR

## Recommended execution order

Run the analysis in the following order after updating input and output paths:

1. `Data_processing/notebooks/reinstatement_scanpy_processing.ipynb`
2. `Data_processing/scanpy_to_seurat.Rmd`
3. `Data_processing/plotting_for_paper.Rmd`
4. Scripts in `DEG/DGE_scripts/`, followed by `DEG/DEG_analysis.Rmd`
5. `MELD/Step01ab_run_MELD_reinstatement.ipynb`, followed by `MELD/Step02_plotting.R`
6. `CellChat/Step01_CellchatV2.R`, followed by `CellChat/Step02_CellchatV2_annotation_plotting.R`
7. `hdWGCNA/hdWGCNA_reinstatement.r`

The DEG, MELD, CellChat, and hdWGCNA analyses are downstream branches and can be run independently once the annotated Seurat or AnnData object has been generated.

## Samples

| GEO samples | Description | Sample ID |
|---|---|---|
| GSM9973884 | Mouse habenula, No Reinstatement | Sample-6
| GSM9973885 | Mouse habenula, Reinstatement | Sample-7
| GSM9973886 | Mouse habenula, No Reinstatement | Sample-10
| GSM9973887 | Mouse habenula, Reinstatement | Sample-15
| GSM9973888 | Mouse habenula, No Reinstatement | Sample-1-B3
| GSM9973889 | Mouse habenula, Reinstatement | Sample-2-B3
| GSM9973890 | Mouse habenula, No Reinstatement | Sample-3-B3
| GSM9973891 | Mouse habenula, Reinstatement | Sample-4-B3
| GSM9973892 | Mouse habenula, Reinstatement | Sample-5-B3
| GSM9973893 | Mouse habenula, Reinstatement | Sample-6-B3
| GSM9973894 | Mouse habenula, No Reinstatement | Sample-7-B3
| GSM9973895 | Mouse habenula, No Reinstatement | Sample-8-B3


## Processed data

The processed snRNA-seq dataset is publicly available from NCBI Gene Expression Omnibus under accession [GSE344264](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE344264). The deposited files provide a sparse count matrix, its corresponding barcode and feature annotations, and per-nucleus metadata.

| File | Size | Format | Description |
|---|---:|---|---|
| `GSE344264_Reinstatement_2022_barcodes.tsv.gz` | 340.3 KB | TSV | Nucleus barcodes corresponding to the count matrix |
| `GSE344264_Reinstatement_2022_cell_metadata.tsv.gz` | 2.5 MB | TSV | Per-nucleus sample, condition, quality-control, clustering, and cell-annotation metadata |
| `GSE344264_Reinstatement_2022_features.tsv.gz` | 226.8 KB | TSV | Gene features corresponding to the count matrix |
| `GSE344264_Reinstatement_2022_raw_counts.mtx.gz` | 356.0 MB | Matrix | Sparse raw gene-count matrix |

Download all four files to reconstruct the deposited expression dataset. Keep the files together and preserve their row and column ordering when matching the count matrix to the feature and barcode tables. The cell metadata can then be joined to the barcodes using the corresponding nucleus identifier.


## Contact

For questions about the analysis, please open a GitHub [issue](https://github.com/swaruplabUCI/snRNAseq_n_Spatial_mouse_habenula_after_reinstatement/issues).

