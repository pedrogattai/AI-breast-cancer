Librarys to be downloaded on R
library(tidyverse)
library(edgeR)        # for cpm normalization
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(pdftools)      # text extraction from PDF
library(stringr)
library(dplyr)
library(purrr)
library(tibble)
library(readr)
# optional for scraped tables:
# install.packages("tabulizer")  # may require Java
# optional OCR if PDF is scanned:
# install.packages("tesseract")
library(GEOquery)



Download all PDFs locally.

-Extract raw text from every page (pdftools). If PDFs are scanned images, fall back to OCR (tesseract).

-Split into sentences and search for (a) words like cell, cells, single-cell, nuclei paired with numeric expressions*,* and (b) accession IDs / repository names (GEO, SRA, ENA, ArrayExpress, Zenodo, Figshare, etc.).

-Normalize all numeric mentions (handle commas, periods, “k”, scientific notation, ranges like “~5k–7k”).

-Link numeric mentions to accession IDs when they appear nearby (same sentence / same paragraph / same figure/table reference).

-For each accession ID found, optionally query the repository APIs (GEOquery, rentrez/SRAdb, ENA REST) or download the supplementary matrix and count cells directly (most robust).

-Output a CSV/TSV with one row per (manuscript, accession, raw matched sentence, normalized cell count, confidence score, notes).



On R
# extract_cells_from_pdfs.R
# Run: Rscript extract_cells_from_pdfs.R
library(pdftools)      # text extraction from PDF
library(stringr)
library(dplyr)
library(purrr)
library(tibble)
library(readr)
# optional for scraped tables:
# install.packages("tabulizer")  # may require Java
# optional OCR if PDF is scanned:
# install.packages("tesseract")

# ----- user inputs -----
pdf_urls <- c(
  "your_path/your_directory/1.pdf",
  "your_path/your_directory/2.pdf"
  # ... add all 44 URLs here
)
out_dir <- "pdfs"
dir.create(out_dir, showWarnings = FALSE)

# ----- helper functions -----
download_pdfs <- function(urls, dir = ".") {
  out_files <- map_chr(urls, function(u){
    fn <- file.path(dir, basename(u))
    tryCatch({
      download.file(u, fn, mode = "wb", quiet = TRUE)
      fn
    }, error = function(e) {
      message("Failed download: ", u, " -> ", e$message)
      NA_character_
    })
  })
  out_files
}

# read pdf -> character vector, one element per page
read_pdf_text_safely <- function(pdf_path) {
  tryCatch({
    txt <- pdftools::pdf_text(pdf_path)
    paste(txt, collapse = "\n\n") # combine pages, keep separators
  }, error = function(e) {
    message("pdftools failed on ", pdf_path, ": ", e$message)
    NA_character_
  })
}

# regex patterns
accession_patterns <- c(
  "GSE\\d+",        # GEO Series
  "GSM\\d+",        # GEO Sample
  "SRP\\d+", "SRA\\d+", "SRX\\d+", "SRR\\d+", # SRA/SRP/SRX/SRR
  "PRJNA\\d+", "PRJEB\\d+",                     # BioProject
  "ERP\\d+", "ERS\\d+",
  "E-MTAB-\\d+", "E-MEXP-\\d+"                  # ArrayExpress
)
accession_regex <- paste0("(", paste(accession_patterns, collapse = "|"), ")")

# numbers close to "cell(s)" or "single cell" etc.
cell_word_regex <- "(?:cell|cells|single-cell|single cell|nuclei|nucleus|droplets)"
# pattern to find numeric-like expressions: 1,234 or 1.2k or 3e4 or 5 × 10^3
number_like <- "\\b\\~?\\d{1,3}(?:[,.]\\d{3})*(?:[.,]\\d+)?\\s*(?:k|K|thousand|m|M|million)?\\b|\\b\\d+(?:\\.\\d+)?[eE][+-]?\\d+\\b|\\b\\d+\\s*×\\s*10\\^\\d+\\b"
# look for numbers within 50 characters of a cell word
cells_with_numbers_regex <- paste0("(?i)(.{0,60}(", number_like, ").{0,60}(", cell_word_regex, ").{0,60})|(?i)(.{0,60}(", cell_word_regex, ").{0,60}(", number_like, ").{0,60})")

# normalize numeric expression to integer
normalize_count <- function(x) {
  x2 <- str_trim(x)
  # remove tildes etc.
  x2 <- str_replace_all(x2, "[~\\s]+", "")
  # lowercase
  x2l <- tolower(x2)
  # scientific (e.g., 1.2e4)
  if (str_detect(x2l, "e[+-]?\\d+")) {
    return(as.numeric(x2l))
  }
  # × 10^n
  if (str_detect(x2l, "×")) {
    # handle "5 × 10^3"
    parts <- str_match(x2l, "(\\d+(?:\\.\\d+)?)\\s*×\\s*10\\^(\\d+)")
    if (!is.na(parts[1,2])) {
      return(as.numeric(parts[1,2]) * 10^as.numeric(parts[1,3]))
    }
  }
  # suffix k/m
  if (str_detect(x2l, "[km]$")) {
    num <- as.numeric(str_replace_all(x2l, "[^0-9\\.]", ""))
    if (str_ends(x2l, "k")) return(round(num * 1e3))
    if (str_ends(x2l, "m")) return(round(num * 1e6))
  }
  # remove commas and dots used as thousands separators (heuristic)
  if (str_detect(x2l, "\\d{1,3}(?:[,.]\\d{3})+")) {
    x_clean <- str_replace_all(x2l, "[,]", "")
    return(as.numeric(x_clean))
  }
  # final numeric parse
  num <- suppressWarnings(as.numeric(str_replace_all(x2l, "[^0-9\\.\\-]", "")))
  if (!is.na(num)) return(num)
  NA_real_
}

# find accession IDs and cell mentions
extract_matches <- function(text, pdf_file) {
  sentences <- unlist(str_split(text, "(?<=[\\.!?;\\n])\\s+", perl = TRUE))
  tibble(sentence = sentences) %>%
    mutate(
      pdf = pdf_file,
      accessions = str_extract_all(sentence, accession_regex),
      has_accession = map_lgl(accessions, ~ length(.x) > 0),
      cell_matches = str_extract_all(sentence, regex(cells_with_numbers_regex, ignore_case = TRUE)),
      has_cell_number = map_lgl(cell_matches, ~ length(.x) > 0)
    ) %>%
    filter(has_accession | has_cell_number) %>%
    mutate(
      numbers_raw = map(cell_matches, ~ str_extract_all(.x, number_like)[[1]]),
      numbers_norm = map(numbers_raw, ~ map_dbl(.x, ~ normalize_count(.x)))
    )
}

# ----- run -----
files <- download_pdfs(pdf_urls, out_dir)
files <- files[!is.na(files)]

results <- map_dfr(files, function(f){
  txt <- read_pdf_text_safely(f)
  if (is.na(txt)) return(tibble())
  extract_matches(txt, basename(f))
})

# convert lists to readable strings
results2 <- results %>%
  mutate(
    accessions = map_chr(accessions, ~ paste(unique(.x), collapse = ";")),
    numbers_raw = map_chr(numbers_raw, ~ if(length(.x)==0) NA_character_ else paste(.x, collapse=";")),
    numbers_norm = map_chr(numbers_norm, ~ if(all(is.na(.x))) NA_character_ else paste(.x, collapse=";"))
  ) %>%
  select(pdf, sentence, accessions, numbers_raw, numbers_norm)

# write results
write_csv(results2, file.path(out_dir, "extracted_cell_mentions.csv"))
message("Wrote ", file.path(out_dir, "extracted_cell_mentions.csv"))

___________________________________________________________________
Put all URLs in pdf_urls.

The script writes extracted_cell_mentions.csv with sentence-level matches, raw number strings, and normalized numeric values (if parsable).

You should manually inspect low-confidence rows (NA numeric or ambiguous ranges).



If you find a GEO accession (GSExxxx or GSMxxxx): use Bioconductor GEOquery to download supplementary files or expression matrices; then count columns (for many scRNA matrices each column = cell). Example:

On R

BiocManager::install("GEOquery")
library(GEOquery)
gse <- "GSE123456"
getGEOSuppFiles(gse, makeDirectory = TRUE)
# find a matrix file (matrix.mtx or counts.csv), read and count columns / cells.

_________________________________________________________________________

For SRA/BioProject (PRJNA...), use rentrez or the ENA REST API to list runs/samples; then inspect whether processed matrices/supplementary files are available on GEO/ENA/figshare/Zenodo. Example using ENA REST (httr/json) — ENA often lists sample/run metadata and links to submitted processed files.

If the paper points to a Supplementary Table giving the cell counts per sample, extract tables from PDF via tabulizer::extract_tables() or check supplemental material (often on the journal website) and parse it.


If only raw FASTQ are available (SRA), the number of cells is not directly stored there — you will need processed matrices or sample metadata from GEO/Zenodo, or re-process raw reads (heavy).

___________________________________________________________________________

how to create heatmaps from specific gene expression extracted from public repositories used in a pool of manuscripts, using the R package?

Load expression matrix (genes × samples or genes × cells).

Preprocess (filter, normalize, optionally log-transform).

Subset to the genes of interest.

(For scRNA) optionally aggregate per cluster/sample or downsample cells.

Create heatmap(s) with annotations, clustering and color control.

Save to PNG/PDF.

Packages to install (one-time)
# CRAN
install.packages(c("tidyverse","pheatmap","ComplexHeatmap","RColorBrewer","circlize"))

# Bioconductor (for GEO, single-cell helpers)
if (!requireNamespace("BiocManager",quiet=TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery","edgeR","SingleCellExperiment","scater"))

# Seurat (for single-cell workflows)
install.packages("Seurat")  # CRAN version
__________________________________________________________________________
Example: bulk / processed matrix (genes × samples)
Assume you have a matrix file counts.csv or a matrix object expr with rownames = gene symbols and columns = samples.

library(tidyverse)
library(edgeR)        # for cpm normalization
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Load matrix (rows = genes, cols = samples)
expr <- read.csv("counts.csv", row.names = 1, check.names = FALSE)
# Or if matrix already loaded: expr <- my_matrix

# Basic QC / filtering
min_counts_per_gene <- 5
keep <- rowSums(expr >= min_counts_per_gene) >= 2   # keep genes with counts in at least 2 samples
expr_f <- expr[keep, ]

# Normalize to CPM and log-transform (robust for visualization)
cpm_mat <- cpm(expr_f, log=FALSE, prior.count=1)   # edgeR::cpm gives normalized counts-per-million
log_cpm <- log2(cpm_mat + 1)

# Genes of interest (vector of gene symbols)
genes_of_interest <- c("TP53","GAPDH","CD3D","CD8A","MYC") # replace with your gene list
genes_present <- intersect(genes_of_interest, rownames(log_cpm))
if (length(genes_present)==0) stop("None of your genes of interest are present in the matrix")

mat_subset <- log_cpm[genes_present, , drop=FALSE]

# Optional: scale rows (z-score) to highlight relative expression
mat_z <- t(scale(t(mat_subset)))   # scale genes (rows)

# Sample annotation (data.frame with rownames matching column names of matrix)
sample_meta <- data.frame(
  condition = c("ctrl","ctrl","treat","treat"),
  row.names = colnames(mat_z)
)

# Color palette and breaks
col_fun <- colorRamp2(c(-2, 0, 2), c("navy","white","firebrick"))

# ComplexHeatmap (more flexible)
Heatmap(mat_z,
        name = "z-score",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        top_annotation = HeatmapAnnotation(df = sample_meta),
        column_title = "Samples",
        row_title = "Genes")

# Or pheatmap (simple)
pheatmap::pheatmap(mat_z,
                   annotation_col = sample_meta,
                   clustering_distance_rows = "euclidean",
                   clustering_method = "complete",
                   show_rownames = TRUE)

__________________________________________________________________________
# PNG
png("heatmap_genes.png", width = 1600, height = 1000, res = 150)
ComplexHeatmap::draw(Heatmap(mat_z, name="z-score", col=col_fun, top_annotation = HeatmapAnnotation(df=sample_meta)))
dev.off()

# PDF
pdf("heatmap_genes.pdf", width = 10, height = 7)
ComplexHeatmap::draw(Heatmap(mat_z, name="z-score", col=col_fun, top_annotation = HeatmapAnnotation(df=sample_meta)))
dev.off()

___________________________________________________________________________
Example: single-cell RNA-seq (cells × genes) — many cells
Single-cell data usually has thousands of cells — you’ll commonly either (A) make heatmaps of averaged expression per cluster/sample, or (B) downsample cells and plot per-cell heatmap (only for moderate numbers).

(A) Best practice: make heatmap of average expression per cluster (Seurat)
library(Seurat)
library(circlize)
library(ComplexHeatmap)

# Suppose you have a Seurat object 'seu' already created and normalized
# seu <- CreateSeuratObject(counts = counts_matrix); seu <- NormalizeData(seu); seu <- FindVariableFeatures(seu); seu <- ScaleData(seu)

# If cluster identity exists:
Idents(seu) <- "seurat_clusters"  # or "cell_type" or sample metadata

# Genes of interest
genes_of_interest <- c("CD3D","CD8A","MS4A1","LYZ","MKI67")
genes_present <- intersect(genes_of_interest, rownames(seu))
avg_expr <- AverageExpression(seu, features = genes_present, assays = "RNA", slot = "data")$RNA
# avg_expr now: genes x clusters (or samples, depending on Idents)

# z-score per gene (row)
avg_z <- t(scale(t(avg_expr)))

Heatmap(avg_z,
        name = "z-score",
        column_title = "Clusters",
        row_title = "Genes",
        cluster_rows = TRUE,
        cluster_columns = TRUE)
(B) Per-cell heatmap (only for a subset of cells)
# Subset to e.g., 500 random cells across conditions
set.seed(42)
cells_sampled <- sample(colnames(seu), size = 500)
mat_cells <- as.matrix(GetAssayData(seu, assay="RNA", slot="data")[genes_present, cells_sampled])

# scale rows
mat_cells_z <- t(scale(t(mat_cells)))

pheatmap(mat_cells_z, show_colnames = FALSE, show_rownames = TRUE)
Tip: Per-cell heatmaps can be huge — if you want to visually inspect patterns across many cells, aggregate by cluster or compute metacells.

__________________________________________________________________________

If your gene list is large (>50), split into biological modules or visualize with hierarchical clustering and annotate modules.

Scale rows (z-score) if you want to emphasize relative expression across samples/clusters. Don’t z-score if you need absolute expression values.

Use ComplexHeatmap for complex annotations (multiple annotation tracks, color bars, dendrogram control).

Use pheatmap for quick exploratory plots.

For scRNA: prefer aggregated (AverageExpression) heatmaps, or violin/dot plots for individual genes (Seurat VlnPlot, DotPlot) when per-cell heatmaps are impractical.

Example complete script (single file) — from matrix to saved heatmap

# full_example_heatmap.R
library(tidyverse); library(edgeR); library(ComplexHeatmap); library(circlize)

# Load
expr <- read.csv("counts.csv", row.names=1, check.names=FALSE)
genes_of_interest <- c("TP53","CD3D","CD8A","GAPDH","MYC")

# Filter & normalize
keep <- rowSums(expr >= 5) >= 2
expr <- expr[keep, ]
log_cpm <- log2(cpm(expr, prior.count = 1) + 1)

# Subset genes
genes_present <- intersect(genes_of_interest, rownames(log_cpm))
mat <- log_cpm[genes_present, , drop=FALSE]
mat_z <- t(scale(t(mat)))

# sample metadata
sample_meta <- data.frame(condition = c("A","A","B","B"), row.names = colnames(mat_z))

# Heatmap
col_fun <- colorRamp2(c(-2,0,2), c("navy","white","firebrick"))
ht <- Heatmap(mat_z, name="z-score", col = col_fun,
              top_annotation = HeatmapAnnotation(df = sample_meta),
              cluster_rows = TRUE, cluster_columns = TRUE)

# Save
pdf("heatmap_genes_example.pdf", width=8, height=6)
draw(ht)
dev.off()