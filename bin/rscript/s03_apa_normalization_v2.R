#!/usr/bin/env Rscript

# load library ------------------------------------------------------------
library(dplyr, quietly = T, warn.conflicts = F)
library(argparse, quietly = T, warn.conflicts = F)

# parameters --------------------------------------------------------------
parser <- ArgumentParser(description='Normalize APA matrix')

parser$add_argument('--sample', 
                    dest='sample_name',
                    action='store',
                    help='the sample name (eg. MOBV1)')
parser$add_argument('--method', 
                    dest='method',
                    action='store',
                    help='the APA identification method (eg. SCAPE)')
parser$add_argument('--APAannotation', 
                    dest='APAannotation',
                    action='store',
                    default='',
                    help='path to the APA annotation file (qs formart)')
parser$add_argument('--APAdir', 
                    dest='APAdir',
                    action='store',
                    default='',
                    help='path to APA output')      
parser$add_argument('--genedir', 
                    dest='genedir',
                    action='store',
                    default='',
                    help='path to single cell gene expression matrix (qs file or output dir of cell/space ranger)')
parser$add_argument('--barcodes',
                    dest='barcodes',
                    action='store',
                    default='barcodes.csv',
                    help='path to the cell barcodes or spatial barcodes (csv format)')
parser$add_argument('--nthreads', 
                    dest='nthreads',
                    action='store',
                    type='integer',
                    default=8,
                    help='the number of threads (default: 8)')
parser$add_argument('--pipelinedir',
                    dest='pipelinedir',
                    action='store',
                    default='.',
                    help='the pipeline folder')

args <- parser$parse_args()

## input load once
core_num <- args[["nthreads"]]
sample_name <- args[["sample_name"]] 
apa_method <- args[["method"]] 
barcode_file <- args[["barcodes"]]
apa_mat_dir <- args[["APAdir"]]
illumina_matrix_dir <- args[["genedir"]]
annot_file <- args[["APAannotation"]]
pipeline_dir <- args[["pipelinedir"]]
scale_factor <- 10000

# function ----------------------------------------------------------------
import::here(
  glue::glue("{pipeline_dir}/apa_analysis_functions.R"),
  "apa_matrix_extract",
  .character_only = TRUE
  )

# data --------------------------------------------------------------------
# Load barcodes
cell_ids <- readr::read_tsv(barcode_file, col_names = FALSE, show_col_types = FALSE)$X1 %>% gsub("-1", "", .)

# Load single cell or spot gene expression matrix for normalization
illumina_matrix_files <- list.files(illumina_matrix_dir, full.names = TRUE)
if(any(grepl("sparseMartix\\.qs$", illumina_matrix_files, ignore.case = TRUE))){
  message("Found .qs file，loading count sparse matrix...")
  illumina_matrix <- qs::qread(glue::glue("{illumina_matrix_dir}/{sample_name}_illumina_counts_sparseMartix.qs"), nthreads = core_num)
  colnames(illumina_matrix) <- gsub("-1", "", colnames(illumina_matrix))
  illumina_matrix <- illumina_matrix[, cell_ids]
} else if (any(grepl("\\.mtx\\.gz$", illumina_matrix_files, ignore.case = TRUE))){
  message("Found .mtx.gz file，loading count sparse matrix...")
  ## file path
  matrix_file <- file.path(illumina_matrix_dir, "matrix.mtx.gz")
  gene_file <- file.path(illumina_matrix_dir, "features.tsv.gz")
  barcode_file <- file.path(illumina_matrix_dir, "barcodes.tsv.gz")

  ## Read in expression matrix
  illumina_matrix <- Matrix::readMM(matrix_file)
  genes <- readr::read_tsv(gene_file, col_names = FALSE, show_col_types = FALSE)$X1
  barcodes <- readr::read_tsv(barcode_file, col_names = FALSE, show_col_types = FALSE)$X1

  ## Make the column names as the cell IDs and the row names as the gene IDs
  rownames(illumina_matrix) <- genes
  colnames(illumina_matrix) <- gsub("-1", "", barcodes)
  illumina_matrix <- illumina_matrix[, cell_ids]
}


## apa expression matrix
feature_list <- c("3UTRs", "Exon", "Intron", "3UTRs_1kb", "3UTRs_2kb", "LastExon1Kb") # this is based on SCAPE annotation strategy
apa_annot <- qs::qread(annot_file, nthreads =  core_num) %>% dplyr::filter(Type %in% feature_list) %>% tidyr::unite(col = "pa_site", c("seqnames", "pa_pos", "strand"), sep = ":", remove=F)

apa_matrix <- apa_matrix_extract(expr_file = apa_mat_dir, annot = apa_annot, method = apa_method, barcodes = cell_ids)

## save raw counts
output_file1 <- glue::glue("{sample_name}_{apa_method}_count_matrix.qs")
qs::qsave(apa_matrix, file = output_file1, nthreads = core_num)

# 2nd strategy
apa_matrix_annot <- apa_matrix %>% 
  dplyr::distinct(peakID, annot.symbol, .keep_all = TRUE) %>% 
  dplyr::rename(transcriptId = annot.tx_name, geneId = annot.symbol, ensemblId = annot.gene_id) %>% 
  dplyr::select(peakID, ensemblId, geneId, transcriptId) # 15032 -- 15027 peakID, ensemblId, geneId, transcriptId

apa_matrix <- apa_matrix %>%
  dplyr::distinct(peakID, annot.symbol, !!!syms(cell_ids), .keep_all = TRUE) %>% 
  dplyr::group_by(peakID) %>% 
  dplyr::summarise(across(all_of(cell_ids), sum)) %>% 
  tibble::column_to_rownames("peakID") 

apa_matrix_rownames <- rownames(apa_matrix)

apa_matrix <- apa_matrix %>% 
  Matrix::as.matrix() %>% 
  as(., "dgCMatrix")

if(isTRUE(all.equal(colnames(apa_matrix), colnames(illumina_matrix)))){

  col_sums <- Matrix::colSums(illumina_matrix)
  diag_col_sums_inv <- Matrix::Diagonal(x = 1 / col_sums)
  normalize_apa_matrix <- apa_matrix %*% diag_col_sums_inv
  normalize_apa_matrix <- normalize_apa_matrix * scale_factor
  colnames(normalize_apa_matrix) <- colnames(apa_matrix)
  
}

normalized_matrix_rowname <- rownames(normalize_apa_matrix)

normalize_apa_matrix <- normalize_apa_matrix %>% 
  as.matrix() %>% 
  dplyr::as_tibble() %>% 
  dplyr::mutate(peakID = normalized_matrix_rowname) %>% 
  dplyr::left_join(apa_matrix_annot, by = "peakID")

## save normalized counts
output_file2 <- glue::glue("{sample_name}_{apa_method}_normalized_by_ranger_matrix.qs")
qs::qsave(normalize_apa_matrix, file = output_file2, nthreads = core_num)
