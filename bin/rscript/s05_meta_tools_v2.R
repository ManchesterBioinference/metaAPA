#!/usr/bin/env Rscript
library(dplyr, warn.conflicts = F)
library(argparse, warn.conflicts = F)

set.seed(42)
# nohup Rscript rscripts/s01_meta_tools.R --sample MOBV1 --nthreads 4 --apacutoff 0 --distance_method norm --cluster_method kmeans-asw &
# parameters --------------------------------------------------------------
parser <- ArgumentParser(description='Meta tools for APA identification')

parser$add_argument('--sample', 
                    dest='sample_name',
                    action='store',
                    help='the sample name (eg. MOBV1)')
parser$add_argument('--apa_count_folder', 
                    dest='apa_count_folder',
                    action='store',
                    default='~/nanopore_validation/predicted_sites',
                    help='path to apa counts with annotation for single sample')
parser$add_argument('--nthreads', 
                    dest='nthreads',
                    action='store',
                    type='integer',
                    default=8,
                    help='the number of threads (default: 8)')
parser$add_argument('--apacutoff', 
                    dest='apa_cutoff',
                    action='store',
                    type='integer',
                    default=0,
                    help='the cutoff of APA method (default: 0)')
parser$add_argument('--gene_list', 
                    dest='gene_list',
                    action='store',
                    default="",
                    help='path to merged BED file.')
parser$add_argument('--barcodes',
                    dest='barcodes',
                    action='store',
                    default='barcodes.csv',
                    help='path to the cell barcodes or spatial barcodes (csv format)')
parser$add_argument('--distance_method', 
                    dest='distance_method',
                    action='store',
                    default="",
                    help='Distance method.')
parser$add_argument('--cluster_method', 
                    dest='cluster_method',
                    action='store',
                    default="",
                    help='Center method.')
parser$add_argument('--pipelinedir',
                    dest='pipelinedir',
                    action='store',
                    default='.',
                    help='the pipeline folder')
args <- parser$parse_args()

## input load once
core_num <- args[["nthreads"]]
sample_name <- args[["sample_name"]] 
apa_count_folder <- args[["apa_count_folder"]]
apa_reads_cutoff <- args[["apa_cutoff"]]
barcode_file <- args[["barcodes"]]
distance_method <- args[["distance_method"]]
cluster_method <- args[["cluster_method"]]
gene_list <- args[["gene_list"]]
pipeline_dir <- args[["pipelinedir"]]
n_starts <- 25
# functions ---------------------------------------------------------------
import::here(
  glue::glue("{pipeline_dir}/metatools_functions.R"), 
  "calculate_cluster_metrics", 
  "extract_pA_sites",
  "calculate_distance",
  "find_elbow_point",
  "find_optimal_center",
  "align_pA_sites",
  .character_only = TRUE
)
# meta-tools --------------------------------------------------------------
gene_list <- vroom::vroom(gene_list, col_names = c("seqnames", "start", "end", "name", "score", "strand")) %>% 
  tidyr::separate(col = "name", into = c("method", "gene_name", "chrom", "pa_site", "strand"), sep = ":") %>% 
  dplyr::pull("gene_name") %>% 
  unique() 
cell_ids <- readr::read_tsv(barcode_file, col_names = FALSE, show_col_types = FALSE)$X1 %>% gsub("-1", "", .)

filename_pattern <- glue::glue("count_matrix.qs")
count_file_fields <- list.files(apa_count_folder, pattern = filename_pattern, full.names = T, recursive = T) %>% 
  basename() %>% 
  stringr::str_split_fixed(pattern="_", n = 4)

n_sample <- count_file_fields[, 1] %>% unique() %>% length()
if(n_sample == 1) {
  apa_method_list <- count_file_fields[, 2] %>% unique()
  } else {
  stop("Please only give one sample!")
  }

sites_data <- parallel::mclapply(
  gene_list, 
  align_pA_sites, 
  cell_ids = cell_ids, 
  apa_method_list = apa_method_list,
  apa_count_folder = apa_count_folder,
  count_cutoff = apa_reads_cutoff, 
  n_starts = n_starts,
  dist_method = distance_method,
  cluster_method = cluster_method,
  core_num = 1L,
  mc.cores = core_num)

qs::qsave(sites_data, file = glue::glue("{sample_name}_{distance_method}_{cluster_method}_{apa_reads_cutoff}filter_site_cluster.qs"), nthreads = core_num)


