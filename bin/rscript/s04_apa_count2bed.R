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
parser$add_argument('--barcodes',
                    dest='barcodes',
                    action='store',
                    default='barcodes.csv',
                    help='path to the cell barcodes or spatial barcodes (csv format)')
parser$add_argument('--apa_count',
                    dest='apa_count_file',
                    action='store',
                    default='~/nanopore_validation/predicted_sites',
                    help='path to apa counts with annotation (folder or qs file)')
parser$add_argument('--nthreads',
                    dest='nthreads',
                    action='store',
                    type='integer',
                    default=8,
                    help='the number of threads (default: 8)')

args <- parser$parse_args()

## input load once
core_num <- args[["nthreads"]]
sample_name <- args[["sample_name"]]
apa_method <- args[["method"]]
barcode_file <- args[["barcodes"]]
apa_count_file <- args[["apa_count_file"]]
# function ----------------------------------------------------------------
count2bed <- function(apa_count_file, sample_name, apa_method, barcode_file, core_num=4) {
  # Load barcodes
  cell_ids <- readr::read_tsv(barcode_file, col_names = FALSE, show_col_types = FALSE)$X1 %>% gsub("-1", "", .)

  # Load count matrix (qs format) and convert into bed format
  peak_bed <- qs::qread(apa_count_file, nthreads = core_num) %>%
    dplyr::mutate(total = rowSums(select(., all_of(cell_ids)))) %>%
    dplyr::select(total, peakID, annot.symbol) %>%
    dplyr::filter(!is.na(annot.symbol)) %>%
    dplyr::filter(total > 0) %>%
    dplyr::mutate(name = glue::glue("{apa_method}:{annot.symbol}:{peakID}")) %>%
    tidyr::separate(col = peakID, into = c("seqnames", "sites", "strand"), sep = ":") %>%
    # tidyr::unite(c("annot.symbol", "peakID"), col = "name" , sep = ":") %>%
    dplyr::mutate(start = dplyr::if_else(strand == "+", as.integer(sites), as.integer(sites) - 1),
                  end = dplyr::if_else(strand == "+", as.integer(sites) + 1, as.integer(sites))) %>%
    dplyr::rename(score = total) %>%
    dplyr::select(seqnames, start, end, name, score, strand) # score is the total count

  return(peak_bed)
}


filename_pattern <- glue::glue("{sample_name}_{apa_method}_count_matrix.qs")
if(grepl(filename_pattern, apa_count_file, ignore.case = TRUE)){
  count_file <- apa_count_file
} else if(dir.exists(apa_count_file)) {
  count_file <- list.files(apa_count_file, pattern = filename_pattern, full.names = T, recursive = T)
} else {
  stop("Wrong input file! Should be folder or qs file.")
}

peak_bed <- count2bed(count_file, sample_name, apa_method, barcode_file, core_num)
vroom::vroom_write(peak_bed, file = glue::glue("{sample_name}_{apa_method}_pA_site.bed"), delim = "\t", col_names = FALSE)
