# extract expression matrix
apa_matrix_extract <- function(expr_file, annot, method, barcodes){
  
  if(method == "Sierra"){
    
    ## file path
    sierra_mat_file <- file.path(expr_file, "Sierra_counts", "matrix.mtx.gz")
    sierra_sit <- file.path(expr_file, "Sierra_counts", "sitenames.tsv.gz")
    sierra_bcd <- file.path(expr_file, "Sierra_counts", "barcodes.tsv.gz")
    
    ## Read in expression matrix
    counts <- Matrix::readMM(sierra_mat_file)
    gene_ids <- readr::read_tsv(sierra_sit, col_names = FALSE, show_col_types = FALSE)$X1
    cell_ids <- readr::read_tsv(sierra_bcd, col_names = FALSE, show_col_types = FALSE)$X1
    
    ## Make the column names as the cell IDs and the row names as the gene IDs
    rownames(counts) <- gene_ids
    colnames(counts) <- gsub("-1", "", cell_ids)
    counts <- as.data.frame(as.matrix(counts))
    
    counts <- counts %>% 
      tibble::rownames_to_column("peakID") %>% 
      tidyr::separate(col = peakID, into = c("gene_name", "chr", "region", "strand"), sep = ":") %>% 
      tidyr::separate(col = region, into = c("start", "end"), sep = "-") %>% 
      dplyr::mutate(strand = dplyr::if_else(strand == 1, "+", "-"),
                    coord = dplyr::if_else(strand == "+", end, start)) %>% 
      tidyr::unite(col = "peakID", c("chr", "coord", "strand"), sep = ":")
  }
  
  if(method == "SCAPE"){
    
    ## file path
    scape_mat_file <- file.path(expr_file, "pasite.csv.gz")
    
    ## rename and reformat
    counts <- vroom::vroom(scape_mat_file, show_col_types = F) %>% 
      tibble::column_to_rownames("...1") %>% 
      dplyr::rename_with(~gsub("-1", "", .x)) %>% 
      tibble::rownames_to_column("peakID") %>% 
      tidyr::separate(col = "peakID", into = c("chr", "coord", "score", "strand"), sep = ":") %>% 
      tidyr::unite(col = "peakID", c("chr", "coord", "strand"), sep = ":") 
    
  }
  
  if(method == "scAPAtrap"){
    ## file path
    scapatrap_expr_file <- file.path(expr_file, "expma.qs")
    
    ## rename
    counts <- qs::qread(scapatrap_expr_file) %>% tidyr::unite(col = "peakID", c("chr", "coord", "strand"), sep = ":")
    
  }
  
  if(method == "polyApipe"){
    
    ## file path
    polyapipe_mat_file <- list.files(expr_file, pattern = "counts.tab.gz", full.names = T)
    
    ## rename and reformat
    counts <- vroom::vroom(polyapipe_mat_file, show_col_types = FALSE) %>% 
      dplyr::filter(cell %in% barcodes) %>% 
      tidyr::separate(col = gene, into = c("chr", "coord", "strand"), sep = "_") %>%
      dplyr::mutate(strand = recode(strand, "r" = "-", "f" = "+")) %>% 
      tidyr::unite(col = "peakID", c("chr", "coord", "strand"), sep = ":") %>% 
      tidyr::pivot_wider(names_from = cell, values_from = count, values_fill = 0) 
    
  }
  
  if(method == "SCAPTURE"){
    
    ## file path
    scapture_mat_file <- list.files(expr_file, pattern = "PASquant.KeepCell.UMIs.tsv.gz", full.names = T)
    scapture_bed_file <- list.files(expr_file, pattern = "PASquant.KeepPAS.bed", full.names = T)
    ## rename and reformat
    
    scapture_bed <- vroom::vroom(scapture_bed_file, col_names = F) %>% 
      dplyr::mutate(pa_site = dplyr::if_else(X6 == "+", X3, X2)) %>% 
      tidyr::unite(col = "peakID", c(X1, pa_site, X6), sep = ":") %>% 
      dplyr::rename(peak_name = X4) %>% 
      dplyr::select(peakID, peak_name)
    
    counts <- vroom::vroom(scapture_mat_file, show_col_types = FALSE) %>% 
      tibble::column_to_rownames("gene") %>% 
      dplyr::rename_with(~gsub("-1", "", .x)) %>% 
      tibble::rownames_to_column("peak_name") %>%
      dplyr::left_join(scapture_bed, by = "peak_name")

  }

  if(method == "scAPA"){
    
    ## file path
    peak_file <- list.files(expr_file, pattern = "Peaks.RDS", full.names = T, recursive = T)
    
    ## rename and reformat
    scAPA_bed <- peak_file %>% 
      readRDS(.) %>% 
      .@row.Data %>% 
      tidyr::separate_rows(Chr, Start, End, Strand, sep = ";") %>%
      dplyr::group_by(GeneID) %>%
      dplyr::summarise(
        Chr = unique(Chr),
        Start = dplyr::first(Start),
        End = dplyr::last(End),
        Length = dplyr::first(Length),
        Strand = unique(Strand)
      ) %>% 
      dplyr::mutate(
        Chr = stringr::str_c("chr", Chr),
        coord = dplyr::if_else(Strand == "+", End, Start)) %>% 
      tidyr::unite(col = "peakID", c("Chr", "coord", "Strand"), sep = ":") %>% 
      dplyr::select(GeneID, peakID)
    
    counts <- peak_file %>% 
      readRDS(.) %>% 
      .@cells.counts %>% 
      tibble::remove_rownames() %>% 
      tibble::column_to_rownames("Peak_ID") %>% 
      dplyr::rename_with(~gsub(".*_", "", .x)) %>% 
      dplyr::rename_with(~gsub("-1", "", .x)) %>% 
      tibble::rownames_to_column("GeneID") %>%
      dplyr::left_join(scAPA_bed, by = "GeneID")
    
  }

  if(method == "MAAPER"){
    peak_file <- list.files(expr_file, pattern = "pas.txt", full.names = T, recursive = T)
    pas_name <- vroom::vroom(peak_file) %>% dplyr::pull(pas)
    ## create a pseudo counts matrix
    counts <- matrix(
      sample(1:10, length(barcodes) * length(pas_name), replace = TRUE), 
      nrow = length(pas_name), 
      dimnames = list(pas_name, barcodes)
    ) %>% 
      dplyr::as_tibble(.) %>% 
      dplyr::mutate(peakID = pas_name)
    
  }

  if(method == "DarPars2"){
    
    file_list <- list.files(expr_file, pattern = "^Dapars2_result_temp.*\\.txt$", recursive = TRUE, full.names = TRUE)
    
    pas_name <- lapply(file_list, function(darpars_output){
      
      pa_site <- vroom::vroom(darpars_output) %>% 
        tidyr::separate(col = Gene, into = c("esembl_id", "symbol", "chr", "strand"), sep = "\\|") %>% 
        tidyr::separate(col = Loci, into = c("chr", "region"), sep = ":") %>% 
        tidyr::separate(col = region, into = c("start", "end"), sep = "-") %>% 
        dplyr::mutate(
          distal_sites = as.integer(dplyr::if_else(strand == "+", end, start)),
          Predicted_Proximal_APA = as.integer(Predicted_Proximal_APA)) %>% 
        tidyr::pivot_longer(col = c("Predicted_Proximal_APA", "distal_sites"), names_to = "type", values_to = "coord") %>% 
        tidyr::unite(col = "pas_name", c("chr", "coord", "strand"), sep = ":") %>% 
        dplyr::distinct(pas_name)
      
      return(pa_site)
    }) %>% 
      dplyr::bind_rows() %>% 
      dplyr::pull(pas_name)
    
    ## create a pseudo counts matrix
    counts <- matrix(
      sample(1:10, length(barcodes) * length(pas_name), replace = TRUE), 
      nrow = length(pas_name), 
      dimnames = list(pas_name, barcodes)
    ) %>% 
      dplyr::as_tibble(.) %>% 
      dplyr::mutate(peakID = pas_name)
  }

  ## add annotation
  counts <- counts %>% 
    dplyr::left_join(annot, by = c("peakID" = "pa_site")) %>% 
    dplyr::filter(!is.na(seqnames))
  
  
  return(counts)
}

# Create annotation database for each genomic feature based on the gtf file.
annotate_from_gtf <- function (
    gtf_file,
    cores = 1L
    ) {
  
  # Create TxDb obj from gene annotation file
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf") 
  
  # Get gene symbol from gtf file
  # Returns dataframe includes 3 columns: 'transcript_id', 'gene_id', 'gene_name'
  id_maps <- rtracklayer::import(gtf_file) |>
    as.data.frame(row.names = NULL) |> 
    dplyr::filter(type == "transcript") |> 
    dplyr::select(gene_id, gene_name, transcript_id) |> 
    dplyr::rename(symbol = gene_name, ensembl = gene_id)
  
  ## ==================== DEFINING REGION ====================
  # Define a function to convert gr list into obj and add meta info, such as tx_name, gene_id, symbol
  add_meta_data <- function(gr_list, id_maps) {
    # Filter out the empty elemets in GRanges object
    gr_list <- gr_list[S4Vectors::elementNROWS(gr_list) > 0]
    
    # Check if filtered grlist is empty
    if (length(gr_list) == 0) {
      stop("All elements in gr_list are empty.")
    }
    
    # Convert gr list into gr obj (set use.names = TRUE to extract transcript id later)
    gr_obj <- unlist(gr_list, use.names = TRUE)
    
    # Check if the length of names is consistent with the length of gr obj
    if (length(names(gr_obj)) != length(gr_obj)) {
      stop("Mismatch between names and GRanges length after unlist. Possibly due to malformed input.")
    }
    
    # Add transcript ID in metadata
    GenomicRanges::mcols(gr_obj)$tx_name <- names(gr_obj)
    # Add gene_id and symbol
    # Add Entrez ID, symbol, and type
    GenomicRanges::mcols(gr_obj)$gene_id <- id_maps[match(GenomicRanges::mcols(gr_obj)$tx_name, id_maps$transcript_id),'ensembl']
    GenomicRanges::mcols(gr_obj)$symbol <- id_maps[match(GenomicRanges::mcols(gr_obj)$gene_id, id_maps$ensembl), 'symbol']
    
    return(gr_obj)
  }
  
  ## ==================
  ## ====== EXON ======
  # Boundary of exon grouping by transcript ID. Return GRangeList object.
  # each entry's key is transcript ID and their values are exons and their boundary positions.
  exons_grl = GenomicFeatures::exonsBy(txdb, by = 'tx', use.names = TRUE)
  exons_gr <- add_meta_data(gr_list = exons_grl, id_maps = id_maps)
  GenomicRanges::mcols(exons_gr)$id = paste0('Exon:ExonRank', exons_gr$exon_rank)
  
  ## =======================
  ## ====== LAST EXON ======
  # Find the last exon
  lastexons_grl <- parallel::mclapply(
    exons_grl[S4Vectors::elementNROWS(exons_grl) != 0], # Remove empty elements
    function(x) {
      if (length(x) == 0) stop("All elements in gr_list are empty.") # Safety Check
      # Get the index of element
      idx <- if (as.character(strand(x)[1]) == "+") length(x) else 1 # Same gene, same strand, so use the first element to simplify
      x[idx]
    },
    mc.cores = cores
  ) 
  # Convert list to gr list and to gr obj with meta info
  lastexons_grl <- as(lastexons_grl, "GRangesList") # If use GRangesList(x), the speed will be slow and the memory will be small. 
  lastexons_gr <- add_meta_data(gr_list = lastexons_grl, id_maps = id_maps)
  GenomicRanges::mcols(lastexons_gr)$id = paste0('LastExon:ExonRank', lastexons_gr$exon_rank)
  
  ## ========================================
  ## ====== LAST EXON WITHIN 1K WINDOW ======
  pos = lastexons_gr[strand(lastexons_gr) == '+', ]
  neg = lastexons_gr[strand(lastexons_gr) == '-', ]
  start(pos) <- end(pos) + 1
  end(pos) <- end(pos) + 1000
  end(neg) <- start(neg) - 1
  start(neg) <- start(neg) - 1000
  
  lastexons1k_gr_ <- c(pos, neg)
  GenomicRanges::mcols(lastexons1k_gr_)$id = paste0('LastExon1Kb:ExonRank', lastexons_gr$exon_rank)
  
  ## =============================
  ## ====== CODING SEQUENCE ======
  cds_grl = GenomicFeatures::cdsBy(txdb, by = 'tx', use.names = TRUE)
  cds_gr <- add_meta_data(gr_list = cds_grl, id_maps = id_maps)
  GenomicRanges::mcols(cds_gr)$id = paste0('CDS:ExonRank', cds_gr$exon_rank)
  
  ## ====================
  ## ====== 5' UTR ======
  fiveUTRs_grl = GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
  fiveUTRs_gr <- add_meta_data(gr_list = fiveUTRs_grl, id_maps = id_maps)
  GenomicRanges::mcols(fiveUTRs_gr)$id = paste0('5UTR:ExonRank', fiveUTRs_gr$exon_rank)
  
  ## ====================
  ## ====== 3' UTR ======
  threeUTRs_grl = GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)
  threeUTRs_gr <- add_meta_data(gr_list = threeUTRs_grl, id_maps = id_maps)
  GenomicRanges::mcols(threeUTRs_gr)$id = sprintf('3UTRs')
  
  ## ======================================
  ## ====== 3' UTR within 1kb window ======
  pos = threeUTRs_gr[strand(threeUTRs_gr) == '+', ]
  neg = threeUTRs_gr[strand(threeUTRs_gr) == '-', ]
  
  start(pos) <- end(pos) + 1
  end(pos) <- end(pos) + 1000
  
  end(neg) <- start(neg) - 1
  start(neg) <- start(neg) - 1000
  
  threeUTRs1Kb <- c(neg, pos)
  threeUTRs1Kb$id <- '3UTRs_1kb'
  
  ## ======================================
  ## ====== 3' UTR within 2kb window ======
  ###for 2kb
  pos = threeUTRs1Kb[strand(threeUTRs1Kb) == '+', ]
  neg = threeUTRs1Kb[strand(threeUTRs1Kb) == '-', ]
  
  start(pos) <- end(pos) + 1
  end(pos) <- end(pos) + 1000
  
  end(neg) <- start(neg) - 1
  start(neg) <- start(neg) - 1000
  
  threeUTRs2Kb <- c(pos, neg)
  threeUTRs2Kb$id <- '3UTRs_2kb'
  
  ## ====================
  ## ====== INTRON ======
  introns_grl = GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)
  introns_grl <- parallel::mclapply(
    introns_grl[S4Vectors::elementNROWS(introns_grl) != 0], 
    function(x) {
      # Check if the gr is empty
      if (length(x) == 0) stop("All elements in gr_list are empty.")
      # Generate the rank number
      ids <- if (as.character(strand(x)[1]) == "+") seq(1:length(x)) else rev(seq(1:length(x)))
      x$id <- paste0('Intron:Rank', ids)
      x
    }, 
    mc.cores = cores
  )
  # Convert list to gr list and then to gr obj adding meta info
  introns_grl <- as(introns_grl, "GRangesList")
  introns_gr <- add_meta_data(gr_list = introns_grl, id_maps = id_maps)
  
  ## ========================
  ## ====== INTERGENIC ======
  genic_gr <- GenomicFeatures::genes(txdb, single.strand.genes.only=FALSE)
  genic_gr <- unlist(genic_gr)
  GenomicRanges::strand(genic_gr) = '*'
  intergenic_gr = GenomicRanges::gaps(genic_gr)
  intergenic_gr = intergenic_gr[GenomicRanges::strand(intergenic_gr) == '*']
  
  GenomicRanges::mcols(intergenic_gr) <- S4Vectors::DataFrame(
    tx_name = NA_character_,
    gene_id = NA_character_,
    symbol = NA_character_,
    id = "Intergenic"
  )
  
  ## ====================
  ## ====== RESULT ======
  use_ind <- c("tx_name", "gene_id", "symbol", "id") # Keep metadata consistent
  
  annotation_db <- c(
    cds_gr[, use_ind],
    threeUTRs_gr[, use_ind],
    threeUTRs1Kb[, use_ind],
    threeUTRs2Kb[, use_ind],
    fiveUTRs_gr[, use_ind],
    intergenic_gr[, use_ind],
    introns_gr[, use_ind],
    exons_gr[, use_ind],
    lastexons1k_gr_[, use_ind]
  )
  
  return(annotation_db)
}

# Create a function to annotate sites by database
annotation_site <- function(
    pAsite,
    gtf_file,
    minoverlap=1L,
    annot_levels = c(
      "3UTRs",
      "Exon",
      "Intron",
      "3UTRs_1kb",
      "3UTRs_2kb",
      "5UTR", ##tien change order
      "CDS",
      "LastExon1Kb",
      "INTERGENIC"
    ),
    cores = 1) {
  #     """
  #     Input:
  #         - pAsite: list of pAsite (chr:position:strand). Example: c(chr1:6275401:+, chr1:6275770:+)
  #     Return:
  #         - annot_res : dataframe of annotation, including seqnames, pa_pos, strand, annot.start, annot.end, annot.tx_name, annot.gene_id, annot.symbol, Type, Rank, distance, group_flag
  #     """
  
  # Convert pAsite into dataframe (tibble) include columns: chr, start, end, strand
  pa_info <-  data.frame(pAsite) |>
    tidyr::separate(
      col = pAsite, 
      into = c("chr", "pa_pos", "strand"),
      sep = ":"
    ) |>
    dplyr::mutate(
      pa_pos = as.integer(pa_pos), 
      start = dplyr::if_else(strand == "+", pa_pos - 1L, pa_pos),
      end = dplyr::if_else(strand == "+", pa_pos, pa_pos + 1L)
    )
  
  # Convert dataframe to GRangers object
  pa_info <- GenomicRanges::makeGRangesFromDataFrame(pa_info, keep.extra.columns = TRUE)
  cat("Done making Granges\n")
  # path directory
  gtf_prefix <- strsplit(basename(gtf_file), split = '[.]')[[1]][1]
  annot_db_file <- file.path(dirname(gtf_file), paste0(gtf_prefix, '_site_annotation.qs')) 
  
  # prepare GRanges of region of interest from GTF (the same as what is used in calculating PAsite)
  if (file.exists(annot_db_file)) {
    annot_db <- qs::qread(annot_db_file, nthreads = cores)
    cat("Read existing annotation database", "\n")
  } else {
    cat("Creating annotation database ")
    start.time <- Sys.time()
    annot_db <- annotate_from_gtf(gtf_file = gtf_file, cores = cores)
    qs::qsave(annot_db, file = annot_db_file, nthreads = cores)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    cat("in ", time.taken, "\n")
  }
  
  # Search each entry of pa_info against GRanges of all regions
  # return Granges object
  annot_res <- annotatr::annotate_regions(
    regions = pa_info,
    annotations = annot_db,
    minoverlap = minoverlap,
    ignore.strand = FALSE,
    quiet = FALSE
  ) # seqnames should have chr, so the example they give have problems.
  cat("Done annotating regions \n")

  inds <- c(
    "seqnames", 
    "pa_pos", 
    "strand", 
    "annot.start", 
    "annot.end", 
    "annot.tx_name", 
    "annot.gene_id", 
    "annot.symbol", 
    "annot.id"
    )
  
  # Separate annot.id into type and rank, and then order by factor
  names(annot_res@elementMetadata@listData$annot) <- NULL # fix row.names error 2025.5.27
  annot_res <- annot_res |> 
    as.data.frame(row.names = NULL) |> 
    tibble::as_tibble() |> 
    dplyr::select(dplyr::all_of(inds)) |>
    tidyr::separate(col = annot.id, into = c("Type", "Rank"), sep = ":", fill = "right") |>
    dplyr::mutate(
      Type = factor(Type, levels = annot_levels),
    )
  cat("Select annotation!")
  inds <- colnames(annot_res)
  # Strategy 1: evaluate the order of annotation by type and distance
  annot_res <- annot_res |>
    dplyr::mutate(
      distance = dplyr::case_when(
        strand == "+" ~ abs(annot.end - pa_pos),
        strand == "-" ~ abs(pa_pos - annot.start),
        TRUE ~ NA_real_
      ),
      peakID = glue::glue("{seqnames}:{pa_pos}:{strand}") 
    ) |>
    dplyr::group_by(peakID) |>
    dplyr::arrange(Type, distance, .by_group = TRUE) |> 
    dplyr::mutate(
      group_flag = if(dplyr::n() > 1) {
        cumsum(c(TRUE, diff(as.integer(Type)) != 0 | diff(distance) != 0))
      } else {
        1
      }
    ) |>
    dplyr::group_by(peakID, group_flag) |>
    dplyr::summarise(
      annot.tx_name = paste(unique(annot.tx_name), collapse = ", "),
      Type = dplyr::first(Type),
      distance = dplyr::first(distance),
      dplyr::across(!c("annot.tx_name", "Type", "distance"), ~paste(unique(.), collapse = ", ")),
      .groups = "drop_last"
    ) |>
    dplyr::group_by(peakID) |>
    dplyr::slice(1)|>
    dplyr::ungroup() |>
    dplyr::select(dplyr::all_of(c(inds, "distance", "group_flag")))
  
  return(annot_res)
}

# annotate pa sites by SCAPE
apa_annotation <- function(apa_dir, method="polyApipe", gtf_file, barcode_file, core_num=8){
  require(magrittr)
  require(GenomicRanges)
  Sys.setenv(VROOM_CONNECTION_SIZE = 500072)
  
  if(method == "SCAPE"){
    peak_file <- list.files(apa_dir, pattern = "pasite.csv.gz", full.names = T, recursive = T)
    pas_name <- vroom::vroom(peak_file, col_select = 1, col_names = c("pas_site"))
    pas_name <- pas_name %>% 
      dplyr::filter(!is.na(pas_site)) %>% 
      tidyr::separate(col = "pas_site", into = c("chr", "coord", "score", "strand"), sep = ":") %>% 
      tidyr::unite(col = "pas_site", c("chr", "coord", "strand"), sep = ":") %>% 
      dplyr::pull(pas_site)
  }
  
  if(method == "Sierra"){
    peak_file <- list.files(apa_dir, pattern = "sitenames.tsv.gz", full.names = T, recursive = T)
    pas_name <- vroom::vroom(peak_file, col_names = c("gene_name", "chr", "region", "strand"))
    pas_name <- pas_name %>% 
      dplyr::mutate(strand = dplyr::if_else(strand=="1", "+", "-")) %>% 
      tidyr::separate(col="region", into=c("start", "end"), sep = "-") %>% 
      dplyr::mutate(sitename = dplyr::if_else(strand=="+", end, start)) %>% 
      tidyr::unite(col = "peakID", c("chr", "sitename", "strand"), sep = ":") %>% 
      dplyr::pull(peakID)
  }
  
  if(method == "scAPAtrap"){
    peak_file <- list.files(apa_dir, pattern = "expma.qs", full.names = T, recursive = T)
    expma <- qs::qread(peak_file, nthreads=8)
    pas_name <- expma %>% 
      dplyr::select(1:6) %>% 
      tidyr::unite(col = "peakID", c("chr", "coord", "strand"), sep = ":") %>% 
      dplyr::pull(peakID)
  }
  
  if(method == "MAAPER"){
    peak_file <- list.files(apa_dir, pattern = "pas.txt", full.names = T, recursive = T)
    pas_name <- vroom::vroom(peak_file) %>% dplyr::pull(pas)
    
  }
  
  if(method == "SCAPTURE"){
    peak_file <- list.files(apa_dir, pattern = "PASquant.KeepPAS.bed", full.names = T, recursive = T)
    pas_name <- vroom::vroom(peak_file, col_names = F)
    pas_name <- pas_name %>% 
      dplyr::mutate(pa_site = dplyr::if_else(X6 == "+", X3, X2)) %>% 
      dplyr::mutate(pas_name = stringr::str_c(X1, pa_site, X6, sep = ":")) %>% 
      dplyr::pull(pas_name)
  }
  
  if(method == "polyApipe"){
    peak_file <- list.files(apa_dir, pattern = "counts.tab.gz", full.names = T, recursive = T)
    barcodes <- readr::read_tsv(barcode_file, col_names = FALSE, show_col_types = FALSE)$X1 %>% gsub("-1", "", .)
    
    pas_name <- vroom::vroom(peak_file) %>% 
      dplyr::filter(cell %in% barcodes) %>% 
      dplyr::group_by(gene) %>% 
      dplyr::summarise(across(count, sum)) %>% 
      dplyr::filter(count > 0) %>% 
      tidyr::separate(col = gene, into = c("chr", "coord", "strand"), sep = "_") %>% 
      dplyr::mutate(strand = recode(strand, "r" = "-", "f" = "+")) %>% 
      tidyr::unite(col = "pas_name", c("chr", "coord", "strand"), sep = ":") %>% 
      dplyr::pull(pas_name)
  }
  
  if(method == "scAPA"){
    require(scAPA)
    peak_file <- list.files(apa_dir, pattern = "Peaks.RDS", full.names = T, recursive = T)
    
    pas_name <- peak_file %>% 
      readRDS(.) %>% 
      .@row.Data %>%
      tidyr::separate_rows(Chr, Start, End, Strand, sep = ";") %>%
      dplyr::group_by(GeneID) %>%
      dplyr::summarise(
        Chr = unique(Chr),
        Start = dplyr::first(Start),
        End = dplyr::last(End),
        Length = dplyr::first(Length),
        Strand = unique(Strand)
      ) %>% 
      dplyr::mutate(
        Chr = stringr::str_c("chr", Chr),
        coord = dplyr::if_else(Strand == "+", End, Start)) %>% 
      tidyr::unite(col = "pas_name", c("Chr", "coord", "Strand"), sep = ":") %>% 
      dplyr::pull(pas_name)
  }

  if(method == "DarPars2"){
    file_list <- list.files(apa_dir, pattern = "^Dapars2_result_temp.*\\.txt$", recursive = TRUE, full.names = TRUE)
    
    pas_name <- lapply(file_list, function(darpars_output){
      pa_site <- vroom::vroom(darpars_output) %>% 
        tidyr::separate(col = Gene, into = c("esembl_id", "symbol", "chr", "strand"), sep = "\\|") %>% 
        tidyr::separate(col = Loci, into = c("chr", "region"), sep = ":") %>% 
        tidyr::separate(col = region, into = c("start", "end"), sep = "-") %>% 
        dplyr::mutate(
          distal_sites = as.integer(dplyr::if_else(strand == "+", end, start)),
          Predicted_Proximal_APA = as.integer(Predicted_Proximal_APA)) %>% 
        tidyr::pivot_longer(col = c("Predicted_Proximal_APA", "distal_sites"), names_to = "type", values_to = "coord") %>% 
        tidyr::unite(col = "pas_name", c("chr", "coord", "strand"), sep = ":") %>% 
        dplyr::distinct(pas_name)

      return(pa_site)
    }) %>% 
      dplyr::bind_rows() %>% 
      dplyr::pull(pas_name)
    
  }

  annot_info <- annotation_site(
      pas_name,
      gtf_file,
      cores = core_num
      ) 
  
  return(annot_info)
}

