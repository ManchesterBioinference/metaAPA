# funtion to decide how many clusters to look at
calculate_cluster_metrics <- function(matrix, n_centers, cluster_method = "kmeans", n_starts = 25) {
  if(cluster_method == "kmeans") {
    kmeans_result <- kmeans(matrix, centers = n_centers, nstart = n_starts)
    wss <- kmeans_result$tot.withinss
    silhouette_scores <- mean(cluster::silhouette(kmeans_result$cluster, dist(matrix))[, 3])
  }
  if(cluster_method == "pam") {
    # https://stackoverflow.com/questions/38306259/compute-within-sum-of-squares-from-pam-cluster-analysis-in-r
    pam_result <- cluster::pam(matrix, k = n_centers, diss = TRUE)
    wss <- pam_result$objective[2] # not within cluster sum of square.
    silhouette_scores <- pam_result$silinfo$avg.width
  }
  
  return(c(n_centers, wss, silhouette_scores))
}


# function to extract pA sites
extract_pA_sites <- function(apa_count_folder, apa_method, cell_ids, gene_symbol, core_num = 4) {
  filename_pattern <- glue::glue("{apa_method}_count_matrix.qs")
  count_file <- list.files(apa_count_folder, pattern = filename_pattern, full.names = T, recursive = T) 
  pA_sites <- glue::glue("{count_file}") %>% 
    qs::qread(., nthreads = core_num) %>% 
    dplyr::filter(annot.symbol == gene_symbol) %>% 
    dplyr::select(all_of(c("peakID", cell_ids))) %>% 
    dplyr::mutate(peakID = glue::glue("{apa_method}:{peakID}"))
  
  return(pA_sites)
}

# function to calculate distance between sites based on different algorithms
calculate_distance <- function(merge_sites, dist_method = "spearman") {
  
  dist_matrix <- NULL
  
  if (dist_method == "spearman") {
    # calculate spearman correlation
    dist_matrix <- merge_sites %>%
      tibble::column_to_rownames("peakID") %>%
      t() %>%
      cor(., method = "spearman")
    dist_matrix <- 1 - dist_matrix
  } else if (dist_method == "pearson") {
    # calculate pearson correlation
    dist_matrix <- merge_sites %>%
      tibble::column_to_rownames("peakID") %>%
      t() %>%
      cor(., method = "pearson")  
    dist_matrix <- 1 - dist_matrix
  } else if (dist_method == "cosine") {
    # calculate cosine distance
    dist_matrix <- merge_sites %>%
      tibble::column_to_rownames("peakID") %>%
      proxy::dist(., method = "cosine") %>% 
      as.matrix()
  } else if (dist_method == "Tanimoto") {
    # calculate Tanimoto distance
    dist_matrix <- merge_sites %>%
      tibble::column_to_rownames("peakID") %>%
      proxy::simil(., method = "Tanimoto") %>% 
      as.matrix()
    dist_matrix <- 1 - dist_matrix
  } else if (dist_method == "Jaccard") {
    # calculate Jaccard distance
    # binary_matrix <- merge_sites %>% tibble::column_to_rownames("peakID")
    # binary_matrix <- binary_matrix > 0
    # binary_matrix <- binary_matrix * 1 
    # dist_matrix <- binary_matrix %>%
    #   proxy::dist(., method = "Jaccard") %>% 
    #   as.matrix()
    dist_matrix <- merge_sites %>%
      tibble::column_to_rownames("peakID") %>%
      proxy::dist(., method = "Jaccard") %>% 
      as.matrix()
  } else if (dist_method == "none") {
    # raw counts
    dist_matrix <- merge_sites %>% tibble::column_to_rownames("peakID")
    
  } else if (dist_method == "norm") {
    dist_matrix <- merge_sites %>% tibble::column_to_rownames("peakID")
    dist_matrix <- sqrt((dist_matrix**2)/rowSums((dist_matrix**2)))
  } else if (dist_method == "euclidean") {
    dist_matrix <- merge_sites %>%
      tibble::column_to_rownames("peakID") %>%
      proxy::dist(., method = "Euclidean") %>% 
      as.matrix()
  } else {
    cat("Distance method can not be recognized!")
    
  }
  
  return(dist_matrix)
}

# find the elbow point
find_elbow_point <- function(x, y) {
  # refer to Jonas anwser in the stackoverflow: https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
  # Calculate the coefficients of the line (y = mx + b) based on the start point and the end point of elbow plot
  fit <- lm(y[c(1, length(y))] ~ x[c(1, length(x))])
  m <- coef(fit)[2]
  b <- coef(fit)[1]
  
  # Calculate perpendicular distance from a point to the line for all points
  distances <-  abs(m * x - y + b) / sqrt(m^2 + 1)
  
  # Find the index of the maximum distance
  max_distance_index <- which.max(distances)
  
  return(max_distance_index)
  
}

# find optimal center based on asw or wss for Kmeans
find_optimal_center <- function(center_metrics, center_method) {
  # get the optimal centers based silhouette or wss --- maybe need to a new metrics or new cluster method
  if (center_method == "asw") {
    optimal_centers <- center_metrics[which.max(center_metrics$silhouette), ]$centers
  } else if (center_method == "wss") {
    optimal_centers_wss_index <- find_elbow_point(center_metrics$centers, center_metrics$wss)
    optimal_centers <- center_metrics[optimal_centers_wss_index, ]$centers
  } else if (center_method == "mix") {
    wss_max <- max(center_metrics$wss)
    center_metrics <- center_metrics %>% dplyr::mutate(mix = silhouette * wss_max - wss)
    optimal_centers <- center_metrics[which.max(center_metrics$mix),]$centers
  } else {
    message("Unknown method to identify optimal center!")
    return(NULL)
  }
  
  return(optimal_centers)
}

# align pA sites
align_pA_sites <- function(
    gene_symbol, 
    cell_ids, 
    apa_method_list = c("Sierra", "SCAPE", "polyApipe"), 
    dist_method = "spearman",
    cluster_method = "kmeans-wss",
    apa_count_folder = "",
    count_cutoff = 0,
    n_starts = 25,
    cell_pct_threshold = 0,
    core_num = 4) {
  
  require(dplyr)
  
  # Check if the argument is valid
  if (!dist_method %in% c("spearman", "cosine", "Jaccard", "pearson", "none", "norm", "euclidean")) {
    stop("Invalid choice. Please select 'spearman', 'cosine', 'pearson', 'none', 'norm', 'euclidean' or 'Jaccard'.")
  }
  if (!cluster_method %in% c("hclust", "pam-wss", "pam-asw", "kmeans-wss", "kmeans-asw", "hdbscan", "kmeans-sparse", "zinbmm")) {
    stop("Invalid choice. Please select 'hclust', 'pam-wss', 'pam-asw' 'kmeans-wss', 'kmeans-asw', 'hdbscan', 'kmeans-sparse', or 'zinbmm'.")
  }

  # The expression matrix for the given genes and methods needs to be extracted. 
  merge_sites <- lapply(
    apa_method_list, 
    extract_pA_sites, 
    cell_ids =  cell_ids, 
    gene_symbol = gene_symbol,
    apa_count_folder = apa_count_folder,
    core_num = core_num) %>% 
    dplyr::bind_rows() %>%
    dplyr::group_by(peakID) %>% 
    dplyr::summarise(across(all_of(cell_ids), sum), .groups = 'drop')
  
  if (nrow(merge_sites) == 0) {
    message(glue::glue("No sites in {gene_symbol}!"))
    return(NULL)
  } else {
    # Sites less than cutoff need to be removed.
    merge_sites <- merge_sites %>%
      dplyr::filter(rowSums(select(., !peakID)) > count_cutoff) %>% 
      dplyr::distinct()
    
    if (nrow(merge_sites) == 0) {
      message(glue::glue("No sites in {gene_symbol}!"))
      return(NULL)
    } 
  }
  
  # Sites expressed less than cell_pct_threshold need to be removed.
  if(cell_pct_threshold > 0) {
    merge_sites <- merge_sites %>% 
      dplyr::mutate(
      cell_pct = rowSums(dplyr::select(., !peakID) != 0)/ncol(dplyr::select(., !peakID))
      ) %>% 
      dplyr::filter(cell_pct >= cell_pct_threshold) %>% 
      dplyr::select(all_of(c("peakID", cell_ids)))
  }
  
  # # The number of unique methods in the data will decide how to deal with it.
  # methods_in_id <- merge_sites %>% 
  #   dplyr::select(peakID) %>% 
  #   tidyr::separate(peakID, into = c("method", "chr", "site", "strand"), sep = ":") %>% 
  #   dplyr::distinct(method) %>% 
  #   dplyr::pull(method)
  # 
  # method_condition <- all(apa_method_list %in% methods_in_id) # may be tested later
  
  # task: divide data into two group based on cell percentage, over 50% to calculate
  
  dist_matrix <- calculate_distance(merge_sites, dist_method = dist_method)
  n_duplicate_rows <- dist_matrix %>% duplicated() %>% sum()
  
  if (n_duplicate_rows == 0) {
    message("There are no duplicate rows.")
  } else {
    message(glue::glue("There are {n_duplicate_rows} duplicated rows."))
    dist_matrix <- dist_matrix %>%
      as.data.frame() %>% 
      tibble::rownames_to_column("peakID") %>% 
      group_by(across(-peakID)) %>%
      summarise(peakID = paste(peakID, collapse = "_"), .groups = 'drop') %>% 
      tibble::column_to_rownames("peakID")
  }
  
  if (nrow(dist_matrix) > length(apa_method_list)) {
    message(glue::glue("No site > No methods!"))
    
    if (cluster_method == "kmeans-wss") {
      # calculate asw and wss
      n_clusters <- nrow(dist_matrix) - 1
      
      cluster_metrics <- lapply(
        2:n_clusters, 
        calculate_cluster_metrics, 
        matrix = dist_matrix, 
        cluster_method = "kmeans",
        n_starts = n_starts
      ) %>% 
        purrr::reduce(rbind) %>% #Reduce(rbind, .)
        dplyr::as_tibble(.name_repair = "unique") %>% 
        setNames(c("centers", "wss", "silhouette"))
      
      # find optimal center
      optimal_centers <- find_optimal_center(cluster_metrics, "wss")
      
      message(glue::glue("The optimal centers of {gene_symbol}: {optimal_centers}"))
      # k means using the optimal center
      kmeans_result <- kmeans(dist_matrix, centers = optimal_centers, nstart = n_starts)
      labels <- kmeans_result$cluster
    } else if (cluster_method == "kmeans-asw") {
      # calculate asw and wss
      n_clusters <- nrow(dist_matrix) - 1
      
      cluster_metrics <- lapply(
        2:n_clusters, 
        calculate_cluster_metrics, 
        matrix = dist_matrix,
        cluster_method = "kmeans",
        n_starts = n_starts
      ) %>% 
        purrr::reduce(rbind) %>% #Reduce(rbind, .)
        dplyr::as_tibble(.name_repair = "unique") %>% 
        setNames(c("centers", "wss", "silhouette"))
      
      # find optimal center
      optimal_centers <- find_optimal_center(cluster_metrics, "asw")
      
      message(glue::glue("The optimal centers of {gene_symbol}: {optimal_centers}"))
      # k means using the optimal center
      kmeans_result <- kmeans(dist_matrix, centers = optimal_centers, nstart = n_starts)
      labels <- kmeans_result$cluster
    } else if (cluster_method == "pam-wss") {
      # calculate asw and wss
      n_clusters <- nrow(dist_matrix) - 1
      
      cluster_metrics <- lapply(
        2:n_clusters, 
        calculate_cluster_metrics, 
        matrix = dist_matrix, 
        cluster_method = "pam",
        n_starts = n_starts
      ) %>% 
        purrr::reduce(rbind) %>% #Reduce(rbind, .)
        dplyr::as_tibble(.name_repair = "unique") %>% 
        setNames(c("centers", "wss", "silhouette"))
      
      # find optimal center
      optimal_centers <- find_optimal_center(cluster_metrics, "wss")
      
      message(glue::glue("The optimal centers of {gene_symbol}: {optimal_centers}"))
      # k means using the optimal center
      pam_result <- cluster::pam(dist_matrix, k = optimal_centers, diss = TRUE)
      labels <- pam_result$clustering
    } else if (cluster_method == "pam-asw") {
      # calculate asw and wss
      n_clusters <- nrow(dist_matrix) - 1
      
      cluster_metrics <- lapply(
        2:n_clusters, 
        calculate_cluster_metrics, 
        matrix = dist_matrix, 
        cluster_method = "pam",
        n_starts = n_starts
      ) %>% 
        purrr::reduce(rbind) %>% #Reduce(rbind, .)
        dplyr::as_tibble(.name_repair = "unique") %>% 
        setNames(c("centers", "wss", "silhouette"))
      
      # find optimal center
      optimal_centers <- find_optimal_center(cluster_metrics, "asw")
      
      message(glue::glue("The optimal centers of {gene_symbol}: {optimal_centers}"))
      # k means using the optimal center
      pam_result <- cluster::pam(dist_matrix, k = optimal_centers, diss = TRUE)
      labels <- pam_result$clustering
    } else if (cluster_method == "hdbscan") {
      hdbscan_result <- dbscan::hdbscan(dist_matrix, minPts = 2) 
      labels <- hdbscan_result$cluster
      names(labels) <- rownames(dist_matrix)
      cluster_metrics <- NULL
    } else if (cluster_method == "hclust") {
      n_clusters <- nrow(dist_matrix) - 1
      hc_results <- hclust(as.dist(dist_matrix), method = "ward.D2")
      cluster_metrics <- purrr::map_dfr(2:n_clusters, function(n_cluster){
        labels <- cutree(hc_results, k = n_cluster)
        metrics <- c(n_cluster, all(table(labels) <= 3), sum(table(labels) == 3))
        names(metrics) <- c("n_cluster", "is_normal_size", "n_high_conf")
        metrics
      })
      
      optimal_centers <- cluster_metrics %>% 
        dplyr::filter(is_normal_size == 1) %>%
        dplyr::filter(n_high_conf == max(n_high_conf)) %>% 
        dplyr::pull(n_cluster) %>% 
        dplyr::first()
      
      labels <- cutree(hc_results, k = optimal_centers)
    } else if (cluster_method == "kmeans-sparse") {
      # choose tuning parameter and optimal K based on Gap statics
      express_mat <- merge_sites %>% tibble::column_to_rownames("peakID") %>% as.matrix()
      n_clusters <- nrow(dist_matrix) - 1
      km_perm_list <- purrr::map_dfr(
        2:n_clusters,
        function(k) {
          km_perm <- sparcl::KMeansSparseCluster.permute(
            express_mat, 
            K = k,
            wbounds = 2:15,
            nperms = 5
          )
          index <- which(km_perm$wbounds==km_perm$bestw)
          gap_statistic <- c(
            bestw = km_perm$bestw,
            ncluster = k,
            nonzero = km_perm$nnonzerows[index],
            bestgap = km_perm$gaps[index],
            bestsd = km_perm$sdgaps[index]
          )
          return(gap_statistic)
        }
      )
      
      bestw <- km_perm_list[which.max(km_perm_list$bestgap),]$bestw
      optimal_centers <- km_perm_list[which.max(km_perm_list$bestgap),]$ncluster
      
      # run sparse k means with the optimal parameters
      km_out <- sparcl::KMeansSparseCluster(
        express_mat,
        K = optimal_centers,
        wbounds = bestw
      )
      
      labels <- km_out[[1]]$Cs
      cluster_metrics <- NULL
    } else if (cluster_method == "zinbmm") {
      source("rscripts/zinbmm.R")
      require(parallel)
      require(stats)
      require(Brobdingnag)
      express_mat <- merge_sites %>% tibble::column_to_rownames("peakID") %>% as.data.frame() %>% as.matrix()
      batch_site <- merge_sites %>%
        dplyr::select(peakID) %>%
        tidyr::separate(peakID, into = c("method", "chr", "site", "strand"), sep = ":") %>%
        dplyr::mutate(
          v1 = as.integer(method == "SCAPE"),
          v2 = as.integer(method == "Sierra"),
          v3 = as.integer(method == "polyApipe")
        ) %>%
        dplyr::select(v1, v2, v3) %>%
        as.matrix()
      
      n_clusters <- nrow(dist_matrix) - 1
      cluster_metrics_list <- parallel::mclapply(
        2:n_clusters,
        function(k){
          results = zinbmm(
            express_mat, 
            batch_site, 
            K=k, 
            ntune=8, 
            tol=1e-02, 
            maxit = 500, 
            ncores=core_num, 
            vrbs=F)
        },
        mc.cores = 1
      )
      
      cluster_metrics <- tibble(
        "k" = 2:n_clusters,
        "modified_bic" = purrr::map_vec(cluster_metrics_list, purrr::pluck, "bic_modified") 
      )
      
      optimal_centers <- cluster_metrics[which.max(cluster_metrics$modified_bic),]$k
      results = zinbmm(
        express_mat, 
        batch_site, 
        K=optimal_centers, 
        ntune=8, 
        tol=1e-02, 
        maxit = 500, 
        ncores=core_num, 
        vrbs=F)
      
      labels <- results$clust
    } else {
      message("Please use right cluster method.")
    }
    
  } else {
    message(glue::glue("No site <= No methods!"))
    labels <- rep("0", nrow(dist_matrix))
    names(labels) <- rownames(dist_matrix)
    cluster_metrics <- NULL
  }
  
  # site labels
  if(nrow(dist_matrix) == 0) {
  labels <- 0
  } else {
  labels <- labels[rownames(dist_matrix)]
  }

  site_labels <- dist_matrix %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("peakID") %>% 
    dplyr::mutate(labels = as.character(labels)) %>% 
    dplyr::select(peakID, labels) %>% 
    tidyr::separate_rows(peakID, sep = "_")
  
  cluster_results <- merge_sites %>% 
    tibble::column_to_rownames("peakID") %>% 
    dplyr::mutate(total = rowSums(.)) %>% 
    tibble::rownames_to_column("peakID") %>% 
    dplyr::select(peakID, total) %>% 
    dplyr::left_join(site_labels, by = "peakID") %>% 
    tidyr::separate(col = "peakID", into = c("method", "chr", "site", "strand"), sep = ":") %>% 
    dplyr::mutate(gene_name = gene_symbol)

  return(list(cluster_results, dist_matrix, cluster_metrics))
}

# calculate base frequency
calculate_base_frequency <- function(feature_detect_data, upstream_width = 50, downstream_width = 50) {
  require(Biostrings)
  require(dplyr)
  # extract DNA sequences
  sequences <- extract_DNA_seq(feature_detect_data, upstream_width = upstream_width, downstream_width = downstream_width)
  # calculate the frequency of base in each position
  afmc <- consensusMatrix(sequences, baseOnly = T, as.prob = T)
  afmc_longer <- afmc %>% 
    t() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("position") %>% 
    tidyr::pivot_longer(cols = !position, names_to = "base", values_to = "percentage") %>% 
    dplyr::mutate(position = as.integer(position)) %>% 
    dplyr::filter(base != "other") # only detect ATCG
  
  return(afmc_longer)
}

# extract DNA sequence
extract_DNA_seq <- function(feature_detect_data, upstream_width = 50, downstream_width = 50) {
  require(Biostrings)
  require(GenomicRanges)
  require(dplyr)
  # load mm10 genome
  genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  
  upstream_width <- as.integer(upstream_width)
  downstream_width <- as.integer(downstream_width)
  
  feature_detect_data <- feature_detect_data %>% 
    dplyr::mutate(site = as.integer(site)) %>% 
    dplyr::mutate(
      start_pos = dplyr::if_else(strand == "+", site - upstream_width, site - downstream_width), 
      end_pos = dplyr::if_else(strand == "+", site + downstream_width, site + upstream_width), 
    )
  # extract position
  chromosome <- feature_detect_data$chr
  strand <- feature_detect_data$strand
  start_pos <- feature_detect_data$start_pos
  end_pos <- feature_detect_data$end_pos
  
  # create GRange object
  region <- GRanges(
    seqnames = chromosome,
    ranges = IRanges(start = start_pos, end = end_pos),
    strand = strand
  )
  
  # get sequence
  sequences <- BSgenome::getSeq(genome, region)
  names(sequences) <- feature_detect_data$site
  
  return(sequences)
}

# convert png to ggplot
convert_png_to_ggplot <- function(imgurl) {
  # load packages
  require(ggplot2)
  require(png)
  require(grid)
  
  # set theme
  image_plot_theme <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
  
  # read png
  png_image <- readPNG(imgurl)
  
  # insert png to a blank plot
  image_plot <- ggplot(d = data.frame(x = c(0, 12), y = c(0, 12)), aes(x, y)) + 
    geom_blank() + # create a blank plot
    annotation_custom(rasterGrob(png_image)) + # insert png to blank plot
    theme_minimal() + 
    image_plot_theme
  
  return(image_plot)
}

# generate GR object based on cluster results
generate_GR <- function(cluster_results, nt_extend_left = 50, nt_extend_right = 50) {
  start_site <- min(cluster_results$site) - nt_extend_left
  end_site <- max(cluster_results$site) + nt_extend_right
  chrom <- unique(cluster_results$chr)
  strand <- unique(cluster_results$strand)
  region <- GenomicRanges::GRanges(seqnames = chrom, strand = strand, ranges = IRanges::IRanges(start = start_site, end = end_site))
  
  return(region)
}

# draw site alignment plot
plot_site_alignment <- function(cluster_results, nt_extend_left = 50, nt_extend_right = 50) {
  require(dplyr)
  require(ggplot2)
  # set range
  left_x <- min(cluster_results$site) - nt_extend_left
  right_x <- max(cluster_results$site) + nt_extend_right
  # plot
  site_align_plot <- cluster_results %>% 
    ggplot(aes(x = site, y = method, col = labels)) + 
    geom_point(size = 3, alpha = 0.7) +
    geom_text(aes(label = total), vjust = -1, hjust = 0.5, size = 3) +
    guides(color = guide_legend(override.aes = list(label = ""))) +
    labs(x = "", y = "", col = "Cluster") +
    coord_cartesian(xlim = c(left_x, right_x)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), panel.grid.minor = element_blank())
  
  return(site_align_plot)
}

# draw transcripts
plot_transcript_structure <- function(gene_of_interest, region, gtf_path) {
  require(ggtranscript)
  require(dplyr)
  require(ggplot2)
  # set genomic range
  left_x <- GenomicRanges::start(region)
  right_x <- GenomicRanges::end(region)
  # strand <- GenomicRanges::strand(region) |> as.character()
  # target_chrom <- GenomicRanges::seqnames(region) |> as.character()
  # filter your gtf for the gene of interest and extract the required annotation columns
  gene_annotation_from_gtf <- rtracklayer::import(gtf_path) %>% 
    dplyr::as_tibble() %>% 
    dplyr::filter(!is.na(gene_name)) %>% 
    dplyr::filter(gene_name == gene_of_interest) %>% 
    dplyr::select(seqnames, start, end, strand, type, gene_name, transcript_name, transcript_type) 
  
  # extract exons
  gene_exons <- gene_annotation_from_gtf %>% dplyr::filter(type == "exon")
  
  # obtain cds
  gene_cds <- gene_annotation_from_gtf %>% dplyr::filter(type == "CDS")
  
  structure_plot <- gene_exons %>%
    ggplot(aes(xstart = start, xend = end, y = transcript_name)) +
    geom_range(fill = "white", height = 0.25) +
    geom_range(data = gene_cds) +
    geom_intron(
      data = to_intron(gene_exons, "transcript_name"),
      aes(strand = strand),
      arrow.min.intron.length = 500,
    ) + 
    coord_cartesian(xlim = c(left_x, right_x)) +
    labs(x="Genomic position",y="") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())
  
  return(structure_plot)
}

# draw coverage plot
plot_coverage <- function(region, bamfile) {
  require(Rsamtools)
  require(ggplot2)
  # extract site info
  strand <- GenomicRanges::strand(region) |> as.character()
  target_chrom <- GenomicRanges::seqnames(region) |> as.character()
  left_x <- GenomicRanges::start(region)
  right_x <- GenomicRanges::end(region)
  # extract the coverage of specific region
  is_minus_strand <- FALSE
  
  if(strand == "-"){
    is_minus_strand <- TRUE
  }
  
  param <- ScanBamParam(which = region, flag = scanBamFlag(isUnmappedQuery = FALSE, isMinusStrand = is_minus_strand)) # save memory
  
  bam_data <- BamFile(bamfile)
  cov_region <- GenomicAlignments::coverage(bam_data, param = param)[[target_chrom]]
  
  # plot
  coverage_plot <- data.frame(Position = left_x:right_x , Coverage = as.numeric(cov_region[left_x:right_x])) %>% 
    ggplot(., aes(x = Position, y = Coverage)) +
    # geom_bar(stat = "identity") +
    geom_area(fill = "lightblue", alpha = 0.7) +
    # geom_line(col = "grey") +
    labs(x = "", y = "", title = "") +
    coord_cartesian(xlim = c(left_x, right_x)) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(coverage_plot)
}

# draw nucleotide composition
plot_nucleotide_composition <- function(
    cluster_data, 
    apa_method = "polyApipe", 
    consistent_site_num = 3, 
    upstream = 100, 
    downstream = 50
    ) {
  require(dplyr)
  require(Biostrings)
  require(GenomicRanges)
  
  if (consistent_site_num %in% c(1, 2, 3)) {
    cluster_data <- cluster_data %>% 
      dplyr::filter(method == apa_method) %>%
      dplyr::mutate(site = as.integer(site)) %>% 
      dplyr::filter(consistent_sites == TRUE & n_site == consistent_site_num)
  } else {
    cluster_data <- cluster_data %>% 
      dplyr::filter(method == apa_method) %>%
      dplyr::mutate(site = as.integer(site)) %>% 
      dplyr::filter(consistent_sites == FALSE)
  }
  
  site_num <- nrow(cluster_data)
  
  # calculate base frequency
  afmc_longer <- calculate_base_frequency(cluster_data, upstream_width = upstream, downstream_width = downstream)
  
  # plot
  composition_plot <- afmc_longer %>% 
    ggplot(aes(x = position - (as.integer(upstream) + 1), y = percentage, col = base)) + 
    geom_line() + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x= "Position", y = "Percentage", col = "Base") +
    ggtitle(glue::glue("{apa_method}: {site_num} sites (n = {consistent_site_num})")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(composition_plot)
}

# draw upset plot
plot_gene_upset <- function(gene_group_list, metadata, sample_name){
  ### this function generate upset plot for genes based on UpSet package
  ### then convert upset base plot into ggplot object based on cowplot package
  ### input:
  ###   - gene_group_list: list. includes gene list from different methods
  ###   - metadata: dataframe. includes sets, Methods, Color
  
  require(cowplot)
  require(UpSetR)
  # Create a deep copy including independent set column (somehow c(), unserialize(serialize()), data.table::copy didn't work).
  # Upset plot needs two independent columns to color the plot.
  metadata <- data.frame(
    "set" = metadata$tools,
    "color" = metadata$color,
    "tools" = metadata$tools
  )
  tool_colors <- setNames(metadata$color, metadata$tools)
  
  p_upset <- UpSetR::upset(
    UpSetR::fromList(gene_group_list), 
    sets = metadata[["set"]],
    keep.order = T,
    order.by = "freq",
    sets.bar.color = metadata[["color"]],
    set.metadata = list(
      data = metadata, 
      plots = list(
        list(type = "matrix_rows", column = "tools", colors = tool_colors, alpha = 0.8)
      )
    )
  )
  
  # Convert base into ggplot: the first solution https://github.com/hms-dbmi/UpSetR/issues/105
  p_overview <- cowplot::plot_grid(
    NULL, NULL, p_upset$Main_bar, 
    NULL, NULL, NULL, 
    p_upset$Sizes, NULL, p_upset$Matrix,
    nrow=3, ncol = 3, align='hv', 
    rel_heights = c(3, -0.3, 1), rel_widths = c(1, -0.1, 3))
  
  # now add the title
  title <- ggdraw() + 
    draw_label(
      sample_name,
      fontface = 'bold',
      # x = 0,
      hjust = 0.5
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 0)
    )
  
  p_overview <- cowplot::plot_grid(
    title, p_overview,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )
  
  return(p_overview)
}

# get legend
get_legend <-function(gplot_obj){
  gtable_obj <- ggplot_gtable(ggplot_build(gplot_obj))
  leg <- which(sapply(gtable_obj$grobs, function(x) x$name) == "guide-box")
  legend <- gtable_obj$grobs[[leg]]
  return(legend)
}
