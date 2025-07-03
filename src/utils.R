gene_copynumber <- function(cn, gene){
  
  message(paste0("Finding states for gene ", gene))
  
  if (!gene %in% signals::gene_locations$hg19$ensembl_gene_symbol){
    stop(paste0(gene, " not in list of genes..."))
  }
  
  message(paste0("Number of cells: ", length(unique(cn$cell_id))))
  
  gene_g <- signals::gene_locations$hg19 %>%
    dplyr::mutate(seqnames = chr) %>%
    dplyr::filter(!stringr::str_detect(seqnames, "_")) %>%
    dplyr::filter(ensembl_gene_symbol == gene) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::mutate(seqnames = chr) %>%
    plyranges::as_granges()
  cn_g <- cn %>% dplyr::rename(seqnames = chr) %>% plyranges::as_granges()
  
  if (!levels(GenomicRanges::seqnames(gene_g)) %in% levels(GenomicRanges::seqnames(cn_g))){
    stop("Chromosome not in data")
  }
  
  overlap <- plyranges::join_overlap_inner(cn_g, gene_g)
  
  if (length(overlap$cell_id) == 0){
    message("No bins overlapping, finding nearest bin...")
    overlap <- plyranges::join_nearest_left(cn_g, gene_g)
  }
  
  gene_g <- as.data.frame(gene_g)
  overlap <- as.data.frame(overlap) %>%
    dplyr::select(-chr) %>%
    dplyr::rename(chr = seqnames) %>%
    dplyr::mutate(start_gene = gene_g$start, end_gene = gene_g$end) %>%
    dplyr::mutate(diff = abs(start_gene - start)) %>%
    dplyr::arrange(diff) %>%
    dplyr::group_by(cell_id) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-diff)
  
  message(paste0("Number of cells with data: ", length(unique(overlap$cell_id))))
  
  if (length(unique(overlap$cell_id)) > length(unique(cn$cell_id))){
    stop("More regions than cells returned")
  }
  
  return(overlap)
}
