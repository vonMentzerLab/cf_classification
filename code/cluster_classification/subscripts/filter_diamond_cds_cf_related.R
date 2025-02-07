library(tidyverse)
library(Biostrings)
library(doParallel)
library(foreach)
library(yaml)

load_fa_seqs <- function(fa_file){
  Biostrings::getSeq(Rsamtools::FaFile(fa_file))
}

myalign <- function(seqs, seq2){
  pairwiseAlignment(seqs, seq2, substitutionMatrix = 'BLOSUM62', scoreOnly = TRUE, gapOpening = 10, gapExtension = 4)
}

myalign2 <- function(seqs, seq2){
  submat = 'BLOSUM62'
  gapopen = 10
  gapext = 4

  alignment <- pairwiseAlignment(seqs, seq2, substitutionMatrix = submat, scoreOnly = FALSE, gapOpening = gapopen, gapExtension = gapext)
  tibble(
    seq1_name = names(seqs),
    seq2_name= names(seq2),
    score = score(alignment),
    max_score = pairwiseAlignment(seq2, seq2, substitutionMatrix = submat, scoreOnly = TRUE, gapOpening = gapopen, gapExtension = gapext),
    pid = pid(alignment))
}



bacterial_gen_code <- getGeneticCode('11')




args <- commandArgs(trailingOnly = TRUE)
config <- read_yaml(args[1])

config$diamond_prot_dir <- paste0(config$gene_search_dir, '/diamond_aligned_amino_cf_db/')
config$cf_filtered_cds_coords_dir  <- paste0(config$gene_search_dir, '/cf_cds_coords/')
config$global_alignment_scores_dir  <- paste0(config$gene_search_dir, '/cds_global_alignment_scores/')

config$cds_coords_dir <- paste0(config$cds_features_dir, '/cds_coords/')
config$nucl_seqs_dir <- paste0(config$cds_features_dir, '/cds_nucl_seqs/')

system(paste('mkdir -p', config$cf_filtered_cds_coords_dir, config$global_alignment_scores_dir))

files <- list.files(config$diamond_prot_dir)

CF_reference_aa_seqs <- load_fa_seqs(config$clustered_genes_reference_db_nucl_seq) %>% 
  translate(genetic.code = bacterial_gen_code)



cores <- config$cores
cl <- makeCluster(cores)
registerDoParallel(cl)
silence <- foreach(file = files, .packages = c('tidyverse', 'Biostrings')) %dopar% {


  dmd <- read_tsv(paste0(config$diamond_prot_dir, '/', file), col_names = FALSE, col_types = 'ccdiiiid')
  names(dmd) <- c('qseqid', 'sseqid', 'pident', 'length', 'slen', 'qlen', 'gaps', 'evalue')
  dmd <- dmd %>% mutate(cov = 100 * (length - gaps) / slen)

  # Filter cf diamond hits with minimum pident and coverage
  # Currently filtering away toxins, but we should just separate the CF and toxin databases
  dmd_qseqid_filtered <- dmd %>% 
    filter(
      cov >= config$min_cov,
      pident >= config$min_pident) %>% 
    select(qseqid) %>% 
    unique()

  # foreach "continue" if there are no CDSs remaining after filtering
  if(nrow(dmd_qseqid_filtered) == 0){
    return(NULL)
  }

  # Get genome_id so that we can load the cds features file
  genome_id <- dmd_qseqid_filtered %>% 
    separate(qseqid, into = c('genome_id', 'cds_id'), sep = '\\.') %>% 
    pull(genome_id) %>% 
    unique()
  cds_features <- read_tsv(paste0(config$cds_coords_dir, '/cds_features_', genome_id, '.tsv.gz'), col_types = 'ciiiciccc')

  # Get first start and last end of diamond hit genes
  cf_cds_start_end <- cds_features %>% 
    filter(my_cds_id %in% dmd_qseqid_filtered$qseqid) %>% 
    group_by(seqnames) %>% 
    summarise(first_start = min(start), last_end = max(end))


  # Get all genes within a maximum distance of the diamond hit genes (in case it missed something, or we want to look at surrounding genes)
  cds_filtered <- cds_features %>% 
    inner_join(cf_cds_start_end) %>% 
    left_join(dmd_qseqid_filtered %>% dplyr::rename(my_cds_id = qseqid) %>% mutate(diamond_cf_candidate = TRUE)) %>% 
    replace_na(list(diamond_cf_candidate = FALSE)) %>% 
    group_by(seqnames) %>% 
    filter(first_start - end <= config$max_dist & start - last_end <= config$max_dist) %>% 
    ungroup() %>% 
    select(-first_start, -last_end)
  cds_filtered %>% write_tsv(paste0(config$cf_filtered_cds_coords_dir, "/cds_around_cf_", genome_id, '.tsv.gz'))

  # Get sequences of relevant cds
  nucl_seqs <- load_fa_seqs(paste0(config$nucl_seqs_dir, '/cds_nucl_seqs_', genome_id, '.fa.gz'))
  nucl_seqs <- nucl_seqs[cds_filtered$my_cds_id]
  aa_seqs <- nucl_seqs %>% translate(genetic.code = bacterial_gen_code, if.fuzzy.codon = 'solve')

  # Perform global alignments against CF ref db
  alignment_scores <- bind_rows(map(1:length(aa_seqs), function(x) myalign2(CF_reference_aa_seqs, aa_seqs[x]))) %>% 
    select(my_cds_id = seq2_name, ref_cf = seq1_name, score, max_score, pid)

  alignment_scores %>% 
    write_tsv(paste0(config$global_alignment_scores_dir, '/global_alignment_scores_cds_around_cf_', genome_id, '.tsv.gz'))

  NULL
}

cat('Done!\n')
