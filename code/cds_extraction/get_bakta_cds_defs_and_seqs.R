library(tidyverse)
library(plyranges)
library(BSgenome)
library(foreach)
library(doParallel)

load_fa_seqs <- function(fa_file){
  Biostrings::getSeq(Rsamtools::FaFile(fa_file))
}

make_bakta_contignames <- function(dnastringset){
  names(dnastringset) <- paste0('contig_', 1:length(dnastringset))
  return(dnastringset)
}

args <- commandArgs(trailingOnly = TRUE)
in_dir <- args[1]
out_dir <- args[2]

nucl_dir <- paste0(out_dir, '/cds_nucl_seqs/')
prot_dir <- paste0(out_dir, '/cds_amino_seqs/')
cds_coords_dir <- paste0(out_dir, '/cds_coords/')

system(paste('mkdir -p', nucl_dir, prot_dir, cds_coords_dir))

genome_ids <- tibble(file = list.files(in_dir)) %>% 
  filter(str_detect(file, '.gff3')) %>% 
  mutate(genome_id = str_replace(file, config$gff3_suffix, '')) %>% 
  pull(genome_id) %>% 
  unique()



cores  <- 30
cl <- makeCluster(cores)
registerDoParallel(cl)


silly_null <- foreach(
  genome_id = genome_ids,
  .inorder = FALSE,
  .packages = c('tidyverse', 'plyranges', 'BSgenome')) %dopar% {

  # Get Bakta CDS features
  cds <- import.gff3(paste0(in_dir, '/', genome_id, config$gff3_suffix)) %>% 
    filter(type == 'CDS') %>% 
    select(phase, ID) %>% 
    mutate(genome_id = genome_id, my_cds_id = paste0(genome_id, '.', ID))

  # The phase column indicates where in a codon the feature starts.
  # Probably 0 in all of them, but other values would be disastrous if not accounted for.
  if(any(cds$phase != 0)){
    cds <- cds %>% 
      anchor_3p() %>% 
      mutate(width = width - phase)
  }

  contig_seqs <- load_fa_seqs(paste0(in_dir, '/', genome_id, config$fa_suffix)) %>% make_bakta_contignames()

  bacterial_gen_code <- getGeneticCode('11')
  nucl_seqs <- BSgenome::getSeq(contig_seqs, cds)
  names(nucl_seqs) <- cds$my_cds_id
  aa_seqs <- translate(nucl_seqs, genetic.code = bacterial_gen_code, if.fuzzy.codon = 'solve')

  # Use compression here if we handle a lot of data. Quite slow though.
  nucl_seqs %>% writeXStringSet(paste0(nucl_dir, '/cds_nucl_seqs_', genome_id, '.fa.gz'), compress = TRUE)
  aa_seqs %>% writeXStringSet(paste0(prot_dir, '/cds_amino_seqs_', genome_id, '.fa.gz'), compress = TRUE)
  cds %>% as_tibble() %>% write_tsv(paste0(cds_coords_dir, '/cds_features_', genome_id, '.tsv.gz'))
  
  NULL # Seems to be impossible to skip outputting something, so let's at least make it something small.
}
