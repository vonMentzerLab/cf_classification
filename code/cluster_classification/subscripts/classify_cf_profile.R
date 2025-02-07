library(tidyverse)
library(cowplot) # This is unnecessary when we don't call the plot function, which isn't done atm.
library(doParallel)
library(foreach)
library(yaml)

read_cf_rules <- function(cf_rules_path){
  read_tsv(cf_rules_path, col_types = 'cc') %>% 
    group_by(cf) %>% 
    mutate(cf_gene_order = row_number()) %>% 
    ungroup() %>% 
    mutate(gene = ifelse(str_detect(gene, ':'), gene, paste0(gene, ':', gene))) %>% 
    separate(gene, sep = ':', into = c('equiv_gene', 'gene')) %>% 
    mutate(gene = str_split(gene, ',')) %>% 
    unnest(gene)
}

get_cf_required_gene_counts <- function(cf_rules){
  cf_rules %>%
    filter(required) %>% 
    select(cf, equiv_gene) %>% 
    unique() %>%
    count(cf, name = 'req_gene_count_in_cf_rules')
}

find_prob_cf_coords <- function(cf_coords){
  operon_zone <- cf_coords %>%
    filter(diamond_cf_candidate) %>% 
    mutate(first_cf_in_operon = (start - lag(end)) > config$max_dist) %>% 
    replace_na(list(first_cf_in_operon = TRUE)) %>% 
    mutate(operon_id = paste0(seqnames, '.', cumsum(first_cf_in_operon))) %>% 
    group_by(operon_id, seqnames) %>% 
    summarise(operon_zone_start = min(start), operon_zone_end = max(end), .groups = 'drop') %>% 
    rowwise() %>% 
    mutate(operon_zone_start = max(operon_zone_start - config$max_dist, 1), operon_zone_end = operon_zone_end + config$max_dist) %>% 
    ungroup()

  # OBS! This definition allows for genes that are not diamond cf candidates to belong to two operons, if the distances are just right
  # This is so that we can capture the surroundings for multiple operons, even if they're close.
  operon_cf_coords <- bind_rows(
    map(
      unique(operon_zone$operon_id),
      function(x) inner_join(cf_coords, filter(operon_zone, operon_id == x), by = 'seqnames') %>% filter(start >= operon_zone_start, end <= operon_zone_end))
  )

  # We currently require probable cf genes (diamond alignment-based) at the ends of the of prob_cf_coords operons
  # This may mean we lose some interesting things, so we should probably go 1 out or so at some point.
  prob_cf_coords <- operon_cf_coords %>% 
    group_by(seqnames, operon_id, diamond_cf_candidate) %>% 
    mutate(
      first_cf_gene_start = ifelse(diamond_cf_candidate & start == min(start), start, 0),
      last_cf_gene_end = ifelse(diamond_cf_candidate & end == max(end), end, 0))  %>% 
    group_by(seqnames, operon_id) %>% 
    mutate(
      first_cf_gene_start = sum(first_cf_gene_start),
      last_cf_gene_end = sum(last_cf_gene_end)
    ) %>% 
    ungroup() %>% 
    mutate(prob_cf_operon_gene = start >= first_cf_gene_start & end <= last_cf_gene_end) %>% 
    filter(prob_cf_operon_gene)

  operon_strands <- prob_cf_coords %>% 
    count(operon_id, strand) %>% 
    arrange(operon_id, strand) %>% 
    group_by(operon_id) %>% 
    slice_max(order_by = n, with_ties = FALSE) %>% # default to plus strand if equal numbers 
    ungroup() %>% 
    select(operon_id, operon_strand = strand)

  prob_cf_coords <- prob_cf_coords %>% 
    left_join(operon_strands, by = 'operon_id') %>% 
    group_by(operon_id) %>% 
    mutate(gene_order = ifelse(operon_strand == '-', dplyr::n()-row_number()+1, row_number()))  %>% 
    ungroup()

  prob_cf_coords

}

make_cf_score_df <- function(prob_cf_coords){
  prob_cf_coords %>% 
    left_join(cf_align_scores, by = 'my_cds_id') %>% 
    mutate(score_fraction = score / max_score) %>% 
    left_join(cf_names %>% dplyr::rename(ref_cf = name, gene = short_name), by = 'ref_cf') %>% 
    left_join(cf_rules, by = 'gene')
}





classify_cfs <- function(cf_score_df, cf_classify_min_score_fraction){

  genome_id <- unique(cf_score_df$genome_id)

  df_max_per_cf <- cf_score_df %>% 
      group_by(operon_id, cf, gene_order) %>%
      slice_max(order_by = score_fraction, with_ties = FALSE) %>% 
      ungroup()

  operon_gene_count <- cf_score_df %>%
    group_by(operon_id) %>%
    slice_max(gene_order, with_ties = FALSE) %>%
    ungroup() %>%
    select(operon_id, operon_gene_count = gene_order)

  # Proper gene presence as opposed to just counting
  df_all_genes_present <- df_max_per_cf  %>%
    filter(
      score_fraction >= cf_classify_min_score_fraction,
      required) %>%
    select(operon_id, cf, equiv_gene) %>% 
    unique() %>% 
    group_by(operon_id, cf) %>%
    arrange(operon_id, cf, equiv_gene) %>% 
    summarise(cf_present_genes =  paste0(equiv_gene, collapse = ' + '), .groups = 'drop') %>% 
    left_join(
      cf_rules %>%
      filter(required) %>% 
      select(cf, equiv_gene) %>% 
      unique() %>% 
      group_by(cf) %>%
      arrange(cf, equiv_gene) %>% 
      summarise(cf_required_genes = paste0(equiv_gene, collapse = ' + ')),
      by = 'cf') %>% 
    mutate(all_genes_present = cf_present_genes == cf_required_genes) %>% 
    select(operon_id, cf, all_genes_present)

# cf_req_gene_counts probably not necessary to keep anymore after adding df_all_genes_present. Delete?
df_cf_gene_count <- df_max_per_cf %>% 
    filter(
      score_fraction >= cf_classify_min_score_fraction) %>%
    count(operon_id, cf, name = 'gene_count') %>% 
    bind_rows(df_max_per_cf %>% select(operon_id, cf) %>% mutate(gene_count = 0)) %>% 
    group_by(operon_id, cf) %>% 
    summarise(gene_count = sum(gene_count), .groups = 'drop') %>% 
    left_join(cf_req_gene_counts, by = 'cf')

  # Currently "ordered" means that:
  # 1: there can't be any other genes in the middle of the gene cluster, even if the CF genes themselves are in the right order
  # 2: For each CF, we check if a stretch of genes in a cluster match what we're supposed to see. Other genes before or after first/last gene matching a CF are fine - this is mostly to allow missing regulators etc
  # 3: only required genes are taken into account.
  df_cf_in_order <- df_max_per_cf %>%
    group_by(operon_id, cf) %>% 
    filter(required) %>% 
    mutate(equiv_gene = ifelse(score_fraction >= cf_classify_min_score_fraction, equiv_gene, 'too_low_score')) %>% 
    summarise(obs_gene_order = paste0(equiv_gene, collapse = '_')) %>%
    ungroup() %>% 
    mutate(
      obs_gene_order_trimmed = str_replace(obs_gene_order, '^(too_low_score_)+', ''),
      obs_gene_order_trimmed = str_replace(obs_gene_order_trimmed, '(_too_low_score)+$', ''),) %>% 
    left_join(
      cf_rules %>%
        filter(required) %>%
        select(cf, equiv_gene, cf_gene_order) %>%
        unique() %>%
        group_by(cf) %>%
        arrange(cf, cf_gene_order) %>%
        summarise(req_gene_order = paste0(equiv_gene, collapse = '_')),
        by = 'cf') %>% 
    mutate(genes_in_order = str_detect(req_gene_order, obs_gene_order_trimmed)) %>%
    select(operon_id, cf, genes_in_order)

  # OBS! We're calculating the mean/min score fraction of all the genes in the detected cluster, and not caring
  # about whether there are any missing required genes here. We could add those to drag the score of incomplete
  # operons down, but this would have pros and cons. A negative aspect would be penalising an excellent match
  # to half of the required genes in favour of a mediocre match to all genes in another CF.
  # Note that CF genes that are not required are included in this if they're detected.
  df_score_fraction_summary <- df_max_per_cf %>% 
    select(operon_id, cf, gene, score_fraction) %>% 
    group_by(operon_id, cf) %>% 
    summarise(mean_score_fraction = mean(score_fraction), .groups = 'drop')

  df_summary <- df_cf_gene_count %>% 
    left_join(df_cf_in_order, by = c('operon_id', 'cf')) %>% 
    left_join(df_score_fraction_summary, by = c('operon_id', 'cf')) %>% 
    left_join(df_all_genes_present, by = c('operon_id', 'cf')) %>% 
    replace_na(list(genes_in_order = FALSE, all_genes_present = FALSE))


  complete_cfs <- df_summary %>% 
    filter(all_genes_present, genes_in_order) %>% 
    group_by(operon_id) %>% 
    slice_max(order_by = mean_score_fraction) %>% # Ties?
    ungroup()

  complete_cf_best_hits <- df_max_per_cf %>% 
    inner_join(complete_cfs %>% select(operon_id, cf), by = c('operon_id', 'cf')) %>% 
    select(seqnames, start, end, strand,  my_cds_id, operon_id, operon_strand, genome_id, cf, gene, score_fraction, gene_order) %>% 
    mutate(classification_type = 'complete')


  incomplete_cfs <- df_summary %>% 
    filter(!(operon_id %in% complete_cfs$operon_id)) %>% 
    group_by(operon_id) %>% 
    slice_max(order_by = mean_score_fraction) %>% 
    ungroup()

  initial_cf_classification <- complete_cfs %>% 
    bind_rows(incomplete_cfs) %>% 
    mutate(
      all_genes_and_ordered = all_genes_present & genes_in_order,
      genome_id = genome_id)


  # We need to remember which operons didn't have any proteins that scored over the threshold in order to be able to have "Unknown" only classifications below
  no_match_cfs <- initial_cf_classification %>%
    anti_join(
      df_max_per_cf %>%
        filter(score_fraction >= cf_classify_min_score_fraction) %>% 
        select(operon_id) %>% 
        unique(),
      by = 'operon_id') 

  ##################################
  # Attempt at classifying hybrids
  #################################

  if(nrow(initial_cf_classification) > 0){
    # Take best-matching cf from incomplete_cfs, and distinguish between genes above and below cf_classify_min_score_fraction
    hybrid_cf_tmp <- initial_cf_classification %>% 
      dplyr::rename(best_operon_cf_match = cf) %>% 
      left_join(df_max_per_cf, by = c('operon_id', 'genome_id')) %>% 
      rowwise() %>% 
      mutate(cf_matching_gene = (score_fraction >= !!cf_classify_min_score_fraction && best_operon_cf_match == cf)) %>% # !! used to make sure we're referencing global variable, as we previously had a column with the same name (though it's renamed now)
      ungroup()
    # Genes to keep matching "main" cf
    hybrid_cf_best_cf_matching_genes <- hybrid_cf_tmp %>%
      filter(cf_matching_gene) 
    # Select best gene from all cfs for genes that we didn't keep above
    hybrid_cf_best_cf_no_pass_matching <- hybrid_cf_tmp %>% 
      anti_join(
        hybrid_cf_best_cf_matching_genes %>% 
          select(operon_id, gene_order),
        by = c('operon_id', 'gene_order')) %>% 
      group_by(operon_id, gene_order) %>% 
      slice_max(order_by = score_fraction, with_ties = FALSE) 

    # Are all genes in the hybrid cf above cf_classify_min_score_fraction?
    hybrid_cf_genes <- hybrid_cf_best_cf_matching_genes %>% 
      bind_rows(hybrid_cf_best_cf_no_pass_matching) %>% 
      arrange(operon_id, gene_order) %>% 
      group_by(operon_id) %>% 
      mutate(hybrid_min_score_qualifying = all(score_fraction >= cf_classify_min_score_fraction)) %>% 
      ungroup()

    # Make classification for both hybrids where all genes pass cf_classify_min_score_fraction, and for those that don't
    hybrid_cfs <- bind_rows(
      hybrid_cf_genes %>%
        filter(hybrid_min_score_qualifying) %>%
        group_by(operon_id, best_operon_cf_match) %>% 
        summarise(
          cf_classification = paste0(unique(cf), collapse = ' / '),
          gene_profile = paste0(gene, collapse = ' + '),
          .groups = 'drop'),
      hybrid_cf_genes %>% 
        filter(!hybrid_min_score_qualifying) %>% 
        group_by(operon_id, best_operon_cf_match) %>% 
        filter(score_fraction >= cf_classify_min_score_fraction) %>% 
        summarise(cf_classification = paste0(c(unique(cf), 'Unknown'), collapse = ' / '), .groups = 'drop') %>% 
        bind_rows(
          no_match_cfs %>%
            select(operon_id, best_operon_cf_match = cf) %>% 
            mutate(cf_classification = 'Unknown')
        ) %>% 
        left_join(hybrid_cf_genes %>% 
          filter(!hybrid_min_score_qualifying) %>% 
          mutate(gene = ifelse(score_fraction >= cf_classify_min_score_fraction, gene, 'Unknown')) %>% 
          group_by(operon_id) %>% 
          summarise(gene_profile = paste0(gene, collapse = ' + ')), by = 'operon_id')
      
    ) %>% mutate(genome_id = genome_id)


    cf_classification <- hybrid_cfs %>% 
      left_join(hybrid_cf_genes %>% group_by(operon_id) %>% summarise(mean_score_fraction = mean(score_fraction)), by = c('operon_id')) %>% 
      left_join(initial_cf_classification %>%
        dplyr::rename(
          best_operon_cf_match = cf,
          best_operon_cf_match_gene_count = gene_count,
          best_operon_cf_match_all_genes_present = all_genes_present,
          best_operon_cf_match_genes_in_order = genes_in_order,
          best_operon_cf_match_mean_score_fraction = mean_score_fraction,
          best_operon_cf_match_all_genes_and_ordered = all_genes_and_ordered),
          by = c('operon_id', 'genome_id', 'best_operon_cf_match')) %>% 
        left_join(operon_gene_count, by = 'operon_id') %>% 
        relocate(
          genome_id,
          operon_id,
          cf_classification,
          gene_profile,
          mean_score_fraction,
          best_operon_cf_match,
          best_operon_cf_match_mean_score_fraction,
          operon_gene_count,
          best_operon_cf_match_gene_count,
          req_gene_count_in_cf_rules,
          best_operon_cf_match_all_genes_present,
          best_operon_cf_match_genes_in_order,
          best_operon_cf_match_all_genes_and_ordered)

      classifications_by_gene <- hybrid_cf_genes %>% 
        select(seqnames, start, end, strand, genome_id, operon_id, operon_strand, my_cds_id, best_operon_cf_match, cf, gene, score_fraction, score, max_score, pid, gene_order)

      operon_coords <- classifications_by_gene %>% 
        select(-strand) %>% 
        group_by(genome_id, operon_id, seqnames, strand = operon_strand) %>% 
        summarise(start = min(start), end = max(end), .groups = 'drop')

      cf_classification <- cf_classification %>% left_join(operon_coords, by = c('genome_id', 'operon_id'))

  } else {
    cf_classification <- tibble()
    classifications_by_gene <- tibble()
  }


  best_matching_gene_any_cf <- cf_score_df %>%
    group_by(operon_id, gene_order) %>%
    slice_max(score_fraction, with_ties = FALSE) %>% 
    ungroup() %>% 
    select(seqnames, start, end, strand , genome_id, operon_id, operon_strand, my_cds_id, cf, gene, score_fraction, score, max_score, pid, gene_order)

  ##################
  # Return list of tibbles. 

  return(list(cf_classification = cf_classification, classifications_by_gene = classifications_by_gene, best_matching_gene_any_cf = best_matching_gene_any_cf))

}



plot_cf_matches <- function(cf_score_df, min_score_fraction){
  genome_id <- df$genome_id %>% unique()
  operon_ids <- unique(df$operon_id)
  plist <- list()
  # Visualise possible CF matches

  genes_not_in_cf_rules <- df %>% filter(is.na(cf)) %>% pull(gene) %>% unique()

  tmp <- df %>%
    group_by(operon_id, cf) %>% 
    mutate(max_score_fraction_for_this_cf = max(score_fraction)) %>% 
    ungroup() %>%
    filter(max_score_fraction_for_this_cf >= min_score_fraction) %>% 
    mutate(
      gene_order = factor(gene_order), # Easy way to avoid axis text displaying e.g. 2.5 instead of individually numbering the genes
      gene = factor(gene, levels = c(cf_rules$gene, genes_not_in_cf_rules))) %>% # We have to use cf_rules to set levels to order the genes in each operon, but this leaves out some genes that are belong to a particular cf, but are not required. Adding these back in here to at least show them in a separate facet (which will be NA)
    rowwise() %>% # Format med digits funkar inte vettigt utan rowwise - blir fel antal
    mutate(label = ifelse(score_fraction > 0.1, as.character(format(score_fraction, scientific = FALSE, digits = 2)), '')) %>% 
    mutate(score_fraction = ifelse(score_fraction < 0, 0, score_fraction))  %>%  # set floor of score_fraction to 0 to get consistent colours between subplots. We care about things that are similar, not exactly how dissimilar they are - 0 is plenty dissimilar.
    ungroup()

  gene_counts <- df %>%
    group_by(factor(operon_id, levels = operon_ids)) %>% # Måste ha samma ordning som plottarna kommer i sen.
    summarise(gene_count = max(gene_order)) %>% 
    pull(gene_count)

  for(i in 1:length(operon_ids)){
    operon_id_iter <- operon_ids[i]

    p <- tmp %>% 
      filter(operon_id == operon_id_iter) %>% 
      ggplot(aes(x = gene_order, y = gene, fill = score_fraction)) +
      geom_tile() +
      geom_text(aes(label = label)) +
      scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c('<= 0', '0.5', '1')) +
      facet_wrap(~cf, ncol = 1, scales = 'free', strip.position = 'right') +
      ggtitle(paste0(operon_id_iter, ' (', genome_id,')')) +
      theme(
        aspect.ratio = 1,
        legend.position = 'bottom',
        legend.title = element_blank()
      )

    plist[[i]] <- p
  }

  plot_grid(plotlist = plist, ncol = length(operon_ids), rel_widths = gene_counts)
}

plot_from_genome_id <- function(genome_id, min_score_fraction){
  file <- paste0('cds_around_cf_', genome_id, '.tsv.gz')

  cf_coords <- read_tsv(paste0(config$cf_filtered_cds_coords_dir, '/', file))

  cf_align_scores <- read_tsv(paste0(config$global_alignment_scores_dir, '/global_alignment_scores_cds_around_cf_', genome_id, '.tsv.gz'))

  prob_cf_coords <- find_prob_cf_coords(cf_coords)

  df <- make_cf_score_df(prob_cf_coords)

  plot_cf_matches(df, min_score_fraction)
}



##############################
# Config

args <- commandArgs(trailingOnly = TRUE)
config <- read_yaml(args[1])

config$cf_filtered_cds_coords_dir <- paste0(config$gene_search_dir, '/cf_cds_coords/')
config$global_alignment_scores_dir  <- paste0(config$gene_search_dir, '/cds_global_alignment_scores/')

cf_names <- read_tsv(config$cf_names)
cf_rules <- read_cf_rules(config$cf_rules)
cf_req_gene_counts <- get_cf_required_gene_counts(cf_rules)


system(paste('mkdir -p', config$out_dir))

# Loop
files <- list.files(config$cf_filtered_cds_coords_dir)

comb <- function(x, y){
  list(
    cf_classification = rbind(x[[1]], y[[1]]),
    classifications_by_gene = rbind(x[[2]], y[[2]]),
    best_matching_gene_any_cf = rbind(x[[3]], y[[3]])
  )
}

cores  <- config$cores
cl <- makeCluster(cores)
registerDoParallel(cl)

results <- foreach(
  i=1:length(files),
  .combine = 'comb',
  .init = list(tibble(), tibble(), tibble()),
  .inorder = FALSE,
  .packages = c('tidyverse')) %dopar% {

  file <- files[i]

  cf_coords <- read_tsv(paste0(config$cf_filtered_cds_coords_dir, '/', file))
  genome_id <- unique(cf_coords$genome_id) # Used to be genome_id here

  cf_align_scores <- read_tsv(paste0(config$global_alignment_scores_dir, '/global_alignment_scores_cds_around_cf_', genome_id, '.tsv.gz'))

  prob_cf_coords <- find_prob_cf_coords(cf_coords)

  df <- make_cf_score_df(prob_cf_coords)

  results <- classify_cfs(df, config$cf_classify_min_score_fraction)

  results
}




results$cf_classification %>% write_tsv(paste0(config$out_dir, '/cf_classification_min_classify_score_fraction_', config$cf_classify_min_score_fraction, '.tsv'))
results$classifications_by_gene %>% write_tsv(paste0(config$out_dir, '/classifications_by_gene_min_classify_score_fraction_', config$cf_classify_min_score_fraction, '.tsv'))
results$best_matching_gene_any_cf %>% write_tsv(paste0(config$out_dir, '/best_matching_gene_any_cf.tsv'))




if(config$output_gene_seqs || config$output_operon_seqs){
  library(plyranges)
  library(BSgenome)

  load_fa_seqs <- function(fa_file){
    Biostrings::getSeq(Rsamtools::FaFile(fa_file))
  }

  make_bakta_contignames <- function(dnastringset){
    names(dnastringset) <- paste0('contig_', 1:length(dnastringset))
    return(dnastringset)
  }


  if(config$output_gene_seqs){
    gene_seq_out_dir <- paste0(config$out_dir, '/seqs/genes/')
    system(paste('mkdir -p', gene_seq_out_dir))
    genome_ids <- results$classifications_by_gene$genome_id %>% unique()
    for(genome_id in genome_ids){
      gr_tmp <- results$classifications_by_gene %>% filter(genome_id == !!genome_id) %>% as_granges() # Filter before converting to get correct seqlevels
      contig_seqs <- load_fa_seqs(paste0(config$genome_seqs_dir, '/', genome_id, config$fa_suffix)) %>% make_bakta_contignames()
      nucl_seqs <- BSgenome::getSeq(contig_seqs, gr_tmp)
      names(nucl_seqs) <- paste0(gr_tmp$operon_id, ',', gr_tmp$my_cds_id, ',', gr_tmp$best_operon_cf_match, ',', gr_tmp$gene)
      nucl_seqs %>% writeXStringSet(paste0(gene_seq_out_dir, genome_id, '.fa.gz'), compress= TRUE)
    }
  }
  if(config$output_operon_seqs){
    operon_seq_out_dir <- paste0(config$out_dir, '/seqs/operons/')
    system(paste('mkdir -p', operon_seq_out_dir))
    operons.df <- results$classifications_by_gene %>%
      group_by(seqnames, strand = operon_strand, genome_id, operon_id, best_operon_cf_match) %>% 
      summarise(start = min(start), end = max(end), .groups = 'drop')
    genome_ids <- results$classifications_by_gene$genome_id %>% unique()
    for(genome_id in genome_ids){
      gr_tmp <- operons.df %>% filter(genome_id == !!genome_id) %>% as_granges() # Filter before converting to get correct seqlevels
      contig_seqs <- load_fa_seqs(paste0(config$genome_seqs_dir, '/', genome_id, config$fa_suffix)) %>% make_bakta_contignames()
      nucl_seqs <- BSgenome::getSeq(contig_seqs, gr_tmp)
      names(nucl_seqs) <- paste0(gr_tmp$operon_id, ',', gr_tmp$best_operon_cf_match)
      nucl_seqs %>% writeXStringSet(paste0(operon_seq_out_dir, genome_id, '.fa.gz'), compress= TRUE)
    }
  }
}



