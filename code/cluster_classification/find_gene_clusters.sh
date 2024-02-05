#!/usr/bin/bash

### Usage ###
# find_gene_clusters.sh config.yaml


# Activate the right conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate adhesins_661k

# shyaml get-value

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# clustered_genes_reference_db_aa_seq=$(cat "$1" | shyaml get-value clustered_genes_reference_db_aa_seq)
# cds_amino_seqs_dir=$(cat "$1" | shyaml get-value cds_amino_seqs_dir)
# diamond_out_dir=$(cat "$1" | shyaml get-value diamond_out_dir)

clustered_genes_reference_db_aa_seq=$(cat "$1" | shyaml get-value clustered_genes_reference_db_aa_seq)
cds_features_dir=$(cat "$1" | shyaml get-value cds_features_dir)
cds_amino_seqs_dir="$cds_features_dir"/cds_amino_seqs/
gene_search_dir=$(cat "$1" | shyaml get-value gene_search_dir)
diamond_out_dir="$gene_search_dir"/diamond_aligned_amino_cf_db/


# echo "$gene_search_dir"
# exit



# Parts that vary based on ref db
# conda activate adhesins_661k
"$SCRIPT_DIR"/subscripts/parallel_diamond_cds_amino_cf_ref.sh "$clustered_genes_reference_db_aa_seq" "$cds_amino_seqs_dir" "$diamond_out_dir"


# Filter diamond
conda activate R
Rscript "$SCRIPT_DIR"/subscripts/filter_diamond_cds_cf_related.R "$1"

# Classify
conda activate R
Rscript "$SCRIPT_DIR"/subscripts/classify_cf_profile.R "$1"
