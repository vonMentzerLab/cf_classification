## Installation

The required packages are detailed in two conda environment yaml files in the `env` directory. To create the environments, run `conda env create --file environment_file_path.yaml -n desired_environment_name`. `miniforge` installations should just require replacing `conda` with `mamba`, and is recommended for faster environment solving.

The `bakta` environment is meant for the `bakta_annotation.sh` script. All other scripts should be run in the `cluster_classification` environment.

## CDS extraction

The pipeline requires gff3 files to extract CDS locations from. This has currently only been tested with Bakta. Other tools may work, but please note that the pipeline expects CDSs not to overlap in a manner where order along the DNA sequence can't be determined easily - essentially, the edges of CDSs may overlap, but no more, for the pipeline to work as expected.


1. Annotate contigs with Bakta. Either use `bakta_annotation.sh` directly, or look at its contents and adapt.
    Usage: `bakta_annotation.sh directory_with_contig_fasta_files bakta_database`
    Note that a bakta database has to be generated if this has not already been done. Please see Bakta's documentation.

2. Run `get_bakta_cds_defs_and_seqs.R` to extract CDS locations as well as nucleotide and amino acid sequences.
    Usage: `get_bakta_cds_defs_and_seqs.R dir_with_gff3_files dir_to_save_coords_in`


The CDS extraction steps have to be performed once per set of genomes, but do not have to be repeated when searching the same genome dataset for different groups of genes with the pipeline.

## Search for gene clusters using predefined cluster definitions

1. Run `find_gene_clusters.sh config.yaml`

Pipeline parameters are defined in `config.yaml`. Please see the example config file in addition to the following explanation.

`cds_features_dir`: The output directory of `get_bakta_cds_defs_and_seqs.R`

`gene_search_dir`: Directory where Diamond alignments and other intermediary files will be stored.

Fasta files of the gene sequences to search for need to be available as both nucleotide (`clustered_genes_reference_db_nucl_seq`) and amino acid (`clustered_genes_reference_db_aa_seq`) sequences. 

`cf_names`: To reduce the length of names in the fasta files for more compact output, a file with columns `short_name` and `name` (i.e. the fasta sequence name) is used. If no name changes are desired, provide a file with the same names in both columns.

`cf_rules`: A table specifying which genes belong to which CF (or other group/cluster), with the following columns:
- `cf`: The CF (or other group/cluster) a gene belongs to
- `gene`: the name from the `short_name` column in the `cf_names` file.
- `required`: Takes values `TRUE` and `FALSE`. Specifies whether the gene is required in order for a cluster to be considerd complete. Regulators are a good example of genes one might want to set `required: FALSE` for.

`min_cov`: Minimum coverage accepted for a Diamond hit to be considered probable clustered gene database hit

`min_pident`: Minimum percent identity accepted for a Diamond hit to be considered probable clustered gene database hit

`max_dist`: Maximum distance between probable genes from clustered gene database to be considered part of the same cluster.

`cf_classify_min_score_fraction`: Minimum alignment score fraction for a gene to be classified as a gene from the clustered gene database, rather than "Unknown".

`out_dir`: Directory to save output files to.

## Output

Three output files will be generated:

- `cf_classification_min_classify_score_fraction_*.tsv`
- `classifications_by_gene_min_classify_score_fraction_*.tsv`
- `best_matching_gene_any_cf.tsv`

The latter two are lists of found genes, with the best match for each gene. `best_matching_gene_any_cf.tsv` lists the best matching genes regardless of whether they belong to the same CF/cluster or not. `classifications_by_gene_min_classify_score_fraction_*.tsv` is nearly the same, but gives preference to genes in the CF/cluster that matches the gene cluster best, provided that they have a score fraction of at least `cf_classify_min_score_fraction` (otherwise the best match regardless of CF is listed). This is done to avoid situations where very similar genes are misclassified as belonging to a different CF, when they're an almost equally good match for the CF that the rest of the genes match.


### Column explanations for the two gene level files

- `seqnames`: Sequentially numbered contig name.
- `start` and `end`: 1-based start and end positions.
- `strand`: The strand the gene is on.
- `operon_strand`: The strand that most genes in a cluster are on.
- `genome_id`: Genome name taken from the contigs fasta file.
- `operon_id`: Name given to the an identified cluster of genes, based on the contig and then sequentially numbered to distinguish between multiple clusters.
- `my_cds_id`: CDS identifier taken from the gff3 file.
- `gene`: Best matching gene (either with preference for the same CF as other genes in the cluster or not, see above).
- `cf`: CF that `gene` belongs to.
- `best_operon_cf_match`: The best matching CF for the cluster of genes taken together.
- `score`: Alignment score from aligning CDS encoded protein to protein in reference database.
- `max_score`: Alignment score from aligning CDS encoded protein to itself.
- `score_fraction`: `score` / `max_score` fraction.
- `pid`: PID calculated as (identical positions) / (aligned positions + internal gap positions).
- `gene_order`: Order index of the gene in the identified cluster.

### Additional column explanations for the gene cluster level file

- `cf_classification`: Overall classification of the identified cluster of genes. Can be a mix of different CFs or contain "Unknown" if there isn't a sufficiently good match for the best matching operon.
- `gene_profile`: String with the genes from `classifications_by_gene_min_classify_score_fraction_*.tsv` concatenated in order, but with "Unknown" listed for genes under the score threshold.
- `mean_score_fraction`: The mean of the score fractions (see above) for each of the genes in the cluster, including those listed as "Unknown".
- `best_operon_cf_match_mean_score_fraction`: Same as `mean_score_fraction`, but taking the best matching genes only from the CF that matches the cluster best as a whole.
- `operon_gene_count`: The number of genes in the identified cluster.
- `best_operon_cf_match_gene_count`: The number of genes in the cluster that match a gene in `best_operon_cf_match`.
- `req_gene_count_in_cf_rules`: The number of genes that are listed as required in `cf_rules`. Please note that direct comparisons with `best_operon_cf_match_gene_count` may be misleading, as the latter includes non-required genes.
- `best_operon_cf_match_all_genes_present`: Whether all of the required genes in `cf_rules` for `best_operon_cf_match` have a match in the cluster of genes.
- `best_operon_cf_match_genes_in_order`: Whether the genes are in the right order compared to `cf_rules`, with no other genes in between.
- `best_operon_cf_match_all_genes_and_ordered`: Whether both of the previous two conditions are met.