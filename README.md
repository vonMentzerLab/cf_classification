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


