#!/usr/bin/python

#Local imports
from process_hmm_results import process
from tblout_to_panadas import parse_tblout_to_df as parse_tblout
from extract_unknown_genes import extract
from pull_and_cluster import run_pullseq, double_cluster
from parse_mmseq_results import mmseq_clusters_info_to_df as parse_mmseq

#Standard imports.
import click
from datetime import datetime
import os
from timeit import default_timer as timer

#Globals
date = datetime.today().strftime('%Y_%m_%d')
outpath = os.getcwd()



@click.group()
def main():
    pass


@main.command('pullseq', short_help = 'Use pullseq to create a fasta file of unknown genes')
@click.argument('unknown_genes_path', type = click.Path())
                # help = "Path to csv/txt file that can be used for pullseq")
@click.argument('seqs_db_path', type = click.Path())
                # help = "Path to the sequences database to be used")
def pull(seqs_db_path, unknown_genes_path):
    '''
    cmd line tool for running pullseq using python. Assuming the pullseq module is loaded.
    Selcting names {unknown_genes_path} from the seqs db in {db_path}.
    cmd is : $ pullseq -i {seqs_db_path} -n {unknown_genes_path}
    :param unknown_genes_path: Path to unknonwn genes csv
    :param seqs_db_path: Path to All sequences database (fasta file).
    :return: 0 upon success.
    '''
    run_pullseq(seqs_db_path, unknown_genes_path)
    return 0



@main.command('cluster', short_help = 'Double cluster sequences hased on mmseqs easy-cluster')
@click.argument("input_path",
                type = click.Path())
@click.option("--cov_rate1", '-c1',
                help = "Coverage rate for mmseqs2 easy-cluster first run, "
                       "this run is used to filterout identical proteins, default value is 0.2",
                default = 0.2,
                type = float
              )
@click.option("--cov_rate2", '-c2',
                help = "Coverage rate for mmseqs2 easy-cluster first run, "
                       "this run is used for the final clustering on the first run's representatives 0.8",
                default = 0.8,
                type = float
              )
@click.option("--min_seq_id1", '-m1',
                help = "minimal sequence id rate for mmseqs2 easy-cluster first run, "
                       "this run is used to filterout identical proteins, default value is 0.99",
                default = 0.99,
                type = float
              )
@click.option("--min_seq_id2", '-m2',
                help = "minimal sequence id rate for mmseqs2 easy-cluster first run, "
                       "this run is used for the final clustering on the first run's representatives 0.3",
                default = 0.3,
                type = float
              )
@click.option("--prefix1", '-p2',
                help = "prefix of clustering output files for 1st run",
                default = 'clusRes_1st_stage',
                type = str
              )
@click.option("--prefix2", '-p2',
                help = "prefix of clustering output files for 2st run",
                default = 'clusRes_final',
                type = str
)
def cluster(input_path: str,
            cov_rate1 : float, cov_rate2: float,
            min_seq_id1 : float, min_seq_id2 : float,
            prefix1 :str, prefix2 :str) -> int:
    '''
    Command line tool for double-clustering operation via mmseqs2 easy-cluster command.
    First cluster with high identity and low coverage to get representatives proteins for each cluster.
    Later cluster with low identity and higher coverage to cluster representatives by similiarity.
    outputs are as declared by mmseqs documentation:
    https://github.com/soedinglab/MMseqs2/tree/master/src/clustering.
    While running these commands a temporary folder is created and the later removed, this is valid.
    output files of the final call will be in the format of {prefix2}_cluster, {prefix2}_all_seqs, {prefix2}_rep_seq
    and can be further analyzed with the parse_mmseqs module in this package.
    cmd line used by module:
    mmseqs easy-cluster -c {cov_rate} --min-seq-id {min_seq_id} {input_path} {prefix} tmp
    the mmseqs2 module must be loaded while running this module.
    :param input_path: of sequences to be clustered.
    :return: 0 upon success. Output files can be located by collecting files starting with {results_prefix}
    and further processed using the parse_mmseq_results module.
    '''
    double_cluster(input_path, cov_rate1, cov_rate2, min_seq_id1, min_seq_id2, prefix1, prefix2)
    return 0



@main.command('process_tblout', short_help = 'Process hmmsearch outputs')
@click.argument("files_path", type = click.Path())
@click.option('--max_gene_distance', '-n',
              default = 5,
              help = 'Maximal distance between two known genes. Genes within this distance that are not recognized will'
                     ' be regarded as potential - included in a possible AMR cassete.',
              type = int)
@click.option('--suffix', '-s',
              default = 'tblout',
              help = 'Suffix of input files. this is used to collect all tblout files from {files_path}',
              type = str)
@click.option('--debug', '-d',
              default = True,
              help = 'When debug mode is on (True), Intermidate printouts and some files are also printed.',
              type = bool)
def process_hmmsearch_outputs(files_path, max_gene_distance, suffix, debug):
    '''
    Command line tool for processing tblout outputs.
    Output file formats:
    1. df_filtered_{date}.csv - A dataframe with only the best hit for each gene in each contig.
    2. final_df_{date}_n{max_gene_distance}.csv - A dataframe describing the AMR cassettes,
    3. unknown_genes_{date}.csv - A list of gene ids, a preperation for pullseq.
    files 1 and 2 will only be printed when debug is set to True (default value).
    :param files_path: path of tblout files.
    :return: 0 opon Success.
    '''
    input_path = files_path
    n = max_gene_distance
    debug_mode = debug

    #Parse into tblout and concatenate results into 1 df.
    df = parse_tblout(input_path, suffix, debug_mode)
    #Build AMR cassetes.
    structured_df = process(df, n, debug_mode)
    #Write to file unknown genes.
    unknown_genes_path = extract(structured_df)

    return 0



@main.command('parse_mmseq', short_help = 'Parse mmseq results')
@click.argument("mmseq_results_path", type = click.Path())
@click.option("--results_prefix", '-pref', '-p',
                help = "Prefix of mmseq results files",
                default = 'clusRes_final',
                type = str
              )
@click.option("--min_cluster_size", "-min_clus_size", "-m",
              help = 'Clusters below this size will not be included in the fimal output',
              default = 10,
              type = int
              )
@click.option('--keep_seq/--dkeep_seq', ' /-d',
              help = 'When flag is True, Final cluster representatives seqs will also be included in the output',
              # type = bool,
              default=False,
              )
def parse_mmseq_results(mmseq_results_path, results_prefix, min_cluster_size, keep_seq):
    '''
    Parse all the mmseq post clustering results.
    Ouput a csv with info about clusters of size bigger then {min_cluster_size}
    :param mmseq_results_path: Path to results (Path to a folder containing the mmseq results files).
    :return: 0 upon succes
    '''
    parse_mmseq(mmseq_results_path, results_prefix, min_cluster_size, keep_seq)
    return 0


@main.command('easy-pipeline', short_help = 'Run whole pipeline at once, most params are set to default.')
@click.argument("tblout_files_path", type = click.Path())
@click.argument('seqs_db_path', type = click.Path())
@click.option('--max_gene_distance', '-n', default = 5, type = int)
@click.option('--suffix', '-s', default = 'tbl', type = str)
@click.option('--debug', '-d', default = True, type = bool)
@click.option('--keep_seq/--dkeep_seq',
              help = 'When flag is True, Final cluster representatives seqs will also be included in the output',
              default=False,
              )
@click.option("--min_cluster_size", "-min_clus_size", "-m",
              help = 'Clusters below this size will not be included in the fimal output',
              default = 100,
              type = int
              )
@click.option("--verbose", "-v",
              help = 'When True, all functions that do not run multiple times (e.g functions that are called within loops)'
                     'will print message when they start and when they end. this can be used for debugging. default is False',
              default=False,
              type = bool
              )
def easy_pipeline(tblout_files_path, seqs_db_path,
                  max_gene_distance, suffix, debug,  #proccess hmm opts
                  keep_seq, min_cluster_size, verbose,
                  ) -> int:
    #TIMER
    start = timer()

    #Proccess tbltout files
    input_path = tblout_files_path
    n = max_gene_distance
    debug_mode = debug
    df = parse_tblout(input_path, suffix, debug)
    structured_df = process(df, n, debug_mode)

    unknown_genes_df, unknown_genes_path = extract(structured_df, verbose)
    #Pull unknown genes
    unknown_genes_fasta_path = run_pullseq(seqs_db_path, unknown_genes_path, verbose)

    #Cluster
    double_cluster(unknown_genes_fasta_path, 0.2, 0.8, 0.99, 0.3, 'clusRes_1st_stage', 'clusRes_final', verbose)

    #Analyze clustering results
    clusters_by_representatives, all_cluster_members = parse_mmseq(unknown_genes_df, 'clusRes_final',
                                                                   min_cluster_size,
                                                                   False, #keep_sequences_flag
                                                                   verbose,
                                                                   )

    # TIMER
    end = timer()
    print(f"\nFinished a full pipeline, runtime :{round(round(end - start,0)//60,1)}, minutes")

    return 0

if __name__ == "__main__":
    main()
