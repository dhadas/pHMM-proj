#!/usr/bin/python

#Local imports
from process_hmm_results import process
from tblout_to_panadas import parse_tblout_to_df as parse_tblout
from extract_unknown_genes import extract
from pull_and_cluster import run_pullseq, double_cluster
from parse_mmseq_results import mmseq_clusters_info_to_df as parse_mmseq
from process_hmm_results import get_min_max_dict
from process_hmm_results import hits_df_to_structural_df

#Standard imports.
import click
from datetime import datetime
import os
import pandas as pd
from timeit import default_timer as timer

#Globals
date = datetime.today().strftime('%Y_%m_%d')
outpath = os.getcwd()
#
#
# def easy_pipeline(tblout_files_path, seqs_db_path,
#                   max_gene_distance, suffix, debug,  #proccess hmm opts
#                   keep_seq, min_cluster_size,
#                   ) -> int:
    #TIMER
    # start = timer()
    #Proccess tbltout files

input_path = '/Users/dror/PycharmProjects/AMR_pipeline/data/small_data_test_folder'
suffix = 'tbl'

df_dict = {}

df = parse_tblout(input_path, suffix)
df = df.sort_values('contig_id')

structured_df = process(df, 5, True)

unknown_genes_df = extract(structured_df)

unknown_genes_path = outpath + 'unknown_genes_' + date

# unknown_genes_df.to_csv('unknown_genes_df.csv')

# Pull unknown genes
# unknown_genes_fasta_path = run_pullseq(seqs_db_path, unknown_genes_path)




    # unknown_genes_path = extract(structured_df)
    #
    # #Pull unknown genes
    # unknown_genes_fasta_path = run_pullseq(seqs_db_path, unknown_genes_path)
    #
    # #Cluster
    # double_cluster(unknown_genes_fasta_path, 0.2, 0.8, 0.99, 0.3, 'clusRes_1st_stage', 'clusRes_final')
    #
    # #Analyze clustering results
    # write_flag = True
    # mmseq_results_filename = parse_mmseq(outpath, 'clusRes_final', min_cluster_size, write_flag, keep_seq)


# df_filtered = pd.read_csv('/Users/dror/PycharmProjects/AMR_pipeline/data/df_filtered_2020_02_21.csv', index_col = 0)
# min_max_dict = get_min_max_dict(df_filtered)
# final_df = hits_df_to_structural_df(df_filtered, min_max_dict, 5)
# final_df.to_csv(f"/Users/dror/PycharmProjects/AMR_pipeline/data/results/final_df_{date}_n5.csv")

