import os
from datetime import datetime
import pandas as pd
date = datetime.today().strftime('%Y_%m_%d')
from general_utils import FUNC_START, FUNC_END
VERBOSE = False
out_path = os.getcwd()

def fetch_unknown_genes(df: pd.DataFrame) -> list:
    '''
    TODO - Improve this so it runs faster, maybe create another column in the original final_df with the genes names.
    Receive a dataframe of format 'final' (as outputted by 'hits_df_to_structural_df2' function.
    for each sample and contig, extract gene numbers from the unknown genes column into a list.
    for each unknown gene in the list, format it's name back to the original fasta contigs file format.
    function does not change the original df, instead it creates a local copy.
    :param df:
    Original dataframe contating data
    :return:
    A list of gene names that can be used to pull genes from a fasta file.
    '''
    if (VERBOSE):
        FUNC_START()

    unknown_genes_list = []
    original_cassetes_list = []
    local_df = df.copy()
    for index, row in local_df.iterrows():
        full_id = row.full_id

        min_gene = min(
            row.known_gene_ids)  # This is also the number of digits that has to be replaced by a an unknwon gene id
        min_gene_len = len(str(min_gene))

        for unknown in row.unknown_gene_ids:
            unknown_gene_len = len(str(unknown))
            if unknown_gene_len == min_gene_len:  # this is usually the case, the known gene and the unknown gene are similiar, like 005 and 006
                curr_id = full_id[0:-min_gene_len]
                curr_id += str(unknown)
            elif unknown_gene_len >= min_gene_len:  # e.g known gene id is 9 and unknown gene is 11, needs to replace another digit
                #             print(f"unknown gene is longer the known gene case\n unknown is : {unknown}, known is {min_gene}")
                curr_id = full_id[0: -unknown_gene_len]
                curr_id += str(unknown)
            #             print(f"known gene id is:\n {full_id}")
            #             print(f"final curr id is:\n {curr_id}")
            else:  # e.g Known gene id is 10 and unknown gene id is 9, replace with 0s
                #             print(f"known gene is longer the unknown gene case\n unknown is : {unknown}, known is {min_gene}")
                curr_id = full_id[0:- min_gene_len]
                num_of_zeros = min_gene_len - unknown_gene_len
                for i in range(min_gene_len - unknown_gene_len):
                    curr_id += '0'
                curr_id += str(unknown)
            #             print(f"known gene id is:\n {full_id}")
            #             print(f"final curr id is:\n {curr_id}")



            unknown_genes_list.append(curr_id)
            original_cassetes_list.append(full_id)

    # cols = ['unknown_gene', 'cassete_representing_gene']
    res_df = pd.DataFrame({'unknown_gene':unknown_genes_list, 'cassette_representing_gene':original_cassetes_list})
    res_df = res_df.drop_duplicates(subset='unknown_gene', keep="first")

    # res_df.columns = cols

    if (VERBOSE):
        FUNC_END()
    return res_df

def write_unknown_genes(unknown_genes_df: pd.DataFrame) -> str:
    '''
    Recive a list of genes to write to ouput file.
    The fucntion assumes nothing about gene list.
    :param unkonwn_genes_list: list of gene names to be written.
    :param output_path: desired path and name of output file, default value is 'txt'
    :return:
    output path for pullseq command.
    Raises 'ValueError' if list is empty.
    other values to be added in the future.
    '''

    if (VERBOSE):
        FUNC_START()

    if unknown_genes_df.empty == True :
        raise ValueError

    out_path1 = out_path + f"/unknown_genes_for_pullseq_{date}.csv"
    unknown_genes_df.to_csv(path_or_buf=out_path1, columns=['unknown_gene'], header=False, index=False)

    out_path2 = out_path + f"/unknown_genes_for_analysis_{date}.csv"
    unknown_genes_df.to_csv(out_path2)

    if (VERBOSE):
        FUNC_END()
    return out_path1

def extract(structured_df : pd.DataFrame, VERBOSE_FLAG: bool = False) -> pd.DataFrame:
    '''
    Prepare output for pullseq make a list of unknown genes - ORF ID per line
    Then write it to a csv file.
    return the path of the csv file
    :param structured_df: the dataframe of the cassettes, this dataframe contains the identified ARGs and potential ARGs.
    A potential ARG is any gene in the cassette that was not identified earlier as an ARG.
    :return: a new dataframe containing all the potential ARGs and the cassette they originate from. The dataframe is
    internally used in the pipeline. During the module's run it also prints the a 1 column csv of the unknown genes that
    is later used for in the run_and_pullseq modulecassette_representing_gene.
    '''

    global VERBOSE
    VERBOSE = VERBOSE_FLAG

    if (VERBOSE):
        FUNC_START()
    # out_path = f"{os.getcwd()}/unknown_genes_{date}.csv"

    print("Extracing uknonwn genes from dataframe.\n"
          f"Output path is: {out_path}\n")

    unknown_genes_df = fetch_unknown_genes(structured_df)
    out_path_for_pullseq = write_unknown_genes(unknown_genes_df)

    print("Finished extracting unknonwm genes.\n"
          f"Writing result to output path: {out_path}, returning a tuple: (unknown_genes_df, output_path of unknown_genes_csv).")

    if (VERBOSE):
        FUNC_END()
    return (unknown_genes_df, out_path_for_pullseq)


