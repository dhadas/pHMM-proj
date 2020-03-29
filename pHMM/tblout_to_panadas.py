#!/usr/bin/python

# from __future__ import annotations
import glob
import os
import numpy as np
import pandas as pd
from datetime import datetime
outpath = os.getcwd()


def get_headers(type : int = 0) -> dict:
    #TODO get rid of this function in the future, it's not really necessary
    headers_dict = {'id_cols': ["target_name", "target_accession", "query_name", "query_accession"],
                    'full_seq_cols': ["eval_full_seq", "score_full_seq", "bias_full_seq"],
                    'best_dom_cols': ["e_val_best_dom", "score_best_dom", "bias_best_dom"],
                    'misc_cols': ["exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "target_description"]
                    }

    # if type is 1, return a list
    if (type == 1):
        return headers_dict['id_cols'] + headers_dict['full_seq_cols'] + headers_dict['best_dom_cols'] + \
               headers_dict['misc_cols']

    # default - return a dictionary
    return headers_dict

def tblout_to_pd(path: str, debug: bool = False) -> pd.DataFrame:
    '''
    function reads a csv tblout file from {path} then turns it into a dataframe.
    function assumed csv is in a tblout format -> output of hmmsearch with -tblout
    function ignores lines that start with '#'
    function separtes cols by 2 or mote spaces
    :param path:
     path to read tblout csv from.
    :return:
    a dataframe, with same columns as csv in {path} and a row for each row in the file
    '''

    # read file to dataframe with a single column containing all data
    if(debug):
        print(f"Now reading tblout file at path:\t{path}\n")
    headers = get_headers(1)
    df = pd.read_csv(path, header=None, index_col=0, sep='\t')

    # drop lines statrting with # from data
    mask_idx = df.index.map(lambda x: x.startswith("#") == False)
    masked_df = df.loc[mask_idx, :].reset_index().copy().rename(columns={0: 'col'})
    if(masked_df.empty):
        return masked_df

    # Seperate cols by spaces, drop orginial 'col'
    reg = r"\s{1,}"  # seprate by 2 or more spaces

    # Create a dictionary for column names
    col_dict = dict(zip(np.arange(len(headers)), headers))
    # split
    masked_df = masked_df.join(masked_df['col'].str.split(reg, expand=True)).rename(columns=col_dict)

    # drop
    masked_df = masked_df.drop('col', axis=1)

    # replace Nones created in the end of each row by extra spaces by empty strings and the rejoin them
    masked_df = masked_df.fillna('')
    masked_df['target_description'] = masked_df['target_description'].str.cat(masked_df.iloc[:, 19:], sep=' ')
    df = masked_df.iloc[:, :19].copy()
    # print(masked_df)

    df = force_data_types(df)

    return df

def multiple_tbl_to_pd(file_path_lst: list, debug : bool = False) -> pd.DataFrame:
    '''
    Function recives a list {file_path_lst} of file paths and returns a concatenated dataframe from all of them
    function reaptedly calls tblout_to_pd to create the returned dataframe then concatenates
    :param file_path_lst:
    a list of paths to files. all files will be in the originial output format of hmmsearch with a -tblout flag.
    :return:
    A concatenated dataframe, with same columns as csv in {path} and a row for each row in the file
    '''
    df_lst = []
    failed_lst = []
    num_of_files = len(file_path_lst)
    skipped_counter = 0

    #for each file path, use tblout_to_pd function to transform it into a dataframe
    #then append to a list of dataframes
    for hmm_out_file in file_path_lst:
        try:
            tmp = tblout_to_pd(hmm_out_file)
        except:
            failed_lst.append(hmm_out_file)
            skipped_counter += 1
            print(f"Failed to parse. skipping file:\n{hmm_out_file}")
            print("Number of skipped files so far : ", skipped_counter)


        if(tmp.empty):
            continue
        df_lst.append(tmp)

    #concatenate all dataframes in the list into one dataframe and return it
    if(df_lst == []):
        raise ValueError("tblout files from the given path list are empty, no dataframe can be created.")
    df = pd.concat(df_lst)
    df = df.reset_index()

    print(f"Finished creating dataframe from tblout files. Number of concatenated files is {num_of_files-skipped_counter} "
          f"out of {num_of_files}. Skipped {skipped_counter} files. returning dataframe.")

    if debug and len(failed_lst) > 0:
        print(f"DEBUG = True. Printing problematic files list. list length: {skipped_counter}")
        for file_name in failed_lst:
            print(file_name)

    return df

def force_data_types(df: pd.DataFrame) -> pd.DataFrame:
    '''
    Transform dataframe datatypes into numeric or string format, based on the data types in the originial
    hmmsearch output format.
    :param df: he dataframe containing the data.
    :return: A copy of the original dataframe with the new data types.
    '''

    headers = get_headers()
    # Force data typesge for the cells of the DataFrame for later manipulations
    df.loc[:, headers['full_seq_cols'] + headers['best_dom_cols']] = df.loc[:, headers['full_seq_cols'] + headers[
        'best_dom_cols']].astype(float)
    df.loc[:, "exp"] = df.loc[:, "exp"].astype(float)
    df.loc[:, headers['misc_cols'][1:-1]] = df.loc[:, headers['misc_cols'][1:-1]].astype(int)
    df.loc[:, headers['misc_cols'][:-1]] = df.loc[:, headers['misc_cols'][:-1]].astype(str)

    return df.copy()

def split_sample_names(tblout_df: pd.DataFrame, debug: bool = False) -> pd.DataFrame:
    '''
    Add new columns to dataframe (for all rows), based on the target names.
    fucntion assumes the target name format is: {sample_name}_ctg_{contig_id}_{gene_id}
    :param tblout_df: Originial dataframe containing all tblout format data.
    :return:
    A copy of the original dataframe with the newly added colums: 'gene_id', 'contig_id', 'sample_contig'
    Original 'target_name' coulmn remains unchanged.
    '''
    cols_to_add = ['contig_id', 'gene_id']
    temp1 = tblout_df.target_name.str.split('_ctg_', expand = True)
    temp2 = temp1[1].str.split('_', expand = True).rename(columns = {0:'contig_id', 1:'gene_id'})


    # Remove rows that contain NaN's. This will be caused by incompatible sample name format's.
    # For example. a sample that has no 'ctg' in the name will not split up well and it will later crash when we
    # try to separate the name into contig and gene.
    result_df = tblout_df.join(temp2[cols_to_add])


    if debug:
        NaN_rows = tblout_df.isnull().any(axis=1).sum().max()
        print(f"In split_sample_names, Debug = True.\n"
              f"Dropping NaN rows. num of dropped rows : {NaN_rows}"
              f"Writing dropped files to {outpath}/faulty_file_names.csv")
        faulty_df = result_df[result_df.isnull().any(axis = 1)]
        faulty_df.to_csv(f'{outpath}/faulty_file_names.csv')

    result_df.dropna(inplace=True)
    result_df['sample_contig'] = result_df.target_name.apply(lambda name: '_'.join(name.split('_')[:-1]))
    result_df['gene_id'] = result_df['gene_id'].astype(int)
    result_df['sample_contig'] = result_df['sample_contig'].astype(str)

    return result_df

def parse_tblout_to_df(files_path: str, file_suffix: str = 'tbl', debug: bool = True) -> pd.DataFrame:
    '''
    recieve the path to a folder containing hmmsearch results, collect them and parse into a dataframe.
    :param files_path: Path to files folder when hmmsearch outputs are, for example : '/davidbur/drohadas'.
    :param file_suffix: suffix of the output files of the hmmsearch command that exist in the folder in {files_path}, default is 'tbl'
    :param debug: defualt True, when false debugging printouts will not be displayed (recommended on).
    :return: a datafrane of the results, a concatenation of all hmm output files.
    the dataframe will contain the original columns of the hmm output table and the additional columns:
        1. 'gene_id' - representing the gene number of the current entry (the last number in the the target name)
        2. 'contig_id' - representing the number of contig the gene belongs too. the the second to last number in the target name.
        3. 'sample_contig' - this is all of the target name apart from the gene, it the future it will replace 'contig_id'

    '''

    if (files_path == '' or files_path == None):
        raise ValueError("Input path is empty or None.")

    '''
    #Set project path to collect sample files from, set suffix, set out folder path,
    #collect all in files {path} ending with {suffix} into {files_list}
    #Set a data object, for later use in file writing
    '''
    path = os.path.join("~", files_path)
    files_lst = glob.glob(os.path.join(path, f"*.{file_suffix}"))

    if(debug):
        print(f"Parsing tblout files:\n" \
            f"Parsing files from path: {files_path}.\n"
            f"With suffix: '.{file_suffix}'")

    if files_lst == []:
        raise ValueError(f"Input folder had no files with suffix: {file_suffix}, \
                              glob.glob from path with *.{file_suffix} returned an empty list")

    # Create a pandas dataframe from all tblout files, force numeric data types
    df = multiple_tbl_to_pd(files_lst, debug)


    # Force numeric types on data
    df = force_data_types(df.copy())

    '''
    Split the sample names to add sample_name, contig_id, gene_id  columns
    Get the best hits for each gene in the current dataframe
    Remove contigs that have only one hit, they have no potential
    Possibly write the filtered dataframe out as a csv for further inspection or use
    '''
    # Enrich data with new columns -> sammple_name, ctg_id, gene_id, , it is possible to supply a diffrent separator. defualt is '_'
    df = split_sample_names(df.copy())



    print(f"Finished parsing tblout files, {len(files_lst)} files with suffix {file_suffix} were found and processed.\n"
          f"parsed dataframe is of shape: {df.shape} - 3 extra columns added.\n"
          f"Returninig dataframe")

    return df

