from __future__ import annotations
import glob
import os
import numpy as np
import pandas as pd
from datetime import datetime


class process_hmm():

    # Class column definitions


    @staticmethod
    def get_headers(type: int = 0) -> dict:
        #TODO get rid of the list output possibility and safley update all code accordingly
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

    def get_min_max_dict(df: pd.DataFrame) -> pd.DataFrame:
        # Create an entry for each contig, determining it's min possible gene and max possible gene
        sample_contig = pd.Series(df.target_name).str.split('_')
        f = lambda x: x[0] + '_' + x[1] + '_' + x[2]
        sample_contig = sample_contig.apply(lambda x: f(x))
        df['sample_contig'] = sample_contig

        temp = df.groupby(['sample_contig']).agg({'gene_id': ['min', 'max']})
        min_max_dict = {}
        for row in temp.iterrows():
            # print(row[0])
            # print("Row max\min are : ",  row[1][0],  row[1][1])
            min_max_dict[row[0]] = dict({'min': row[1][0], 'max': row[1][1]})
        return min_max_dict

    @staticmethod
    def tblout_to_pd(path: str) -> pd.DataFrame:
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
        headers = process_hmm.get_headers(1)
        # read file to dataframe with a single column containing all data
        df = pd.read_csv(path, header=None, index_col=0, sep='\t')

        # drop lines statrting with # from data
        mask_idx = df.index.map(lambda x: x.startswith("#") == False)
        masked_df = df.loc[mask_idx, :].reset_index().copy().rename(columns={0: 'col'})

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

        df = process_hmm.force_data_types(df)

        return df


    @staticmethod
    def multiple_tbl_to_pd(file_path_lst: list) -> pd.DataFrame:
        '''
        Function recives a list {file_path_lst} of file paths and returns a concatenated dataframe from all of them
        function reaptedly calls tblout_to_pd to create the returned dataframe then concatenates
        :param file_path_lst:
        a list of paths to files. all files will be in the originial output format of hmmsearch with a -tblout flag.
        :return:
        A concatenated dataframe, with same columns as csv in {path} and a row for each row in the file
        '''
        df_lst = []

        #for each file path, use tblout_to_pd function to transform it into a dataframe
        #then append to a list of dataframes
        for hmm_out_file in file_path_lst:
            tmp = process_hmm.tblout_to_pd(hmm_out_file)
            df_lst.append(tmp)

        #concatenate all dataframes in the list into one dataframe and return it
        df = pd.concat(df_lst)
        df = df.reset_index()

        return df
    @staticmethod
    def force_data_types(df: pd.DataFrame) -> pd.DataFrame:
        '''
        Transform dataframe datatypes into numeric or string format, based on the data types in the originial
        hmmsearch output format.
        :param df:
        The dataframe containing the data
        :return:
        A copy of the original dataframe with the new data types
        '''

        #TODO copy df into df_local right on the begninng
        #TODO add the newly forced data types to the documentation
        headers = process_hmm.get_headers()
        # Force data typesge for the cells of the DataFrame for later manipulations
        df.loc[:, headers['full_seq_cols'] + headers['best_dom_cols']] = df.loc[:, headers['full_seq_cols'] + headers[
            'best_dom_cols']].astype(float)
        df.loc[:, "exp"] = df.loc[:, "exp"].astype(float)
        df.loc[:, headers['misc_cols'][1:-1]] = df.loc[:, headers['misc_cols'][1:-1]].astype(int)
        df.loc[:, headers['misc_cols'][:-1]] = df.loc[:, headers['misc_cols'][:-1]].astype(str)

        return df.copy()

    @staticmethod
    def split_sample_name(df: pd.DataFrame, separator: str = '_') -> pd.DataFrame:
        '''
        Add new columns to dataframe (for all rows), based on the target names.
        fucntion assumes the target name format is: {sample_name}_{contig_id}_{gene_id}
        {separator} between the target name parts can be exchanged. defualt parameter is '_'
        :param df: Originial dataframe containing all data
        :param separator:
        An optional parameter. represents the seperatot between the parts of the targetname/
        defaulat values is '_'.
        :return:
        A copy of the original dataframe with the newly added colums: 'gene_id', 'contig_id'
        Original 'target_name' coulmn remains unchanged.
        '''
        #TODO consider adding all necessary columns already.
        cols = {0: 'sample_names', 1: 'ctg', 2: 'contig_id', 3: 'gene_id'}
        local_df = df.copy().join(df['target_name'].str.split(separator, expand=True)).rename(columns=cols)
        local_df['gene_id'] = local_df['gene_id'].astype(int)

        local_df['contig_id'] = 'ctg_' + local_df.contig_id
        local_df.contig_id = local_df.contig_id.astype(str)
        local_df = local_df.drop(['ctg', 'index'], axis=1)

        return local_df


    @staticmethod
    def get_best_hits(df: pd.DataFrame) -> pd.DataFrame:
        '''
        Recives a dataframe {df} of tblout outputs and returns a dataframe with the best (1) hit for each gene.
        Best hit is determined by the 'score_full_seq' coulmn value as appears in the hmmsearch tblout output.
        :param df:
        Pandas dataframe containing original data
        :return:
        A new dataframe containg *only* the best hit for each gene.
        '''
        local_df = df.copy()
        local_df = local_df.sort_values('score_full_seq', ascending=False).drop_duplicates('target_name')
        local_df.reset_index(inplace=True)  # create a best hits dataframe

        return local_df

    ''' Improve later'''
    @staticmethod
    def get_best_hits_old(df: pd.DataFrame) -> pd.DataFrame:
        '''
        DEPRECATED - OLDCODE
        :param df:
        :return:
        '''
        tmp = df.copy()
        r = [x for x in range(3)]
        grouped = tmp.groupby(['target_name'])
        grouped = grouped.agg({'score_full_seq': 'max'})
        grouped = grouped.rename(columns={'score_full_seq': 'max_score'})
        tmp = pd.merge(tmp, grouped, how='left', on=['target_name'])
        tmp = tmp[tmp['score_full_seq'] == tmp['max_score']]

        return tmp.copy()

    @staticmethod
    def get_sample_and_contig_margins_dict(df: pd.DataFrame)-> dict:
        df.groupby(temp = df.groupby(['sample_names', 'contig_id']).agg({'gene_id': ['min', 'max']})
    )

    @staticmethod
    def exract_cassetes_from_group(grp: pd.DataFrame, threshold: int = 1) -> list:
        # contig contains only 1 gene with high score, hence it does not contain a cassete
        # grp.to_csv("group.csv")
        # if (grp.shape[0] <= 1):
        #     print("Contig with 1 hit found  ")
        #     return []

        # Get if here if the contig has more than 1 gene, hence it should be checked
        contig_id = grp.contig_id
        gene_list = grp.gene_id.values
        gene_dict = {int(g): [g] for g in gene_list}

        # Remove lonley genes
        gene_dict = process_hmm.remove_lonley_genes(gene_dict, gene_list, threshold)

        # Build cassetes
        gene_dict = process_hmm.build_cassete(gene_dict, gene_list, threshold)

        return [g for g in gene_dict.values() if len(g) > 1]

    @staticmethod
    def remove_lonley_genes(gene_dict: dict, gene_list: list, threshold: int = 1) -> dict:
        for gene in gene_list:
            gene = int(gene)
            flag = 0
            i = 1

            # Check each gene's surroundings
            while (flag == 0 and i < threshold + 1):

                if ((gene + i) in gene_dict) or ((gene - i) in gene_dict):
                    flag = 1 #not lonley, turn flag on
                i += 1

            # Flag is on meaning gene is not lonley, so don't pop it
            if (flag):
                continue
            gene_dict.pop(gene)

        return gene_dict

    @staticmethod
    def build_cassete(gene_dict: dict, gene_list: list, threshold: int = 1) -> dict:

        #TODO check if the list is sorted
        #Recive a list of known genes (scored high enough to stay)
        #Recieve a dictionary, contianing each gene and itself only
        #transfrom genes to int
        for gene in gene_list:
            g = int(gene)
            flag = 0
            i = g + 1

            #check for the gene's surrounding
            #if there is another known gene under range of {thrshold} add it to the dictionary and then pop it out
            #items never pop out of the list, only from the dictionary, thus we can still check the surroundings of a gene

            while (flag == 0 and i < g + threshold + 1):
                if i in gene_dict:
                    gene_dict[i] = gene_dict[i] + gene_dict[g]
                    gene_dict.pop(g)
                    flag = 1
                i += 1

        return gene_dict

    def get_architecture(df: pd.DataFrame, contig: str, genes_list: list) -> dict:
        genes_list.sort()
        dict = {}
        for gene in genes_list:
            mask = (df['contig_id'] == contig) & (df['gene_id'] == gene)
            # & (df['gene_id'] == gene)
            tmp = df.loc[mask, :]
            profile_name = tmp['query_name'].to_string()
            dict[gene] = profile_name
            # print(type(tmp.query_name.astype(str)))

        return dict

    @staticmethod
    def structure_to_str(struct: list) -> str:
        res = ""
        for word in struct:
            res += f"__({word})"

        res += '__'
        return res


    '''
    index = 1: return members
    index = 2: return architecture
    index = 3: return specific score
    index = 4: return specific eval
    '''
    def extract_architectures(df: pd.DataFrame, contig_min_max:dict, sample_contig: str, contig_genes_lst: list,
                              index: int = 1, k: int = 0) -> str:
        df_local = df.copy()
        # index = str(index)
        members = []
        architecture = []
        scores_specific = []
        e_val_specific = []

        min_gene, max_gene = min(contig_genes_lst), max(contig_genes_lst)
        for j in range(max(1, min_gene - k), max_gene + 1 + k):

            # print("contig min/max dict is : \n", contig_min_max)
            # margin_cond = j > = 1 #Todo take care of maximum of each contig
            # if not (margin_cond):
            #     continue

            # If no such gene was identified, add '.?.' marker to the specific structure list
            if (j not in contig_genes_lst):
                architecture.append('??')
                scores_specific.append('??')
                e_val_specific.append('??')

            else:
                '''Gene number k was identified and scored as a hit, 
                extract it's query name from the dataframe and append it to both lists '''
                mask = (df_local.sample_contig == sample_contig) & (df_local['gene_id'] == j)
                tmp = df_local.loc[mask, :].reset_index()

                # tmp = df_local.loc[df_local.contig_id == contig_id]


                profile = tmp.at[0, 'query_name']
                members.append(profile)
                architecture.append(profile)
                eval = tmp.at[0, 'eval_full_seq']
                score = tmp.at[0, 'score_full_seq']
                scores_specific.append(str(score))
                e_val_specific.append(str(eval))

        '''sort the general structure of the cassete to get a general desciption of the included identified proteins'''
        members.sort()

        architecture = process_hmm.structure_to_str(architecture)
        scores_specific = process_hmm.structure_to_str(scores_specific)
        e_val_specific = process_hmm.structure_to_str(e_val_specific)

        res_dict = {}
        res_dict[1] = members
        res_dict[2] = architecture
        res_dict[3] = scores_specific
        res_dict[4] = e_val_specific

        return res_dict[index]

    @staticmethod
    def hits_df_to_structural_df(df: pd.DataFrame, min_max_dict:dict
                                 , inner_spacing: int = 3, margin_spacing: int = 1) -> pd.DataFrame:

        grp = df.groupby(['sample_names', 'contig_id'])

        # n - maximal space between inner genes
        inner_spacing = 3
        # k - maximal space before and after first and last genes
        margin_spacing = 1

        contigs_dict = {}
        curr_cassete = []

        for name, group in grp:
            curr_cassete = process_hmm.exract_cassetes_from_group(group, inner_spacing)
            if (curr_cassete == None or curr_cassete == []):
                continue
            contigs_dict[name] = curr_cassete

        result_df = pd.DataFrame()
        # Get all general structures

        data = []
        # iterate throughout all groups
        for name, group in grp:
            sample_id, contig_id = name[0], name[1]
            if (name not in contigs_dict):
                continue
            # Iterate through all intresting chunks in the contig
            for chunk in contigs_dict[name]:
                chunk.sort(reverse=False)
                left_margin, right_margin = max(min(chunk) - margin_spacing, 1), min(max(chunk) + 1 + margin_spacing, df.shape[0] + 1)

                # Get extract_architecture for the current chunk of current contig in a specific sample
                members = process_hmm.extract_architectures(df, min_max_dict, contig_id, chunk, 1, margin_spacing)
                arch = process_hmm.extract_architectures(df, min_max_dict, contig_id, chunk, 2, margin_spacing)
                score = process_hmm.extract_architectures(df, min_max_dict, contig_id, chunk, 3, margin_spacing)
                eval = process_hmm.extract_architectures(df, min_max_dict, contig_id, chunk, 4, margin_spacing)
                cassete_length = len(np.arange(left_margin, right_margin))

                mask = (df['sample_names'] == sample_id) & (df['gene_id'] == min(chunk)) & (
                        df['contig_id'] == contig_id)
                full_id = df.loc[mask, :].reset_index().at[0, 'target_name']

                row = [full_id,  # full hit id
                       contig_id,
                       # contig_id of the first gene of the first hit, currently we assume they all belong to one contig
                       len(chunk),  # num of known genes
                       chunk,  # gene numbers known members of the hit
                       [x for x in range(left_margin, right_margin) if x not in chunk],
                       # gene numbers of unknown members of the hit
                       cassete_length,
                       ','.join(members),  # Known group members of the hit, by name, transformed to string
                       arch,  # Specific architecture of the group, by name
                       score,  # scores of known genes, adjust to architecture's shape
                       eval  # scores of known genes, adjust to architecture's shape
                       ]

                data.append(row)

        final_df = pd.DataFrame(data=data,
                                columns=['full_id', 'contig_id', '#members', 'known_gene_ids', 'unknown_gene_ids',
                                         'span',
                                         'members', 'architecture', 'member_scores', 'members_evals'])

        return final_df.copy()
#

    def remove_single_hit_contigs(df: pd.DataFrame) -> pd.DataFrame:

        local_df = df.copy()
        #group dataframe by samples and contigs
        #each group is a specific contig in a specific sample


        #remove groups with 1 gene in the contig, as they are not interesting
        local_df = df.copy()
        sample_contig = pd.Series(local_df.target_name).str.split('_')
        def f(x): return x[0] + '_' + x[1] + '_' + x[2]
        sample_contig = sample_contig.apply(lambda x: f(x))
        local_df['sample_contig'] = sample_contig

        grp = local_df.groupby(['sample_contig'])

        group_size_mask = local_df['sample_contig'].groupby(local_df['sample_contig']).transform('size')
        return local_df[group_size_mask > 1]




    def hits_df_to_structural_df2(df: pd.DataFrame, min_max_dict: dict
                                  , inner_spacing: int = 3, margin_spacing: int = 1) -> pd.DataFrame:


        local_df = df.copy()

        # TODO rewrite function in a more simple manner by creating a pseudo column sample_name+contig_id and grouping by it
        grp = local_df.groupby(['sample_contig'])

        contigs_dict = {}
        curr_cassete = []

        for name, group in grp:
            curr_cassete = process_hmm.exract_cassetes_from_group(group, inner_spacing)
            if (curr_cassete == None or curr_cassete == []):
                continue
            contigs_dict[name] = curr_cassete

        # Get all general structures

        data = []

        # iterate throughout all groups
        for name, group in grp:
            sample_contig = name
            # sample_id, contig_id = name[0], name[1]
            contig_margins = min_max_dict[name]

            #TODO take care of minimum and maximum of the contig
            #TODO


            if (name not in contigs_dict):
                continue

            # Iterate through all intresting chunks in the contig
            for chunk in contigs_dict[name]:
                chunk.sort(reverse=False)
                left_margin, right_margin = max(min(chunk) - margin_spacing, 1), min(max(chunk) + 1 + margin_spacing,
                                                                                     local_df.shape[0] + 1)


                # Get extract_architecture for the current chunk of current contig in a specific sample
                members = process_hmm.extract_architectures(local_df, contig_margins, sample_contig, chunk, 1, margin_spacing)
                arch = process_hmm.extract_architectures(local_df, contig_margins, sample_contig, chunk, 2, margin_spacing)
                score = process_hmm.extract_architectures(local_df, contig_margins, sample_contig, chunk, 3, margin_spacing)
                eval = process_hmm.extract_architectures(local_df, contig_margins, sample_contig, chunk, 4, margin_spacing)
                cassete_length = len(np.arange(left_margin, right_margin))

                mask = (local_df.sample_contig == sample_contig) & (local_df['gene_id'] == min(chunk))
                full_id = local_df.loc[mask, :].reset_index().at[0, 'target_name']


                row = [full_id,  # full hit id
                       sample_contig.split('_')[2],
                       # contig_id of the first gene of the first hit, currently we assume they all belong to one contig
                       len(chunk),  # num of known genes
                       chunk,  # gene numbers known members of the hit
                       [x for x in range(left_margin, right_margin) if x not in chunk],
                       # gene numbers of unknown members of the hit
                       cassete_length,
                       ','.join(members),  # Known group members of the hit, by name, transformed to string
                       arch,  # Specific architecture of the group, by name
                       score,  # scores of known genes, adjust to architecture's shape
                       eval, # scores of known genes, adjust to architecture's shape
                       sample_contig #add a sample_contig data for further rebuilding of unknonwn genes for pull seq
                       ]

                data.append(row)

        final_df = pd.DataFrame(data=data,
                                columns=['full_id', 'contig_id', '#members', 'known_gene_ids', 'unknown_gene_ids',
                                         'span',
                                         'members', 'architecture', 'member_scores', 'members_evals', 'sample_contig'])

        return final_df


    @staticmethod
    def fetch_unknown_genes(df: pd.DataFrame) -> list:
        '''
        Receive a dataframe of format 'final' (as outputted by 'hits_df_to_structural_df2' function.
        for each sample and contig, extract gene numbers from the unknown genes column into a list.
        for each unknown gene in the list, format it's name back to the original fasta contigs file format.
        function does not change the original df, instead it creates a local copy.
        :param df:
        Original dataframe contating data
        :return:
        A list of gene names that can be used to pull genes from a fasta file.
        '''
        #TODO consider returnining a dataframe
        unknown_genes_list = []
        local_df = df.copy()
        for index, row in local_df.iterrows():
            full_id = row.full_id

            min_gene = min(row.known_gene_ids)  # This is also the number of digits that has to be replaced by a an unknwon gene id
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

        return unknown_genes_list

    @staticmethod
    def write_unknown_genes(unknown_genes_list: list, output_path : str = 'txt') -> int:
        '''
        Recive a list of genes to write to ouput file.
        The fucntion assumes nothing about gene list.
        :param unkonwn_genes_list: list of gene names to be written.
        :param output_path: desired path and name of output file, default value is 'txt'
        :return:
        0 upon success
        Raises 'ValueError' if list is empty.
        other values to be added in the future.
        '''
        #TODO handle writer exceptions.

        if unknown_genes_list == [] or unknown_genes_list == None:
            raise ValueError

        with open(output_path, 'w') as f:
            for gene in unknown_genes_list:
                f.write("%s\n" % gene)

        return 0

    @staticmethod
    def process(files_path: str, out_path : str, file_suffix: str = 'hmmertbl', n : int = 5 ) -> int:
        '''

        :param files_path: path of folder where 'hmmer search tblout'  where files are located.
        :param out_path: path of output folder, files will be written to this folder.
        :param file_suffix: optional, input files suffix. default value is hmmertbl, this suffix will help to
        exclude all the hmmer search output files from the other files in the folder.
        :param n: maximum distance between two knonwn genes. all genes within this might be included in the output file.
        default value is 5.
        :return:
        0 upon success.
        Raises ValueError if parameters are missing.
        '''
        #TODO take care of empty out_path, in this case write output files to source directory
        #TODO take care of empty in_path

        if( files_path == '' or files_path == None):
            raise ValueError("Input path is empty or None.")

        if (out_path == '' or out_path == None):
            raise ValueError("Output path is empty or None.")

        '''
        #Set project path to collect sample files from, set suffix, set out folder path,
        #collect all in files {path} ending with {suffix} into {files_list}
        #Set a data object, for later use in file writing
        '''

        date = datetime.today().strftime('%Y_%m_%d')

        path = os.path.join("~", files_path)
        files_lst = glob.glob(os.path.join(path, f"*.{file_suffix}"))


        #Create a pandas dataframe from all tblout files, forece numeric data types
        df = process_hmm.multiple_tbl_to_pd(files_lst)
        # Force numertic types on data
        df = process_hmm.force_data_types(df.copy())
        # display(df.shape)


        '''
        Split the sample names to add sample_name, contig_id, gene_id  columns
        Get the best hits for each gene in the current dataframe
        Remove contigs that have only one hit, they have no potential
        Possibly write the filtered dataframe out as a csv for further inspection or use
        '''

        # Enrich data with new columns -> sammple_name, ctg_id, gene_id, , it is possible to supply a diffrent separator. defualt is '_'
        df = process_hmm.split_sample_name(df.copy())

        # get_best_hits for each gene
        df_filtered = process_hmm.get_best_hits(df)

        # remove contigs with only one hit, they are not significant enough
        df_filtered = process_hmm.remove_single_hit_contigs(df_filtered)

        # write hits df to csv for further inspection if needed
        df_filtered.to_csv(f"{out_path}/df_filtered_{date}.csv")

        '''
        Create a pandas dataframe from all tblout files, forece numeric data types 
        '''
        df = process_hmm.multiple_tbl_to_pd(files_lst)
        # Force numertic types on data
        df = process_hmm.force_data_types(df.copy())
        # display(df.shape)

        '''
        Split the sample names to add sample_name, contig_id, gene_id  columns
        Get the best hits for each gene in the current dataframe
        Remove contigs that have only one hit, they have no potential
        Possibly write the filtered dataframe out as a csv for further inspection or use
        '''

        # Enrich data with new columns -> sammple_name, ctg_id, gene_id, , it is possible to supply a diffrent separator. defualt is '_'
        df = process_hmm.split_sample_name(df.copy())

        # get_best_hits for each gene
        df_filtered = process_hmm.get_best_hits(df)

        # remove contigs with only one hit, they are not significant enough
        df_filtered = process_hmm.remove_single_hit_contigs(df_filtered)

        # write hits df to csv for further inspection if needed
        df_filtered.to_csv(f"{out_path}/df_filtered_{date}.csv")

        '''
        Construct the final dataframe, set {n} to the maximal space between two genes
        possibly write it out as a csv for further inspection
        '''

        #TODO Take this out of the functions as well
        min_max_dict = process_hmm.get_min_max_dict(df)

        final_df = process_hmm.hits_df_to_structural_df2(df_filtered, min_max_dict, n)

        # write df to csv
        final_df.to_csv(f"{out_path}/final_df_{date}_n{n}.csv")

        '''
        Prepare output for pullseq make a list of unknown genes - ORF ID per line
        Then write it to a csv file.
        '''

        unknown_genes_list = process_hmm.fetch_unknown_genes(final_df)
        status = process_hmm.write_unknown_genes(unknown_genes_list, f"{out_path}/unknown_genes_ids_{date}.csv")

        return status
