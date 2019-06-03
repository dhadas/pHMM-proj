from __future__ import annotations
import glob
import os
import numpy as np
import pandas as pd


class process_hmm():

    # Class column definitions
    @staticmethod
    def get_headers(type: int = 0) -> dict:
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
        headers = process_hmm.get_headers(1)
        # read file to dataframe with a single column containing all data
        df = pd.read_csv(path, header=None, index_col=0, sep='\t')

        # drop lines statrting with # from data
        mask_idx = df.index.map(lambda x: x.startswith("#") == False)
        masked_df = df.loc[mask_idx, :].reset_index().copy().rename(columns={0: 'col'})



        # Seperate rows by spaces, drop orginial 'col'
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


    def multiple_tbl_to_pd(files_lst: list ) -> pd.DataFrame:


        df_lst = []

        for hmm_out_file in files_lst:
            tmp = process_hmm.tblout_to_pd(hmm_out_file)
            df_lst.append(tmp)

        df = pd.concat(df_lst)
        df = df.reset_index()

        return df

    def force_data_types(df: pd.DataFrame) -> pd.DataFrame:
        headers = process_hmm.get_headers()
        # Force data typesge for the cells of the DataFrame for later manipulations
        df.loc[:, headers['full_seq_cols'] + headers['best_dom_cols']] = df.loc[:, headers['full_seq_cols'] + headers[
            'best_dom_cols']].astype(float)
        df.loc[:, "exp"] = df.loc[:, "exp"].astype(float)
        df.loc[:, headers['misc_cols'][1:-1]] = df.loc[:, headers['misc_cols'][1:-1]].astype(int)
        df.loc[:, headers['misc_cols'][:-1]] = df.loc[:, headers['misc_cols'][:-1]].astype(str)

        return df.copy()

    def split_sample_name(df: pd.DataFrame, separator: str = '_') -> pd.DataFrame:


        cols = {0: 'sample_names', 1: 'ctg', 2: 'contig_id', 3: 'gene_id'}
        local_df = df.copy().join(df['target_name'].str.split(separator, expand=True)).rename(columns=cols)
        local_df['gene_id'] = local_df['gene_id'].astype(int)

        local_df['contig_id'] = 'ctg_' + local_df.contig_id
        local_df.contig_id = local_df.contig_id.astype(str)
        local_df = local_df.drop(['ctg', 'index'], axis=1)

        return local_df

    ''' Works well but Imporve later'''

    def get_best_hits(df: pd.DataFrame) -> pd.DataFrame:
        local_df = df.copy()
        local_df = local_df.sort_values('score_full_seq', ascending=False).drop_duplicates('target_name')
        local_df.reset_index(inplace=True)  # create a best hits dataframe

        return local_df

    ''' Improve later'''

    def get_best_hits2(df: pd.DataFrame) -> pd.DataFrame:
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


class assembler():

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
        gene_dict = assembler.remove_lonley_genes(gene_dict, gene_list, threshold)

        # Build cassetes
        gene_dict = assembler.build_cassete(gene_dict, gene_list, threshold)

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

        architecture = assembler.structure_to_str(architecture)
        scores_specific = assembler.structure_to_str(scores_specific)
        e_val_specific = assembler.structure_to_str(e_val_specific)

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
            curr_cassete = assembler.exract_cassetes_from_group(group, inner_spacing)
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
                members = assembler.extract_architectures(df, min_max_dict, contig_id, chunk, 1, margin_spacing)
                arch = assembler.extract_architectures(df, min_max_dict, contig_id, chunk, 2, margin_spacing)
                score = assembler.extract_architectures(df, min_max_dict, contig_id, chunk, 3, margin_spacing)
                eval = assembler.extract_architectures(df, min_max_dict, contig_id, chunk, 4, margin_spacing)
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
            curr_cassete = assembler.exract_cassetes_from_group(group, inner_spacing)
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
                members = assembler.extract_architectures(local_df, contig_margins, sample_contig, chunk, 1, margin_spacing)
                arch = assembler.extract_architectures(local_df, contig_margins, sample_contig, chunk, 2, margin_spacing)
                score = assembler.extract_architectures(local_df, contig_margins, sample_contig, chunk, 3, margin_spacing)
                eval = assembler.extract_architectures(local_df, contig_margins, sample_contig, chunk, 4, margin_spacing)
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

    import regex

    def get_max_gene(path : str) -> int:
        #Get a path of contig file and return the id of it's last gene
        with open("/Users/dror/PycharmProjects/Lab/power/SRR5976181.faa") as f:
            contig = f.read()
            separated_contig = contig.split('>N')[1:-1]
            genes_list = []
            for gene_block in separated_contig:

                header = gene_block.split("\n")[0]
                try:
                    gene_id2 = header.split()[0]
                    gene_id = gene_id2.split('_')[1]

                except IndexError:
                    print(gene_id2)
                #         display(int(gene_id))
                genes_list.append(int(gene_id))
            n = len(genes_list)

        return genes_list[n]

# suffix = "faa.tbl"
# proj_path = "/Users/dror/PycharmProjects/Lab/power/resFamHighSpec_vs_pigs_130219/"
#
# path = os.path.join("~",proj_path)
# files_lst = glob.glob(os.path.join(path, f"*.{suffix}"))
# df = process_hmm.multiple_tbl_to_pd(files_lst)
#
#
#
# suffix = "faa.tbl"
# proj_path = "/Users/dror/PycharmProjects/Lab/power/resFamHighSpec_vs_pigs_130219/"
#
# path = os.path.join("~",proj_path)
# files_lst = glob.glob(os.path.join(path, f"*.{suffix}"))
# df = process_hmm.multiple_tbl_to_pd(files_lst)
# df = process_hmm.force_data_types(df.copy())
# df = process_hmm.split_sample_name(df.copy())
#
# min_max_dict = process_hmm.get_min_max_dict(df.copy())
#
#
# df_filtered = pd.read_csv("data/cassetes_dbs/df_filtered_2019_04_20.csv")
# # df = assembler.hits_df_to_structural_df2(df_filtered, min_max_dict)
