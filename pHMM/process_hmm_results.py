import pandas as pd
import numpy as np
import os
from datetime import datetime
from general_utils import *

out_path = os.getcwd()
date = datetime.today().strftime('%Y_%m_%d')
VERBOSE = True

def get_best_hits(df: pd.DataFrame) -> pd.DataFrame:
    """
    Recives a dataframe {df} of tblout outputs and returns a dataframe with the best (1) hit for each gene.
    Best hit is determined by the 'score_full_seq' coulmn value as appears in the hmmsearch tblout output.
    :param df:
    Pandas dataframe containing original data
    :return:
    A new dataframe containg *only* the best hit for each gene.
    """
    if(VERBOSE):
        FUNC_START()

    local_df = df.copy()
    local_df = local_df.sort_values('score_full_seq', ascending=False).drop_duplicates('target_name')
    local_df.reset_index(inplace=True, drop=True)  # create a best hits dataframe

    if (VERBOSE):
        FUNC_END()
    return local_df


def exract_cassetes_from_group(grp: pd.DataFrame, spacing: int = 5) -> list:
    """
    Recieve a small dataframe which is a result of a groupby('sample_contig') and build a dictionary that contains genes
    that are the minimal genes in an ARG cassette. E.G dict['gene_001']: ['gene_001', 'gene_004', 'gene_007'].
    Genes that have no other identified ARGs in distance of less than {n} from them, will not be included in the
    dictionary.
    :param grp: the dataframe of the current contig.
    :param spacing: maximal distance between two identified ARGs, if a gene has no identified ARGs in distance of less
    than {n} it will be defined as lonely and later removed.
    :return: a dictionary of genes and their cassettes of other identified genes.
    """
    # When we get here we are already grouped by sample_contig which includdes both the sample name and the contig id
    # Get if here if the contig has more than 1 gene, hence it should be checked
    contig_id = grp.contig_id

    gene_list = grp.gene_id.values
    gene_dict = {int(g): [g] for g in gene_list}

    # Remove lonley genes
    gene_dict = remove_lonley_genes(gene_dict, gene_list, spacing)

    #reset list
    gene_list = list(gene_dict.keys())
    # Build cassetes
    gene_dict = build_cassete(gene_dict, gene_list, spacing)

    return [g for g in gene_dict.values()]

def remove_lonley_genes(gene_dict: dict, gene_list: list, spacing: int = 5) -> dict:
    """
    recieve a gene list (nums of genes, e.g [1, 5, 25] and a gene dictionary of the format {x : [x]} for each gene number
    in the gene list. iterate through the gene list and remove genes that are lonley, meaning that there are no other
    genes in the list within distance smaller than {spacing} from them.
    :param gene_dict: a dictionary of the format {x : [x]} where 'x' is a number of a gene in gene_list.
    :param gene_list: a list containing numbers (integers) of genes.
    :param spacing: maximal distance between two genes in the list, otherwise they are treated as lonley and removed.
    :return: the same dictionary, but without the numbers of the lonely genes.
    """
    gene_list.sort()
    n = len(gene_list)

    for j in range(1, n-1):
        gene = gene_list[j]

        if(gene_list[j] - gene_list[j-1] > spacing and gene_list[j+1]-gene_list[j]) > spacing:
            gene_dict.pop(gene)

    #Take care of special cases of the last and first genes
    if gene_list[1]-gene_list[0] > spacing:
        gene_dict.pop(gene_list[0])

    if gene_list[-1] - gene_list[-2] > spacing:
        gene_dict.pop(gene_list[-1])

    return gene_dict



def build_cassete(gene_dict: dict, gene_list: list, spacing) -> dict:
    """
    build a cassette of genes in the format of {x: [x, y, z, w...]} where y, z, w follow the rule: z + spacing > w
    meaning that each gene in the cassette is no more than {spacing} genes away from the genes that are before and it.
    :param gene_list: a list of genes, e.g [1, 5, 8 ..], the list is sorted.
    :param gene_dict: a dictionary of gene numbers in the format of {x:[x]} for each gene in gene list. it is assumed
    that there are no lonley genes in the cassette (they were previously removed by the remove_lonley_genes func.
    :param spacing: maximal distance between two genes in the list, this is used to create the cassettes, a gene could be
    more than {spacing} away from gene1 but if it's in the list above then it belongs to the cassette of the next gene
    in  the list because all lonley genes were already removed and the list is sorted.
    :return: a dictionary representing cassettes of genes where the minimal gene in the cassette is they key and
    the rest of the genes (including the key) are the values, wrapped as a list.
    e.g  {x: [x, y, z, w...]}
    """
    # Recive a list of known genes (scored high enough to stay)
    # Recieve a dictionary, contianing each gene and itself only
    # transfrom genes to int
    gene_list.sort()
    n = len(gene_list)

    for i in range(n-1):

        curr_gene = gene_list[i]
        next_gene = gene_list[i + 1]

        # check if the next gene is less than {sapcing} away, if it is add it to the cassette of the current gene and pop
        # so it will not have a cassette of itself.
        if(next_gene-curr_gene <= spacing):
            gene_dict[next_gene] = gene_dict[next_gene]+gene_dict[curr_gene]
            gene_dict.pop(curr_gene)

    return gene_dict


def structure_to_str(struct: list) -> str:
    '''
    A helper fucntion to create the visual representation of the ARG cassettes.
    Used for diffrent data types, both description of a gene, or it's e-value etc..
    :param struct: a list of the known members in the cassette
    :return: a string representation of the format __(gene1_data)__(gene2_data)__...
    '''
    res = ""
    for word in struct:
        res += f"__({word})"

    res += '__'
    return res


def get_cassete_structure(contig_group: pd.DataFrame, chunk: list, outer_spacing: int = 2):
    '''

    :param contig_group: a pandas dataframe
    :param chunk:
    :param outer_spacing:
    :return:
    '''
    res_dict = {}
    members = []
    architecture = []
    scores_specific = []
    e_val_specific = []

    #TODO take better care of min/max genes
    left_margin, right_margin = max(min(chunk) - outer_spacing, 1), max(chunk) + 1 + outer_spacing


    chunk_min, chunk_max = min(chunk), max(chunk)
    for j in range(max(1, chunk_min - outer_spacing),
                   chunk_max + 1 + outer_spacing):

        if (j not in chunk):
            architecture.append('??')
            scores_specific.append('??')
            e_val_specific.append('??')

        else:
            '''Gene number k was identified and scored as a hit, 
            extract it's query name from the dataframe and append it to both lists '''
            mask = contig_group.gene_id == j
            tmp = contig_group[mask].iloc[0,:]

            profile = tmp.query_name
            members.append(profile)
            architecture.append(profile)
            eval = tmp.eval_full_seq
            score = tmp.score_full_seq
            scores_specific.append(str(score))
            e_val_specific.append(str(eval))


    members.sort()

    architecture = structure_to_str(architecture)
    scores_specific = structure_to_str(scores_specific)
    e_val_specific = structure_to_str(e_val_specific)

    res_dict['members'] = ','.join(members)
    res_dict['#members'] = len(chunk)
    res_dict['architecture'] = architecture
    res_dict['member_scores'] = scores_specific
    res_dict['members_evals'] = e_val_specific
    res_dict['span'] = len(np.arange(left_margin, right_margin))
    res_dict['contig_id'] = tmp.contig_id
    res_dict['sample_contig'] = tmp.sample_contig
    res_dict['known_gene_ids'] = chunk
    res_dict['unknown_gene_ids'] = [x for x in range(left_margin, right_margin) if x not in chunk]

    return res_dict

def remove_single_hit_contigs(df: pd.DataFrame) -> pd.DataFrame:
    """
    Receive a dataframe group it by sample_contigs (accounts for grouping by contig numbers) and remove contigs that
    only have 1 hit. contigs with only 1 hit cannot have cassette in them because it requires at list two identified
    ARGs to make a cassette.
    :param df: containing all genes, must have a column of 'sample_contig' and 'gene_id'.
    :return: a copy of the same dataframe, all contigs in the new dataframe have at list two identified genes in them.
    """

    if (VERBOSE):
        FUNC_START()
    # group dataframe by samples and contigs
    # each group is a specific contig in a specific sample
    # remove groups with 1 gene in the contig, as they are not interesting
    df_filtered = df.groupby('sample_contig').filter(lambda subgroup: len(subgroup) > 1)
    df_filtered = df_filtered.sort_values(['sample_contig', 'gene_id'])
    df_filtered = df_filtered.reset_index(drop=True)
    if (VERBOSE):
        FUNC_END()
    return df_filtered



def hits_df_to_structural_df(df: pd.DataFrame,
                             inner_spacing: int = 5) -> pd.DataFrame:
    """
    receive a dataframe of single hits, that has all the single hit contigs removed from it.
    construct a dataframe in which is entry is a gene representing a whole ARG cassette including the knonwn and unknown
    genes in the cassette.
    :param df: a dataframe containing all genes, must have a column of 'sample_contig' and 'gene_id'.
    :param inner_spacing: maximal distance between two known genes for them to belong to the same ARG cassette.
    :return: a dataframe of ARG cassettes.
    with the columns:
        1. 'full_id' - the 'target_name' of the cassette represnting gene
        (the minimal gene in the cassette, this is arbitrary).
        2. 'contig_id' - the contig number in which all the cassette members are located in.
        3. '#members' - number of identified ARGs, this is the number of genes in the cassette that had hit in the
        hmmsearch run.
        4. 'known_gene_ids' - gene numbers (int) of the known genes (genes that scored a hit in the hmmsearch).
        5. 'unknown_gene_ids' - gene numbers (int) of the unknown genes that are found between two known genes in the
        same cassette, or in the cassette's margins.
        6. 'span' - number of genes (known and unknown) in the cassette.
        7. 'members' - description of the cassette's known genes as it appears in the raw data.
        8. 'architecture' - a visual representation of the cassette architecture of the format
        __(known_gene1_description)__(??)__(??)__(known_gen4_description)__(??)__(??)__
        the number of (??) is no more than {spacing}//2 on ther margins, and no more than {spacing} in between 2 known
         ARGs.
        9. 'member_scores' - visual representation of the architecture, in the same format as above
        only the instead of description, what appears is  __(known_gene1_hmmer-seq-score)__
        10. 'members_evals' - same as above only with known member (identified ARG) e-value from the hmmsearch
        11. 'sample_contig' - a unique identifier of the contigs of all samples in the current cassette representing
        the contig they all belong to.
    """

    if (VERBOSE):
        FUNC_START()

    grp = df.groupby(['sample_contig'])
    contigs_dict = {}

    outer_spacing = inner_spacing // 2

    for name, group in grp:
        cassetes_list = exract_cassetes_from_group(group, inner_spacing)
        # curr_contig = grp.contig_id
        if (cassetes_list == None or cassetes_list == []):
            continue
        contigs_dict[name] = cassetes_list

    # Get all general structures

    csts_dict = {}

    # iterate throughout all groups
    for name, group in grp:
        if (name not in contigs_dict):
            continue

        # Iterate through all intresting chunks in the contig
        for chunk in contigs_dict[name]:
            chunk.sort(reverse=False)
            mask = group.gene_id == min(chunk)

            full_id = group[mask].iloc[0,:].target_name
            csts_dict[full_id] = get_cassete_structure(group, chunk, outer_spacing)

            # print(f"Appended chunk {chunk} from the row: {full_id}\n")

    cols_order = ['full_id', 'contig_id', '#members', 'known_gene_ids', 'unknown_gene_ids', 'span',
               'members', 'architecture', 'member_scores', 'members_evals', 'sample_contig']

    final_df = pd.DataFrame.from_dict(csts_dict, orient='index')
    final_df = final_df.reset_index().rename(columns = {'index':'full_id'})
    final_df = final_df[cols_order]

    if (VERBOSE):
        FUNC_END()
    return final_df




def assign_min_max_genes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Receive a dataframe and group it by sample_contigs (accounts for grouping by contig numbers) add two new columns for
    for each entry - the minimal and maximal genes in the contig the gene (of the entry) belongs to.
    :param df: dataframe with all samples in it, must have a sample_contig column to groupby
    :return: dataframe with two new columns for each sample that contain an integer representing the minimal and
    maximal gene number in this contig
    """

    df['min_gene_for_contig'] = df.groupby('sample_contig')['gene_id'].transform(min)
    df['max_gene_for_contig'] = df.groupby('sample_contig')['gene_id'].transform(max)

    return df

def get_min_max_dict(df: pd.DataFrame) -> pd.DataFrame:
    """
    Receive a dataframe and group it by sample_contigs (accounts for grouping by contig numbers) add two new columns for
    for each entry - the minimal and maximal genes in the contig the gene (of the entry) belongs to.
    :param df: dataframe with all samples in it, must have a sample_contig column to groupby
    :return: a dictionary where the key is the sample name ('target_name') and the value is a another dictionary of the
    format {'min': x, 'max': y} where x and y are integers.
    """

    temp = df.groupby(['sample_contig']).agg({'gene_id': ['min', 'max']})
    min_max_dict = {}
    for row in temp.iterrows():
        # print(row[0])
        # print("Row max\min are : ",  row[1][0],  row[1][1])
        min_max_dict[row[0]] = dict({'min': row[1][0], 'max': row[1][1]})
    return min_max_dict

def process(results_df : pd.DataFrame, n: int = 5, debug: bool = True, VERBOSE_FLAG = False) -> pd.DataFrame:
    '''

        Process a group (>=1) of files of the tblout format. Stages:
        1. Collect all files into a Pandas (pd) dataframe.
        2. Get best hits for each gene, get the best hit.
        3. form cassetes of size 'n' (determined by user, default 5) of known genes.
        4. Prepare output in a pullseq compatible format.

        notes:
        1. in stage (2), Contigs with a single hits are removed as they cannot have a cassete.
        2. All input files must have the same suffix (weather default or user determined).

        output:
        if debug = 1
        cwd/df_filtered_{date}.csv will be written to cwd (current working directory) after stage 2.
        cwd/final_df_{date}_n{n}.csv will be written to cwd after stage 3.

        else - only this is outputted to cwd:
        cwd/unknown_genes_ids_{date}.csv is written after stage 4.


        :param files_path: path of folder where 'hmmer search tblout'  where files are located.
        :param file_suffix: optional, input files suffix. default value is hmmertbl, this suffix will help to
        exclude all the hmmer search output files from the other files in the folder.
        :param n: maximum distance between two knonwn genes. all genes within this might be included in the output file.
        default value is 5.
        :param debug: A boolean parameter, by defulat is True, if set to false, no intermidate tables will be written
        to cwd.

        :return: 0 upon success.



        :raise ValueError: if parameters are missing.
        '''
    global VERBOSE
    VERBOSE = VERBOSE_FLAG

    if (VERBOSE):
        FUNC_START()

    print("Processing hmm results dataframe.\n"\
          f"Output folder is: {out_path}.")

    df_filtered = get_best_hits(results_df)
    df_filtered = remove_single_hit_contigs(df_filtered)
    df_filtered = assign_min_max_genes(df_filtered)

    # remove contigs with only one hit, they are not significant enough
    # df_filtered = remove_single_hit_contigs(df_filtered).drop(['level_0', 'index'], axis=1).reset_index(drop=True)

    # write hits df to csv for further inspection if needed
    if debug:
        df_filtered.to_csv(f"{out_path}/df_filtered_{date}.csv")

    '''
    Construct the final dataframe, set {n} to the maximal space between two genes
    possibly write it out as a csv for further inspection
    '''

    # min_max_dict = get_min_max_dict(df_filtered)

    final_df = hits_df_to_structural_df(df_filtered, n)

    # write df to csv
    if (debug):
        final_df.to_csv(f"{out_path}/final_df_{date}_n{n}.csv")

    '''
    Prepare output for pullseq make a list of unknown genes - ORF ID per line
    Then write it to a csv file.
    '''
    print(f"Finished processing results. files written:\n" 
          f"df_filtered_{date}.csv - Contains best results but unstructured cassetes.\n"
          f"final_df_{date}_n{n}.csv - Contains sturctured cassetes.\n"
          f"Retuning final_df")

    if(VERBOSE):
        FUNC_END
    return final_df