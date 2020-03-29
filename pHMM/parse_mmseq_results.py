import pandas as pd
from Bio import SeqIO
from datetime import datetime
import os
from general_utils import FUNC_START, FUNC_END

out_path = os.getcwd()
date = datetime.today().strftime('%Y_%m_%d')

# out_path = '/Users/dror/PycharmProjects/AMR_pipeline/data/pHMM_out'
VERBOSE = False


def mmseq_cluster_tsv_to_df(path : str) -> pd.DataFrame:
    '''
    the **_cluster.tsv file is an output of the mmseqs2 easy-cluster command.
    the data is made up of two columns, tab delimited:
    1. id of the representative of the cluster
    2. id of a member of a cluster
    the function reads the input and processes it to a pandas dataframe.
    :param path: absolute path of .tsv file to be parsed
    :return: a dataframe. each entry has 2 columns:
    1. Representative column (left column of tsv file)
    2. member column (right column of tsv file)
    '''

    df = pd.read_csv(path,
                     sep='\t',
                     header=None,
                     names=['rep_id', 'member_id'],
                     )

    return df

def get_cluster_sizes(all_clusts_df: pd.DataFrame) -> pd.DataFrame:
    '''
    receive a datafrane of all sequences and the representatives of their cluster (symboling for which cluster they belong to)
    infer the cluster sizes by grouping the dataframe by the representatives id. return the inferred sizes as part of a new dataframe.
    :param all_clusts_df: a dataframe containing all the clusters and their members. cluster representatives id column must
    be names 'rep_id'.
    :return: a new dataframe, each entry has a representative id ('rep_id') and the rep's cluster size 'cluster size'.
    '''
    if(VERBOSE):
        FUNC_START()
    grp = all_clusts_df.groupby(['rep_id'])
    df_filtered = grp.count()
    temp = df_filtered.reset_index()
    temp.columns = ['rep_id', 'cluster_size']
    if (VERBOSE):
        FUNC_END()
    return temp

def filter_big_clusters(all_clusts_df: pd.DataFrame, n: int) -> pd.DataFrame:
    '''
    receive a clusters info dataframe and filter out clusters of size smaller the 'n'.
    :param cluster_df: Parsed dataframe representing the tsv file.
    :param n: Minimal cluster size. Clusters of size smaller then n will be filtered out.
    :return: a filtered dataframe, left column for cluster representative name, right column - cluster size.
    '''
    grp = all_clusts_df.groupby(['rep_id'])
    temp = grp.filter(lambda x: len(x) >= n)
    return temp

def fasta_to_df(fasta_input_path) -> pd.DataFrame:
    '''
    Use the SeqIO module to extract the data from the fasta formatted file and create a pandas dataframe out of it.
    :param fasta_input_path: Absolute path of the fasta file.
    :return: a dataframe with the data from the fasta file, columns are: 'id', 'description', 'seq_length', 'sequence'
    '''

    # For debugging purposes
    if (VERBOSE):
        FUNC_START()

    # Lists for holding the data
    ids = []
    descriptions = []
    seqs = []
    lens = []

    data = SeqIO.parse(open(fasta_input_path), 'fasta')

    # Iterate through all genes in the fasta file, extract their data and add it to the relevant list.
    for gene in data:
        desc = (' ').join(gene.description.split('_')[-1].split(' ')[1:])
        id = gene.id
        seq = gene.seq
        
        descriptions.append(desc)
        ids.append(id)
        seqs.append(str(seq))
        lens.append(len(seq))

    # Create pandas series from the lists
    ids = pd.Series(ids)
    seqs = pd.Series(seqs)
    descriptions = pd.Series(descriptions)
    lens = pd.Series(lens)

    # Create a dataframe from the series' and rename the columns, then return it.
    cols = ['id', 'description', 'seq_length', 'sequence', ]
    reps_df = pd.concat([ids, descriptions, lens, seqs], axis = 1)
    reps_df.columns = cols

    if (VERBOSE):
        FUNC_END()

    return reps_df


def add_cluster_members_info(all_clusters: pd.DataFrame, unknown_genes_faa_path:str)-> pd.DataFrame:
    '''
    Enrich the all_clusters dataframe which only has the 'member_id' and the 'rep_id' columns with all the data from the
    unknown_genes fasta file. this includes - description, sequence, sequence length.
    :param all_clusters: a pandas dataframe containing all the genes. must have a 'member_id' column to merge by.
    :param unknown_genes_faa_path: Absolute path to unknown genes fasta file. if originates from this pipeline, it is
     created by the pull_and_cluster module.
    :return: A dataframe off all genes and their cluster representatives, each entry will have the member's original annotation
    and it's sequence.
    '''
    if (VERBOSE):
        FUNC_START()

    # All the genes we deal were originally found as unknown genes and therefore they still exist on the original
    # unknown genes fasta file that was previously created and should be located on the working directory
    # Parse this into a pd dataframe so we can take the additional features that are provided with each sequence such as length etc. and use it.

    all_unknown_seqs = fasta_to_df(unknown_genes_faa_path)
    # all_unknown_seqs = fasta_to_df(f'{out_path}/unknown_genes_2020_02_28.faa') #For local use
    all_unknown_seqs = all_unknown_seqs.rename(columns = {'id':'member_id'})

    # merge the parsed dataframe to the 'all_clusters' dataframe, resulting that for each cluster member we can now know it's description, length etc.
    all_clusters = all_clusters.merge(all_unknown_seqs, on='member_id')
    all_clusters = all_clusters.sort_values('rep_id')

    if (VERBOSE):
        FUNC_END()
    return all_clusters

def add_annotation_frequencies_to_clusters(clusts_by_reps: pd.DataFrame,
                                           all_clusters: pd.DataFrame)-> pd.DataFrame:
    '''
    Add to the cluster representations dataframe a column describing the distribution of the annotation of the cluster members.
     for example, a cluster could be made of 100% (1.0) members that had the annotation of a 'hypothetical protein'.
    the new column format is 'annotation_1 frequency_1, annotation_2 frequency_2, ...'.
    :param clusts_by_reps: A dataframe with an entry for each cluster represented by it's representing member.
    :param all_clusters:  pandas dataframe containing all the genes. each entry has a rep_id column that indicates which
    cluster it belongs to.
    :return:
    '''
    if (VERBOSE):
        FUNC_START()
    annotation_freqs_vec = []

    # For representative in the rep's dataframe use it's id to mask out all the members of it's cluster, then add the
    # ... frequency of the diffrerent annotations in the cluster to the reps dataframe.
    for cluster_index, cluster_data in clusts_by_reps.iterrows():
        mask = (all_clusters.rep_id == cluster_data.rep_id)
        tmp = all_clusters[mask]

        annotation_freqs = [f'{anno} ({round(freq, 3)})' for anno, freq in
                            dict(tmp.description.value_counts(normalize=True)).items()]
        annotation_freqs_vec.append(', '.join(annotation_freqs))

    # Assign the annotation frequencies vector to the reps dataframe.
    clusts_by_reps = clusts_by_reps.assign(annotations_freqs=annotation_freqs_vec)

    if (VERBOSE):
        FUNC_END()
    return clusts_by_reps

def add_cassete_source(all_clusters: pd.DataFrame, unknown_genes_df: pd.DataFrame) -> pd.DataFrame:
    '''
    Add new data to the all_clusters df. for each gene, add to it's entry a new column stating what cassette it comes
    from. The cassette will be represented by the ID of the cassette representing gene as it appears in
    the final_df_{date}_n5 file that is produced by the proccess_hmm_results module.
    :param unknown_genes_df: A dataframe containing all unknown genes with their original fasta name.
    :param all_clusters: A dataframe off all genes and their cluster representatives.
    :return: A new dataframe, a copy of all_clusters but with the new column - 'cassette_representing_gene'.
    '''
    if (VERBOSE):
        FUNC_START()
    #Rename the unknown dataframe columns to suit the merge command.
    unknown_genes_df.rename(columns={'unknown_gene': 'member_id'}, inplace=True)

    # Add this data to the members table
    all_clusters = pd.merge(all_clusters, unknown_genes_df, on='member_id', how='inner')
    if (VERBOSE):
        FUNC_END()
    return all_clusters


def add_ResFam_members_from_cassete(all_clusters: pd.DataFrame, csts_df: pd.DataFrame) -> pd.DataFrame:
    '''
    Use the cassttes dataframe (csts_df) to add to the all clusters datafrane (all_clusters) a new column.
    For each gene in all clusters, add the known and identified members of it's cassette, and the architecture of the
    cassette.
    :param csts_df: A dataframe contaitning all the cassettes, this is the same dataframe as 'final_df_{date}_n5'
     that is produced by the proccess_hmm_results module.
    :param all_clusters:  A dataframe off all genes and their cluster representatives.
    :return: A new dataframe, a copy of all_clusters but with the new columns -  'members', 'architecture'
    '''
    if (VERBOSE):
        FUNC_START()

    # For all unknown genes find out which ResFam members are located nearby
    # Extract original cassete for each unknonwn gene so it's possible to trace back what comprehends that cassete

    # Subset the cassetes df then merge it the members table.
    # Resulting a new columm for each member, containng the ResFam gene's it located near to.
    # Read final cassetes dataframe. iterate throughout it's contigs and look for their members in the
    # csts_df = pd.read_csv('/Users/dror/PycharmProjects/AMR_pipeline/data/pHMM_out/final_df_2020_02_28_n5.csv', index_col = 0) #For local use


    cols = ['full_id', 'members', 'architecture']
    tmp = csts_df[cols].rename(columns={'full_id': 'cassette_representing_gene'})
    all_clusters = all_clusters.merge(tmp, on='cassette_representing_gene', how='inner')

    if (VERBOSE):
        FUNC_END()
    return all_clusters


def add_ResFam_members_counts(clusts_by_reps, all_clusters):
    '''
    For each cluster in the clusts_by_reps dataframe, add the distribution of the known ARGs in the different cassettes
    the the cluster members come from.
    :param clusts_by_reps: A dataframe of all the representatives of the clusters. this dataframe represents the data
    about the different clusters.
    :param all_clusters: A dataframe off all genes and their cluster representatives.
    :return: A new dataframe, a copy of the clusts_by_reps dataframe, with a new column representing the distribution of
    the identified members in the cassettes the cluster members originate from.
    '''
    if (VERBOSE):
        FUNC_START()

    counts_vector = []
    for cluster_index, cluster_data in clusts_by_reps.iterrows():
        mask = (all_clusters.rep_id == cluster_data.rep_id)
        tmp = all_clusters[mask]
        cluster_size = cluster_data.cluster_size

        counts_dict = tmp.members.str.get_dummies(',').apply(lambda col: col.sum(), axis=0).to_dict()
        counts = [f'{resFam_member}; {count}; {round(count / cluster_size, 3)}' for resFam_member, count in
                  sorted(counts_dict.items(), reverse=True, key=lambda item: item[1])]


        counts_vector.append(', '.join(counts))

    clusts_by_reps = clusts_by_reps.assign(ResFam_neighbours=counts_vector)

    if (VERBOSE):
        FUNC_END()

    return clusts_by_reps


def mmseq_clusters_info_to_df(unknown_genes_df: pd.DataFrame, results_prefix: str = 'clusRes_final',
                              minimal_cluster_size: int = 100,
                              keep_rep_seqs_flag: bool = False,
                              VERBOSE_FLAG: bool = False) -> pd.DataFrame:
    '''
    Recive cluster file path and suffix, and parse results from the mmseqs module into a pandas dataframe for further
    processing.
    A minimum cluster size can be determined by the user.
    output can be also written into file by using the write_to_csv flag.

    :param unknown_genes_df: A dataframe containing the unknonwn genes that were previously extacted, and the cassette
    that they were extracted from.
    :param results_path: a path to the tsv files
    :param results_prefix: default 'clusRes_final' results prefix as decided by the user
    :param minimal_cluster_size: default 30. clusters below this size will not be included in the report.
    :param write_to_csv: default False. if True, a file describing the clusters will be written to the work directory.
    :param keep_rep_seqs_flag: default False. if True, for each cluster, it's representative sequence will be also included
    in the data.
    :return: a pd.DataFrame containing one entry per each cluster who has more than {minimal_cluster_size} members.
    additional data - cluster id, it's representative description, it's size and it's representative length.
    '''
    global VERBOSE
    VERBOSE = VERBOSE_FLAG

    #Read mmseqs output files
    print("In parse_mmseq_results module. Parsing mmseq output files.")
    #The table is a result of the parse-mmseq command. This is a description of the clusters using their representatives.
    cluster_file_path = out_path + f'/{results_prefix}_cluster.tsv'
    rep_seqs_path = out_path + f'/{results_prefix}_rep_seq.fasta'

    #Read all clusters output of mmseqs
    all_clusters = mmseq_cluster_tsv_to_df(cluster_file_path)

    #Add to the representors df of each cluster some data: 'rep_id', 'description', 'rep_seq_length', 'rep_seq'
    clusts_by_reps = fasta_to_df(rep_seqs_path)
    clusts_by_reps = clusts_by_reps.rename(columns = {'id':'rep_id', 'seq_length':'rep_seq_length', 'sequence':'rep_seq'})
    
    print("Built cluster-reps dataframe, and member-rep dataframe")
    print(f"Clusters data: {clusts_by_reps.shape[0]} different clusters, made of {all_clusters.shape[0]} different members.")

    #Filter small clusters, they are not intreseting
    filtered = filter_big_clusters(all_clusters, minimal_cluster_size)
    sized_clusters = get_cluster_sizes(filtered)
    clusts_by_reps = clusts_by_reps.merge(sized_clusters, how='inner', on='rep_id')
    clusts_by_reps = clusts_by_reps.sort_values(by='cluster_size', ascending=False).reset_index(drop=True)

    print(f"Clusters data: clusters under the size of {minimal_cluster_size} were filtered out."
          f" biggest cluster size is: {clusts_by_reps.cluster_size[0]}.")
    print(f"Updated number of clusters is: {clusts_by_reps.shape[0]}.")

    new_cols_order = ['rep_id', 'description', 'cluster_size', 'rep_seq_length', 'rep_seq']
    clusts_by_reps = clusts_by_reps[new_cols_order]

    if(keep_rep_seqs_flag == False):
        print("keep_rep_seq_flag is set to False, dropping cluster representatives seqs from final output (Seq source is the pullseq output).")
        clusts_by_reps.drop('rep_seq', axis=1, inplace= True)

    # if (write_to_csv):
    #     mmseq_results_to_csv(clusts_by_reps)

    ### Add cluster members info based on the cassete the cluster members originate from
    unknown_genes_faa_path = f'{out_path}/unknown_genes_{date}.faa'
    all_clusters = add_cluster_members_info(all_clusters, unknown_genes_faa_path)


    #Add annotation freuencies of the different cluster members for each cluster
    print("Adding annotation frequencies to clusters.")
    clusts_by_reps = add_annotation_frequencies_to_clusters(clusts_by_reps, all_clusters)

    #Add the original cassete and it's members to each cluster member (of all clusters)
    print("Adding identified ResFam neighbour counts to clusters.")
    all_clusters = add_cassete_source(all_clusters, unknown_genes_df)

    #For each cluster members (of all clusters), add data about which ResFam members are in it's original cassete (from final_df previously produced)
    csts_df = pd.read_csv(f'{out_path}/final_df_{date}_n5.csv', index_col=0)
    all_clusters = add_ResFam_members_from_cassete(all_clusters, csts_df)
    clusts_by_reps = add_ResFam_members_counts(clusts_by_reps, all_clusters)

    new_order = ['rep_id', 'description', 'cluster_size', 'annotations_freqs', 'ResFam_neighbours', 'rep_seq_length']
    new = 'cluster_members_annotations_freq', 'cluster_members_resfam_neighbours',
    clusts_by_reps = clusts_by_reps[new_order].rename(columns={'annotations_freqs': 'members_annotations_freq',
                                                               'ResFam_neighbours': 'members_ResFam_neighbours'}
                                                      )

    file_name = f'{out_path}/final_clusters_profiling_{date}.csv'
    print(f"Writing cluster profiling results to {file_name}.")
    clusts_by_reps.to_csv(file_name)
    file_name = f'{out_path}/all_clusters_{date}.csv'
    print(f"Writing all_clusters dataframe to {file_name}.")
    all_clusters.to_csv(file_name)

    print("Returning a tuple of: (clusters_by_reps_df, all_clusters_df)")
    return (clusts_by_reps, all_clusters)




