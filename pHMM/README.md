# Usage instructions for the pHMM python package

### Modules installed and loaded
1. Python 3.7 or above  (e.g $ module load python/python-anaconda3.7-itaym )   
2. mmseqs2 (e.g $ module load mmseq2_new)  
3. pullseq (e.g $ module load pullseq)   
4. hmmer  

#### Non trivial Python packages in use:
1. pandas  
2. numpy  
3. click  
4. shutil  
5. subproccess  
6. datetime.datetime  
7. Bio.Seqio (Bio python)  

#### Abstract
This package is meant to ease the process of analyzing contigs from the stage of runnning an HMM search against a 
known HMM profiles DB and till it's final output - a list of clustered, unknown genes that were generated according to
the context of the contigs.
A full pipeline is consisted the following stages:
1. Preform an hmmsearch of sequences of the arranged as contigs. This stage is preformed outside of the scope of the 
provided code modules.  
2. Parse hmmsearch outputs into a pandas dataframe.  
3. Process hmmsearch outputs, Build AMR cassettes.  
4. Fetch unknown genes, generate a fasta file containing all these genes.  
4. Cluster unknown genes in two stages both using mmseqs package:  
    1. Cluster together similiar sequences by using high identity and low coverage filter.  
    2. Cluster the representatives of each cluster by using low identity and higher coverage.  
5. Further investigate the clustering results and add more data about each cluster.  

### File composition
The package has several modules, each one has a different role in the process. see further details in the code documentation.  
1. **tblout_to_pandas.py** - receive a path of a folder where the hmmsearch results are, concatenate them and parse 
the results into a pandas dataframe while conserving the original fields of the hmmsearch output.  
2. **process_hmm_results.py** - process the parsed tblout results dataframe and construct from them a new dataframe 
of ARG cassettes and potential unknown genes. This module has two outputs:  
    * _df_filtered_{date}.csv_ - A final dataframe off all the hmm results after removing contigs with 1 hit and after 
    choosing only the best hit for each gene.  
    * _final_df_{date}_n{x}.csv_ - A dataframe containing all the cassettes that were found, each entry is a cassette and 
    includes all the details about the cassette members, and it's architecture (visualization, scores and e-values).  
3. **extract_unknown_genes.py** - This module is all about creating tables for the pullseq module to run later and 
for creating a dataframe (and write it out as csv) of uknown gens and the cassette they originated from. it has two
 outputs:  
    * _unknown_genes_for_analysis_{date}.csv_ - This csv has two columns: the unkonwn gene, the cassette it came from.
    It is printed in the case one want to preform the anlysis by self (running the parse mmseq module separatley and not as part of the pipeline).
    In general, when running a full pipeline, this file can be ignored as a dataframe with exactly the same contents is
    returned when the module finishes it's work.  
    * _unknown_genes_for_pullseq_{date}.csv_ - A csv of one column containing the gene names of all the unknown genes
    that were located inside the ARG cassettes. This csv must be printed in this format for the pullseq module to run.  
4. **pull_and_cluster.py** - This module runs the pullseq and mmsqes modules from outside of the python shell, the outputs
files are in accordance with the packages doumentation. Three of the outputs are being used later:  
    * unknown_genes_{date}.faa - This file is the outcome of the pullseq module. It contains all the unknown genes
     that were extracted from the ARG cassettes in a fasta format. This file is later used for both clustering and 
     adding data to the clusters (parse_mmseq_results module).  
    * _clusRes_final_rep_seq.faa_ - used to infer the data about the representatives of the clusters.  
    * _clusRes_final_cluster.tsv_ - used to infer the clusters.  
5. **parse_mmseq_results.py** - This final module is used to analyze the clustering results. It creates two main dataframes
inside it, one is the all_clusters that has the data about all the genes and the clusters the belong to, the 
second dataframe is clusts_by_reps and is used to save data about the clusters themselves.
This module has two outputs, one for each dataframe. These outputs are:     
    * _all_clusters_{date}.csv_ - containing all the sequences of the unknown genes and which cluster they belong to, 
    and also which cassette they came from.   
    * _final_clusters_profiling_{date}.csv_ - Containing all the data about the clusters, and the distribution 
    of their members.   
These two outputs are considered the final outputs of the pipeline.  
6. **general_utils.py** - This file only contains some lambda functions for debugging purposes (when the verbose flag is 
True, they are used for printouts.)  


### Usage example with a qsub job
1. locate the pHMM_jobs.sh in a directory youre comfortable with.  
2. edit the pHMM_jobs with updated file locations  
3. $qsub pHMM_jobs.sh  
4. cat {user_name}{job_id}.o - to see the runtime log of the modules.  

Note that it is not necessary to load all the modules stated above when running this job as it does so itself.  

### Usage example - easy-pipeline -> RECOMMENDED
pHMM easy-pipeline -m 100 -s tbl -v False /davidb/drorhadas/Jan2020/all_prots_vs_resFamHighSpec/ /davidb/drorhadas/data/all_proteins_concatenated/all_proteins_concatenated.faa  


### Usage example (e.g on the Phoenix server)
* $ module load python/python-anaconda3.7-itaym
* $ module load pullseq
* $ module load mmseq2_new
* $ ln -s {package_path} pHMM 
* $ python pHMM --help
* $ python phMM pprocess_tblout --help
* $ python pHMM process_tblout /davidb/drorhadas/temp/resFam_HighSpec_vs_Chinese_Pigfarm_hmmsearch_results --max_gene_distance 5 -s tbl
* $ python pHMM pullseq /davidb/drorhadas/data/Chinese_pigs/pigs_DB/China_pigfarms_concatenated.faa unknown_genes_2019_09_24.csv
* $ head -n 20 unknown_genes_2019_09_24.faa
* $ python pHMM cluster ./unknown_genes_ids_2019_09_24.faa
* $ python pHMM parse-mmseq ./



