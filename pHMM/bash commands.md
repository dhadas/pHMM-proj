#BASH tips

##HMMER
### General HMMER hmmsearch workflow:
$ module load hmmr/hmmr-3.1b2  
$ hmmsearch [options] <hmmfile> <seqdb> 
or 
$ hmmsearch -h

usage example:
hmmsearch --tblout desired_output_name.tbl -E 0.001 hmm_db_file_path protein_target_file_path


### Running hmmer on all files in a certrain directory:
$for file in 
    ../all_prots_folder/*.faa ; do
     hmmsearch --tblout "./resFamHighSpec_vs_$(basename $file).tbl" -E 0.001 "../resFam_high_spec_classfication" "$file" 
     ; done
     

more about basename : https://stackoverflow.com/questions/3362920/get-just-the-filename-from-a-path-in-a-bash-script  
more about for :  https://stackoverflow.com/questions/10523415/execute-command-on-all-files-in-a-directory

 ### Counting the lines in a file
 #wc -l <file name>
  

     