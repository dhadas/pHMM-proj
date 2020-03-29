#!/bin/bash

#submit to dudus queue
#PBS -q dudu24h

# use this to name your job
#PBS -N dhadas

# number of CPUs o use, here 1 (you can also use -l mem=10gb to specify you need a specific ammount of RAM).
#PBS -l ncpus=4

# use environemnt variable from your login
#PBS -V

# Join standard output and error to a single file
#PBS -j oe

### here is where the execution starts
### in the case of hmmer this will also be the output folder for the .tbl files.
cd /davidb/drorhadas/Jan2020/pHMM_out

##### Load below the modules you need #####
module load python/python-anaconda3.7-itaym
module load hmmer/hmmer-3.2.1
module load mmseq2_new
module load pullseq

##### below enter your commands #####
python pHMM easy-pipeline -m 100 -s tbl -v False /davidb/drorhadas/Jan2020/all_prots_vs_resFamHighSpec/ /davidb/drorhadas/data/all_proteins_concatenated/all_proteins_concatenated.faa


#### run th command using qsub <file_name>
#### view the command status using qstat -u drorhadas
#### output file is dhadas.o(number of job) if you use print commands in python they will appear here
#### Use: $ cat dhadas.o(number of job) to see the output

