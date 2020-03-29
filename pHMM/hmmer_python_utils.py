import subprocess
import os
from datetime import datetime

out_path = os.getcwd()
date = datetime.day

def run_hmmer_search(seqdb_path, hmm_file_path, options, output_filename) -> int:
    '''
    Use python to run the equal of : hmmsearch [options] <hmmfile> <seqdb>
    Assume hmmer module is already loaded
    :param seqdb_path: path to sequences db
    :param hmm_file_path: hmm profile file path
    :param options: list of the format ['-opt', '{value}']
    :return: 0 upon success. See hmmer documentation for output file description  - http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf
    '''

    print("Running hmmersearch .\n"
          f"input Seqs db path is: {seqdb_path}\n"
          f"HMM profiles path is:: {hmm_file_path}\n"
          f"with options: {str(options)}")

    output_filename = out_path + f'/{output_filename}_{date}.tblout'

    cmd_line = ['hmmearch', options, hmm_file_path, seqdb_path]
    subprocess.run(cmd_line, stdout=open(output_filename, 'w'), stderr=subprocess.STDOUT)
    print("Finshed rynning pullseq, unknown genes file:\n"
          f"{output_filename}")

    output_filename = out_path + f'/unknown_genes_ids_{date}.faa'
    return 0


print(out_path)