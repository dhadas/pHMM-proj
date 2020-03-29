import subprocess
from datetime import datetime
import os
import shutil
date = datetime.today().strftime('%Y_%m_%d')
out_path = os.getcwd()
from general_utils import FUNC_START, FUNC_END
VERBOSE = False


def run_pullseq(db_path, names_path, VERBOSE_FLAG) -> str:
    '''
    run the package pullseq from within the pHMM pipeline
    :param db_path: absolute path to the directory where the sequences database is (the database to pull from)
    :param names_path: absolute path to the directory where the the names file (the file containing the genes that
    are to be pulled) is located.
    :return: a path (string) to the fasta file of the all the pulled genes (sequences and data). file name format is:
    unknown_genes_{date}.faa.
    '''
    global VERBOSE
    VERBOSE = VERBOSE_FLAG

    if (VERBOSE):
        FUNC_START()
    print("Running pullseq.\n"
          f"Pulling sequences from names path: {names_path}\n"
          f"Against db path: {db_path}\n"
          f"Output path is: {out_path}")

    output_filename = out_path+f'/unknown_genes_{date}.faa'

    cmd_line = ['pullseq', f'-i', f'{db_path}', f'-n', f'{names_path}']
    subprocess.run(cmd_line, stdout=open(output_filename, 'w'), stderr=subprocess.STDOUT)
    print("Finshed running pullseq, unknown genes file:\n"
          f"{output_filename}")

    if (VERBOSE):
        FUNC_END()
    return output_filename

def run_mmseqs_easy_clust(coverage, min_seq_id, input_path, results_prefix, temporary_output_dump: str = 'tmp'):
    '''
    Python wrapper function for running cmd package mmseqs2:
    cmd is : $mmseqs easy-cluster -c {coverage} --min-seq-id {min_seq_id} {input_path} {results_prefix} {temporary_output_dump}
    Assumes the the module mmseqs2 is already loaded to the terminal session.
    :param coverage: Minimal coverage of the substring in each sequence. for the -c flag.
    :param min_seq_id: Minimal identity between sequences. for the --min-seq-id option.
    :param input_path: Path to sequence database to be clustered.
    :param results_prefix: Prefix of all output files outputted by operation.
    :param temporary_output_dump: Temporary location to create temporary runtime files. this is later deleted by the
    function after completion of of clustering operation.
    :return: 0 upon success. Output files can be located by collecting files starting with {results_prefix}
    and further processed using the parse_mmseq_results module.
    '''
    if (VERBOSE):
        FUNC_START()
    cmd_line = ['mmseqs',
                'easy-cluster',
                f"-c",
                f"{coverage}",
                f"--min-seq-id",
                f"{min_seq_id}",
                input_path,
                results_prefix,
                temporary_output_dump]

    subprocess.run(cmd_line, stderr=subprocess.STDOUT)
    shutil.rmtree('tmp')
    if (VERBOSE):
        FUNC_END()
    return 0

def double_cluster(input_path, cov_rate1, cov_rate2, min_id1, min_id2, pref1, pref2, VERBOSE_FLAG):
    global VERBOSE
    VERBOSE = VERBOSE_FLAG

    if (VERBOSE):
        FUNC_START()
    run_mmseqs_easy_clust(cov_rate1, min_id1, input_path, pref1, 'tmp')
    reps_filepath = out_path + f'/{pref1}_rep_seq.fasta'
    run_mmseqs_easy_clust(cov_rate2, min_id2, reps_filepath, pref2, 'tmp')
    if (VERBOSE):
        FUNC_END()
    return 0