import os
import subprocess

import pandas as pd
from Bio import SeqIO

from constants import *


def divide_fasta(input_file, output_dir):
    with open(input_file, 'r') as in_file:
        for seq_record in SeqIO.parse(input_file, "fasta"):
            with open(output_dir + str(seq_record.id).replace(':', '_'), "w") as out_file:
                out_file.write('>' + str(seq_record.id) + '\n' + str(seq_record.seq) + '\n')


def parse_pssm(filename):
    colnames = ['res'] + AA_ORDER
    usecols = [1] + list(range(22, 42))
    df = pd.read_table(filename, sep=r"\s*", header=None, skiprows=3, skipfooter=5, usecols=usecols)
    df.columns = colnames
    df.set_index('res', inplace=True)
    df = df / 100
    return df


divide_fasta(PROFILES_BLIND_INPUT_FASTA_FILE, PROFILES_BLIND_INPUT_FASTA_DIR)

subprocess.check_call([SCRIPT_MAKEBLASTDB, PROFILES_BLIND_DB])
subprocess.check_call([SCRIPT_PSIBLAST, PROFILES_BLIND_DB, PROFILES_BLIND_INPUT_FASTA_DIR,
                       PROFILES_BLIND_OUTPUT_PSSM_DIR, PROFILES_BLIND_OUTPUT_BLAST_DIR])

for filename in os.listdir(PROFILES_BLIND_OUTPUT_PSSM_DIR):
    df = parse_pssm(PROFILES_BLIND_OUTPUT_PSSM_DIR + filename)
    df.to_csv(PROFILES_BLIND_OUTPUT_CSV_DIR + os.path.splitext(filename)[0] + '.csv')
