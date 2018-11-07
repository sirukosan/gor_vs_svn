import os
import subprocess

import pandas as pd

from constants import *

subprocess.check_call([SCRIPT_MAKEBLASTDB, PROFILES_DB])
subprocess.check_call([SCRIPT_PSIBLAST, PROFILES_DB, INPUT_FASTA_FOLDER,
                       PROFILES_OUTPUT_PSSM_DIR, PROFILES_OUTPUT_BLAST_DIR])


def parse_pssm(filename):
    colnames = ['res'] + AA_ORDER
    usecols = [1] + list(range(22, 42))
    df = pd.read_table(filename, sep=r"\s*", header=None, skiprows=3, skipfooter=5, usecols=usecols)
    df.columns = colnames
    df.set_index('res', inplace=True)
    df = df / 100
    return df


for filename in os.listdir(PROFILES_OUTPUT_PSSM_DIR):
    df = parse_pssm(PROFILES_OUTPUT_PSSM_DIR + filename)
    df.to_csv(PROFILES_OUTPUT_CSV_DIR + os.path.splitext(filename)[0] + '.csv')
