import subprocess

from utils import *

subprocess.check_call([SCRIPT_MAKEBLASTDB, PROFILES_DB])
subprocess.check_call([SCRIPT_PSIBLAST, PROFILES_DB, INPUT_FASTA_FOLDER,
                       PROFILES_OUTPUT_PSSM_DIR, PROFILES_OUTPUT_BLAST_DIR])

for filename in os.listdir(PROFILES_OUTPUT_PSSM_DIR):
    df = parse_pssm_file(PROFILES_OUTPUT_PSSM_DIR + filename)
    df.to_csv(PROFILES_OUTPUT_CSV_DIR + os.path.splitext(filename)[0] + '.csv')
