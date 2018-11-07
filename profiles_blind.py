import subprocess

from utils.utils import *

subprocess.check_call([SCRIPT_MAKEBLASTDB, PROFILES_DB])
subprocess.check_call([SCRIPT_PSIBLAST, PROFILES_DB, BLIND_FASTA_DIR,
                       PROFILES_BLIND_OUTPUT_PSSM_DIR, PROFILES_BLIND_OUTPUT_BLAST_DIR])

for filename in os.listdir(PROFILES_BLIND_OUTPUT_PSSM_DIR):
    df = parse_pssm_file(PROFILES_BLIND_OUTPUT_PSSM_DIR + filename)
    df.to_csv(PROFILES_BLIND_OUTPUT_CSV_DIR + os.path.splitext(filename)[0] + '.csv')
