import subprocess
import glob
from utils import *

combine_fasta_files(glob.glob(INPUT_FASTA_FOLDER + '*'), TRAIN_FASTA_FILE_NAME)
subprocess.check_call([SCRIPT_MAKEBLASTDB, TRAIN_FASTA_FILE_NAME])
prepare_fasta_for_blastclust(ROW_BLIND_FASTA_FILE_NAME, BLIND_FASTA_FILE_NAME)
subprocess.check_call([SCRIPT_BLASTCLUST, BLIND_FASTA_FILE_NAME, BLASTCLUST_OUT_FILE_NAME])
unprepare_fasta_after_blastclust(BLIND_FASTA_FILE_NAME)
remained_ids = get_blastclustout_first_ids(BLASTCLUST_OUT_FILE_NAME)
remove_from_fasta(BLIND_FASTA_FILE_NAME, remained_ids, True)
subprocess.check_call([SCRIPT_BLASTP, BLIND_FASTA_FILE_NAME, TRAIN_FASTA_FILE_NAME, BLASTP_OUT_FILE_NAME])
remained_ids = get_blastp_ids(BLASTP_OUT_FILE_NAME, 30)
remove_from_fasta(BLIND_FASTA_FILE_NAME, remained_ids, True)

with open(BLIND_DSSP_FILE_NAME, 'w') as blind_dssp_file:
    for seq_record in SeqIO.parse(BLIND_FASTA_FILE_NAME, "fasta"):
        pdb_id = seq_record.id.split(':')[0]
        chain = seq_record.id.split(':')[1]
        blind_dssp_file.write(make_dssp(pdb_id, chain) + '\n')
