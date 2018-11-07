import os
import subprocess
import urllib.request

from Bio import SeqIO
from Bio.PDB import DSSP
from Bio.PDB import PDBParser

from constants import *


def remove_from_fasta(fasta_file_name, remained_ids):
    tmp_fasta_file_name = fasta_file_name + TMP
    with open(tmp_fasta_file_name, 'w') as tmp_fasta_file:
        for seq_record in SeqIO.parse(fasta_file_name, "fasta"):
            if seq_record.id in remained_ids:
                tmp_fasta_file.write('>' + str(seq_record.id) + '\n' + str(seq_record.seq) + '\n')

    os.remove(fasta_file_name)
    os.rename(tmp_fasta_file_name, fasta_file_name)


def make_dssp(pdb_id, chain):
    url = PDB_URL + pdb_id + '.pdb'
    urllib.request.urlretrieve(url, pdb_id)

    parser = PDBParser()
    structure = parser.get_structure(pdb_id, pdb_id)
    dssp = DSSP(structure[0], pdb_id, dssp='mkdssp')

    aa = ''
    ss = ''
    for key in dssp.keys():
        if key[0] == chain:
            aa += dssp[key][1]
            ss += SS_MAP[dssp[key][2]]
    result = ">" + pdb_id + ':' + chain + '\n' + aa + '\n' + ss
    os.remove(pdb_id)
    return result


# Make train fasta file with all sequences an blast DB
with open(TRAIN_FASTA_FILE_NAME, 'w') as train_fasta_file:
    for filename in os.listdir(INPUT_FASTA_FOLDER):
        for seq_record in SeqIO.parse(INPUT_FASTA_FOLDER + filename, "fasta"):
            train_fasta_file.write('>' + str(seq_record.id) + '\n' + str(seq_record.seq) + '\n')

subprocess.check_call([SCRIPT_MAKEBLASTDB, TRAIN_FASTA_FILE_NAME])

with open(BLIND_FASTA_FILE_NAME, 'w') as blind_fasta_file:
    i = 0
    for seq_record in SeqIO.parse(ROW_BLIND_FASTA_FILE_NAME, "fasta"):
        if len(seq_record.seq) > 5 and 'XXXXX' not in seq_record.seq and 'UUUUU' not in seq_record.seq:
            blind_fasta_file.write('>' + seq_record.id.split('|')[0] + '_' + str(i) + '\n' + str(seq_record.seq) + '\n')
            i += 1

subprocess.check_call([SCRIPT_BLASTCLUST, BLIND_FASTA_FILE_NAME, BLASTCLUST_OUT_FILE_NAME])


tmp_blind_fasta_file_name = BLIND_FASTA_FILE_NAME + TMP

with open(tmp_blind_fasta_file_name, 'w') as tmp_blind_fasta_file:
    for seq_record in SeqIO.parse(BLIND_FASTA_FILE_NAME, "fasta"):
        tmp_blind_fasta_file.write('>' + seq_record.id.split('_')[0] + '\n' + str(seq_record.seq) + '\n')

os.remove(BLIND_FASTA_FILE_NAME)
os.rename(tmp_blind_fasta_file_name, BLIND_FASTA_FILE_NAME)

remained_ids = []
with open(BLASTCLUST_OUT_FILE_NAME, 'r') as blastclust_out_file:
    for line in blastclust_out_file:
        remained_ids.append(line.split(' ')[0].split('_')[0])

remove_from_fasta(BLIND_FASTA_FILE_NAME, remained_ids)

subprocess.check_call([SCRIPT_BLASTP, BLIND_FASTA_FILE_NAME, TRAIN_FASTA_FILE_NAME, BLASTP_OUT_FILE_NAME])

id_map = {}
remained_ids = []
with open(BLASTP_OUT_FILE_NAME, 'r') as blastp_out_file:
    for line in blastp_out_file:
        line_arr = line.split('\t')
        if line_arr[0] in id_map:
            if float(line_arr[2]) > id_map[line_arr[0]]:
                id_map[line_arr[0]] = float(line_arr[2])
        else:
            id_map[line_arr[0]] = float(line_arr[2])

for key, value in id_map.items():
    if value < 30:
        remained_ids.append(key)
# remained_ids = random.sample(remained_ids, 100)

remove_from_fasta(BLIND_FASTA_FILE_NAME, remained_ids)

with open(BLIND_DSSP_FILE_NAME, 'w') as blind_dssp_file:
    for seq_record in SeqIO.parse(BLIND_FASTA_FILE_NAME, "fasta"):
        pdb_id = seq_record.id.split(':')[0]
        chain = seq_record.id.split(':')[1]
        blind_dssp_file.write(make_dssp(pdb_id,chain) + '\n')