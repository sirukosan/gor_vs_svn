import os
from constants import *
from Bio import SeqIO
from Bio.PDB import DSSP
from Bio.PDB import PDBParser
import urllib.request


def make_dssp(pdb_id, chain, with_aa=True, chain_delimiter=':'):
    """
    Retrieve dssp string from PDB database

    :param pdb_id: pdb id
    :param chain: pdb chain id
    :param with_aa: add amino acids to result
    :param chain_delimiter: delimiter between id and chain id
    :return: fasta-like string with dssp string (and aa string in case with_aa = True)
    """
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

    body = '\n' + ss
    if with_aa:
        body = '\n' + aa + body

    result = ">" + pdb_id + chain_delimiter + chain + body
    os.remove(pdb_id)
    return result


def remove_from_fasta(fasta_file_name, ids, leave=False):
    """
    Removing ids from fasta file

    :param fasta_file_name: path to fasta file
    :param ids: ids to remove/leave
    :param leave: leave this ids if True
    """
    tmp_fasta_file_name = fasta_file_name + TMP
    with open(tmp_fasta_file_name, 'w') as tmp_fasta_file:
        for seq_record in SeqIO.parse(fasta_file_name, "fasta"):
            if (seq_record.id not in ids) == (not leave):
                tmp_fasta_file.write('>' + str(seq_record.id) + '\n' + str(seq_record.seq) + '\n')

    os.remove(fasta_file_name)
    os.rename(tmp_fasta_file_name, fasta_file_name)


def combine_fasta_files(fastas_paths, out_file):
    """
    Combine fasta files into one

    :param fastas_paths: list of fasta files
    :param out_file: fasta file to write in
    """
    with open(out_file, 'w') as out:
        for filename in fastas_paths:
            for seq_record in SeqIO.parse(filename, "fasta"):
                out.write('>' + str(seq_record.id) + '\n' + str(seq_record.seq) + '\n')


def prepare_fasta_for_blastclust(in_fasta, out_fasta):
    """
    Remove from fasta sequences which fail blastclust and adding indices
    prepare/unprepare needed because of bug (?) in blastclust
    :param in_fasta: in fasta file
    :param out_fasta: out asta file
    """
    with open(out_fasta, 'w') as out:
        i = 0
        for seq_record in SeqIO.parse(in_fasta, "fasta"):
            if len(seq_record.seq) > 5 and 'XXXXX' not in seq_record.seq and 'UUUUU' not in seq_record.seq:
                out.write(
                    '>' + seq_record.id.split('|')[0] + '_' + str(i) + '\n' + str(seq_record.seq) + '\n')
                i += 1


def unprepare_fasta_after_blastclust(in_fasta):
    tmp_file_name = in_fasta + TMP

    with open(tmp_file_name, 'w') as tmp_file:
        for seq_record in SeqIO.parse(in_fasta, "fasta"):
            tmp_file.write('>' + seq_record.id.split('_')[0] + '\n' + str(seq_record.seq) + '\n')

    os.remove(in_fasta)
    os.rename(tmp_file_name, in_fasta)


def get_blastclustout_first_ids(blastclust_out_file):
    """
    Get first id in each cluster from blastclust output file

    :param blastclust_out_file:  blastclust output file path
    :return: list of ids
    """
    remained_ids = []
    with open(blastclust_out_file, 'r') as blastclust_out:
        for line in blastclust_out:
            remained_ids.append(line.split(' ')[0].split('_')[0])
    return remained_ids


def get_blastp_ids(blastp_out_file, identity_threshold=1000):
    """
    Get from blastp output ids with identity threshold

    :param blastp_out_file: blastp output file
    :param identity_threshold: identity threshold
    :return: list of ids
    """
    id_map = {}
    remained_ids = []
    with open(blastp_out_file, 'r') as blastp_out:
        for line in blastp_out:
            line_arr = line.split('\t')
            if line_arr[0] in id_map:
                if float(line_arr[2]) > id_map[line_arr[0]]:
                    id_map[line_arr[0]] = float(line_arr[2])
            else:
                id_map[line_arr[0]] = float(line_arr[2])

    for key, value in id_map.items():
        if value < identity_threshold:
            remained_ids.append(key)

    return remained_ids
