import os
import urllib.request
import pandas as pd

from constants import *
from Bio import SeqIO
from Bio.PDB import DSSP
from Bio.PDB import PDBParser


def make_dssp(pdb_id, chain, chain_delimiter=':'):
    """
    Retrieve dssp string from PDB database

    :param pdb_id: pdb id
    :param chain: pdb chain id
    :param chain_delimiter: delimiter between id and chain id
    :return: id, aa sequence, ss sequence
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

    os.remove(pdb_id)
    return pdb_id + chain_delimiter + chain, aa, ss


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


def divide_fasta_file(input_file, output_dir, ext='', is_fasta=True):
    """
    Divide fasta file by single files

    :param ext: extension of new files
    :param is_fasta: is it fasta file or just fasta-like file (like dssp)
    :param input_file: input fasta file
    :param output_dir: output directory
    """
    if not is_fasta:
        divide_fasta_like_file(input_file, output_dir, ext)
        return

    with open(input_file, 'r'):
        for seq_record in SeqIO.parse(input_file, "fasta"):
            with open(output_dir + str(seq_record.id).replace(':', '_') + '.' + ext, "w") as out_file:
                out_file.write('>' + str(seq_record.id) + '\n' + str(seq_record.seq) + '\n')


def divide_fasta_like_file(input_file, output_dir, ext=''):
    """
    Divide file like fasta file by single files

    :param ext: extension of new files
    :param input_file: input fasta file
    :param output_dir: output directory
    """
    with open(input_file, 'r') as file:
        body = ''
        p_id = ''
        for line in file:
            if line[0] == '>':
                if len(p_id) > 0:
                    with open(output_dir + p_id.replace(':', '_') + '.' + ext, "w") as out_file:
                        out_file.write('>' + p_id.replace(':', '_') + '\n' + body + '\n')
                    body = ''
                p_id = line.strip()[1:]
            else:
                body += line.strip()
        with open(output_dir + p_id.replace(':', '_') + '.' + ext, "w") as out_file:
            out_file.write('>' + p_id.replace(':', '_') + '\n' + body + '\n')


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
    """
    remove indecies from fasta after blastclust
    :param in_fasta: input fasta file
    """
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


def parse_pssm_file(pssm_file):
    """
    Parsing pssm file

    :param pssm_file: input pssm file
    :return: pandas dataframe
    """
    colnames = ['res'] + AA_ORDER
    usecols = [1] + list(range(22, 42))
    df = pd.read_table(pssm_file, sep=r"\s*", header=None, skiprows=3, skipfooter=5, usecols=usecols)
    df.columns = colnames
    df.set_index('res', inplace=True)
    df = df / 100
    return df


def restore_profile_from_csv(csv_file):
    """
    Restoring profile from csv file

    :param csv_file: csv file
    :return: numpy 2d array
    """
    with open(csv_file) as profile_csv:
        profile = pd.read_csv(profile_csv).values[:, 1:].astype(float)
        return profile


def get_relative_file(in_file, directory, ext):
    """
    get path to file in other directory with other extention
    :param in_file: input path
    :param directory: other directory
    :param ext: new extention
    :return:
    """
    filename_w_ext = os.path.basename(in_file)
    filename, file_extension = os.path.splitext(filename_w_ext)
    return os.path.join(directory, filename + '.' + ext)


def get_dssp_from_file(in_file):
    """
    Get dssp string from file with single inline dssp entity

    :param in_file: input dssp file
    :return: dssp string
    """
    with open(in_file) as file:
        return file.readlines()[1].strip()


def is_zero_profile(in_file):
    profile = restore_profile_from_csv(in_file)
    for i in range(0, profile.shape[0]):
        for j in range(0, profile.shape[1]):
            if profile[i,j] != 0:
                return False
    return True


def remove_zero_profiles(profile_files):
    for profile_file in profile_files:
        if is_zero_profile(profile_file):
            os.remove(profile_file)
