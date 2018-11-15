import gzip
import pickle
from sys import stdout

from sklearn import svm

from utils.utils import *


def build_vectors_for_profile(profile_file, dssp_dir):
    """
    Builds vectors Y and X for profile
    :param profile_file: path to the profile csv file
    :param dssp_dir: fath to the dssp folder
    :return: Y and X vectors
    """
    profile = restore_profile_from_csv(profile_file)
    dssp = get_dssp_from_file(get_relative_file(profile_file, dssp_dir, 'dssp'))
    Y = transform_dssp_to_vector(dssp)
    X = np.empty((0, window * 20))
    for i in range(0, profile.shape[0]):
        row = np.empty((0, 0))
        for j in range(wind_l, wind_r + 1):
            if i + j < 0 or i + j >= profile.shape[0]:
                row = np.append(row, np.zeros((1, 20)))
            else:
                row = np.append(row, profile[i + j, :])
        X = np.vstack((X, row))
    return Y, X


def transform_dssp_to_vector(dssp):
    """
    Converts dssp string to number array
    :param dssp: dssp string
    :return: numbers array
    """
    res = []
    for s in dssp:
        if s.casefold() == 'h':
            res.append(0)
        elif s.casefold() == 'e':
            res.append(1)
        elif s.casefold() == '-' or s.casefold() == 'c':
            res.append(2)
    return res


def build_vectors(profile_files, out_x_file, out_y_file, dssp_dir, log_out=stdout, excep=[], exclude=True):
    """
    Builds vectors for number of profiles
    :param profile_files: array of profile files
    :param out_x_file: out file for x vectors
    :param out_y_file: out file for y vectors
    :param dssp_dir: dir with dssp sequences
    :param log_out: log file
    :param excep: ids of exception sequences
    :param exclude: True - exclude exception sequences, False - leave only exception sequences
    """
    i = 0
    with open(out_x_file, 'ab') as x_file, open(out_y_file, 'ab') as y_file:
        for profile_file in profile_files:
            if (os.path.splitext(os.path.basename(profile_file))[0] not in excep) == (exclude is True):
                log_out.write(str(i) + ' processing ' + profile_file + '\n')
                y, x = build_vectors_for_profile(profile_file, dssp_dir)
                np.savetxt(x_file, x, delimiter=",", fmt='%1.3f')
                np.savetxt(y_file, y, delimiter=",", fmt='%i')
                i += 1


def loop(i, test_dir, test_file, par_tup, x_train, y_train):
    """
    Train CSV model
    :param i: number of test
    :param test_dir: directory ov test
    :param test_file: file containes excluded ids for cross-validation
    :param par_tup: C and gamma parameters tuple
    :param x_train: X vector for training
    :param y_train: Y vector for training
    """
    with open(os.path.join(test_dir, 'log_' + str(i)), 'a') as curr_log_file:
        curr_log_file.write('building model for ' + test_file + '\nwith parameters: C = ' + str(
            par_tup[0]) + '\ngamma = ' + str(par_tup[1]) + '\n')
        model = svm.SVC(C=par_tup[0], kernel='rbf', gamma=par_tup[1], verbose=True)
        curr_log_file.write('SVC inited' + '\n')
        curr_log_file.write('fitting' + '\n')
        curr_log_file.flush()

        model.fit(x_train, y_train)
        curr_log_file.write('dumping' + '\n')
        pickle.dump(model, gzip.open(os.path.join(test_dir, 'model_' + str(i) + '.pkl.gz'), 'w'))
        curr_log_file.flush()


def mock_loop(i, test_dir, test_file, par_tup, x_train, y_train):
    """

    Imitation of model training (crates all files and folders)
    :param i: number of test
    :param test_dir: directory ov test
    :param test_file: file containes excluded ids for cross-validation
    :param par_tup: C and gamma parameters tuple
    :param x_train: X vector for training
    :param y_train: Y vector for training
    """
    with open(os.path.join(test_dir, 'log_' + str(i)), 'a') as curr_log_file:
        curr_log_file.write('building model for ' + test_file + '\nwith parameters: C = ' + str(
            par_tup[0]) + '\ngamma = ' + str(par_tup[1]) + '\n')
        model = svm.SVC(C=par_tup[0], kernel='rbf', gamma=par_tup[1], verbose=True)
        curr_log_file.write('SVC inited' + '\n')
        curr_log_file.write('fitting' + '\n')
        curr_log_file.flush()

        #model.fit(x_train, y_train)
        curr_log_file.write('dumping' + '\n')
        #pickle.dump(model, gzip.open(os.path.join(test_dir, 'model_' + str(i) + '.pkl.gz'), 'w'))
        curr_log_file.flush()
