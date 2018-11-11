from sys import stdout

import numpy as np
from sklearn import svm
import pickle
import gzip
import time
from utils.utils import *


def build_vectors_for_profile(profile_file, dssp_dir):
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
    res = []
    for s in dssp:
        if s.casefold() == 'h':
            res.append(0)
        elif s.casefold() == 'e':
            res.append(1)
        elif s.casefold() == '-' or s.casefold() == 'c':
            res.append(2)
    return res


def build_vectors(profile_files, out_x_file, out_y_file, dssp_dir, log_out=stdout, exclude=[]):
    i = 0
    with open(out_x_file, 'ab') as x_file, open(out_y_file, 'ab') as y_file:
        for profile_file in profile_files:
            if os.path.splitext(os.path.basename(profile_file))[0] not in exclude:
                log_out.write(str(i) + ' processing ' + profile_file + '\n')
                y, x = build_vectors_for_profile(profile_file, dssp_dir)
                np.savetxt(x_file, x, delimiter=",", fmt='%1.3f')
                np.savetxt(y_file, y, delimiter=",", fmt='%i')
                i += 1


def loop(i, test_dir, test_file, par_tup, x_train, y_train):
    with open(os.path.join(test_dir, 'log_' + str(i)), 'a') as curr_log_file:
        curr_log_file.write('building model for ' + test_file + '\nwith parameters: C = ' + str(
            par_tup[0]) + '\ngamma = ' + str(par_tup[1]) + '\n')
        model = svm.SVC(C=par_tup[0], kernel='rbf', gamma=par_tup[1], verbose=True)
        curr_log_file.write('SVC inited' + '\n')
        curr_log_file.write('fitting' + '\n')
        model.fit(x_train, y_train)
        curr_log_file.write('dumping' + '\n')
        pickle.dump(model, gzip.open(os.path.join(test_dir, 'model_' + str(i) + '.pkl.gz'), 'w'))
