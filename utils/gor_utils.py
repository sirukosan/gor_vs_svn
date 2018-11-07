import numpy as np

from utils.utils import *


def gor_train(profile_files, dssp_dir):
    res = None
    i = 0
    for profile_file in profile_files:
        i += 1
        processed_profile = process_profile(profile_file, get_relative_file(profile_file, dssp_dir, 'dssp'))
        if res is None:
            res = processed_profile
        else:
            res += processed_profile
    return res


def window_index(ss, i):
    if ss == 'h':
        return window // 2 + i
    elif ss == 'e':
        return window + window // 2 + i
    elif ss == 'c':
        return 2 * window + window // 2 + i
    elif ss == 't':
        return 3 * window + window // 2 + i


def process_profile(profile_csv_path, dssp_path):
    with open(dssp_path) as dssp:
        profile = restore_profile_from_csv(profile_csv_path)
        ss = (dssp.readlines()[1].lower().strip()).replace('-', 'c')
        res = np.zeros(shape=(window * len(SS_ORDER), len(AA_ORDER)))

        for i in range(0, profile.shape[0] - 1):
            for j in range(wind_l, wind_r + 1):
                if i + j < 0 or i + j >= profile.shape[0]:
                    continue
                curr_row = profile[i]
                res[window_index(ss[i], j)] += curr_row
                res[window_index('t', j)] += curr_row
        res = res / profile.shape[0]
        return res


def predict_profile(profile_csv_path, gor_csv_path):
    profile = restore_profile_from_csv(profile_csv_path)
    gor = restore_gor_from_csv(gor_csv_path)
    res = ''

    for i in range(0, profile.shape[0]):
        pred_score = 0
        pred_ss = SS_ORDER[0]
        for s in SS_ORDER[:-1]:
            new_pred_score = 0
            for j in range(wind_l, wind_r + 1):
                if i + j < 0 or i + j >= profile.shape[0]:
                    continue
                curr_prow = profile[i]
                curr_grow = gor[window_index(s, j)]
                new_pred_score += np.sum(curr_prow * curr_grow)
            if new_pred_score > pred_score:
                pred_score = new_pred_score
                pred_ss = s
        res += pred_ss
    return res


def save_gor_as_csv(numpy_gor_matrix):
    iterables = [SS_ORDER, list(range(wind_l, wind_r + 1))]
    m_index = pd.MultiIndex.from_product(iterables, names=['ss', 'pos'])
    pd_result = pd.DataFrame(data=numpy_gor_matrix, index=m_index, columns=AA_ORDER)
    pd_result.to_csv(GOR_TRAINED_FILE)


def restore_gor_from_csv(csv_file):
    with open(csv_file) as gor_csv:
        gor = pd.read_csv(gor_csv).values[:, 2:].astype(float)
        return gor
