import os
import numpy as np
import pandas as pd
from constants import *


def gor_train(profiles_dir, dssp_dir):
    res = None
    i = 0
    for filename in os.listdir(profiles_dir):
        i += 1
        if res is None:
            res = process_profile(profiles_dir + filename, os.path.splitext(dssp_dir + filename)[0] + '.dssp')
        else:
            res += process_profile(profiles_dir + filename, os.path.splitext(dssp_dir + filename)[0] + '.dssp')
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
    with open(profile_csv_path) as profile_csv, open(dssp_path) as dssp:
        df = pd.read_csv(profile_csv).values[:, 1:].astype(float)
        ss = (dssp.readlines()[1].lower().strip()).replace('-', 'c')
        res = np.zeros(shape=(window * len(SS_ORDER), len(AA_ORDER)))

        for i in range(0, df.shape[0] - 1):
            for j in range(wind_l, wind_r + 1):
                if i + j < 0 or i + j >= df.shape[0]:
                    continue
                curr_row = df[i]
                res[window_index(ss[i], j)] += curr_row
                res[window_index('t', j)] += curr_row
        res = res / df.shape[0]
        return res


def predict_profile(profile_csv_path, gor_csv_path):
    with open(profile_csv_path) as profile_csv, open(gor_csv_path) as gor_csv:
        profile = pd.read_csv(profile_csv).values[:, 1:].astype(float)
        gor = pd.read_csv(gor_csv).values[:, 2:].astype(float)
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


# np_result = gor_train(PROFILES_OUTPUT_CSV_DIR, INPUT_DSSP_FOLDER)
#
# iterables = [SS_ORDER, list(range(wind_l, wind_r + 1))]
# m_index = pd.MultiIndex.from_product(iterables, names=['ss', 'pos'])
# pd_result = pd.DataFrame(data=np_result, index=m_index, columns=AA_ORDER)
#
# pd_result.to_csv(GOR_TRAINED_FILE)

print(predict_profile('data/profiles_blind/out/csv_profiles/5D14_A.csv', GOR_TRAINED_FILE))
