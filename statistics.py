import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from constants import *

sns.set_style("darkgrid")
sns.set(font_scale=2)
dssp_folder = INPUT_DSSP_FOLDER
fasta_folder = INPUT_FASTA_FOLDER


def init_stat():
    res = pd.DataFrame({'index': AA_ORDER,
                        'residue': AA_ORDER,
                        'helix': [0] * 20,
                        'strand': [0] * 20,
                        'coil': [0] * 20})
    res.set_index('index', inplace=True)
    return res


def init_heat_stat():
    return pd.DataFrame(0, index=AA_ORDER, columns=np.arange(-8, 9))


def fill_heat_stat(stat, fasta_str, i):
    start = i - 8 if (i - 8 >= 0) else 0
    end = i + 9 if (i + 9 <= len(fasta_str)) else len(fasta_str)
    for j in range(start, end):
        if fasta_str[j] != 'x':
            stat.at[fasta_str[j], j - i] += 1


def plot_pie(stat):
    sums = (stat.helix.sum(), stat.strand.sum(), stat.coil.sum())
    labels = 'Helix', 'Strand', 'Coil'
    explode = (0.0, 0.1, 0.0)
    fig, ax = plt.subplots(figsize=(15, 10))
    plt.axis('equal')
    plt.pie(sums, explode=explode, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)


def plot_bar(stat):
    fig, ax = plt.subplots(figsize=(15, 10))
    sns.barplot(data=stat, x='residue', y='perc', hue='ss', ax=ax)


def plot_heat(stat):
    fig, ax = plt.subplots(figsize=(15, 10))
    sns.heatmap(stat, cmap="YlGnBu", ax=ax)


def collect_stats():
    stat = init_stat()
    heat_stat_e = init_heat_stat()
    heat_stat_h = init_heat_stat()

    for file_name in os.listdir(dssp_folder):
        with open(dssp_folder + file_name) as dssp_file, \
                open(fasta_folder + file_name.replace('.dssp', '.fasta')) as fasta_file:
            dssp_str = dssp_file.readlines()[1].lower().strip()
            fasta_str = fasta_file.readlines()[1].lower().strip()
            for i in range(len(dssp_str)):
                if dssp_str[i] == '-' and fasta_str[i] != 'x':
                    stat.at[fasta_str[i], 'coil'] += 1
                elif dssp_str[i] == 'e':
                    fill_heat_stat(heat_stat_e, fasta_str, i)
                    stat.at[fasta_str[i], 'strand'] += 1
                elif dssp_str[i] == 'h':
                    fill_heat_stat(heat_stat_h, fasta_str, i)
                    stat.at[fasta_str[i], 'helix'] += 1

    stat.insert(1, 'total', stat.sum(axis=1))

    perc_stat = stat.copy()
    perc_stat[['total', 'helix', 'strand', 'coil']] = perc_stat[['total', 'helix', 'strand', 'coil']] * 100 / perc_stat[
        ['total', 'helix', 'strand', 'coil']].sum()

    colls = np.concatenate((np.arange(-8, 9)), axis=None)
    heat_stat_e[colls] = heat_stat_e[colls] * 100 / heat_stat_e[colls].sum()
    heat_stat_h[colls] = heat_stat_h[colls] * 100 / heat_stat_h[colls].sum()

    pie_stat = stat.copy()

    bar_stat = pd.melt(perc_stat, id_vars='residue', var_name="ss", value_name="perc")

    return pie_stat, bar_stat, heat_stat_e, heat_stat_h


pie_stat, bar_stat, heat_stat_e, heat_stat_h = collect_stats()

print('plotting pie')
plot_pie(pie_stat)
print('plotting bar')
plot_bar(bar_stat)
print('plotting heat e')
plot_heat(heat_stat_e)
print('plotting heat h')
plot_heat(heat_stat_h)
plt.show()

