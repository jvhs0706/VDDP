# use seaborn
import seaborn as sns

import matplotlib.pyplot as plt
import pandas as pd

import os, sys

sns.set_palette('bright')  # or 'deep', 'bright', etc.
sns.set_style('whitegrid')  # or 'darkgrid', 'white', 'dark', etc.

# ---- Font Size Settings ----
plt.rcParams.update({
    'font.size': 15,
    'axes.labelsize': 20,
    'axes.titlesize': 20,
    'xtick.labelsize': 15,
    'ytick.labelsize': 15,
    'legend.fontsize': 20,
    'figure.titlesize': 25,
    'lines.markersize': 10,   # default marker size
    'lines.linewidth': 3,     # default line width
})

if __name__ == '__main__':
    log_file = f'logs/{sys.argv[1]}.csv'
    df = pd.read_csv(log_file)

    df['prec'] = 2 ** df['log_prec']

    # plot (setup, committing, computing, proving, verifying) |vs inv_prec with error bar
    plt.figure(figsize=(15, 5))
    sns.lineplot(x = 'prec', y = 'setup', data = df, marker='o', label='Setup', markersize=10)
    sns.lineplot(x = 'prec', y = 'committing', data = df, marker='X', label='Committing', markersize=10)
    sns.lineplot(x = 'prec', y = 'computing', data = df, marker='s', label='Computing', markersize=10)
    sns.lineplot(x = 'prec', y = 'proving', data = df, marker='P', label='Proving', markersize=10)
    sns.lineplot(x = 'prec', y = 'verifying', data = df, marker='D', label='Verifying', markersize=10)
    plt.legend()

    # move the legend to the bottom

    plt.legend(loc='lower center', bbox_to_anchor=(0.475, -0.35), ncol=5)

    # y axis in log scale
    plt.xscale('log')
    plt.gca().set_xticks([64, 128, 256, 512, 1024, 2048, 4096])
    plt.gca().get_xaxis().set_major_formatter(plt.ScalarFormatter())
    plt.yscale('log')

    # x axis legend: precision bits
    plt.xlabel(r'$|\Omega|$')

    # y axis legend: time in seconds
    plt.ylabel('Time (s)')
    plt.savefig('plots/vrr-decomp.pdf', bbox_inches='tight')
