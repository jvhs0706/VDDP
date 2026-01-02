import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys
import pandas as pd

from matplotlib.colors import LogNorm
from matplotlib.cm import get_cmap

sns.set_palette('bright')  # or 'deep', 'bright', etc.
sns.set_style('whitegrid')  # or 'darkgrid', 'white', 'dark', etc.

# ---- Font Size Settings ----
plt.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 25,
    'axes.titlesize': 25,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'legend.fontsize': 25,
    'figure.titlesize': 30,
    'lines.markersize': 10,   # default marker size
    'lines.linewidth': 3,     # default line width
    'mathtext.fontset': 'stix',
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
})




def regather_data(df_vddlm, df_vdbm, item: str):
    eps_concatenated = np.concatenate([df_vddlm['eps'].values, df_vdbm['orig_eps'].values])
    dim_concatenated = np.concatenate([df_vddlm['dim'].values, df_vdbm['dim'].values])
    algo_concatenated = ['VDDLM'] * len(df_vddlm) + ['VDBM'] * len(df_vdbm)
    item_concatenated = np.concatenate([df_vddlm[f'vddlm_{item}'].values, df_vdbm[f'vdbm_{item}'].values])

    df = pd.DataFrame({
        'eps': eps_concatenated,
        '$d$': dim_concatenated,
        'Algorithm': algo_concatenated,
        item: item_concatenated
    })

    return df


if __name__ == '__main__':
    vddlm_log, vdbm_log = sys.argv[1], sys.argv[2]
    vddlm_df = pd.read_csv(f'logs/{vddlm_log}.csv')
    vdbm_df = pd.read_csv(f'logs/{vdbm_log}.csv')

    vddlm_df['vddlm_prover'] = vddlm_df['computing'] + vddlm_df['proving']
    vddlm_df['vddlm_comm'] = vddlm_df['comm'] / (1 << 20)
    vddlm_df['vddlm_verifier'] = vddlm_df['verifying']

    vdbm_df.rename(columns={
        'delta': 'vdbm_delta',
        'prover': 'vdbm_prover',
        'verifier': 'vdbm_verifier',
        'comm': 'vdbm_comm'
    }, inplace=True)

    xticks = [1e-3, 1e-2, 1e-1, 1]
    xtick_labels = [f'{x:.3f}' for x in xticks]

    # --- Plot Settings ---
    fig, axs = plt.subplots(1, 3, figsize=(25, 5), sharex = True)
    xticks = [1e-3, 1e-2, 1e-1, 1]
    xtick_labels = [f'{x:.3f}' for x in xticks]
    linestyles = ['-', ':']  # Solid for VDDLM, dashed for VDBM
    markers = ['o', 's', '^', 'P', 'X', 'D', '*', 'v', 'H']  # Enough for up to 9 unique $d$ values
    cmap = plt.cm.tab10

    # --- Helper to plot one metric ---
    
    def plot_metric(ax, df, metric, ylabel, logabove=None):
        dims = sorted(df['$d$'].unique())
        algos = ['VDDLM', 'VDBM']
        norm = LogNorm(vmin=min(dims), vmax=max(dims))

        for i, d in enumerate(dims):
            color = cmap(norm(d))
            marker = markers[i % len(markers)]
            for j, algo in enumerate(algos):
                linestyle = linestyles[j % len(linestyles)]
                
                subset = df[(df['$d$'] == d) & (df['Algorithm'] == algo)]
                sns.lineplot(x='eps', y=metric, data=subset,
                        ax=ax, label=f'$d$={d}, {algo}',
                        color=color, marker=marker, linestyle=linestyle, errorbar=('sd', 0.95))
        ax.legend_.remove()
                

        ax.set_xscale('log')
        if logabove is not None:
            ax.set_yscale('symlog', linthresh=logabove, linscale=10)
            ax.set_ylim(bottom=0)  # prevent showing negatives
        else:
            ax.set_yscale('log')
        ax.set_xticks(xticks)
        ax.set_xticklabels(xtick_labels)
        ax.set_xlabel('$\epsilon$')
        ax.set_ylabel(ylabel)

    # --- Plot each metric ---
    for ax, metric, ylabel, logabove in zip(axs, ['prover', 'comm', 'verifier'],
                                ["Servers' RT (s)", 'Communication (MB)', "Verifier's RT (s)"], [None, 5, 30]):
        ax_data = regather_data(vddlm_df, vdbm_df, metric)
        plot_metric(ax, ax_data, metric, ylabel, logabove)

    # --- Legend ---
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=5, bbox_to_anchor=(0.5, -0.2))

    plt.tight_layout()
    plt.savefig('plots/vddlm.pdf', format='pdf', bbox_inches='tight')


    # Another plot for L1 error
    fig, ax = plt.subplots(figsize=(12, 4))
    # plot three lines
    vddlm_df['orig_l1'] /= vddlm_df['dim']
    vddlm_df['vdbm_l1'] /= vddlm_df['dim']
    vddlm_df['l1'] /= vddlm_df['dim']
    sns.lineplot(data=vddlm_df, x='orig_eps', y='orig_l1', label='DDLM (expectation)', marker='o', linestyle=':', errorbar = None) # no error bars
    sns.lineplot(data=vddlm_df, x='orig_eps', y='vdbm_l1', label='DBM/VDBM (expectation)', marker='s', linestyle='--', errorbar = None)
    sns.lineplot(data=vddlm_df, x='eps', y='l1', label='VDDLM', marker='^', linestyle='-')
    ax.set_xscale('log')
    # ax.set_yscale('symlog', linthresh=10**1.5, linscale=5)  # use symlog for better visibility of small values
    ax.set_yscale('log')
    ax.set_ylim(bottom=0)  # prevent showing negatives
    ax.set_xlabel('$\epsilon$')
    ax.set_ylabel('L1 Error\n(Per-Dimension)')
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)
    ax.legend(loc='lower left', fontsize=20)
    plt.tight_layout()
    plt.savefig('plots/vddlm-l1.pdf', format='pdf', bbox_inches='tight')


    