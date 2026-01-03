import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys
from matplotlib.colors import LogNorm

sns.set_palette('bright')
sns.set_style('whitegrid')
# sns.set_theme(style="whitegrid")

plt.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 25,
    'axes.titlesize': 25,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'legend.fontsize': 25,
    'figure.titlesize': 30,
    'lines.markersize': 10,
    'lines.linewidth': 3,
    'mathtext.fontset': 'stix',
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
})

if __name__ == "__main__":
    df = pd.read_csv(f'logs/{sys.argv[1]}.csv')
    
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=False, tight_layout=True)

    axs = axs.flatten()

    df['prover'] = df['computing'] + df['proving']
    df['verifier'] = df['verifying']
    df['communication'] = df['comm'] / (1024 * 1024)  # Convert to MB
    xticks = [1, 4, 16, 64]
    xtick_labels = [str(x) for x in xticks]

    sns.lineplot(x='sigma', y='prover', data=df, style='dim', hue='dim', hue_norm=LogNorm(),
                 ax=axs[0], markers=True, dashes=True)
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_xticks(xticks)
    axs[0].set_xticklabels(xtick_labels)
    axs[0].set_xlabel('$\sigma$')
    axs[0].set_ylabel('Server RT (s)')

    

    sns.lineplot(x='sigma', y='communication', data=df, style='dim', hue='dim', hue_norm=LogNorm(),
                 ax=axs[1], markers=True, dashes=True)
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].set_xticks(xticks)
    axs[1].set_xticklabels(xtick_labels)
    axs[1].set_xlabel('$\sigma$')
    axs[1].set_ylabel('Communication (MB)')

    sns.lineplot(x='sigma', y='verifier', data=df, style='dim', hue='dim', hue_norm=LogNorm(),
                 ax=axs[2], markers=True, dashes=True)
    axs[2].set_xscale('log')
    axs[2].set_yscale('log')
    axs[2].set_xticks(xticks)
    axs[2].set_xticklabels(xtick_labels)
    axs[2].set_xlabel('$\sigma$')
    axs[2].set_ylabel('Verifier RT (s)')

    sns.lineplot(x='sigma', y='l1', data=df, style='dim', hue='dim', hue_norm=LogNorm(),
                 ax=axs[3], markers=True, dashes=True)
    axs[3].set_xscale('log')
    axs[3].set_yscale('log')
    axs[3].set_xticks(xticks)
    axs[3].set_xticklabels(xtick_labels)
    axs[3].set_xlabel('$\sigma$')
    axs[3].set_ylabel('L1 Error')

    for ax in axs:
        ax.get_legend().remove()

    handles, labels = axs[0].get_legend_handles_labels()
    labels = [f'$d={label}$' for label in labels]
    fig.legend(handles, labels, loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.06))

    plt.savefig('plots/vddgm.pdf', bbox_inches='tight')
