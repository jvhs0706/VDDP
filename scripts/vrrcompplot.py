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
    'font.size': 15,
    'axes.labelsize': 20,
    'axes.titlesize': 20,
    'xtick.labelsize': 15,
    'ytick.labelsize': 15,
    'legend.fontsize': 20,
    'figure.titlesize': 25,
    'lines.markersize': 10,   # default marker size
    'lines.linewidth': 3     # default line width
})



def regather_data(df_ours, df_kcy21, item: str):
    prec_concatenated = np.concatenate([1 << df_ours['log_prec'].values, 1 << df_kcy21['log_prec'].values])
    num_class_concatenated = np.concatenate([df_ours['num_class'].values, df_kcy21['num_class'].values])
    algo_concatenated = np.array(['Ours'] * len(df_ours) + ['KCY21'] * len(df_kcy21))
    def get_category(nclass, algo):
        return f'{algo}, $K={nclass}$'
    get_category_vectorized = np.vectorize(get_category)
    catetories = get_category_vectorized(num_class_concatenated, algo_concatenated)
    
    item_concatenated = np.concatenate([df_ours[item].values, df_kcy21[item].values])

    df = pd.DataFrame({
        'prec': prec_concatenated,
        '$K$': num_class_concatenated,
        'Algorithm': algo_concatenated,
        'cat': catetories,
        item: item_concatenated
    })

    return df


if __name__ == '__main__':
    vrr_ours_log, vrr_kcy21_log = sys.argv[1], sys.argv[2]
    ours_df = pd.read_csv(f'logs/{vrr_ours_log}.csv')
    kcy21_df = pd.read_csv(f'logs/{vrr_kcy21_log}.csv')

    valid_precs = [64, 256, 1024]  # Valid precisions for Ours
    valid_num_classes = [8, 32, 128]

    ours_df = ours_df[ours_df['log_prec'].isin(np.log2(valid_precs)) &
                      ours_df['num_class'].isin(valid_num_classes)]
    kcy21_df = kcy21_df[kcy21_df['log_prec'].isin(np.log2(valid_precs)) &
                        kcy21_df['num_class'].isin(valid_num_classes)]

    ours_df['client'] = ours_df['committing'] + ours_df['computing'] + ours_df['proving']
    ours_df['comm'] /= (1 << 20)  # Convert to MB
    ours_df['verifier'] = ours_df['verifying']

    # --- Plot Settings ---
    fig, axs = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    hatch_patterns = ['o', '+', '*', '/o', '/+', '/*']

    # --- Helper to plot one metric ---
    
    def plot_metric(ax, df, metric, ylabel):
        precs = sorted(df['prec'].unique())
        nclasses = sorted(df['$K$'].unique())
        hue_order = [f'Ours, $K={K}$' for K in nclasses] + [f'KCY21, $K={K}$' for K in nclasses]

        # Define unique hatches for each category
        
        hatch_map = {cat: hatch_patterns[i % len(hatch_patterns)] for i, cat in enumerate(hue_order)}

        # Plot with seaborn
        barplot = sns.barplot(x='prec', y=metric, hue='cat', data=df, ax=ax,
                            palette='tab10', hue_order=hue_order)

        # Apply hatch to each bar
        for container, cat in zip(barplot.containers, hue_order):
            hatch = hatch_map[cat]
            for bar in container:
                bar.set_hatch(hatch)

        ax.set_yscale('log')
        ax.set_xlabel(r'$|\Omega|$')
        ax.set_ylabel(ylabel)
        ax.legend().remove()


    # --- Plot each metric ---
    for ax, metric, ylabel in zip(axs, ['client', 'comm', 'verifier'],
                                ["Client's RT (s)", 'Communication (MB)', "Verifier's RT (s)"]):
        ax_data = regather_data(ours_df, kcy21_df, metric)
        plot_metric(ax, ax_data, metric, ylabel)

    import matplotlib.patches as mpatches

    # Create custom legend patches with hatches
    hue_order = [f'Ours, $K={K}$' for K in sorted(ours_df['num_class'].unique())] + \
                [f'KCY21, $K={K}$' for K in sorted(kcy21_df['num_class'].unique())]

    hatch_map = {cat: hatch_patterns[i % len(hatch_patterns)] for i, cat in enumerate(hue_order)}

    legend_patches = []
    cmap = plt.cm.tab10
    for i, cat in enumerate(hue_order):
        patch = mpatches.Patch(
            facecolor=cmap(i % 10),
            hatch=hatch_map[cat],
            label=cat,
            edgecolor='white'
        )
        legend_patches.append(patch)

    fig.legend(handles=legend_patches, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.08))


    plt.tight_layout()
    plt.savefig('plots/vrr.pdf', format='pdf', bbox_inches='tight')

    