import seaborn as sns

import pandas as pd
import matplotlib.pyplot as plt

import os, sys

if __name__ == "__main__":
    df = pd.read_csv(f'logs/{sys.argv[1]}.csv')

    sns.set_theme(style="whitegrid")
    plt.rcParams.update({'mathtext.fontset': 'stix',
        'font.family': 'serif',
        'font.serif': ['Times New Roman']})
    
    # 1 by 2 grid of plots
    fig, axs = plt.subplots(1, 2, figsize=(8, 3), sharex=True, sharey=False, tight_layout=True)
    # change names of algorithms
    df['algorithm'] = df['algorithm'].replace({'bc23': 'VDBM', 'ours': 'VDDLM'})

    # plot eps vs num_bit for bc23 and ours
    sns.scatterplot(x = 'eps', y = 'num_coin', data = df, hue = 'algorithm', style = 'algorithm', ax = axs[1], markers = ['o', 'P'], s=100)
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].set_xlabel('$\epsilon$', fontsize=15)
    axs[1].set_ylabel('Number of Coins', fontsize=15)
    axs[1].legend()
    
    # plot eps vs l1 for bc23 and ours
    sns.scatterplot(x = 'eps', y = 'l1', data = df, hue = 'algorithm', style = 'algorithm', ax = axs[0], markers = ['o', 'P'], s=100)
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_xlabel('$\epsilon$', fontsize=15)
    axs[0].set_ylabel('L1 Error', fontsize=15)
    axs[0].legend()

    # # plot eps vs l1 for bc23 and ours
    # sns.scatterplot(x = 'eps', y = 'l2', data = df, hue = 'algorithm', ax = axs[2])
    # axs[2].set_xscale('log')
    # axs[2].set_yscale('log')
    # axs[2].set_xlabel('$\epsilon$', fontsize=15)
    # axs[2].set_ylabel('L2 error', fontsize=15)
    # axs[2].legend()

    plt.savefig('plots/vddlm-anal-plot.pdf', bbox_inches='tight')


    

