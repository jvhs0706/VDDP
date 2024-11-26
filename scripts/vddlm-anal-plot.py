import seaborn as sns

import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
    df = pd.read_csv('logs/vddlm-anal.csv')

    sns.set_theme(style="darkgrid")
    # plt.rcParams["font.family"] = "Liberation Serif"
    
    # 1 by 2 grid of plots
    fig, axs = plt.subplots(2, 1, figsize=(6, 6), sharex=True, sharey=False, tight_layout=True)
    # change names of algorithms
    df['algorithm'] = df['algorithm'].replace({'bc23': 'BC23', 'ours': 'Ours'})

    # plot eps vs num_bit for bc23 and ours
    sns.scatterplot(x = 'eps', y = 'num_coin', data = df, hue = 'algorithm', ax = axs[0])
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_xlabel('$\epsilon$', fontsize=15)
    axs[0].set_ylabel('Number of coins', fontsize=15)
    axs[0].legend()
    
    # plot eps vs l1 for bc23 and ours
    sns.scatterplot(x = 'eps', y = 'l1', data = df, hue = 'algorithm', ax = axs[1])
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].set_xlabel('$\epsilon$', fontsize=15)
    axs[1].set_ylabel('L1 error', fontsize=15)
    axs[1].legend()

    plt.savefig('logs/vddlm-anal-plot.pdf', bbox_inches='tight')


    

