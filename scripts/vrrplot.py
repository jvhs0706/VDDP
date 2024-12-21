# use seaborn
import seaborn as sns

import matplotlib.pyplot as plt
import pandas as pd

if __name__ == '__main__':
    df = pd.read_csv('logs/vrr.csv')
    sns.set_theme(style="darkgrid")

    df['prec'] = 2 ** df['log_prec']

    # plot (setup, committing, computing, proving, verifying) |vs inv_prec with error bar
    plt.figure(figsize=(6, 4))
    sns.lineplot(x = 'prec', y = 'setup', data = df, marker='o', label='Setup')
    sns.lineplot(x = 'prec', y = 'committing', data = df, marker='o', label='Committing')
    sns.lineplot(x = 'prec', y = 'computing', data = df, marker='o', label='Computing')
    sns.lineplot(x = 'prec', y = 'proving', data = df, marker='o', label='Proving')
    sns.lineplot(x = 'prec', y = 'verifying', data = df, marker='o', label='Verifying')
    plt.legend()

    # move the legend to the bottom

    plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.35), ncol=3)

    # y axis in log scale
    plt.xscale('log')
    plt.gca().set_xticks([64, 128, 256, 512, 1024, 2048, 4096])
    plt.gca().get_xaxis().set_major_formatter(plt.ScalarFormatter())
    plt.yscale('log')

    # x axis legend: precision bits
    plt.xlabel('$|\Omega|$')

    # y axis legend: time in seconds
    plt.ylabel('Time (s)')
    plt.savefig('logs/vrr-running-times.pdf', bbox_inches='tight')

    # make a 2 * 4 grid of plots
    fig, axs = plt.subplots(2, 4, figsize=(12, 6), sharex=True, sharey=True, tight_layout=True)

    # plot (setup, committing, computing, proving, verifying) vs inv_prec with error bar

    for ax, NUM_CLASS in zip(axs.flatten(), df['num_class'].unique()):
        sns.lineplot(x = 'log_prec', y = 'actual_eps', hue = 'eps', data = df[df['num_class'] == NUM_CLASS], marker='o', legend='full', ax = ax)

        ax.set_title(f'num_class = {NUM_CLASS}', fontsize='medium')
        # x, y axis in log2 scale
        ax.set_yscale('log')
        ax.set_yticks([0.125, 0.25, 0.5, 1, 2, 4, 8])
        # also display these numbers
        ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        ax.get_yaxis().set_major_formatter(plt.ScalarFormatter())

        # add a title to the legend
        ax.legend(title='$\epsilon$')


        # x axis legend: precision bits
        ax.set_xlabel('Precision bits', fontsize='medium')

        # y axis legend: time in seconds
        ax.set_ylabel('Actual $\epsilon$', fontsize='medium')
        # grasp one legend and put the legend at the buttom, detele the legends in the subplots
        handles, labels = ax.get_legend_handles_labels()
        ax.legend().remove()

    # put the legend at the buttom
    # move a little bit lower so as not to block the figure
    fig.legend(handles, labels, bbox_to_anchor=(0.5, -0.1), loc='lower center', ncol=8, title='Desired $\epsilon$', fontsize='medium')

    plt.show()
    plt.savefig('logs/vrr-precision.pdf', bbox_inches='tight')
