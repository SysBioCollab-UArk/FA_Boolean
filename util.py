import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def create_heatmap(data, xtick_labels):

    fig, ax = plt.subplots(constrained_layout=True, figsize=(6.4*1.5, 4.8*1.5))  # (6.4, 4.8)

    # sns.set(font_scale=1.4)
    res = sns.heatmap(data, annot=False, vmin=0, vmax=1, fmt='.2f', linewidths=0.1, linecolor='grey', cmap='Greys')

    res.set_xticks(np.arange(len(xtick_labels))+0.5)
    res.set_xticklabels(xtick_labels, rotation=90)
    ax.xaxis.tick_top()

    res.set_yticklabels(np.arange(1, len(data)+1), rotation=0)
    ax.set_ylabel('Iteration')

    plt.show()
