
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import seaborn as sns


def create_heatmap(data, xticks_labels):
    print(xticks_labels)
    print (data)

    fig, ax = plt.subplots(figsize=(30, 20))

    ax.xaxis.tick_top()

    sns.set(font_scale=1.4)
    res = sns.heatmap(data, annot=False, vmin=0, vmax=1, fmt='.2f', linewidths=0.1, linecolor='grey', cmap='Greys')

    res.set_xticklabels(xticks_labels, fontsize=14, rotation=90)
    # print([(i + 0.5) for i in range(len(data))])
    # print(['%d' % (i + 1) for i in range(len(data))])
    # print(len(data))
    # print(len(data[0]))
    # print(len(xticks_labels))
    # quit()


    # plt.yticks(ticks=[(i + 0.5) for i in range(len(data))], labels=['%d' % (i + 1) for i in range(len(data))],
    #            rotation=0, fontsize=14)

    plt.ylabel('Iteration', fontsize=14)

    plt.show()
