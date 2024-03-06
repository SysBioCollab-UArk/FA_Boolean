import boolean2
from boolean2 import util
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def run_FA_Boolean_synch(model_text, n_sim_steps):
    model = boolean2.Model(model_text, mode='sync')
    model.initialize()
    model.iterate(steps=n_sim_steps)

    # get the names of the nodes and get the indices so we can order them the same as in Rodriguez et al. (2012)
    nodes = model.states[0].keys()
    nodes_order = ['ICL', 'FANCM', 'FAcore', 'FANCD2I', 'MUS81', 'FANCJBRCA1', 'XPF', 'FAN1', 'ADD', 'DSB', 'PCNATLS',
                   'MRN', 'BRCA1', 'ssDNARPA', 'FANCD1N', 'RAD51', 'HRR', 'USP1', 'KU', 'DNAPK', 'NHEJ', 'ATR', 'ATM',
                   'p53', 'CHK1', 'CHK2', 'H2AX', 'CHKREC']
    index = []
    for node in nodes_order:
        index.append(nodes.index(node))

    # collect the state values for each iteration
    states = []
    for state in model.states:
        # states.append([int(x) for x in state.values()])
        states.append([int(state.values()[index[i]]) for i in range(len(nodes))])  # ordered as Rodriguez et al. (2012)

    # create the heatmap for this initial condition
    # create_heatmap(states, nodes)
    create_heatmap(states, nodes_order)


def run_FA_Boolean_asynch(model_text, n_sim_steps, n_runs):
    # for storing trajectories
    coll = util.Collector()

    # run the simulations
    for i in range(n_runs):
        print(i)

        model = boolean2.Model(model_text, mode='async')
        model.initialize()
        model.iterate(steps=n_sim_steps)
        coll.collect(states=model.states, nodes=model.nodes)

    # get average node values
    avgs = coll.get_averages(normalize=True)

    # plots
    plt.figure(figsize=(6.4*1.5, 4.8*1.5))  # (6.4, 4.8)
    for sp in sorted(avgs.keys()):
        plt.plot(avgs[sp], lw=2, label=sp)
    plt.xticks(ticks=np.arange(0, n_sim_steps+1, 2))
    plt.ylim((-0.1, 1.1))
    plt.xlabel('iteration')
    plt.ylabel('probability ON')
    plt.legend(ncol=1, bbox_to_anchor=[1, 1], loc=0, labelspacing=0.2)
    # defaults: labelspacing=0.5, columnspacing=2)
    plt.tight_layout()


def create_heatmap(data, xtick_labels):

    fig, ax = plt.subplots(constrained_layout=True, figsize=(6.4*1.5, 4.8*1.5))  # (6.4, 4.8)

    # sns.set(font_scale=1.4)
    res = sns.heatmap(data, annot=False, vmin=0, vmax=1, fmt='.2f', linewidths=0.1, linecolor='grey', cmap='Greys')

    res.set_xticks(np.arange(len(xtick_labels))+0.5)
    res.set_xticklabels(xtick_labels, rotation=90)
    ax.xaxis.tick_top()

    res.set_yticklabels(np.arange(1, len(data)+1), rotation=0)
    ax.set_ylabel('Iteration')
