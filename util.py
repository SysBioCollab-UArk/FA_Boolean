import boolean2
from boolean2 import util
import pysb
from pysb.simulator import BngSimulator
from pysb.importers.boolean import model_from_boolean
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def run_FA_Boolean_synch(model_text, n_sim_steps, outfile=None):
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
    create_heatmap(states, nodes_order, outfile=outfile)


def run_FA_Boolean_asynch(model_text, n_sim_steps, n_runs, outfile=None):
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
    plot_timecourses(np.arange(n_sim_steps+1), avgs, outfile=outfile)


def run_FA_Boolean_pysb(in_model, t_end, n_runs, verbose=True, outfile=None):
    if isinstance(in_model, pysb.core.Model):
        model = in_model
    else:
        model = model_from_boolean(in_model, mode='GSP')
    sim = BngSimulator(model, verbose=verbose)
    tspan = np.linspace(0, t_end, int(round(t_end) * 10 + 1))
    output = sim.run(tspan=tspan, n_runs=n_runs)
    avgs = {}
    for obs in model.observables:
        if '_True_' in obs.name:
            avgs[obs.name[:-9]] = np.mean(np.array(output.observables)[obs.name], axis=0)
    plot_timecourses(tspan, avgs, outfile=outfile)


def create_heatmap(data, xtick_labels, outfile=None):

    fig, ax = plt.subplots(constrained_layout=True, figsize=(6.4*1.5, 4.8*1.5))  # (6.4, 4.8)

    # sns.set(font_scale=1.4)
    res = sns.heatmap(data, annot=False, vmin=0, vmax=1, fmt='.2f', linewidths=0.1, linecolor='grey', cmap='Greys')

    res.set_xticks(np.arange(len(xtick_labels))+0.5)
    res.set_xticklabels(xtick_labels, rotation=90)
    ax.xaxis.tick_top()

    res.set_yticklabels(np.arange(1, len(data)+1), rotation=0)
    ax.set_ylabel('Iteration')

    if outfile is not None:
        plt.savefig(outfile, format='pdf')


def plot_timecourses(tspan, data, xlabel='iteration', outfile=None):

    groups = ['Mutations', 'Upstream FA/BRCA', 'Downstream FA/BRCA-HRR', 'NHEJ', 'Checkpoint Control']

    species_to_plot = [
        ['ICL', 'ADD', 'DSB'],
        ['FANCM', 'FAcore', 'FANCD2I', 'MUS81', 'FANCJBRCA1', 'XPF', 'FAN1'],
        ['PCNATLS', 'MRN', 'BRCA1', 'ssDNARPA', 'FANCD1N', 'RAD51', 'HRR', 'USP1'],
        ['KU', 'DNAPK', 'NHEJ'],
        ['ATR', 'ATM', 'p53', 'CHK1', 'CHK2', 'H2AX', 'CHKREC']
    ]

    plt.figure(figsize=(6.4 * 1.5, 4.8 * 2.25))  # (6.4, 4.8)
    ax1 = plt.subplot(8, 1, (1, 2))
    ax2 = plt.subplot(8, 2, (5, 9))
    ax3 = plt.subplot(8, 2, (6, 10))
    ax4 = plt.subplot(8, 2, (11, 15))
    ax5 = plt.subplot(8, 2, (12, 16))
    axes = [ax1, ax2, ax3, ax4, ax5]
    for ax, group, species in zip(axes, groups, species_to_plot):
        for sp in species:
            ax.plot(tspan, data[sp], lw=2, label=sp)
        ax.set_title(group)
        # ax.set_xticks(ticks=np.arange(0, n_sim_steps + 1, 2))
        ax.set_ylim((-0.1, 1.1))
        ax.set_xlabel(xlabel)
        ax.set_ylabel('probability ON')
        ncol = 2 if group == 'Checkpoint Control' else 1
        ax.legend(loc='upper right', labelspacing=0.2, ncol=ncol)
    plt.tight_layout()

    if outfile is not None:
        plt.savefig(outfile, format='pdf')
