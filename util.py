import boolean2
from boolean2 import util
import pysb
from pysb.simulator import BngSimulator
from pysb.importers.boolean import model_from_boolean
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import cycle
import pandas as pd


def run_FA_Boolean_synch(model_text, n_sim_steps, detect_cycles=False, outfile=None):
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
    start, end = util.detect_cycles(model.states) if detect_cycles else (0, len(model.states))
    if detect_cycles and (start, end) == (0, 0):  # make sure a steady state is reached
        raise Exception('No cycles detected. Try increasing the number of synchronous updating iterations.')
    states = []
    for state in model.states[start:start+end]:
        states.append([int(state.values()[index[i]]) for i in range(len(nodes))])  # ordered as Rodriguez et al. (2012)

    # create the heatmap for this initial condition
    create_heatmap(states, nodes_order, outfile=outfile)


def run_FA_Boolean_asynch(model_text, n_sim_steps, n_runs, species_list, verbose=True, outfile=None, **kwargs):
    # for storing trajectories
    coll = util.Collector()

    # run the simulations
    for i in range(n_runs):
        if verbose:
            print(i)

        model = boolean2.Model(model_text, mode='async')
        model.initialize()
        model.iterate(steps=n_sim_steps)
        coll.collect(states=model.states, nodes=model.nodes)

    # get average node values
    avgs = coll.get_averages(normalize=True)

    # plots
    plot_timecourses(np.arange(n_sim_steps+1), avgs, outfile=outfile)


def run_FA_Boolean_pysb(in_model, t_end, n_runs, species_list, param_values=None, verbose=True, outfile=None, **kwargs):

    # run simulations
    if isinstance(in_model, pysb.core.Model):
        model = in_model
    else:
        model = model_from_boolean(in_model, mode='GSP')
    sim = BngSimulator(model, verbose=verbose)
    tspan = np.linspace(0, t_end, int(round(t_end) * 10 + 1))
    output = sim.run(tspan=tspan, n_runs=n_runs, param_values=param_values)
    avgs = {}
    for obs in model.observables:
        if '_True_' in obs.name:
            avgs[obs.name[:-9]] = np.mean(np.array(output.observables)[obs.name], axis=0)
    if outfile is not None:
        fig  = plot_timecourses(tspan, avgs, species_list, outfile=outfile)

    # calculate times that species thresholds are reached
    species = kwargs.get('species', [])
    if isinstance(species, str):
        species = [species]
    if len(species) > 0:
        if outfile is not None:
            axs = fig.get_axes()
            color_iter = [cycle(plt.rcParams['axes.prop_cycle'].by_key()['color']) for _ in species_list]
            color = [None] * len(species_list)
        sp_thresh = kwargs.get('threshold', 0.01)
        t_thresh = []
        for sp in species:
            ax_idx = next(i for i in range(len(species_list)) if sp in species_list[i])  # which plot
            max_idx = np.argmax(avgs[sp])  # index the species is at its max
            # threshold time
            t_thresh.append(next((t for t, val in zip(tspan[max_idx:], avgs[sp][max_idx:]) if val < sp_thresh), np.inf))
            if outfile is not None:
                if color[ax_idx] is None:  # plot horizontal dashed line at threshold
                    axs[ax_idx].axhline(sp_thresh, ls='--', color='k')
                color[ax_idx] = next(color_iter[ax_idx])
                if t_thresh[-1] < np.inf:
                    axs[ax_idx].axvline(t_thresh[-1], ls='--', color=color[ax_idx])  # plot vertical dashed lines
        return t_thresh

    return None


def create_heatmap(data, xtick_labels, outfile=None):

    fig, ax = plt.subplots(constrained_layout=True, figsize=(6.4*1.5, 1+0.25*len(data)))  # 4.8*1.5))  # (6.4, 4.8)

    # sns.set(font_scale=1.4)
    res = sns.heatmap(data, annot=False, vmin=0, vmax=1, fmt='.2f', linewidths=0.1, linecolor='grey', cmap='Greys',
                      cbar=None)

    res.set_xticks(np.arange(len(xtick_labels))+0.5)
    res.set_xticklabels(xtick_labels, rotation=90)
    ax.xaxis.tick_top()

    res.set_yticks(np.arange(1, len(data) + 1)-0.5)
    res.set_yticklabels(np.arange(1, len(data)+1), rotation=0)
    # ax.set_ylabel('Iteration')

    if outfile is not None:
        plt.savefig(outfile, format='pdf')


def plot_timecourses(tspan, data, species_to_plot, xlabel='iteration', outfile=None):

    groups = ['Mutations', 'Upstream FA/BRCA', 'Downstream FA/BRCA-HRR', 'NHEJ', 'Checkpoint Control']

    fig = plt.figure(figsize=(6.4 * 1.5, 4.8 * 2.25))  # (6.4, 4.8)
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

    return fig


def plot_repair_time_distributions(csv_file):
    df = pd.read_csv(csv_file)

    value_cols = df.columns[-3:]   # e.g. ICL, ADD, DSB
    group_col = df.columns[1]      # initials
    mut_col = df.columns[0]        # mutations

    df[mut_col] = df[mut_col].fillna("WT")

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # Ensure numeric, preserving inf and nan
    for col in value_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    for group_value, subdf in df.groupby(group_col):
        subdf = subdf.copy()
        title_str = str(group_value).replace("_", " = ", 1).replace("_", " ")

        # ------------------------------------------------------------
        # Determine mutation order from the FIRST repair-time column
        # using the median of finite values only.
        # Categories with no finite values are placed at the end.
        # ------------------------------------------------------------
        base_col = value_cols[0]
        ordering_info = []

        for mut_value, mut_df in subdf.groupby(mut_col):
            vals = mut_df[base_col].to_numpy(dtype=float)
            finite_vals = vals[np.isfinite(vals)]
            inf_count = np.isinf(vals).sum()
            nan_count = np.isnan(vals).sum()

            if len(finite_vals) > 0:
                sort_key = np.median(finite_vals)
                status_rank = 0
            elif inf_count > 0:
                sort_key = np.inf
                status_rank = 1
            else:
                # all NaN (or effectively no usable values)
                sort_key = np.inf
                status_rank = 2

            ordering_info.append((
                0 if mut_value == "WT" else status_rank + 1,  # WT always first
                sort_key,
                mut_value
            ))

        ordering_info.sort(key=lambda x: (x[0], x[1], str(x[2])))
        mutation_order = [x[2] for x in ordering_info]

        # ------------------------------------------------------------
        # Create one figure with 3 subplots
        # ------------------------------------------------------------
        fig, axes = plt.subplots(
            1, len(value_cols),
            figsize=(18, 13),
            constrained_layout=True
        )

        if len(value_cols) == 1:
            axes = [axes]

        for j, col in enumerate(value_cols):
            ax = axes[j]

            grouped_finite = []
            inf_counts = []
            nan_counts = []
            labels = []

            for mut_value in mutation_order:
                mut_df = subdf[subdf[mut_col] == mut_value]
                vals = mut_df[col].to_numpy(dtype=float)

                finite_vals = vals[np.isfinite(vals)]
                inf_count = np.isinf(vals).sum()
                nan_count = np.isnan(vals).sum()

                grouped_finite.append(finite_vals)
                inf_counts.append(inf_count)
                nan_counts.append(nan_count)
                labels.append(
                    r"$\bf{WT}$" if mut_value == "WT"
                    else str(mut_value).replace("_", " = ")
                )

            positions = np.arange(len(labels))

            # Finite values for scaling
            finite_arrays = [vals for vals in grouped_finite if len(vals) > 0]
            all_finite = np.concatenate(finite_arrays) if len(finite_arrays) > 0 else np.array([])

            if len(all_finite) == 0:
                # No finite data at all for this subplot
                ax.set_yticks(positions)
                ax.set_yticklabels(labels)
                ax.tick_params(labelsize=14)
                ax.set_ylim(-0.5, len(labels) - 0.5)
                ax.invert_yaxis()
                ax.set_xlabel("Repair time", fontsize=16)
                ax.set_title(col, fontsize=16, fontweight='bold', color=colors[j % len(colors)])
                ax.text(
                    0.5, 0.5, "No finite values",
                    ha="center", va="center",
                    transform=ax.transAxes
                )
                continue

            max_finite = np.max(all_finite)
            annot_x = max_finite * 1.18
            x_max = max_finite * 1.32

            # Boxplot only for categories with finite values
            bp_data = []
            bp_pos = []
            for pos, vals in enumerate(grouped_finite):
                if len(vals) > 0:
                    bp_data.append(vals)
                    bp_pos.append(pos)

            ax.boxplot(
                bp_data,
                positions=bp_pos,
                vert=False,
                widths=0.6,
                patch_artist=True,
                boxprops=dict(facecolor=colors[j % len(colors)], alpha=0.5),
                medianprops=dict(color='black'),
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                flierprops=dict(
                    marker='o',
                    markersize=4,
                    markerfacecolor=colors[j % len(colors)],
                    alpha=0.5
                )
            )

            # Overlay finite replicate values
            rng = np.random.default_rng(0)
            for pos, vals in enumerate(grouped_finite):
                if len(vals) > 0:
                    jitter = rng.uniform(-0.12, 0.12, size=len(vals))
                    ax.scatter(
                        vals,
                        np.full(len(vals), pos) + jitter,
                        color=colors[j % len(colors)],
                        alpha=0.8,
                        s=25,
                        edgecolors='none',
                        label=None
                    )

            # Add dashed line for WT mean
            if "WT" in mutation_order:
                wt_index = mutation_order.index("WT")
                wt_vals = grouped_finite[wt_index]

                if len(wt_vals) > 0:
                    wt_mean = np.mean(wt_vals)

                    ax.axvline(
                        wt_mean,
                        linestyle='--',
                        linewidth=2,
                        color='black',
                        alpha=0.7
                    )

            # Annotate inf / nan counts at right edge
            text_x = x_max - max_finite * 0.02

            for pos, (inf_count, nan_count) in enumerate(zip(inf_counts, nan_counts)):
                text_parts = []
                if inf_count > 0:
                    text_parts.append(rf"$\infty \times {inf_count}$")
                if nan_count > 0:
                    text_parts.append(f"N/A × {nan_count}")

                if text_parts:
                    ax.text(
                        text_x,
                        pos,
                        "   ".join(text_parts),
                        va='center',
                        ha='right',
                        fontsize=16,
                        color='red'
                    )

            ax.set_yticks(positions)
            ax.set_yticklabels(labels)
            ax.tick_params(labelsize=14)
            ax.set_ylim(-1, len(labels))
            ax.invert_yaxis()
            ax.set_xlim(0, x_max)
            ax.set_xlabel("Repair time", fontsize=16)
            ax.set_title(col, fontsize=16, fontweight='bold', color=colors[j % len(colors)])

        fig.suptitle("Initial: " + title_str, fontsize=16, fontweight='bold')
        plt.show()


if __name__ == '__main__':
    csvfile = open('REPAIR_TIMES_pysb.csv', 'r')
    # plot_separate_sorted_bars(csvfile)
    plot_repair_time_distributions(csvfile)
