import numpy as np
from typing import Dict
import matplotlib.pyplot as plt
import scipy.stats as ss
import seaborn as sns
import math
from matplotlib.lines import Line2D
from statsmodels.stats.multitest import multipletests

def mutation_comparison(
    sub_0_counts: np.ndarray,
    sub_1_counts: np.ndarray,
    mut2idx: Dict[str, int],
    nmer4norm=None,
    title=r"$log_{2}$" +" ratio of singleton fractions\n in strains with D vs. B haplotypes at QTL",
    outname='heatmap.png',
    plot_type='heatmap',
    vmin: float = -1,
    vmax: float = 1,
):
    """
    plot a comparison of mutation spectra in one subset of strains
    vs. another, using either a heatmap or scatter plot
    sub_0_counts: np.array() of shape (n_strains, n_muts)
    sub_1_counts: np.array() of shape (n_strains, n_muts)
    mut2idx: dictionary mapping mutation types to indices
    """

    plt.rc('font', size=14)

    # map "base" mutation types to output indices
    # in the final plot. the heatmap will have 6
    # "blocks" of 16 mutation types stacked on top
    # of one another, with each block of 16 mutation
    # types corresponding to a single "base" mutation
    mut_out_idx = dict(
        zip(['C>A', 'C>G', 'C>T', 'A>G', 'A>C', 'A>T'],
            np.arange(6) * 4))

    five_pr_idx = dict(zip(['T', 'G', 'C', 'A'], range(4)))
    three_pr_idx = dict(zip(['A', 'C', 'G', 'T'], range(4)))

    # get a list of all 96 possible mutation types
    muts = np.array(list(mut2idx.keys()))

    # convert the counts of each 3-mer mutation type in each
    # subset to fractions
    sub_0_fracs = sub_0_counts / np.sum(sub_0_counts)
    sub_1_fracs = sub_1_counts / np.sum(sub_1_counts)

    n_muts = 96

    # the output heatmap is shape (16, 4)
    out_shape = (int(n_muts / 4), 4)

    # make an array where we'll store log2 ratios
    # of kmer counts for each mutation type
    out_array = np.zeros(out_shape)

    # make an array where we'll store the raw counts of
    # each mutation type in each subset, formatted to match
    # the shape of the 16 * 4 mutation matrix
    sub_0_counts_grid = np.zeros(out_shape)
    sub_1_counts_grid = np.zeros(out_shape)

    # do the same, but for storing p-values of each
    # sub_0 vs. sub_1 comparison
    pvals = np.zeros(out_shape)

    muts_out = np.zeros(out_shape, dtype=object)
    base_muts_out = np.zeros(out_shape, dtype=object)

    # if we want to normalize counts of mutations by the
    # sum of 3-mer nucleotides of each type across the strains
    # in subset 0 or subset 1, we query the nmer4norm dataframe
    if nmer4norm is not None:

        # loop over mutation and its total in sub_0
        for mut_idx, x in enumerate(sub_0_counts):

            orig_mut = muts[mut_idx]
            mut_nmer = orig_mut.split('>')[0]

            # get that mutation's total in sub_1
            y = sub_1_counts[mut_idx]

            # get the sum of 3-mer nucleotides of the base mutation
            # type we're looking at in subset_0 and subset_1
            nmer_comp = nmer4norm[nmer4norm['nmer'] == mut_nmer]

            nmer_sub0_comp = nmer_comp[nmer_comp['haplotype'] ==
                                       0]['count_sum'].values[0]
            nmer_sub1_comp = nmer_comp[nmer_comp['haplotype'] ==
                                       1]['count_sum'].values[0]

            # if the sum of 3-mer nucleotides is higher in subset_1,
            # we adjust the counts of mutations in subset_1 *down*, and
            # vice versa
            if nmer_sub0_comp > nmer_sub1_comp:
                scaling_f = nmer_sub1_comp / nmer_sub0_comp
                sub_0_counts[mut_idx] = int(x * scaling_f)
                sub_1_counts[mut_idx] = y
            elif nmer_sub1_comp > nmer_sub0_comp:
                scaling_f = nmer_sub0_comp / nmer_sub1_comp
                sub_1_counts[mut_idx] = int(y * scaling_f)
                sub_0_counts[mut_idx] = x

    for mut_idx, x in enumerate(sub_0_counts):

        # get that mutation's total in sub_1
        y = sub_1_counts[mut_idx]

        if any([c == 0 for c in [x, np.sum(sub_0_counts), y, np.sum(sub_1_counts)]]): 
            p = 1
        else:
            # chi2 contingency table to compare mutation
            # frequencies in either subset
            _, p, _, _ = ss.chi2_contingency(
                [[x, np.sum(sub_0_counts)], [y, np.sum(sub_1_counts)]], )

        orig_mut = muts[mut_idx]

        # access the "base" 1-mer mutation and its 5'
        # and 3' flanking nucleotides
        base_mut = orig_mut.split('>')[0][1] + '>' + orig_mut.split('>')[1][1]

        five_pr, three_pr = orig_mut.split('>')[0][0], orig_mut.split(
            '>')[0][-1]

        # calculate the index of this particular 3-mer
        # in the output (16 * 4) array
        out_idx_x = mut_out_idx[base_mut] + five_pr_idx[five_pr]
        out_idx_y = three_pr_idx[three_pr]

        # store the raw counts of each mutation type in the
        # formatted (16 * 4) grid
        sub_0_counts_grid[out_idx_x, out_idx_y] = x
        sub_1_counts_grid[out_idx_x, out_idx_y] = y

        # get fractions rather than counts for plotting ratios
        x_frac = sub_0_fracs[mut_idx]
        y_frac = sub_1_fracs[mut_idx]

        if x_frac == 0 or y_frac == 0: ratio = 0
        else: ratio = math.log(x_frac / y_frac, 2)

        out_array[out_idx_x, out_idx_y] = ratio
        muts_out[out_idx_x, out_idx_y] = orig_mut
        base_muts_out[out_idx_x, out_idx_y] = base_mut
        pvals[out_idx_x, out_idx_y] = p

    # find indices where ratio of sub_0:sub_1 is significant
    # at Bonferonni-corrected p-value
    flat_pvals = pvals.flatten()
    reject, adj_pvals, _, _ =  multipletests(flat_pvals, alpha=0.05, method="fdr_bh")
    pvals = adj_pvals.reshape(out_array.shape)
    sig_pvals = np.where(pvals < 0.05)# / flat_pvals.shape[0]))
    non_sig_pvals = np.where(pvals >= 0.05)

    sns.set_style('ticks')

    if plot_type == "scatter":

        f, ax = plt.subplots()

        sns.set_style('ticks')

        # map base mutation types to colors to aesthetically
        # group 3-mer mutation in plot
        mut_colors = dict(
            zip(["A>C", "A>G", "A>T", "C>A", "C>G", "C>T"],
                sns.color_palette('colorblind', 6)[::-1]))

        x_max, y_max = 0, 0

        # plot data, and record maximum x and y values
        # so that we can plot an abline later

        # plot significant values first, with a black stroke
        # on outside of circles
        for (y, x) in zip(sig_pvals[0], sig_pvals[1]):
            x_frac = sub_0_counts_grid[y, x] / np.sum(sub_0_counts_grid)
            y_frac = sub_1_counts_grid[y, x] / np.sum(sub_1_counts_grid)

            if y_frac > y_max: y_max = y_frac
            if x_frac > x_max: x_max = x_frac

            base_mut = base_muts_out[y, x]
            ax.scatter(x_frac,
                       y_frac,
                       edgecolor='k',
                       color=mut_colors[base_mut],
                       s=75)

        # then plot non-significant values, with a white stroke
        # on outside of circles
        for (y, x) in zip(non_sig_pvals[0], non_sig_pvals[1]):
            x_frac = sub_0_counts_grid[y, x] / np.sum(sub_0_counts_grid)
            y_frac = sub_1_counts_grid[y, x] / np.sum(sub_1_counts_grid)

            if y_frac > y_max: y_max = y_frac
            if x_frac > x_max: x_max = x_frac

            base_mut = base_muts_out[y, x]

            ax.scatter(
                x_frac,
                y_frac,
                edgecolor='w',
                color=mut_colors[base_mut],
                s=75,
            )

        ax.axline((0, 0), (x_max, y_max), color='k', zorder=0)

        # generate a custom legend, describing the above
        # mappings of mutations to colors
        custom_legend = [
            Line2D(
                [0],
                [0],
                marker='o',
                color='w',
                label=mut,
                markerfacecolor=mut_colors[mut],
                markersize=10,
            ) for mut in mut_colors
        ]

        ax.legend(handles=custom_legend, frameon=False)

        ax.set_xlabel("Fraction of singletons in CC\nstrains with IAP element")
        ax.set_ylabel("Fraction of singletons in CC\nstrains without IAP element")

        sns.despine(ax=ax, top=True, right=True)

        f.tight_layout()

        f.savefig(outname, bbox_inches='tight', dpi=300)

    elif plot_type == "heatmap":
        f, ax = plt.subplots(figsize=(4, 8))

        cmap = sns.diverging_palette(230, 15, as_cmap=True)

        sns.heatmap(
            out_array,
            cmap="coolwarm",
            edgecolor='w',
            vmin=vmin,
            vmax=vmax,
        )

        # plot "dots" in heatmap where the ratio is significant
        for (y, x) in zip(sig_pvals[0], sig_pvals[1]):
            
            ax.scatter(x + 0.5, y + 0.5, c='w', edgecolor='k')

        # add boundary lines between each block of 16 mutations
        for x in np.arange(0, out_shape[0], 4):
            ax.axhline(y=x, ls=':', c='k')

        # add in y-axis labels
        ylabs = []
        for i, m in enumerate(muts_out[:, 0]):
            m_split = m.split('>')
            ylab = None
            if i == 0: ylab = "5'-" + m_split[0][0]
            else:
                if i in (2, 6, 10, 14, 18, 22):
                    ylab = m_split[0][1] + r'$\to$' + m_split[1][
                        1] + r'  ' + m_split[0][0]
                else:
                    ylab = m_split[0][0]
            ylabs.append(ylab)

        # and x-axis labels
        xlab = [
            "3'-" + m[-1] if i == 0 else m[-1]
            for i, m in enumerate(muts_out[-1])
        ]

        ax.set_xticklabels(xlab)
        ax.set_yticklabels(ylabs, rotation=0)

        ax.set_title(title, fontsize=12)
        f.tight_layout()
        f.savefig(outname, bbox_inches='tight')