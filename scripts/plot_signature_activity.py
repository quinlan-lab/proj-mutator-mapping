import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

plt.rc('font', size=20)

def main(args):

    palette = dict(
            zip(
                ["B-B", "B-D", "D-B", "D-D"],
                ["#398D84", "#E67F3A", "#EBBC2C", "#2F294A"],
            ))
    
    sig_dict = {"SBS1": r"$\bf{SBS1}$" + "\n(Spontaneous deamination of 5-mC)",
                "SBS5": r"$\bf{SBS5}$" + "\n(Unknown clock-like)",
                "SBS18": r"$\bf{SBS18}$" + "\n(Damage by ROS)",
                "SBS30": r"$\bf{SBS30}$" + "\n(Defective BER due to NTHL1)"}

    mutations = pd.read_csv(args.spectra)[["sample", "Haplotypes", "Generations"]]
    activities = pd.read_csv(args.activities, sep="\t")

    # convert activities to a tidy dataframe
    activities_tidy = activities.melt(id_vars="Samples", var_name="Signature", value_name="Count")
    activities_tidy = activities_tidy#.query('Signature == "SBS18"')
    activities_tidy["Signature"] = activities_tidy["Signature"].apply(lambda s: sig_dict[s])
    # figure out the total sum of mutations in each sample that was assigned
    # to a specific SBS signatures
    activity_sums = np.sum(activities.values[:, 1:], axis=1)
    smp2total = dict(zip(activities["Samples"], activity_sums))
    activities_tidy["Total"] = activities_tidy["Samples"].apply(lambda s: smp2total[s])
    activities_tidy["Fraction"] = activities_tidy["Count"] / activities_tidy["Total"]
    # annotate each sample's SBS activities with relevant metadata
    annotated_activities = activities_tidy.merge(mutations, left_on="Samples", right_on="sample").sort_values("Haplotypes")
    f, ax = plt.subplots(figsize=(18, 6))
    sns.barplot(
        data=annotated_activities,
        y="Signature",
        x="Fraction",
        hue="Haplotypes",
        palette=palette,
        capsize=0.05,
        errcolor="k",
        ec="k",
        linewidth=2,
        ax=ax,
    )
    sns.despine(ax=ax, top=True, right=True)
    ax.set_xlabel("Fraction of mutations attributed to signature")
    ax.set_ylabel(None)
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(
        handles,
        labels,
        title="Genotypes at chr4 and chr6 peaks",
        frameon=False,
    )

    # change all spine widths
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)

    # increase tick width
    ax.tick_params(width=1.5)
    f.tight_layout()

    f.savefig(args.out, dpi=300)




if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--spectra")
    p.add_argument("--activities")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
