import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter


def create_fig(
    x,
    y,
    chrmd_df,
    genesmd_df,
    chr_string,
    title_dict_plt,
    plot_background,
    plot_border,
    packed,
    legend,
):
    """Generate the figure and axes fitting the data."""

    fig = plt.figure(figsize=(x, y))
    gs = gridspec.GridSpec(
        len(chrmd_df), 1, height_ratios=chrmd_df.y_height
    )  # size of chromosome subplot according to number of gene rows

    # one plot per chromosome
    axes = []
    for i in range(len(chrmd_df)):
        chrom = chrmd_df.index[i]
        axes.append(plt.subplot(gs[i]))
        ax = axes[i]
        # Adjust plot display
        ax.set_title(chr_string.format(**locals()), fontdict=title_dict_plt)
        ax.set_facecolor(plot_background)
        plt.setp(ax.spines.values(), color=plot_border)
        plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=plot_border)
        ax.xaxis.set_tick_params(bottom=False)
        ax.yaxis.set_tick_params(left=False)

        # set x axis limits
        x_min, x_max = chrmd_df.iloc[i]["min_max"]
        x_rang = x_max - x_min
        ax.set_xlim(
            x_min - 0.05 * x_rang, x_max + 0.05 * x_rang
        )  # add 5% to limit coordinates range
        plt.ticklabel_format(style="plain")
        ax.grid(visible=True, axis="x", linestyle=":")  # , zorder = -1)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.get_major_formatter().set_scientific(False)  # not scientific notation
        ax.xaxis.get_major_formatter().set_useOffset(False)  # not offset notation

        # set y axis limits
        y_min = 0
        y_max = chrmd_df.iloc[i].y_height
        ax.set_ylim(y_min, y_max)
        # gene name as y labels
        y_ticks_val = []
        y_ticks_name = []
        if not packed:
            y_ticks_val = [i + 0.5 for i in range(int(y_max))]
            y_ticks_name = genesmd_df.groupby("chrix", group_keys=False).groups[chrom]
        ax.set_yticks(y_ticks_val)
        ax.set_yticklabels(y_ticks_name)

    plt.subplots_adjust(hspace=0.7)
    # Create legend
    if legend:
        handles = genesmd_df["legend_item"].tolist()
        labels = genesmd_df.index.tolist()
        fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(1, 1))

    return fig, axes
