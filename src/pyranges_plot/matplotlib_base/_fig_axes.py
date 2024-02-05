import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
from matplotlib.patches import Rectangle


def create_fig(
    x,
    y,
    chrmd_df,
    genesmd_df,
    ts_data,
    chr_string,
    title_dict_plt,
    plot_background,
    plot_border,
    packed,
    legend,
    tick_pos_d,
    ori_tick_pos_d,
    tick_val_d,
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
        # consider introns off
        if tick_pos_d:
            original_ticks = [
                int(tick.get_text().replace("−", "-")) for tick in ax.get_xticklabels()
            ]
            jump = original_ticks[1] - original_ticks[0]
            to_add = []
            for ii in range(1, len(tick_pos_d[chrom]) - 1, 2):
                not_shr0 = tick_pos_d[chrom][ii]
                not_shr1 = tick_pos_d[chrom][ii + 1]
                if (not_shr1 - not_shr0) > jump:
                    add_lines0 = jump - (not_shr0 % jump) + not_shr0
                    add_lines1 = jump - (not_shr1 % jump) + not_shr1
                    to_add += [num for num in range(add_lines0, add_lines1, jump)]

            to_add_val = to_add.copy()

            for itick in range(len(to_add)):
                for ix, row in ts_data[chrom].iterrows():
                    if row["End"] - row["cumdelta"] <= to_add[itick]:
                        cdel = row["cumdelta"]
                    else:
                        break
                to_add[itick] -= cdel
            # print(to_add)
            # print(to_add_val)
            # print(tick_pos_d[chrom] + to_add)
            # print(tick_val_d[chrom] + to_add_val)
            ax.set_xticks(sorted(tick_pos_d[chrom] + to_add_val))
            ax.set_xticklabels(sorted(tick_val_d[chrom] + to_add))

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

        # Add shrink rectangles
        if ts_data:
            rects_df = ts_data[chrom]
            rects_df["cumdelta_end"] = rects_df["cumdelta"]
            rects_df["cumdelta_start"] = rects_df["cumdelta"].shift(
                periods=1, fill_value=0
            )
            rects_df["Start"] -= rects_df["cumdelta_start"]
            rects_df["End"] -= rects_df["cumdelta_end"]

            for a, b in zip(rects_df["Start"], rects_df["End"]):
                ts_range = Rectangle(
                    (a, y_min - 1),
                    b - a,
                    y_max + 1,
                    edgecolor="grey",
                    facecolor="grey",
                    alpha=0.2,
                    # fill=True, # option1
                    fill=None,
                    hatch="///",
                    linewidth=0,  # option2
                )
                ax.add_patch(ts_range)

    plt.subplots_adjust(hspace=0.7)
    # Create legend
    if legend:
        handles = genesmd_df["legend_item"].tolist()
        labels = genesmd_df.index.tolist()
        fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(1, 1))

    return fig, axes
