import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
from matplotlib.patches import Rectangle
from pyranges_plot.core import cumdelting
from ._core import make_annotation


def create_fig(
    x,
    y,
    chrmd_df,
    genesmd_df,
    ts_data,
    title_chr,
    title_dict_plt,
    plot_background,
    plot_border,
    packed,
    legend,
    tick_pos_d,
    ori_tick_pos_d,
    tag_background,
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
        ax.set_title(title_chr.format(**locals()), fontdict=title_dict_plt)
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
            # get previous default ticks
            original_ticks = [
                int(tick.get_text().replace("−", "-")) for tick in ax.get_xticklabels()
            ][1:]

            jump = original_ticks[1] - original_ticks[0]

            # find previous ticks that should be conserved
            to_add_val = []
            # there is data to shrink
            if ori_tick_pos_d[chrom]:
                to_add_val += [
                    tick
                    for tick in original_ticks
                    if tick < min(ori_tick_pos_d[chrom])
                    or tick > max(ori_tick_pos_d[chrom])
                ]
                for ii in range(1, len(ori_tick_pos_d[chrom]) - 1, 2):
                    not_shr0 = ori_tick_pos_d[chrom][ii]
                    not_shr1 = ori_tick_pos_d[chrom][ii + 1]
                    to_add_val += [
                        i for i in original_ticks if not_shr0 < i <= not_shr1
                    ]

            # nothing to shrink
            else:
                to_add_val += original_ticks

            # compute new coordinates of conserved previous ticks
            to_add = to_add_val.copy()
            to_add = cumdelting(to_add, ts_data, chrom)

            # set new ticks
            # x_ticks_val = sorted(tick_pos_d[chrom] + to_add)
            x_ticks_val = sorted(to_add)
            # do not add ticks beyond adjusted limits
            x_ticks_val = [
                num for num in x_ticks_val if num <= chrmd_df.loc[chrom]["min_max"][1]
            ]
            # x_ticks_name = sorted(ori_tick_pos_d[chrom] + to_add_val)[
            #     : len(x_ticks_val)
            # ]
            x_ticks_name = to_add_val[: len(x_ticks_val)]

            # adjust names
            ax.set_xticks(x_ticks_val)
            ax.set_xticklabels(x_ticks_name)

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

            for a, b, c, d in zip(
                rects_df["Start"],
                rects_df["End"],
                rects_df["cumdelta_start"],
                rects_df["cumdelta_end"],
            ):
                ts_range = Rectangle(
                    (a, y_min - 1),
                    b - a,
                    y_max + 1,
                    edgecolor="grey",
                    facecolor="lightyellow",
                    alpha=0.7,
                    fill=True,
                    linewidth=0,  # option2
                )
                ax.add_patch(ts_range)
                make_annotation(
                    ts_range,
                    fig,
                    ax,
                    f"Shrinked region:\n[{a+c} - {b+d}]",
                    tag_background,
                )

    plt.subplots_adjust(hspace=0.7)
    # Create legend
    if legend:
        handles = genesmd_df["legend_item"].tolist()
        labels = genesmd_df.index.tolist()
        fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(1, 1))

    return fig, axes
