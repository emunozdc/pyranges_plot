from .core import coord2inches, inches2coord
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter


############ FIGURE AND AXES


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
    """xxx"""

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
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.get_major_formatter().set_scientific(
            False
        )  # turn off scientific notation
        ax.xaxis.get_major_formatter().set_useOffset(False)  # turn off offset notation

        # set y axis limits
        y_min = 0
        y_max = chrmd_df.iloc[i].y_height
        ax.set_ylim(y_min, y_max)
        # gene name as y labels
        y_ticks_val = []
        y_ticks_name = []
        if not packed:
            y_ticks_val = [i + 0.5 for i in range(int(y_max))]
            y_ticks_name = genesmd_df.groupby(genesmd_df["chrix"]).groups[chrom]
        ax.set_yticks(y_ticks_val)
        ax.set_yticklabels(y_ticks_name)

    plt.subplots_adjust(hspace=0.7)
    # Create legend
    if legend:
        handles = genesmd_df["legend_item"].tolist()
        labels = genesmd_df.index.tolist()
        fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(1, 1))

    return fig, axes


############ PLOTTING ITEMS


def make_annotation(item, fig, ax, geneinfo, tag_background):
    """xxx"""

    # create annotation and make it not visible
    annotation = ax.annotate(
        "",
        xy=(0, 0),
        xytext=(20, 20),
        textcoords="offset points",
        bbox=dict(
            boxstyle="round",
            edgecolor=tag_background,
            facecolor=tag_background,
        ),
        arrowprops=dict(arrowstyle="->"),
        color="white",
    )
    annotation.set_visible(False)

    # make annotation visible when over the gene line
    def on_hover(event):
        visible = annotation.get_visible()
        contains_item, _ = item.contains(event)  # Check if mouse is over the gene line
        if contains_item:
            annotation.set_text(geneinfo)
            annotation.xy = (event.xdata, event.ydata)
            annotation.set_visible(True)
            fig.canvas.draw()
        elif visible:
            annotation.set_visible(False)
            fig.canvas.draw()

    fig.canvas.mpl_connect("motion_notify_event", on_hover)


def plot_direction(
    ax,
    strand,
    item_size,
    item_threshold,
    start,
    stop,
    incl,
    gene_ix,
    exon_width,
    arrow_color,
    arrow_style,
    arrow_width,
):
    """xxx"""

    if strand:
        # create and plot direction lines
        if item_size > item_threshold:
            ##diagonal_line = OX arrow extension(item middle point +- incl), OY arrow extension (item middle point + half of exon width)
            top_plus = (
                [(start + stop) / 2 + incl, (start + stop) / 2 - incl],
                [gene_ix, gene_ix + exon_width / 2 - 0.01],
            )
            bot_plus = (
                [(start + stop) / 2 - incl, (start + stop) / 2 + incl],
                [gene_ix - exon_width / 2 + 0.01, gene_ix],
            )
            top_minus = (
                [(start + stop) / 2 + incl, (start + stop) / 2 - incl],
                [gene_ix - exon_width / 2 + 0.01, gene_ix],
            )
            bot_minus = (
                [(start + stop) / 2 - incl, (start + stop) / 2 + incl],
                [gene_ix, gene_ix + exon_width / 2 - 0.01],
            )

            if strand == "+":
                ax.plot(
                    bot_plus[0],
                    bot_plus[1],
                    color=arrow_color,
                    linewidth=arrow_width,
                    solid_capstyle=arrow_style,
                )

                ax.plot(
                    top_plus[0],
                    top_plus[1],
                    color=arrow_color,
                    linewidth=arrow_width,
                    solid_capstyle=arrow_style,
                )

            elif strand == "-":
                ax.plot(
                    bot_minus[0],
                    bot_minus[1],
                    color=arrow_color,
                    linewidth=arrow_width,
                    solid_capstyle=arrow_style,
                )

                ax.plot(
                    top_minus[0],
                    top_minus[1],
                    color=arrow_color,
                    linewidth=arrow_width,
                    solid_capstyle=arrow_style,
                )


def _apply_gene(
    transcript_str,
    df,
    fig,
    ax,
    strand,
    gene_ix,
    exon_color,
    tag_background,
    geneinfo,
    exon_width,
    transcript_utr_width,
    arrow_size_min,
    arrow_color,
    arrow_style,
    arrow_width,
):
    """xxx"""

    # NOT transcript strucutre
    if not transcript_str:
        df.apply(
            _plot_row,
            args=(
                fig,
                ax,
                strand,
                gene_ix,
                exon_color,
                tag_background,
                geneinfo,
                exon_width,
                arrow_size_min,
                arrow_color,
                arrow_style,
                arrow_width,
            ),
            axis=1,
        )

    # WITH transcript structure
    else:
        # transcript has only CDS and exon
        if (
            df.Feature.str.contains("CDS").any()
            and df.Feature.str.contains("exon").any()
        ):
            # get coordinates for utr and cds
            tr_start, cds_start = df.groupby("Feature").Start.apply(min)[
                ["exon", "CDS"]
            ]
            tr_end, cds_end = df.groupby("Feature").End.apply(max)[["exon", "CDS"]]

            # create utr
            start_utr = Rectangle(
                (tr_start, gene_ix - transcript_utr_width / 2),
                cds_start - tr_start,
                transcript_utr_width,
                edgecolor=exon_color,
                facecolor=exon_color,
                fill=True,
            )
            end_utr = Rectangle(
                (cds_end, gene_ix - transcript_utr_width / 2),
                tr_end - cds_end,
                transcript_utr_width,
                edgecolor=exon_color,
                facecolor=exon_color,
                fill=True,
            )
            ax.add_patch(start_utr)
            ax.add_patch(end_utr)

            # make annotation for utr
            make_annotation(start_utr, fig, ax, geneinfo, tag_background)
            make_annotation(end_utr, fig, ax, geneinfo, tag_background)

            # keep CDS data and plot it
            df = df.groupby("Feature").get_group("CDS")
            df.apply(
                _plot_row,
                args=(
                    fig,
                    ax,
                    strand,
                    gene_ix,
                    exon_color,
                    tag_background,
                    geneinfo,
                    exon_width,
                    arrow_size_min,
                    arrow_color,
                    arrow_style,
                    arrow_width,
                ),
                axis=1,
            )

        # transcript only has exon
        elif (
            not df.Feature.str.contains("CDS").any()
            and df.Feature.str.contains("exon").any()
        ):
            # plot just as utr
            df.apply(
                _plot_row,
                args=(
                    fig,
                    ax,
                    strand,
                    gene_ix,
                    exon_color,
                    tag_background,
                    geneinfo,
                    transcript_utr_width,
                    arrow_size_min,
                    arrow_color,
                    arrow_style,
                    arrow_width,
                ),
                axis=1,
            )

        # transcript only has CDS
        elif (
            df.Feature.str.contains("CDS").any()
            and not df.Feature.str.contains("exon").any()
        ):
            df.apply(
                _plot_row,
                args=(
                    fig,
                    ax,
                    strand,
                    gene_ix,
                    exon_color,
                    tag_background,
                    geneinfo,
                    exon_width,
                    arrow_size_min,
                    arrow_color,
                    arrow_style,
                    arrow_width,
                ),
                axis=1,
            )

        # transcript has neither, skip gene
        else:
            return


def _plot_row(
    row,
    fig,
    ax,
    strand,
    gene_ix,
    exon_color,
    tag_background,
    geneinfo,
    exon_width,
    arrow_size_min,
    arrow_color,
    arrow_style,
    arrow_width,
):
    """Plot elements corresponding to one row of one gene."""

    # Exon start and stop
    start = int(row["Start"])
    stop = int(row["End"])

    # Plot EXON as rectangle
    exon_rect = Rectangle(
        (start, gene_ix - exon_width / 2),
        stop - start,
        exon_width,
        edgecolor=exon_color,
        facecolor=exon_color,
        fill=True,
    )
    ax.add_patch(exon_rect)

    # create annotation for exon
    make_annotation(exon_rect, fig, ax, geneinfo, tag_background)

    # Plot DIRECTION ARROW in EXON
    # decide about placing a direction arrow
    arrow_size = coord2inches(fig, ax, 0.05 * start, 0.05 * stop, 0, 0)
    incl = inches2coord(fig, ax, 0.15)  # how long in the plot (OX)

    # create and plot lines
    plot_direction(
        ax,
        strand,
        arrow_size,
        arrow_size_min,
        start,
        stop,
        incl,
        gene_ix,
        exon_width,
        arrow_color,
        arrow_style,
        arrow_width,
    )
