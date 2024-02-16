from ._core import coord2percent, percent2coord, make_annotation
from matplotlib.patches import Rectangle
import pandas as pd


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
    """Plot the direction arrow in the given item if it proceeds."""

    dir_flag = 0

    if strand:
        # create and plot direction lines
        if item_size > item_threshold:
            dir_flag = 1
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

    return dir_flag


def _apply_gene_bridge(
    transcript_str,
    df,
    fig,
    ax,
    strand,
    gene_ix,
    exon_color,
    tag_background,
    genename,
    showinfo,
    exon_width,
    transcript_utr_width,
    arrow_size_min,
    arrow_color,
    arrow_style,
    arrow_width,
    dir_flag,
    arrow_size,
):
    """Evaluate data and provide _plot_row with right parameters."""
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
                genename,
                showinfo,
                exon_width,
                arrow_size_min,
                arrow_color,
                arrow_style,
                arrow_width,
                dir_flag,
                arrow_size,
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
            tr_start, cds_start = df.groupby("Feature", group_keys=False).Start.apply(
                min
            )[["exon", "CDS"]]
            tr_end, cds_end = df.groupby("Feature", group_keys=False).End.apply(max)[
                ["exon", "CDS"]
            ]

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
            # get the gene information to print on hover
            # default
            if strand:
                geneinfo_start = f"[{strand}] ({tr_start}, {cds_start})\nID: {genename}"  # default with strand
                geneinfo_end = f"[{strand}] ({cds_end}, {tr_end})\nID: {genename}"  # default with strand
            else:
                geneinfo_start = f"({tr_start}, {cds_start})\nID: {genename}"  # default without strand
                geneinfo_end = (
                    f"({cds_end}, {tr_end})\nID: {genename}"  # default without strand
                )

            # customized
            # showinfo_dict = row.to_dict()  # first element of gene rows
            # if showinfo:
            #     geneinfo += "\n" + showinfo.format(**showinfo_dict)

            make_annotation(start_utr, fig, ax, geneinfo_start, tag_background)
            make_annotation(end_utr, fig, ax, geneinfo_end, tag_background)

            # keep CDS data and plot it
            df = df.groupby("Feature", group_keys=False).get_group("CDS")
            df.apply(
                _plot_row,
                args=(
                    fig,
                    ax,
                    strand,
                    gene_ix,
                    exon_color,
                    tag_background,
                    genename,
                    showinfo,
                    exon_width,
                    arrow_size_min,
                    arrow_color,
                    arrow_style,
                    arrow_width,
                    dir_flag,
                    arrow_size,
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
                    genename,
                    showinfo,
                    transcript_utr_width,
                    arrow_size_min,
                    arrow_color,
                    arrow_style,
                    arrow_width,
                    dir_flag,
                    arrow_size,
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
                    genename,
                    showinfo,
                    exon_width,
                    arrow_size_min,
                    arrow_color,
                    arrow_style,
                    arrow_width,
                    dir_flag,
                    arrow_size,
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
    genename,
    showinfo,
    exon_width,
    arrow_size_min,
    arrow_color,
    arrow_style,
    arrow_width,
    dir_flag,
    arrow_size,
):
    """Plot elements corresponding to one row of one gene."""

    # Make gene annotation
    # get the gene information to print on hover
    # default
    if strand:
        geneinfo = f"[{strand}] ({row.oriStart}, {row.oriEnd})\nID: {genename}"  # default with strand
    else:
        geneinfo = (
            f"({row.oriStart}, {row.oriEnd})\nID: {genename}"  # default without strand
        )

    # customized
    showinfo_dict = row.to_dict()  # first element of gene rows
    if showinfo:
        geneinfo += "\n" + showinfo.format(**showinfo_dict)

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
    arrow_size = coord2percent(ax, 0.05 * start, 0.05 * stop)
    incl = percent2coord(ax, arrow_size / 2)  # how long in the plot (OX)

    # create and plot lines
    if not dir_flag:
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


def plot_introns(
    sorted_exons,
    ts_chrom,
    fig,
    ax,
    geneinfo,
    tag_background,
    gene_ix,
    exon_color,
    strand,
    intron_threshold,
    exon_width,
    arrow_color,
    arrow_style,
    arrow_width,
    arrow_size,
):
    """Plot intron lines as needed."""

    dir_flag = []

    for i in range(len(sorted_exons) - 1):
        # define intron
        start = sorted_exons["End"].iloc[i]
        stop = sorted_exons["Start"].iloc[i + 1]

        # NOT introns off
        if ts_chrom.empty:
            ts_intron = pd.DataFrame()

        # INTRONS OFF
        else:
            ts_intron = ts_chrom[
                (ts_chrom["Start_adj"] >= start) & (ts_chrom["Start_adj"] < stop)
            ].reset_index()

        # Plot LINES binding exons
        # No to-shrink regions in intron
        if ts_intron.empty:
            # create continuous line
            intron_line = ax.plot(
                [start, stop],
                [gene_ix, gene_ix],
                color=exon_color,
                linewidth=1,
                zorder=1,
            )
            # add to plot with annotation
            make_annotation(intron_line[0], fig, ax, geneinfo, tag_background)

        # Intron has to-shrink regions
        else:
            # alternating fixed and to-shrink regions or vice versa
            prev_tsend = None

            # Iterate over to-shrink regions
            for ix, row in ts_intron.iterrows():
                # (1) Add previous fixed region if needed
                # consider intron starts with fixed region
                if not prev_tsend and row["Start_adj"] != start:
                    prev_tsend = start

                # create continuous line
                intron_line = ax.plot(
                    [prev_tsend, row["Start_adj"]],
                    [gene_ix, gene_ix],
                    color=exon_color,
                    linewidth=1,
                    zorder=1,
                )
                # add to plot with annotation
                make_annotation(intron_line[0], fig, ax, geneinfo, tag_background)

                # (2) Add to-shrink region
                intron_line = ax.plot(
                    [row["Start_adj"], row["End_adj"]],
                    [gene_ix, gene_ix],
                    color=exon_color,
                    linewidth=0.5,
                    linestyle="--",
                    zorder=1,
                )
                # add to plot with annotation
                make_annotation(intron_line[0], fig, ax, geneinfo, tag_background)

                # (3) Add final fixed region if needed
                if (ix == len(ts_intron) - 1) and (row["End_adj"] != stop):
                    # add last fixed region
                    # create continuous line
                    intron_line = ax.plot(
                        [row["End_adj"], stop],
                        [gene_ix, gene_ix],
                        color=exon_color,
                        linewidth=1,
                        zorder=1,
                    )
                    # add to plot with annotation
                    make_annotation(intron_line[0], fig, ax, geneinfo, tag_background)

                # store interval end for next iteration
                prev_tsend = row["End_adj"]

        intron_size = coord2percent(ax, start, stop)
        incl = percent2coord(
            ax, arrow_size / 2
        )  # how long is the arrow in the plot (OX)

        # Plot DIRECTION ARROW in INTRONS if strand is known
        dir_flag.append(
            plot_direction(
                ax,
                strand,
                intron_size,
                intron_threshold,
                start,
                stop,
                incl,
                gene_ix,
                exon_width,
                arrow_color,
                arrow_style,
                arrow_width,
            )
        )

    if 1 in dir_flag:
        return 1
    else:
        return 0
