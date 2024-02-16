from ._core import coord2percent, percent2coord
import plotly.graph_objects as go
import pandas as pd


def plot_direction(
    fig,
    strand,
    genename,
    item_size,
    item_threshold,
    start,
    stop,
    incl,
    gene_ix,
    chrom_ix,
    exon_width,
    arrow_color,
    arrow_line_width,
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
                arrow_bot = go.Scatter(
                    x=bot_plus[0],
                    y=bot_plus[1],
                    mode="lines",
                    line=go.scatter.Line(color=arrow_color, width=arrow_line_width),
                    showlegend=False,
                    name=genename,
                    hoverinfo="skip",
                )
                arrow_top = go.Scatter(
                    x=top_plus[0],
                    y=top_plus[1],
                    mode="lines",
                    line=go.scatter.Line(color=arrow_color, width=arrow_line_width),
                    showlegend=False,
                    name=genename,
                    hoverinfo="skip",
                )
                fig.add_trace(arrow_bot, row=chrom_ix + 1, col=1)
                fig.add_trace(arrow_top, row=chrom_ix + 1, col=1)

            elif strand == "-":
                arrow_bot = go.Scatter(
                    x=bot_minus[0],
                    y=bot_minus[1],
                    mode="lines",
                    line=go.scatter.Line(color=arrow_color, width=arrow_line_width),
                    showlegend=False,
                    name=genename,
                    hoverinfo="skip",
                )
                arrow_top = go.Scatter(
                    x=top_minus[0],
                    y=top_minus[1],
                    mode="lines",
                    line=go.scatter.Line(color=arrow_color, width=arrow_line_width),
                    showlegend=False,
                    name=genename,
                    hoverinfo="skip",
                )
                fig.add_trace(arrow_bot, row=chrom_ix + 1, col=1)
                fig.add_trace(arrow_top, row=chrom_ix + 1, col=1)
    return dir_flag


def _apply_gene_bridge(
    transcript_str,
    df,
    fig,
    strand,
    genename,
    gene_ix,
    exon_color,
    chrom_ix,
    geneinfo,
    showinfo,
    exon_width,
    transcript_utr_width,
    legend,
    arrow_size_min,
    arrow_color,
    arrow_line_width,
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
                strand,
                genename,
                gene_ix,
                exon_color,
                chrom_ix,
                showinfo,
                exon_width,
                legend,
                arrow_size_min,
                arrow_color,
                arrow_line_width,
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
            utr_start, cds_start = df.groupby("Feature", group_keys=False).Start.apply(
                min
            )[["exon", "CDS"]]
            utr_end, cds_end = df.groupby("Feature", group_keys=False).End.apply(max)[
                ["exon", "CDS"]
            ]

            # create start utr
            x0, x1 = utr_start, cds_start
            y0, y1 = (
                gene_ix - transcript_utr_width / 2,
                gene_ix + transcript_utr_width / 2,
            )
            fig.add_trace(
                go.Scatter(
                    x=[x0, x1, x1, x0, x0],
                    y=[y0, y0, y1, y1, y0],
                    fill="toself",
                    fillcolor=exon_color,
                    mode="lines",
                    line=dict(color=exon_color),
                    text=geneinfo,
                    hoverinfo="text",
                    name=genename,
                    showlegend=legend,
                ),
                row=chrom_ix + 1,
                col=1,
            )

            # create end utr
            x0, x1 = cds_end, utr_end
            y0, y1 = (
                gene_ix - transcript_utr_width / 2,
                gene_ix + transcript_utr_width / 2,
            )
            fig.add_trace(
                go.Scatter(
                    x=[x0, x1, x1, x0, x0],
                    y=[y0, y0, y1, y1, y0],
                    fill="toself",
                    fillcolor=exon_color,
                    mode="lines",
                    line=dict(color=exon_color),
                    text=geneinfo,
                    hoverinfo="text",
                    name=genename,
                    showlegend=legend,
                ),
                row=chrom_ix + 1,
                col=1,
            )

            # keep CDS data and plot it
            df = df.groupby("Feature", group_keys=False).get_group("CDS")
            df.apply(
                _plot_row,
                args=(
                    fig,
                    strand,
                    genename,
                    gene_ix,
                    exon_color,
                    chrom_ix,
                    showinfo,
                    exon_width,
                    legend,
                    arrow_size_min,
                    arrow_color,
                    arrow_line_width,
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
                    strand,
                    genename,
                    gene_ix,
                    exon_color,
                    chrom_ix,
                    showinfo,
                    exon_width,
                    legend,
                    arrow_size_min,
                    arrow_color,
                    arrow_line_width,
                    dir_flag,
                    arrow_size,
                ),
                axis=1,
            )

        # trancript only has exon
        elif (
            not df.Feature.str.contains("CDS").any()
            and df.Feature.str.contains("exon").any()
        ):
            # plot just as utr
            df.apply(
                _plot_row,
                args=(
                    fig,
                    strand,
                    genename,
                    gene_ix,
                    exon_color,
                    chrom_ix,
                    showinfo,
                    transcript_utr_width,
                    legend,
                    arrow_size_min,
                    arrow_color,
                    arrow_line_width,
                    dir_flag,
                    arrow_size,
                ),
                axis=1,
            )

        # transcript has neither, skip it
        else:
            return


def _plot_row(
    row,
    fig,
    strand,
    genename,
    gene_ix,
    exon_color,
    chrom_ix,
    showinfo,
    exon_width,
    legend,
    arrow_size_min,
    arrow_color,
    arrow_line_width,
    dir_flag,
    arrow_size,
):
    """Plot elements corresponding to one row of one gene."""

    # Get the gene information to print on hover
    # default
    if strand:
        geneinfo = f"[{strand}] ({row.oriStart}, {row.oriEnd})<br>ID: {genename}"  # default with strand
    else:
        geneinfo = f"({row.oriStart}, {row.oriEnd})<br>ID: {genename}"  # default without strand

    # customized
    showinfo_dict = row.to_dict()  # first element of gene rows
    if showinfo:
        showinfo = showinfo.replace("\n", "<br>")
        geneinfo += "<br>" + showinfo.format(**showinfo_dict)

    # consider legend
    if legend:
        legend = bool(row["legend_tag"])

    # Exon start and stop
    start = int(row["Start"])
    stop = int(row["End"])
    # convert to coordinates for rectangle
    x0, x1 = start, stop
    y0, y1 = (
        gene_ix - exon_width / 2,
        gene_ix + exon_width / 2,
    )  ##gene middle point -+ half of exon size

    # Plot EXON as rectangle
    fig.add_trace(
        go.Scatter(
            x=[x0, x1, x1, x0, x0],
            y=[y0, y0, y1, y1, y0],
            fill="toself",
            fillcolor=exon_color,
            mode="lines",
            line=dict(color=exon_color),
            text=geneinfo,
            hoverinfo="text",
            name=row["legend_tag"],
            showlegend=legend,
        ),
        row=chrom_ix + 1,
        col=1,
    )

    # Plot DIRECTION ARROW in EXON
    # decide about placing a direction arrow
    arrow_size = coord2percent(fig, chrom_ix + 1, 0.05 * start, 0.05 * stop)
    incl = percent2coord(fig, chrom_ix + 1, arrow_size / 2)  # how long in the plot (OX)

    # create and plot lines
    if not dir_flag:
        plot_direction(
            fig,
            strand,
            genename,
            arrow_size,
            arrow_size_min,
            start,
            stop,
            incl,
            gene_ix,
            chrom_ix,
            exon_width,
            arrow_color,
            arrow_line_width,
        )


def plot_introns(
    sorted_exons,
    ts_chrom,
    fig,
    gene_ix,
    exon_color,
    chrom_ix,
    strand,
    genename,
    intron_threshold,
    exon_width,
    arrow_color,
    arrow_line_width,
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
            x0, x1 = start, stop
            y0, y1 = gene_ix, gene_ix
            intron_line = go.Scatter(
                x=[x0, x1],
                y=[y0, y1],
                mode="lines",
                line=dict(color=exon_color, width=0.7, dash="solid"),
                hoverinfo="skip",
                showlegend=False,
            )
            fig.add_trace(intron_line, row=chrom_ix + 1, col=1)

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
                x0, x1 = prev_tsend, row["Start_adj"]
                y0, y1 = gene_ix, gene_ix
                intron_line = go.Scatter(
                    x=[x0, x1],
                    y=[y0, y1],
                    mode="lines",
                    line=dict(color=exon_color, width=0.7, dash="solid"),
                    hoverinfo="skip",
                    showlegend=False,
                )
                fig.add_trace(intron_line, row=chrom_ix + 1, col=1)

                # (2) Add to-shrink region
                x0, x1 = row["Start_adj"], row["End_adj"]
                y0, y1 = gene_ix, gene_ix
                intron_line = go.Scatter(
                    x=[x0, x1],
                    y=[y0, y1],
                    mode="lines",
                    line=dict(color=exon_color, width=0.7, dash="dash"),
                    hoverinfo="skip",
                    showlegend=False,
                )
                fig.add_trace(intron_line, row=chrom_ix + 1, col=1)

                # (3) Add final fixed region if needed
                if (ix == len(ts_intron) - 1) and (row["End_adj"] != stop):
                    # add last fixed region
                    # create continuous line
                    x0, x1 = row["End_adj"], stop
                    y0, y1 = gene_ix, gene_ix
                    intron_line = go.Scatter(
                        x=[x0, x1],
                        y=[y0, y1],
                        mode="lines",
                        line=dict(color=exon_color, width=0.7, dash="solid"),
                        hoverinfo="skip",
                        showlegend=False,
                    )
                    fig.add_trace(intron_line, row=chrom_ix + 1, col=1)

                # store interval end for next iteration
                prev_tsend = row["End_adj"]

        # Plot DIRECTION ARROW in INTRONS if strand is known
        intron_size = coord2percent(fig, chrom_ix + 1, start, stop)
        incl = percent2coord(
            fig, chrom_ix + 1, arrow_size / 2
        )  # how long in the plot (OX)

        dir_flag.append(
            plot_direction(
                fig,
                strand,
                genename,
                intron_size,
                intron_threshold,
                start,
                stop,
                incl,
                gene_ix,
                chrom_ix,
                exon_width,
                arrow_color,
                arrow_line_width,
            )
        )

    if 1 in dir_flag:
        return 1
    else:
        return 0
