from .core import coord2percent, percent2coord
import plotly.subplots as sp
import plotly.graph_objects as go


####### CREATE FIGURE AND AXES
def create_fig(chrmd_df, genesmd_df, chr_string, title_dict_ply, packed):
    """xxx"""

    titles = [chr_string.format(**locals()) for chrom in chrmd_df.index]
    fig = sp.make_subplots(
        rows=len(chrmd_df),
        cols=1,
        row_heights=chrmd_df.y_height.to_list(),
        subplot_titles=titles,
    )

    # one subplot per chromosome
    for i in range(len(chrmd_df)):
        chrom = chrmd_df.index[i]
        fig.add_trace(go.Scatter(x=[], y=[]), row=i + 1, col=1)

        # set title format
        fig.layout.annotations[i].update(font=title_dict_ply)

        # set x axis limits
        x_min, x_max = chrmd_df.iloc[i]["min_max"]
        x_rang = x_max - x_min
        fig.update_xaxes(
            range=[x_min - 0.05 * x_rang, x_max + 0.05 * x_rang],
            tickformat="d",
            showgrid=False,
            zeroline=False,
            row=i + 1,
            col=1,
        )  # add 5% to limit coordinates range
        # fig.update_xaxes(title_text="Chromosome %s" % chrom, row=i + 1, col=1)

        # set y axis limits
        y_min = 0
        y_max = chrmd_df.iloc[i].y_height
        y_ticks_val = []
        y_ticks_name = []
        if not packed:
            y_ticks_val = [i + 0.5 for i in range(int(y_max))]
            y_ticks_name = genesmd_df.groupby(genesmd_df["chrix"]).groups[chrom]
        fig.update_yaxes(
            range=[y_min, y_max],
            tickvals=y_ticks_val,
            ticktext=y_ticks_name,
            showgrid=False,
            row=i + 1,
            col=1,
        )

    return fig


###### PLOTTING ELEEMNTS


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
                arrow_bot = go.Scatter(
                    x=bot_plus[0],
                    y=bot_plus[1],
                    mode="lines",
                    line=go.scatter.Line(color=arrow_color, width=1),
                    showlegend=False,
                    name=genename,
                    hoverinfo="skip",
                )
                arrow_top = go.Scatter(
                    x=top_plus[0],
                    y=top_plus[1],
                    mode="lines",
                    line=go.scatter.Line(color=arrow_color, width=1),
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
                    line=go.scatter.Line(color=arrow_color, width=1),
                    showlegend=False,
                    name=genename,
                    hoverinfo="skip",
                )
                arrow_top = go.Scatter(
                    x=top_minus[0],
                    y=top_minus[1],
                    mode="lines",
                    line=go.scatter.Line(color=arrow_color, width=1),
                    showlegend=False,
                    name=genename,
                    hoverinfo="skip",
                )
                fig.add_trace(arrow_bot, row=chrom_ix + 1, col=1)
                fig.add_trace(arrow_top, row=chrom_ix + 1, col=1)


def _apply_gene(
    transcript_str,
    df,
    fig,
    strand,
    genename,
    gene_ix,
    exon_color,
    chrom_ix,
    geneinfo,
    exon_width,
    transcript_utr_width,
    legend,
    arrow_size_min,
    arrow_color,
):
    """xxx"""

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
                geneinfo,
                exon_width,
                legend,
                arrow_size_min,
                arrow_color,
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
            utr_start, cds_start = df.groupby("Feature").Start.apply(min)[
                ["exon", "CDS"]
            ]
            utr_end, cds_end = df.groupby("Feature").End.apply(max)[["exon", "CDS"]]

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
            df = df.groupby("Feature").get_group("CDS")
            df.apply(
                _plot_row,
                args=(
                    fig,
                    strand,
                    genename,
                    gene_ix,
                    exon_color,
                    chrom_ix,
                    geneinfo,
                    exon_width,
                    legend,
                    arrow_size_min,
                    arrow_color,
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
                    geneinfo,
                    exon_width,
                    legend,
                    arrow_size_min,
                    arrow_color,
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
                    geneinfo,
                    exon_width,
                    legend,
                    arrow_size_min,
                    arrow_color,
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
    geneinfo,
    exon_width,
    legend,
    arrow_size_min,
    arrow_color,
):
    """Plot elements corresponding to one row of one gene."""

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
            name=genename,
            showlegend=legend,
        ),
        row=chrom_ix + 1,
        col=1,
    )

    # Plot DIRECTION ARROW in EXON
    # decide about placing a direction arrow
    arrow_size = coord2percent(fig, chrom_ix + 1, 0.05 * start, 0.05 * stop)
    incl = percent2coord(fig, chrom_ix + 1, 0.003)  # how long in the plot (OX)

    # create and plot lines
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
    )
