import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
import numpy as np
from pyranges.core.names import CHROM_COL, START_COL, END_COL, STRAND_COL

from .core import initialize_dash_app, coord2percent
from .fig_axes import create_fig
from .data2plot import plot_introns, apply_gene_bridge
from ..names import PR_INDEX_COL, BORDER_COLOR_COL


def plot_exons_ply(
    subdf,
    feat_dict,
    genesmd_df,
    chrmd_df,
    chrmd_df_grouped,
    ts_data,
    id_col,
    max_shown=25,
    transcript_str=False,
    tooltip=None,
    legend=False,
    y_labels=False,
    text=True,
    title_chr=None,
    packed=True,
    to_file=None,
    file_size=None,
    warnings=None,
    tick_pos_d=None,
    ori_tick_pos_d=None,
):
    """Create Plotly plot."""

    # Get default plot features
    # tag_background = feat_dict['tag_background']
    fig_bkg = feat_dict["fig_bkg"]
    plot_bkg = feat_dict["plot_bkg"]
    plot_border = feat_dict["plot_border"]
    title_dict_ply = feat_dict["title_dict_ply"]
    grid_color = feat_dict["grid_color"]
    exon_border = feat_dict["exon_border"]
    exon_height = feat_dict["exon_height"]
    transcript_utr_width = feat_dict["transcript_utr_width"]
    v_spacer = feat_dict["v_spacer"]
    text_size = feat_dict["text_size"]
    plotly_port = feat_dict["plotly_port"]
    arrow_line_width = feat_dict["arrow_line_width"]
    arrow_color = feat_dict["arrow_color"]
    arrow_size_min = feat_dict["arrow_size_min"]
    arrow_size = feat_dict["arrow_size"]
    arrow_intron_threshold = feat_dict["arrow_intron_threshold"]
    shrinked_bkg = feat_dict["shrinked_bkg"]
    shrinked_alpha = feat_dict["shrinked_alpha"]

    # Create figure and chromosome plots
    fig = create_fig(
        subdf,
        chrmd_df,
        chrmd_df_grouped,
        genesmd_df,
        ts_data,
        title_chr,
        title_dict_ply,
        grid_color,
        packed,
        y_labels,
        tick_pos_d,
        ori_tick_pos_d,
        shrinked_bkg,
        shrinked_alpha,
        v_spacer,
        exon_height,
    )

    # Plot genes
    subdf.groupby(id_col + [PR_INDEX_COL], group_keys=False, observed=True).apply(
        lambda subdf: gby_plot_exons(
            subdf,
            fig,
            chrmd_df_grouped,
            genesmd_df,
            ts_data,
            tooltip,
            legend,
            transcript_str,
            text,
            text_size,
            exon_height,
            exon_border,
            transcript_utr_width,
            plot_bkg,
            arrow_line_width,
            arrow_color,
            arrow_size_min,
            arrow_size,
            arrow_intron_threshold,
        )
    )  # .reset_index(level=PR_INDEX_COL)

    # Adjust plot display
    fig.update_layout(
        plot_bgcolor=plot_bkg,
        font_color=plot_border,
        showlegend=legend,
        paper_bgcolor=fig_bkg,
    )
    fig.update_xaxes(
        showline=True,
        linewidth=1,
        linecolor=plot_border,
        mirror=True,
        color=plot_border,
    )
    fig.update_yaxes(
        showline=True,
        linewidth=1,
        linecolor=plot_border,
        mirror=True,
        color=plot_border,
    )

    # Provide output
    # insert silent information for warnings
    if warnings:
        fig.data[0].customdata = np.array([0, 0, 0])  # [tot_ngenes_l, 0, 0])
        if (
            "_blackwarning!" in genesmd_df.columns
            and "_iterwarning!" in genesmd_df.columns
        ):
            fig.data[0].customdata = np.array(
                [0, 91124, 91321]
            )  # [tot_ngenes_l, 91124, 91321])
        elif (
            "_blackwarning!" in genesmd_df.columns
            and "_iterwarning!" not in genesmd_df.columns
        ):
            fig.data[0].customdata = np.array(
                [0, 91124, 0]
            )  # [tot_ngenes_l, 91124, 0])
        elif (
            "_blackwarning!" not in genesmd_df.columns
            and "_iterwarning!" in genesmd_df.columns
        ):
            fig.data[0].customdata = np.array(
                [0, 0, 91321]
            )  # [tot_ngenes_l, 0, 91321])
    else:
        fig.data[0].customdata = np.array(["no warnings"])

    if to_file is None:
        app_instance = initialize_dash_app(fig, max_shown)
        app_instance.run(port=plotly_port)

    else:
        fig.update_layout(width=file_size[0], height=file_size[1])
        pio.write_image(fig, to_file)


def gby_plot_exons(
    df,
    fig,
    chrmd_df_grouped,
    genesmd_df,
    ts_data,
    showinfo,
    legend,
    transcript_str,
    text,
    text_size,
    exon_height,
    exon_border,
    transcript_utr_width,
    plot_background,
    arrow_line_width,
    arrow_color,
    arrow_size_min,
    arrow_size,
    arrow_intron_threshold,
):
    """Plot elements corresponding to the df rows of one gene."""

    # Gene parameters
    chrom = df[CHROM_COL].iloc[0]
    pr_ix = df[PR_INDEX_COL].iloc[0]
    genename = df["__id_col_2count__"].iloc[0]
    genemd = genesmd_df.loc[genename]  # store data for the gene
    df["legend_tag"] = [1] + [0] * (
        len(df) - 1
    )  # only one legend entry/linked intervals

    # in case same gene in +1 pr
    if not isinstance(genemd, pd.Series):
        genemd = genemd[genemd[PR_INDEX_COL] == pr_ix]  # in case same gene in +1 pr
        genemd = pd.Series(genemd.iloc[0])
    gene_ix = genemd["ycoord"] + 0.5
    # color of first interval will be used as intron color and utr color for simplicity
    if exon_border is None:
        exon_border = df[BORDER_COLOR_COL].iloc[0]

    chrom_ix = chrmd_df_grouped.loc[chrom]["chrom_ix"]

    if STRAND_COL in df.columns:
        strand = df[STRAND_COL].unique()[0]
    else:
        strand = ""

    # Get the gene information to print on hover
    # default
    if strand:
        geneinfo = f"[{strand}] ({min(df.__oriStart__)}, {max(df.__oriEnd__)})<br>ID: {genename}"  # default with strand
    else:
        geneinfo = f"({min(df.__oriStart__)}, {max(df.__oriEnd__)})<br>ID: {genename}"  # default without strand

    # add annotation for introns to plot
    x0, x1 = min(df[START_COL]), max(df[END_COL])
    y0, y1 = gene_ix - exon_height / 160, gene_ix + exon_height / 160

    fig.add_trace(
        go.Scatter(
            x=[x0, x1, x1, x0, x0],
            y=[y0, y0, y1, y1, y0],
            fill="toself",
            fillcolor=plot_background,
            mode="lines",
            line=dict(color=exon_border, width=0),
            hoverinfo="text",
            text=geneinfo,
            showlegend=False,
        ),
        row=chrom_ix + 1,
        col=1,
    )

    # Plot INTRON lines
    sorted_exons = df[[START_COL, END_COL]].sort_values(by=START_COL)
    if ts_data:
        ts_chrom = ts_data[chrom]
    else:
        ts_chrom = pd.DataFrame()

    if isinstance(arrow_intron_threshold, int):
        arrow_intron_threshold = coord2percent(
            fig, chrom_ix + 1, 0, arrow_intron_threshold
        )
    if isinstance(arrow_size, int):
        arrow_size = coord2percent(fig, chrom_ix + 1, 0, arrow_size)

    dir_flag = plot_introns(
        sorted_exons,
        ts_chrom,
        fig,
        gene_ix,
        exon_border,
        chrom_ix,
        strand,
        genename,
        arrow_intron_threshold,
        exon_height,
        arrow_color,
        arrow_line_width,
        arrow_size,
    )

    # Plot the gene rows
    apply_gene_bridge(
        transcript_str,
        text,
        text_size,
        df,
        fig,
        strand,
        genename,
        gene_ix,
        exon_border,  # this works as "exon_color" used for utr (not interval)
        exon_border,
        chrom_ix,
        geneinfo,
        showinfo,
        exon_height,
        transcript_utr_width,
        legend,
        arrow_size_min,
        arrow_color,
        arrow_line_width,
        dir_flag,
    )
