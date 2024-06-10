import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
import numpy as np
from pyranges.core.names import CHROM_COL, START_COL, END_COL, STRAND_COL

from .core import initialize_dash_app
from .fig_axes import create_fig
from .data2plot import plot_introns, apply_gene_bridge
from ..names import PR_INDEX_COL


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
    exon_width = feat_dict["exon_width"]
    transcript_utr_width = feat_dict["transcript_utr_width"]
    text_size = feat_dict["text_size"]
    # text_pad = feat_dict["text_pad"]
    v_space = feat_dict["v_space"]
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
        plot_bkg,
        plot_border,
        tick_pos_d,
        ori_tick_pos_d,
        shrinked_bkg,
        shrinked_alpha,
        v_space,
    )

    # Plot genes
    subdf.groupby(id_col + [PR_INDEX_COL], group_keys=False, observed=True).apply(
        lambda subdf: gby_plot_exons(
            subdf,
            fig,
            chrmd_df,
            chrmd_df_grouped,
            genesmd_df,
            ts_data,
            id_col,
            tooltip,
            legend,
            transcript_str,
            text,
            text_size,
            exon_width,
            exon_border,
            transcript_utr_width,
            plot_bkg,
            arrow_line_width,
            arrow_color,
            arrow_size_min,
            arrow_size,
            arrow_intron_threshold,
            v_space,
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
    chrmd_df,
    chrmd_df_grouped,
    genesmd_df,
    ts_data,
    id_col,
    showinfo,
    legend,
    transcript_str,
    text,
    text_size,
    exon_width,
    exon_border,
    transcript_utr_width,
    plot_background,
    arrow_line_width,
    arrow_color,
    arrow_size_min,
    arrow_size,
    arrow_intron_threshold,
    v_space,
):
    """Plot elements corresponding to the df rows of one gene."""

    # Gene parameters
    chrom = df[CHROM_COL].iloc[0]
    pr_ix = df[PR_INDEX_COL].iloc[0]
    genename = df["__id_col_2count__"].iloc[0]
    genemd = genesmd_df.loc[genename]  # store data for the gene
    df["legend_tag"] = [genename] + [""] * (len(df) - 1)

    # in case same gene in +1 pr
    if not isinstance(genemd, pd.Series):
        genemd = genemd[genemd[PR_INDEX_COL] == pr_ix]  # in case same gene in +1 pr
        genemd = pd.Series(genemd.iloc[0])
    gene_ix = genemd["ycoord"] + 0.5 * v_space
    exon_color = genemd["color"]

    if exon_border is None:
        exon_border = exon_color

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
    y0, y1 = gene_ix - exon_width / 160, gene_ix + exon_width / 160

    fig.add_trace(
        go.Scatter(
            x=[x0, x1, x1, x0, x0],
            y=[y0, y0, y1, y1, y0],
            fill="toself",
            fillcolor=plot_background,
            mode="lines",
            line=dict(color=exon_color, width=0),
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

    dir_flag = plot_introns(
        sorted_exons,
        ts_chrom,
        fig,
        gene_ix,
        exon_color,
        chrom_ix,
        strand,
        genename,
        arrow_intron_threshold,
        exon_width,
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
        exon_color,
        exon_border,
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
        genemd,
        id_col,
    )
