import matplotlib.pyplot as plt
import pandas as pd

from ._core import plt_popup_warning
from ._fig_axes import create_fig
from ._data2plot import (
    _apply_gene_bridge,
    plot_introns,
)

arrow_style = "round"


def plot_exons_plt(
    subdf,
    tot_ngenes,
    feat_dict,
    genesmd_df,
    chrmd_df,
    ts_data,
    max_ngenes=25,
    id_col="gene_id",
    transcript_str=False,
    showinfo=None,
    legend=False,
    title_chr=None,
    packed=True,
    to_file=None,
    file_size=None,
    warnings=None,
    tick_pos_d=None,
    ori_tick_pos_d=None,
):
    """xxx"""

    # Get default plot features
    tag_background = feat_dict["tag_background"]
    plot_background = feat_dict["plot_background"]
    plot_border = feat_dict["plot_border"]
    title_dict_plt = feat_dict["title_dict_plt"]
    exon_width = feat_dict["exon_width"]
    transcript_utr_width = feat_dict["transcript_utr_width"]
    arrow_line_width = feat_dict["arrow_line_width"]
    arrow_color = feat_dict["arrow_color"]
    arrow_size_min = feat_dict["arrow_size_min"]
    arrow_size = feat_dict["arrow_size"]
    arrow_intron_threshold = feat_dict["arrow_intron_threshold"]

    # Create figure and axes
    if file_size:
        x = file_size[0]
        y = file_size[1]
    else:
        x = 20
        y = (
            sum(chrmd_df.y_height) + 2 * len(chrmd_df)
        ) / 2  # height according to genes and add 2 per each chromosome

    fig, axes = create_fig(
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
    )

    # Plot genes
    subdf.groupby(id_col, group_keys=False).apply(
        lambda subdf: _gby_plot_exons(
            subdf,
            axes,
            fig,
            chrmd_df,
            genesmd_df,
            ts_data,
            id_col,
            showinfo,
            tag_background,
            transcript_str,
            exon_width,
            transcript_utr_width,
            arrow_intron_threshold,
            arrow_line_width,
            arrow_color,
            arrow_size_min,
            arrow_size,
        )
    )

    # Provide output
    if to_file is None:
        # evaluate warning
        if tot_ngenes > max_ngenes and warnings:
            plt_popup_warning(
                "The provided data contains more genes than the ones plotted."
            )
        plt.show()
    else:
        plt.savefig(to_file, format=to_file[-3:])


def _gby_plot_exons(
    df,
    axes,
    fig,
    chrmd_df,
    genesmd_df,
    ts_data,
    id_col,
    showinfo,
    tag_background,
    transcript_str,
    exon_width,
    transcript_utr_width,
    arrow_intron_threshold,
    arrow_line_width,
    arrow_color,
    arrow_size_min,
    arrow_size,
):
    """Plot elements corresponding to the df rows of one gene."""

    # Gene parameters
    genename = df[id_col].iloc[0]
    gene_ix = genesmd_df.loc[genename]["ycoord"] + 0.5
    exon_color = genesmd_df.loc[genename].color
    chrom = genesmd_df.loc[genename].chrix
    chrom_ix = chrmd_df.index.get_loc(chrom)
    ax = axes[chrom_ix]
    if "Strand" in df.columns:
        strand = df["Strand"].unique()[0]
    else:
        strand = ""

    # Make gene annotation
    # get the gene information to print on hover
    # default
    if strand:
        geneinfo = f"[{strand}] ({min(df.oriStart)}, {max(df.oriEnd)})\nID: {genename}"  # default with strand
    else:
        geneinfo = f"({min(df.oriStart)}, {max(df.oriEnd)})\nID: {genename}"  # default without strand

    # Plot INTRON lines
    sorted_exons = df[["Start", "End"]].sort_values(by="Start")
    sorted_exons["intron_dir_flag"] = [0] * len(sorted_exons)
    if ts_data:
        ts_chrom = ts_data[chrom]
    else:
        ts_chrom = pd.DataFrame()

    dir_flag = plot_introns(
        sorted_exons,
        ts_chrom,
        fig,
        ax,
        geneinfo,
        tag_background,
        gene_ix,
        exon_color,
        strand,
        arrow_intron_threshold,
        exon_width,
        arrow_color,
        arrow_style,
        arrow_line_width,
        arrow_size,
    )

    # Plot the gene rows as EXONS
    _apply_gene_bridge(
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
        arrow_line_width,
        dir_flag,
        arrow_size,
    )
