import matplotlib.pyplot as plt
import pandas as pd

from ._core import plt_popup_warning
from ._fig_axes import create_fig
from .plot_vcf_plt import _gby_plot_vcf
from ._data2plot import (
    _apply_gene_bridge,
    plot_introns,
)

arrow_style = "round"


def plot_exons_plt(
    subdf,
    vcf,
    tot_ngenes_l,
    feat_dict,
    genesmd_df,
    chrmd_df,
    chrmd_df_grouped,
    ts_data,
    id_col,
    max_shown=25,
    transcript_str=False,
    showinfo=None,
    legend=False,
    id_ann=True,
    title_chr=None,
    packed=True,
    to_file=None,
    file_size=None,
    warnings=None,
    tick_pos_d=None,
    ori_tick_pos_d=None,
):
    """Create Matplotlib plot."""

    # Get default plot features
    tag_bkg = feat_dict["tag_bkg"]
    fig_bkg = feat_dict["fig_bkg"]
    plot_bkg = feat_dict["plot_bkg"]
    plot_border = feat_dict["plot_border"]
    title_dict_plt = feat_dict["title_dict_plt"]
    grid_color = feat_dict["grid_color"]
    exon_border = feat_dict["exon_border"]
    exon_width = feat_dict["exon_width"]
    transcript_utr_width = feat_dict["transcript_utr_width"]
    v_space = feat_dict["v_space"]
    arrow_line_width = feat_dict["arrow_line_width"]
    arrow_color = feat_dict["arrow_color"]
    arrow_size_min = feat_dict["arrow_size_min"]
    arrow_size = feat_dict["arrow_size"]
    arrow_intron_threshold = feat_dict["arrow_intron_threshold"]
    shrinked_bkg = feat_dict["shrinked_bkg"]
    shrinked_alpha = feat_dict["shrinked_alpha"]

    # Create figure and axes
    if file_size:
        x = file_size[0]
        y = file_size[1]
    else:
        x = 20
        if vcf is None:
            y = (
                sum(chrmd_df_grouped["y_height"]) + len(chrmd_df_grouped.index) * 2
            ) / 2  # height according to genes and add 2 per each chromosome
        else:
            y = (
                sum(chrmd_df_grouped["y_height"])
                + len(chrmd_df.index.columns.get_loc("Chromosome").drop_duplicates())
                * 3
            ) / 2  # increase 1 per chromosome to show variants plot

    fig, axes = create_fig(
        x,
        y,
        chrmd_df,
        chrmd_df_grouped,
        genesmd_df,
        ts_data,
        title_chr,
        title_dict_plt,
        plot_bkg,
        plot_border,
        grid_color,
        packed,
        legend,
        tick_pos_d,
        ori_tick_pos_d,
        tag_bkg,
        fig_bkg,
        shrinked_bkg,
        shrinked_alpha,
        v_space,
    )

    # Plot genes
    subdf.groupby([id_col, "pr_ix"], group_keys=False, observed=True).apply(
        lambda subdf: _gby_plot_exons(
            subdf,
            axes,
            fig,
            chrmd_df,
            chrmd_df_grouped,
            genesmd_df,
            ts_data,
            id_col,
            showinfo,
            tag_bkg,
            plot_border,
            transcript_str,
            id_ann,
            exon_width,
            exon_border,
            transcript_utr_width,
            arrow_intron_threshold,
            arrow_line_width,
            arrow_color,
            arrow_size_min,
            arrow_size,
            v_space,
        )
    )

    # Prevent zoom in y axis
    # for ax in axes:
    #     initial_ylim = ax.get_ylim()
    #     # event handler function to fix y-axis range
    #     def on_xlims_change(axes):
    #         axes.set_ylim(initial_ylim)
    #     ax.callbacks.connect('xlim_changed', on_xlims_change)

    # Provide output
    if to_file is None:
        # evaluate warning
        one_warn = 0
        for tot_ngenes in tot_ngenes_l:
            if tot_ngenes > max_shown and warnings and not one_warn:
                one_warn = 1
                plt_popup_warning(
                    "The provided data contains more genes than the ones plotted."
                )
        plt.show()
    else:
        plt.savefig(to_file, format=to_file[-3:], dpi=400)


def _gby_plot_exons(
    df,
    axes,
    fig,
    chrmd_df,
    chrmd_df_grouped,
    genesmd_df,
    ts_data,
    id_col,
    showinfo,
    tag_bkg,
    plot_border,
    transcript_str,
    id_ann,
    exon_width,
    exon_border,
    transcript_utr_width,
    arrow_intron_threshold,
    arrow_line_width,
    arrow_color,
    arrow_size_min,
    arrow_size,
    v_space,
):
    """Plot elements corresponding to the df rows of one gene."""

    # Gene parameters
    chrom = df["Chromosome"].iloc[0]
    chr_ix = chrmd_df_grouped.loc[chrom]["chrom_ix"]
    if isinstance(chr_ix, pd.Series):
        chr_ix = chr_ix.iloc[0]
    pr_ix = df["pr_ix"].iloc[0]
    genename = df[id_col].iloc[0]
    genesmd_df = genesmd_df.loc[genename]  # store data for the gene
    # in case same gene in +1 pr
    if not isinstance(genesmd_df, pd.Series):
        genesmd_df = genesmd_df[
            genesmd_df["pr_ix"] == pr_ix
        ]  # in case same gene in +1 pr
        gene_ix = genesmd_df["ycoord"].loc[genename] + 0.5 * v_space
        exon_color = genesmd_df["color"].iloc[0]
    # in case gene in single pr
    else:
        gene_ix = genesmd_df["ycoord"] + 0.5 * v_space
        exon_color = genesmd_df["color"]

    if exon_border is None:
        exon_border = exon_color

    ax = axes[chr_ix]
    if "Strand" in df.columns:
        strand = df["Strand"].unique()[0]
    else:
        strand = ""

    # Make gene annotation
    # get the gene information to print on hover
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
        tag_bkg,
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
        id_ann,
        df,
        fig,
        ax,
        strand,
        gene_ix,
        exon_color,
        exon_border,
        tag_bkg,
        plot_border,
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
