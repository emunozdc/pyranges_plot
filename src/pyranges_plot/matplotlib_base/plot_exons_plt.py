import matplotlib.pyplot as plt
from ..core import get_warnings
from ._core import coord2inches, inches2coord, plt_popup_warning
from ._fig_axes import create_fig
from ._data2plot import (
    make_annotation,
    plot_direction,
    _apply_gene_bridge,
)

# plot parameters
arrow_width = 1
arrow_color = "grey"
arrow_style = "round"
arrow_size_max = 0.3
arrow_size_min = 0.1
intron_threshold = 0.3


# PLOT_EXONS FUNCTIONS


def plot_exons_plt(
    subdf,
    tot_ngenes,
    feat_dict,
    genesmd_df,
    chrmd_df,
    max_ngenes=25,
    id_col="gene_id",
    transcript_str=False,
    showinfo=None,
    legend=False,
    chr_string=None,
    packed=True,
    to_file=None,
    file_size=None,
):
    """xxx"""

    # Get default plot features
    tag_background = feat_dict["tag_background"]
    plot_background = feat_dict["plot_background"]
    plot_border = feat_dict["plot_border"]
    title_dict_plt = feat_dict["title_dict_plt"]
    exon_width = feat_dict["exon_width"]
    transcript_utr_width = feat_dict["transcript_utr_width"]

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
        chr_string,
        title_dict_plt,
        plot_background,
        plot_border,
        packed,
        legend,
    )

    # Plot genes
    subdf.groupby(id_col).apply(
        lambda subdf: _gby_plot_exons(
            subdf,
            axes,
            fig,
            chrmd_df,
            genesmd_df,
            id_col,
            showinfo,
            tag_background,
            transcript_str,
            exon_width,
            transcript_utr_width,
        )
    )

    # Provide output
    if to_file is None:
        # evaluate warning
        warnings = get_warnings()
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
    id_col,
    showinfo,
    tag_background,
    transcript_str,
    exon_width,
    transcript_utr_width,
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
        geneinfo = f"[{strand}] ({min(df.Start)}, {max(df.End)})\nID: {genename}"  # default with strand
    else:
        geneinfo = f"({min(df.Start)}, {max(df.End)})\nID: {genename}"  # default without strand

    # customized
    showinfo_dict = df.iloc[0].to_dict()  # first element of gene rows
    if showinfo:
        geneinfo += "\n" + showinfo.format(**showinfo_dict)

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
        geneinfo,
        exon_width,
        transcript_utr_width,
        arrow_size_min,
        arrow_color,
        arrow_style,
        arrow_width,
    )

    # Evaluate each intron
    sorted_exons = df[["Start", "End"]].sort_values(by="Start")

    for i in range(len(sorted_exons) - 1):
        start = sorted_exons["End"].iloc[i]
        stop = sorted_exons["Start"].iloc[i + 1]
        intron_size = coord2inches(fig, ax, start, stop, 0, 0)
        incl = inches2coord(fig, ax, 0.15)  # how long is the arrow in the plot (OX)

        # Plot LINE binding exons
        intron_line = ax.plot(
            [start, stop], [gene_ix, gene_ix], color=exon_color, linewidth=1, zorder=1
        )

        # Create annotation for intron
        make_annotation(intron_line[0], fig, ax, geneinfo, tag_background)

        # Plot DIRECTION ARROW in INTRONS if strand is known
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
