import plotly.graph_objects as go
import plotly.colors
import plotly.io as pio
import pandas as pd
import numpy as np
from ._core import initialize_dash_app
from ._fig_axes import create_fig
from ._data2plot import plot_introns, _apply_gene_bridge


# plot parameters
colormap = plotly.colors.sequential.thermal
arrow_width = 1
arrow_color = "grey"
arrow_size_max = 0.02
arrow_size_min = 0.005
intron_threshold = 0.03


# PLOT_EXONS FUNCTIONS


def plot_exons_ply(
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
    chr_string=None,
    packed=True,
    to_file=None,
    file_size=None,
    warnings=None,
):
    """xxx"""

    # Get default plot features
    # tag_background = feat_dict['tag_background']
    plot_background = feat_dict["plot_background"]
    plot_border = feat_dict["plot_border"]
    title_dict_ply = feat_dict["title_dict_ply"]
    exon_width = feat_dict["exon_width"]
    transcript_utr_width = feat_dict["transcript_utr_width"]

    # Create figure and chromosome plots
    fig = create_fig(
        chrmd_df,
        genesmd_df,
        ts_data,
        chr_string,
        title_dict_ply,
        packed,
        plot_background,
    )

    # Plot genes
    subdf.groupby(id_col, group_keys=False).apply(
        lambda subdf: _gby_plot_exons(
            subdf,
            fig,
            chrmd_df,
            genesmd_df,
            ts_data,
            id_col,
            showinfo,
            legend,
            transcript_str,
            exon_width,
            transcript_utr_width,
            plot_background,
        )
    )

    # Adjust plot display
    fig.update_layout(
        plot_bgcolor=plot_background, font_color=plot_border, showlegend=legend
    )
    fig.update_xaxes(showline=True, linewidth=1, linecolor=plot_border, mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor=plot_border, mirror=True)

    # Provide output
    # insert silent information for warnings
    if warnings:
        fig.data[0].customdata = np.array([tot_ngenes, 0, 0])
        if (
            "_blackwarning!" in genesmd_df.columns
            and "_iterwarning!" in genesmd_df.columns
        ):
            fig.data[0].customdata = np.array([tot_ngenes, 91124, 91321])
        elif (
            "_blackwarning!" in genesmd_df.columns
            and not "_iterwarning!" in genesmd_df.columns
        ):
            fig.data[0].customdata = np.array([tot_ngenes, 91124, 0])
        elif (
            not "_blackwarning!" in genesmd_df.columns
            and "_iterwarning!" in genesmd_df.columns
        ):
            fig.data[0].customdata = np.array([tot_ngenes, 0, 91321])
    else:
        fig.data[0].customdata = np.array(["no warnings"])

    if to_file is None:
        app_instance = initialize_dash_app(fig, max_ngenes)
        app_instance.run_server()

    else:
        if not file_size:
            fig.update_layout(width=1600, height=800)
        else:
            fig.update_layout(width=file_size[0], height=file_size[1])

        pio.write_image(fig, to_file)


def _gby_plot_exons(
    df,
    fig,
    chrmd_df,
    genesmd_df,
    ts_data,
    id_col,
    showinfo,
    legend,
    transcript_str,
    exon_width,
    transcript_utr_width,
    plot_background,
):
    """Plot elements corresponding to the df rows of one gene."""

    # Gene parameters
    genename = df[id_col].iloc[0]
    gene_ix = genesmd_df.loc[genename]["ycoord"] + 0.5
    exon_color = genesmd_df.loc[genename].color
    chrom = genesmd_df.loc[genename].chrix
    chrom_ix = chrmd_df.index.get_loc(chrom)
    if "Strand" in df.columns:
        strand = df["Strand"].unique()[0]
    else:
        strand = ""

    # Get the gene information to print on hover
    # default
    if strand:
        geneinfo = f"[{strand}] ({min(df.Start)}, {max(df.End)})<br>ID: {genename}"  # default with strand
    else:
        geneinfo = f"({min(df.Start)}, {max(df.End)})<br>ID: {genename}"  # default without strand

    # customized
    showinfo_dict = df.iloc[0].to_dict()  # first element of gene rows
    if showinfo:
        showinfo = showinfo.replace("\n", "<br>")
        geneinfo += "<br>" + showinfo.format(**showinfo_dict)

    # add annotation for introns to plot
    x0, x1 = min(df["Start"]), max(df["End"])
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
        ),
        row=chrom_ix + 1,
        col=1,
    )

    # Plot INTRON lines
    sorted_exons = df[["Start", "End", "cumdelta"]].sort_values(by="Start")
    if ts_data:
        ts_chrom = ts_data[chrom]
    else:
        ts_chrom = pd.DataFrame()

    plot_introns(
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
    )

    # Plot the gene rows
    _apply_gene_bridge(
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
    )
