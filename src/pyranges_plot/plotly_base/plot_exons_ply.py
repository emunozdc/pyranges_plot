import plotly.graph_objects as go
import plotly.colors
import plotly.io as pio
import pandas as pd
import numpy as np
from ._core import initialize_dash_app
from ._fig_axes import create_fig, create_fig_with_vfc
from ._data2plot import plot_introns, _apply_gene_bridge


def plot_exons_ply(
    subdf,
    tot_ngenes_l,
    feat_dict,
    genesmd_df,
    vcf,
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
    if vcf is None:
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
            plot_bkg,
            plot_border,
            tick_pos_d,
            ori_tick_pos_d,
            shrinked_bkg,
            shrinked_alpha,
            v_space,
        )

    else:
        fig = create_fig_with_vfc(
            subdf,
            vcf,
            chrmd_df,
            genesmd_df,
            ts_data,
            title_chr,
            title_dict_ply,
            packed,
            plot_bkg,
            tick_pos_d,
            ori_tick_pos_d,
            shrinked_bkg,
            shrinked_alpha,
        )

    # Plot genes
    subdf.groupby([id_col, "pr_ix"], group_keys=False, observed=True).apply(
        lambda subdf: _gby_plot_exons(
            subdf,
            fig,
            chrmd_df,
            chrmd_df_grouped,
            genesmd_df,
            ts_data,
            id_col,
            showinfo,
            legend,
            transcript_str,
            id_ann,
            exon_width,
            exon_border,
            transcript_utr_width,
            plot_bkg,
            arrow_line_width,
            arrow_color,
            arrow_size_min,
            arrow_size,
            arrow_intron_threshold,
            vcf,
            v_space,
        )
    )  # .reset_index(level="pr_ix")

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
            and not "_iterwarning!" in genesmd_df.columns
        ):
            fig.data[0].customdata = np.array(
                [0, 91124, 0]
            )  # [tot_ngenes_l, 91124, 0])
        elif (
            not "_blackwarning!" in genesmd_df.columns
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
        if not file_size:
            fig.update_layout(width=1600, height=800)
        else:
            fig.update_layout(width=file_size[0], height=file_size[1])

        pio.write_image(fig, to_file)


def _gby_plot_exons(
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
    id_ann,
    exon_width,
    exon_border,
    transcript_utr_width,
    plot_background,
    arrow_line_width,
    arrow_color,
    arrow_size_min,
    arrow_size,
    arrow_intron_threshold,
    vcf,
    v_space,
):
    """Plot elements corresponding to the df rows of one gene."""

    # Gene parameters
    chrom = df["Chromosome"].iloc[0]
    pr_ix = df["pr_ix"].iloc[0]
    genename = df[id_col].iloc[0]
    df["legend_tag"] = [genename] + [""] * (len(df) - 1)

    # in case same gene in +1 pr
    genemd = genesmd_df.loc[genename]
    if not isinstance(genemd, pd.Series):
        genemd = genemd[genemd["pr_ix"] == pr_ix]  # in case same gene in +1 pr
        gene_ix = genemd["ycoord"].loc[genename] + 0.5 * v_space
        exon_color = genemd["color"].iloc[0]
    # in case different genes in different pr
    else:
        gene_ix = genemd["ycoord"] + 0.5 * v_space
        exon_color = genemd["color"]

    if exon_border is None:
        exon_border = exon_color

    chrom_ix = chrmd_df_grouped.loc[chrom]["chrom_ix"]

    if vcf is not None:
        chrom_ix = 2 * chrom_ix + 1
    if "Strand" in df.columns:
        strand = df["Strand"].unique()[0]
    else:
        strand = ""

    # Get the gene information to print on hover
    # default
    if strand:
        geneinfo = f"[{strand}] ({min(df.oriStart)}, {max(df.oriEnd)})<br>ID: {genename}"  # default with strand
    else:
        geneinfo = f"({min(df.oriStart)}, {max(df.oriEnd)})<br>ID: {genename}"  # default without strand

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
            showlegend=False,
        ),
        row=chrom_ix + 1,
        col=1,
    )

    # Plot INTRON lines
    sorted_exons = df[["Start", "End"]].sort_values(by="Start")
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
    _apply_gene_bridge(
        transcript_str,
        id_ann,
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
        arrow_size,
    )
