import plotly.graph_objects as go
import plotly.colors
import plotly.io as pio
import numpy as np
from ..core import (
    coord2percent,
    percent2coord,
    print_default,
    get_default,
)
from ..data_preparation import (
    make_subset,
    get_genes_metadata,
    get_chromosome_metadata,
)
from ..ply_func import (
    create_fig,
    plot_direction,
    _apply_gene,
)


# plot parameters
colormap = plotly.colors.sequential.thermal
arrow_width = 1
arrow_color = "grey"
arrow_size_max = 0.02
arrow_size_min = 0.005
intron_threshold = 0.03


# PLOT_EXONS FUNCTIONS


def plot_exons_ply(
    df,
    max_ngenes=25,
    id_col="gene_id",
    transcript_str=False,
    color_col=None,
    colormap=colormap,
    limits=None,
    showinfo=None,
    legend=False,
    chr_string=None,
    packed=True,
    to_file=None,
    file_size=None,
    **kargs,
):
    """
    Create genes plot from PyRanges object DataFrame

    Parameters
    ----------
    df: pandas.DataFrame

        Dataframe with genes.

    max_ngenes: int, default 20

        Maximum number of genes plotted in the dataframe order.

    id_col: str, default 'gene_id'

        Name of the column containing gene ID.

    transcript_str: bool, default False

        Display differentially transcript regions belonging and not belonging to CDS. The CDS/exon information
        must be stored in the 'Feature' column of the PyRanges object or the dataframe.

    color_col: str, default None

        Name of the column used to color the genes.

    colormap: {matplotlib.colors.ListedColormap, list, str, dict}, default plotly.colors.sequential.thermal

        Sequence of colors for the genes, it can be provided as a Matplotlib colormap,
        a Plotly color sequence (built as lists), a string naming the previously mentioned
        color objects from Matplotlib and Plotly, or a dictionary with the following
        structure {color_column_value1: color1, color_column_value2: color2, ...}

    limits: {None, dict, tuple}, default None

        Customization of coordinates for the chromosome plots.
        - None: minimun and maximum exon coordinate plotted plus a 5% of the range on each side.
        - dict: {chr_name1: (min_coord, max coord), chr_name2: (min_coord, max_coord), ...}. Not
        all the plotted chromosomes need to be specified in the dictionary and some coordinates
        can be indicated as None, both cases lead to the use of the default value.
        - tuple: the coordinate limits of all chromosomes will be defined as indicated.
        - pyranges.pyranges_main.PyRanges: for each matching chromosome between the plotted data
        and the limits data, the limits will be defined by the minimum and maximum coordinates
        in the pyranges object defined as limits. If some plotted chromosomes are not present they
        will be left as default.

    showinfo: str, default None

        Dataframe information to show in a tooltip when placing the mouse over a gene. This must be
        provided as a string containing the column names of the values to be shown within curly brackets.
        For example if you want to show the value of the pointed gene for the column "col1" a valid showinfo
        string could be: "Value of col1: {col1}". Note that the values in the curly brackets are not
        strings. If you want to introduce a newline you can use "\n". By default, it shows the ID of the
        gene followed by its start and end position.

    legend: bool, default False

        Whether or not the legend should appear in the plot.

    chr_string: str, default "Chromosome {chrom}"

        String providing the desired titile for the chromosome plots. It should be given in a way where
        the chromosome value in the data is indicated as {chrom}.

    packed: bool, default True

        Disposition of the genes in the plot. Use True for a packed disposition (genes in the same line if
        they do not overlap) and False for unpacked (one row per gene).

    to_file: str, default None

        Name of the file to export specifying the desired extension. The supported extensions are '.png' and '.pdf'.

    file_size: {list, tuple}, default None

        Size of the plot to export defined by a sequence object like: (height, width). The default values
        make the height according to the number of genes and the width as 1600.

    **kargs:

        Customizable plot features can be defined using kargs. Use print_default() function to check the variables'
        nomenclature, description and default values.



    Examples
    --------

    >>> plot_exons_ply(df, max_ngenes=25, colormap='Set3')

    >>> plot_exons_ply(df, color_col='Strand', colormap={'+': 'green', '-': 'red'})

    >>> plot_exons_ply(df, limits = {'1': (1000, 50000), '2': None, '3': (10000, None)})


    """

    # Deal with plot features as kargs
    wrong_keys = [k for k in kargs if k not in print_default(return_keys=True)]
    if len(wrong_keys):
        raise Exception(
            "The following keys do not match any customizable features: {wrong_keys}.\nCheck the print_default function to see the customizable names"
        )

    def getvalue(key):
        if key in kargs:
            value = kargs[key]
            return value  ## add invalid data type??
        else:
            return get_default(key)

    # Get default plot features
    # tag_background = getvalue("tag_background")
    plot_background = getvalue("plot_background")
    plot_border = getvalue("plot_border")
    title_dict_ply = {
        "family": "Arial",
        "color": getvalue("title_color"),
        "size": int(getvalue("title_size")),
    }
    exon_width = getvalue("exon_width")
    transcript_utr_width = 0.3 * exon_width

    # Make DataFrame subset if needed
    subdf, tot_ngenes = make_subset(df, id_col, max_ngenes)

    # Create genes metadata DataFrame
    if color_col is None:
        color_col = id_col
    genesmd_df = get_genes_metadata(subdf, id_col, color_col, packed, colormap)

    # Create chromosome metadata DataFrame
    chrmd_df = get_chromosome_metadata(subdf, id_col, limits, genesmd_df)

    # Create figure and chromosome plots
    fig = create_fig(chrmd_df, genesmd_df, chr_string, title_dict_ply, packed)

    # Plot genes
    subdf.groupby(id_col).apply(
        lambda subdf: _gby_plot_exons(
            subdf,
            fig,
            chrmd_df,
            genesmd_df,
            id_col,
            showinfo,
            legend,
            transcript_str,
            exon_width,
            transcript_utr_width,
        )
    )

    # Adjust plot display
    fig.update_layout(
        plot_bgcolor=plot_background, font_color=plot_border, showlegend=legend
    )
    fig.update_xaxes(showline=True, linewidth=1, linecolor=plot_border, mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor=plot_border, mirror=True)

    # Provide output
    # insert silent information for subset warning
    fig.data[0].customdata = np.array(tot_ngenes)
    # insert legend position
    # fig.update_layout(legend = dict(x=1, y=1))

    if to_file is None:
        return fig
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
    id_col,
    showinfo,
    legend,
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

    # Evaluate each intron
    sorted_exons = df[["Start", "End"]].sort_values(by="Start")

    for i in range(len(sorted_exons) - 1):
        start = sorted_exons["End"].iloc[i]
        stop = sorted_exons["Start"].iloc[i + 1]
        intron_size = coord2percent(fig, chrom_ix + 1, start, stop)
        incl = percent2coord(fig, chrom_ix + 1, 0.003)  # how long in the plot (OX)

        # Plot LINE binding exons
        # line as rectangle to have annotation
        x0, x1 = min(df.Start), max(df.End)
        y0, y1 = gene_ix - exon_width / 150, gene_ix + exon_width / 150
        exon_line = go.Scatter(
            x=[x0, x1, x1, x0, x0],
            y=[y0, y0, y1, y1, y0],
            fill="toself",
            fillcolor=exon_color,
            mode="lines",
            line=dict(color=exon_color, width=0.5),
            text=geneinfo,
        )
        fig.add_trace(exon_line, row=chrom_ix + 1, col=1)

        # Plot DIRECTION ARROW in INTRONS if strand is known
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
        )

    # Plot the gene rows
    _apply_gene(
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
