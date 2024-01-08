import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from ..core import (
    coord2inches,
    inches2coord,
    print_default,
    get_default,
    plt_popup_warning,
)
from ..data_preparation import (
    make_subset,
    get_genes_metadata,
    get_chromosome_metadata,
)
from ..plt_func import (
    create_fig,
    make_annotation,
    plot_direction,
    _apply_gene,
)

# plot parameters
exon_width = 0.4
transcript_utr_width = 0.3 * exon_width
arrow_width = 1
arrow_color = "grey"
arrow_style = "round"
arrow_size_max = 0.3
arrow_size_min = 0.1
intron_threshold = 0.3


# PLOT_EXONS FUNCTIONS


def plot_exons_plt(
    df,
    max_ngenes=25,
    id_col="gene_id",
    transcript_str=False,
    color_col=None,
    colormap=None,
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
        make the height according to the number of genes and the width as 20.

    **kargs:

        Customizable plot features can be defined using kargs. Use print_default() function to check the variables'
        nomenclature, description and default values.



    Examples
    --------

    >>> plot_exons_plt(df, max_ngenes=25, colormap='Set3')

    >>> plot_exons_plt(df, color_col='Strand', colormap={'+': 'green', '-': 'red'})

    >>> plot_exons_plt(df, limits = {'1': (1000, 50000), '2': None, '3': (10000, None)})


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
    tag_background = getvalue("tag_background")
    plot_background = getvalue("plot_background")
    plot_border = getvalue("plot_border")
    title_dict_plt = {
        "family": "sans-serif",
        "color": getvalue("title_color"),
        "size": int(getvalue("title_size")) - 5,
    }

    # Make DataFrame subset if needed
    subdf, tot_ngenes = make_subset(df, id_col, max_ngenes)

    # Create genes metadata DataFrame
    if color_col is None:
        color_col = id_col
    genesmd_df = get_genes_metadata(subdf, id_col, color_col, packed, colormap)

    # Create chromosome metadata DataFrame
    chrmd_df = get_chromosome_metadata(subdf, id_col, limits, genesmd_df)

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
        )
    )

    # Provide output
    if to_file is None:
        # evaluate warning
        if tot_ngenes > max_ngenes:
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
    #default
    if strand:
        geneinfo = f"[{strand}] ({min(df.Start)}, {max(df.End)})\nID: {genename}"  # default with strand
    else:
        geneinfo = f"({min(df.Start)}, {max(df.End)})\nID: {genename}"  # default without strand

    #customized
    showinfo_dict = df.iloc[0].to_dict() # first element of gene rows
    if showinfo:
        geneinfo += '\n' + showinfo.format(**showinfo_dict)

    # Plot the gene rows as EXONS
    _apply_gene(
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


def _plot_row(
    row,
    fig,
    ax,
    strand,
    gene_ix,
    exon_color,
    tag_background,
    geneinfo,
    exon_width,
):
    """Plot elements corresponding to one row of one gene."""

    # Exon start and stop
    start = int(row["Start"])
    stop = int(row["End"])

    # Plot EXON as rectangle
    exon_rect = Rectangle(
        (start, gene_ix - exon_width / 2),
        stop - start,
        exon_width,
        edgecolor=exon_color,
        facecolor=exon_color,
        fill=True,
    )
    ax.add_patch(exon_rect)

    # create annotation for exon
    make_annotation(exon_rect, fig, ax, geneinfo, tag_background)

    # Plot DIRECTION ARROW in EXON
    # decide about placing a direction arrow
    arrow_size = coord2inches(fig, ax, 0.05 * start, 0.05 * stop, 0, 0)
    incl = inches2coord(fig, ax, 0.15)  # how long in the plot (OX)

    # create and plot lines
    plot_direction(
        ax,
        strand,
        arrow_size,
        arrow_size_min,
        start,
        stop,
        incl,
        gene_ix,
        exon_width,
        arrow_color,
        arrow_style,
        arrow_width,
    )
