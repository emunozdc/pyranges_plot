import pandas as pd
import pyranges as pr
import plotly.colors
from .core import (
    get_engine,
    get_id_col,
    print_options,
    get_options,
    get_warnings,
    set_theme,
    get_theme,
    set_options,
)
from .data_preparation import (
    make_subset,
    get_genes_metadata,
    get_chromosome_metadata,
    compute_thresh,
)
from .introns_off import introns_resize, recalc_axis
from .matplotlib_base.plot_exons_plt import plot_exons_plt
from .plotly_base.plot_exons_ply import plot_exons_ply


def plot(
    data,
    *,
    engine=None,
    id_col=None,
    warnings=None,
    max_shown=25,
    packed=True,
    color_col=None,
    shrink=False,
    limits=None,
    thick_cds=False,
    text=False,
    legend=False,
    title_chr="Chromosome {chrom}",
    y_labels=None,
    tooltip=None,
    to_file=None,
    theme=None,
    **kargs,
):
    """
    Create genes plot from 1/+ PyRanges objects.

    Parameters
    ----------
    data: {pyranges.PyRanges or list of pyranges.PyRanges}
        Pyranges, derived dataframe or list of them with annotation data.

    engine: str, default None
        Library in which the plot should be built, it accepts either Matplotlib ['matplotlib'/'plt'] or
        Plotly ['ply'/'plotly'].

    id_col: str, default None
        Name of the column containing gene ID.

    warnings: bool, default True
        Whether the warnings should be shown or not.

    max_shown: int, default 20
        Maximum number of genes plotted in the dataframe order.

    packed: bool, default True
        Disposition of the genes in the plot. Use True for a packed disposition (genes in the same line if
        they do not overlap) and False for unpacked (one row per gene).

    color_col: str, default None
        Name of the column used to color the genes.

    shrink: bool, default False
        Whether to compress the intron ranges to facilitate visualization or not.

    limits: {None, dict, tuple, pyranges.pyranges_main.PyRanges}, default None
        Customization of coordinates for the chromosome plots.
        - None: minimum and maximum exon coordinate plotted plus a 5% of the range on each side.
        - dict: {chr_name1: (min_coord, max coord), chr_name2: (min_coord, max_coord), ...}. Not
        all the plotted chromosomes need to be specified in the dictionary and some coordinates
        can be indicated as None, both cases lead to the use of the default value.
        - tuple: the coordinate limits of all chromosomes will be defined as indicated.
        - pyranges.pyranges_main.PyRanges: for each matching chromosome between the plotted data
        and the limits data, the limits will be defined by the minimum and maximum coordinates
        in the pyranges object defined as limits. If some plotted chromosomes are not present they
        will be left as default.

    thick_cds: bool, default False
        Display differentially transcript regions belonging and not belonging to CDS. The CDS/exon information
        must be stored in the 'Feature' column of the PyRanges object or the dataframe.

    text: {bool, '{string}'}, default False
        Whether ann annotation should appear beside the gene in the plot. If True, the id/index will be used. To
        customize the annotation use the '{string}' option to choose another data column.

    legend: bool, default False
        Whether the legend should appear in the plot.

    title_chr: str, default "Chromosome {chrom}"
        String providing the desired title for the chromosome plots. It should be given in a way where
        the chromosome value in the data is indicated as {chrom}.

    y_labels: list, default None
        Name to identify the PyRanges object/s in the plot.

    tooltip: str, default None
        Dataframe information to show in a tooltip when placing the mouse over a gene, the given
        information will be added to the default: strand, start-end coordinates and id. This must be
        provided as a string containing the column names of the values to be shown within curly brackets.
        For example if you want to show the value of the pointed gene for the column "col1" a valid tooltip
        string could be: "Value of col1: {col1}". Note that the values in the curly brackets are not
        strings. If you want to introduce a newline you can use "\n".

    to_file: {str, tuple}, default None
        Name of the file to export specifying the desired extension. The supported extensions are '.png' and '.pdf'.
        Optionally, a tuple can be privided where the file name is specified as a str in the first position and in the
        second position there is a tuple specifying the height and width of the figure in px.

    theme: str, default "light"
        General color appearance of the plot. Available modes: "light", "dark".

    **kargs
        Customizable plot features can be defined using kargs. Use print_options() function to check the variables'
        nomenclature, description and default values.



    Examples
    --------

    >>> import pyranges as pr, pyranges_plot as prp

    >>> p = pr.PyRanges({"Chromosome": [1]*5, "Strand": ["+"]*3 + ["-"]*2, "Start": [10,20,30,25,40], "End": [15,25,35,30,50], "transcript_id": ["t1"]*3 + ["t2"]*2}, "feature1": ["A", "B", "C", "A", "B"])

    >>> plot(p, engine='plt', id_col="transcript_id",  max_shown=25, colormap='Set3')

    >>> plot(p, engine='matplotlib', id_col="transcript_id", color_col='Strand', colormap={'+': 'green', '-': 'red'})

    >>> plot(p, engine='ply', id_col="transcript_id", limits = {'1': (1000, 50000), '2': None, '3': (10000, None)})

    >>> plot(p, engine='plotly', id_col="transcript_id", shrink=True, tooltip = "Feature1: {feature1}")

    >>> plot(data, engine='plt', id_col="transcript_id", color_col='Strand', packed=False, to_file='my_plot.pdf')
    """

    if not isinstance(data, list):
        data = [data]
    # for df_i in df:

    # Deal with export
    if to_file:
        if isinstance(to_file, str):
            ext = to_file[-4:]
            if ext not in [".pdf", ".png"]:
                raise Exception(
                    "Please specify the desired format to export the file including either '.png' or '.pdf' as an extension."
                )
            file_size = (1600, 800)
        else:
            ext = to_file[0][-4:]
            if ext not in [".pdf", ".png"]:
                raise Exception(
                    "Please specify the desired format to export the file including either '.png' or '.pdf' as an extension."
                )
            file_size = to_file[1]
            to_file = to_file[0]
    else:
        file_size = (1600, 800)

    # Deal with id column
    if id_col is None:
        id_col = get_id_col()

    for df_item in data:
        if id_col is not None and id_col not in df_item.columns:
            raise Exception(
                "Please define a correct name of the ID column using either set_id_col() function or plot_generic parameter as plot_generic(..., id_col = 'your_id_col')"
            )

    # Deal with transcript structure
    if thick_cds:
        for df_item in data:
            if "Feature" not in df_item.columns:
                raise Exception(
                    "The transcript structure information must be stored in 'Feature' column of the data."
                )

    # Deal with warnings
    if warnings is None:
        warnings = get_warnings()

    # Deal with engine
    if engine is None:
        engine = get_engine()

    # PREPARE DATA for plot
    # Deal with plot features as kargs
    wrong_keys = [k for k in kargs if k not in print_options(return_keys=True)]
    if wrong_keys:
        raise Exception(
            f"The following keys do not match any customizable features: {wrong_keys}.\nCheck the customizable variable names using the print_options function."
        )

    def getvalue(key):
        if key in kargs:
            value = kargs[key]
            return value  ## add invalid data type??
        else:
            return get_options(key)

    # Get default plot features
    # store old options to reset them after the plot
    oldtheme = get_theme()
    oldfeat_dict = get_options("values")

    # check option modifications in params
    if theme is None:  # not specified in params
        theme = get_theme()
        if theme is None:  # not specified with set_theme
            theme = "light"
    set_theme(theme)

    feat_dict = {
        "colormap": getvalue("colormap"),
        "tag_bkg": getvalue("tag_bkg"),
        "fig_bkg": getvalue("fig_bkg"),
        "plot_bkg": getvalue("plot_bkg"),
        "plot_border": getvalue("plot_border"),
        "title_dict_plt": {
            "family": "sans-serif",
            "color": getvalue("title_color"),
            "size": int(getvalue("title_size")) - 5,
        },
        "title_dict_ply": {
            "family": "Arial",
            "color": getvalue("title_color"),
            "size": int(getvalue("title_size")),
        },
        "grid_color": getvalue("grid_color"),
        "exon_border": getvalue("exon_border"),
        "exon_width": float(getvalue("exon_width")),
        "transcript_utr_width": 0.3 * float(getvalue("exon_width")),
        "text_pad": float(getvalue("text_pad")),
        "v_space": float(getvalue("v_space")),
        "plotly_port": getvalue("plotly_port"),
        "arrow_line_width": float(getvalue("arrow_line_width")),
        "arrow_color": getvalue("arrow_color"),
        "arrow_size_min": float(getvalue("arrow_size_min")),
        "arrow_size": float(getvalue("arrow_size")),
        "arrow_intron_threshold": getvalue("arrow_intron_threshold"),
        "shrink_threshold": getvalue("shrink_threshold"),
        "shrinked_bkg": getvalue("shrinked_bkg"),
        "shrinked_alpha": float(getvalue("shrinked_alpha")),
    }
    shrink_threshold = feat_dict["shrink_threshold"]
    colormap = feat_dict["colormap"]

    # restore options set before plot is called
    set_theme(oldtheme)
    set_options(oldfeat_dict)

    # Make DataFrame subset if needed
    df_d = {}
    tot_ngenes_l = []
    for pr_ix, df_item in enumerate(data):
        # deal with empty PyRanges
        if df_item.empty:
            continue
        df_item = df_item.copy()

        # consider not known id_col, plot each interval individually
        if id_col is None:
            df_item["id_col"] = [str(i) for i in range(len(df_item))]
            df_d[pr_ix], tot_ngenes = make_subset(df_item, "id_col", max_shown)
            tot_ngenes_l.append(tot_ngenes)

        # known id_col
        else:
            df_d[pr_ix], tot_ngenes = make_subset(df_item, id_col, max_shown)
            tot_ngenes_l.append(tot_ngenes)

    # set not known id_col
    if id_col is None:
        id_col = "id_col"

    # concat subset dataframes and create new column with input list index
    if not df_d:
        raise Exception("The provided PyRanges object/s are empty.")
    subdf = pd.concat(df_d, names=["pr_ix"]).reset_index(
        level="pr_ix"
    )  ### change to pr but doesn't work yet!!

    # Create genes metadata DataFrame
    if color_col is None:
        color_col = id_col
    genesmd_df = get_genes_metadata(
        subdf, id_col, color_col, packed, colormap, feat_dict["v_space"]
    )

    # Create chromosome metadata DataFrame
    chrmd_df, chrmd_df_grouped = get_chromosome_metadata(
        subdf, id_col, limits, genesmd_df, packed, feat_dict["v_space"]
    )

    # Deal with introns off
    # adapt coordinates to shrinked
    ts_data = {}
    subdf["oriStart"] = subdf["Start"]
    subdf["oriEnd"] = subdf["End"]
    tick_pos_d = {}
    ori_tick_pos_d = {}

    if shrink:
        # compute threshold
        if isinstance(shrink_threshold, int):
            subdf["shrink_threshold"] = [shrink_threshold] * len(subdf)
        elif isinstance(shrink_threshold, float):
            subdf["shrink_threshold"] = [shrink_threshold] * len(subdf)
            subdf = subdf.groupby("Chromosome", group_keys=False, observed=True).apply(
                lambda x: compute_thresh(x, chrmd_df_grouped) if not x.empty else None
            )

        subdf = subdf.groupby("Chromosome", group_keys=False, observed=True).apply(
            lambda x: introns_resize(x, ts_data, id_col)  # if not x.empty else None
        )  # empty rows when subset
        subdf["Start"] = subdf["Start_adj"]
        subdf["End"] = subdf["End_adj"]

        # recompute limits
        chrmd_df, chrmd_df_grouped = get_chromosome_metadata(
            subdf,
            id_col,
            limits,
            genesmd_df,
            packed,
            feat_dict["v_space"],
            ts_data=ts_data,
        )

        # compute new axis values and positions if needed
        if ts_data:
            tick_pos_d, ori_tick_pos_d = recalc_axis(
                ts_data, tick_pos_d, ori_tick_pos_d
            )

    else:
        subdf["cumdelta"] = [0] * len(subdf)

    # Sort data to plot chromosomes and pr objects in order
    subdf.sort_values(["Chromosome", "pr_ix", id_col, "Start"], inplace=True)
    chrmd_df.sort_values(["Chromosome", "pr_ix"], inplace=True)
    subdf["exon_ix"] = subdf.groupby(
        ["Chromosome", "pr_ix", id_col], group_keys=False, observed=True
    ).cumcount()

    # Adjust vertical space
    genesmd_df["ycoord"] = genesmd_df["ycoord"] * feat_dict["v_space"]
    chrmd_df["pr_line"] = chrmd_df["pr_line"] * feat_dict["v_space"]

    if engine == "plt" or engine == "matplotlib":
        plot_exons_plt(
            subdf=subdf,
            tot_ngenes_l=tot_ngenes_l,
            feat_dict=feat_dict,
            genesmd_df=genesmd_df,
            chrmd_df=chrmd_df,
            chrmd_df_grouped=chrmd_df_grouped,
            ts_data=ts_data,
            max_shown=max_shown,
            id_col=id_col,
            transcript_str=thick_cds,
            tooltip=tooltip,
            legend=legend,
            y_labels=y_labels,
            text=text,
            title_chr=title_chr,
            packed=packed,
            to_file=to_file,
            file_size=file_size,
            warnings=warnings,
            tick_pos_d=tick_pos_d,
            ori_tick_pos_d=ori_tick_pos_d,
        )

    elif engine == "ply" or engine == "plotly":
        plot_exons_ply(
            subdf=subdf,
            tot_ngenes_l=tot_ngenes_l,
            feat_dict=feat_dict,
            genesmd_df=genesmd_df,
            chrmd_df=chrmd_df,
            chrmd_df_grouped=chrmd_df_grouped,
            ts_data=ts_data,
            max_shown=max_shown,
            id_col=id_col,
            transcript_str=thick_cds,
            tooltip=tooltip,
            legend=legend,
            y_labels=y_labels,
            text=text,
            title_chr=title_chr,
            packed=packed,
            to_file=to_file,
            file_size=file_size,
            warnings=warnings,
            tick_pos_d=tick_pos_d,
            ori_tick_pos_d=ori_tick_pos_d,
        )

    else:
        raise Exception(
            "Please define engine with set_engine() or specifying it with the 'engine' parameter."
        )
