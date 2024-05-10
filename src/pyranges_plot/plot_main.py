import pandas as pd
import pyranges as pr
import plotly.colors
from .core import (
    get_engine,
    get_idcol,
    print_default,
    get_default,
    get_warnings,
)
from .data_preparation import (
    make_subset,
    get_genes_metadata,
    get_chromosome_metadata,
    compute_thresh,
)
from ._introns_off import introns_resize, recalc_axis
from .matplotlib_base.plot_exons_plt import plot_exons_plt
from .plotly_base.plot_exons_ply import plot_exons_ply


def plot(
    df,
    vcf=None,
    engine=None,
    id_col=None,
    warnings=None,
    max_shown=25,
    introns_off=False,
    transcript_str=False,
    color_col=None,
    colormap=None,
    limits=None,
    showinfo=None,
    legend=False,
    id_ann=True,
    title_chr="Chromosome {chrom}",
    packed=True,
    to_file=None,
    file_size=None,
    mode="light",
    **kargs,
):
    """
    Create genes plot from 1/+ PyRanges objects.

    Parameters
    ----------
    df: {pyranges.PyRanges, pandas.DataFrame or list of pyranges.PyRanges}
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

    introns_off: bool, default False
        Whether to compress the intron ranges to facilitate visualization or not.

    transcript_str: bool, default False
        Display differentially transcript regions belonging and not belonging to CDS. The CDS/exon information
        must be stored in the 'Feature' column of the PyRanges object or the dataframe.

    color_col: str, default None
        Name of the column used to color the genes.

    colormap: {matplotlib.colors.ListedColormap, list, str, dict}, default plotly.colors.qualitative.Alphabet
        Sequence of colors for the genes, it can be provided as a Matplotlib colormap,
        a Plotly color sequence (built as lists), a string naming the previously mentioned
        color objects from Matplotlib and Plotly, or a dictionary with the following
        structure {color_column_value1: color1, color_column_value2: color2, ...}. When a specific
        color_col value is not specified in the dictionary it will be colored in black.

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

    showinfo: str, default None
        Dataframe information to show in a tooltip when placing the mouse over a gene, the given
        information will be added to the default: strand, start-end coordinates and id. This must be
        provided as a string containing the column names of the values to be shown within curly brackets.
        For example if you want to show the value of the pointed gene for the column "col1" a valid showinfo
        string could be: "Value of col1: {col1}". Note that the values in the curly brackets are not
        strings. If you want to introduce a newline you can use "\n".

    legend: bool, default False
        Whether the legend should appear in the plot.

    id_ann: bool, default True
        Whether the id annotation should appear beside the gene in the plot.

    title_chr: str, default "Chromosome {chrom}"
        String providing the desired title for the chromosome plots. It should be given in a way where
        the chromosome value in the data is indicated as {chrom}.

    packed: bool, default True
        Disposition of the genes in the plot. Use True for a packed disposition (genes in the same line if
        they do not overlap) and False for unpacked (one row per gene).

    to_file: str, default None
        Name of the file to export specifying the desired extension. The supported extensions are '.png' and '.pdf'.

    file_size: list or tuple, default None
        Size of the plot to export defined by a sequence object like: (height, width). The default values
        make the height according to the number of genes and the width as 20 in Matplotlib and 1600 in Plotly.

    mode: str, default "light"
        General color appearance of the plot. Available modes: "light", "dark".

    **kargs
        Customizable plot features can be defined using kargs. Use print_default() function to check the variables'
        nomenclature, description and default values.



    Examples
    --------

    >>> import pyranges as pr, pyranges_plot as prp

    >>> p = pr.PyRanges({"Chromosome": [1]*5, "Strand": ["+"]*3 + ["-"]*2, "Start": [10,20,30,25,40], "End": [15,25,35,30,50], "transcript_id": ["t1"]*3 + ["t2"]*2}, "feature1": ["A", "B", "C", "A", "B"])

    >>> plot(p, engine='plt', id_col="transcript_id",  max_shown=25, colormap='Set3')

    >>> plot(p, engine='matplotlib', id_col="transcript_id", color_col='Strand', colormap={'+': 'green', '-': 'red'})

    >>> plot(p, engine='ply', id_col="transcript_id", limits = {'1': (1000, 50000), '2': None, '3': (10000, None)})

    >>> plot(p, engine='plotly', id_col="transcript_id", introns_off=True, showinfo = "Feature1: {feature1}")

    >>> plot(df, engine='plt', id_col="transcript_id", color_col='Strand', packed=False, to_file='my_plot.pdf')
    """

    if not isinstance(df, list):
        df = [df]
    # for df_i in df:

    # Deal with export
    if to_file:
        ext = to_file[-4:]
        try:
            if ext not in [".pdf", ".png"]:
                raise Exception(
                    "Please specify the desired format to export the file including either '.png' or '.pdf' as an extension."
                )
        except SystemExit as e:
            print("An error occured:", e)

    # Deal with id column
    if id_col is None:
        id_col = get_idcol()

    try:
        for df_item in df:
            if id_col is not None and id_col not in df_item.columns:
                raise Exception(
                    "Please define a correct name of the ID column using either set_idcol() function or plot_generic parameter as plot_generic(..., id_col = 'your_id_col')"
                )
    except SystemExit as e:
        print("An error occured:", e)

    # Deal with transcript structure
    if transcript_str:
        try:
            for df_item in df:
                if "Feature" not in df_item.columns:
                    raise Exception(
                        "The transcript structure information must be stored in 'Feature' column of the data."
                    )
        except SystemExit as e:
            print("An error occured:", e)

    # Deal with warnings
    if warnings is None:
        warnings = get_warnings()

    # Deal with engine
    if engine is None:
        engine = get_engine()

    try:
        # PREPARE DATA for plot
        # Deal with plot features as kargs
        wrong_keys = [k for k in kargs if k not in print_default(return_keys=True)]
        if wrong_keys:
            raise Exception(
                f"The following keys do not match any customizable features: {wrong_keys}.\nCheck the customizable variable names using the print_default function."
            )

        def getvalue(key, mode):
            if key in kargs:
                value = kargs[key]
                return value  ## add invalid data type??
            else:
                return get_default(key, mode)

        # Get default plot features
        if colormap is None:
            if mode == "light":
                colormap = plotly.colors.qualitative.Alphabet
            elif mode == "dark":
                colormap = "G10"

        feat_dict = {
            "tag_bkg": getvalue("tag_bkg", mode),
            "fig_bkg": getvalue("fig_bkg", mode),
            "plot_bkg": getvalue("plot_bkg", mode),
            "plot_border": getvalue("plot_border", mode),
            "title_dict_plt": {
                "family": "sans-serif",
                "color": getvalue("title_color", mode),
                "size": int(getvalue("title_size", mode)) - 5,
            },
            "title_dict_ply": {
                "family": "Arial",
                "color": getvalue("title_color", mode),
                "size": int(getvalue("title_size", mode)),
            },
            "grid_color": getvalue("grid_color", mode),
            "exon_border": getvalue("exon_border", mode),
            "exon_width": float(getvalue("exon_width", mode)),
            "transcript_utr_width": 0.3 * float(getvalue("exon_width", mode)),
            "v_space": float(getvalue("v_space", mode)),
            "plotly_port": getvalue("plotly_port", mode),
            "arrow_line_width": float(getvalue("arrow_line_width", mode)),
            "arrow_color": getvalue("arrow_color", mode),
            "arrow_size_min": float(getvalue("arrow_size_min", mode)),
            "arrow_size": float(getvalue("arrow_size", mode)),
            "arrow_intron_threshold": getvalue("arrow_intron_threshold", mode),
            "shrink_threshold": getvalue("shrink_threshold", mode),
            "shrinked_bkg": getvalue("shrinked_bkg", mode),
            "shrinked_alpha": float(getvalue("shrinked_alpha", mode)),
        }
        shrink_threshold = feat_dict["shrink_threshold"]

        # Make DataFrame subset if needed
        df_d = {}
        tot_ngenes_l = []
        for pr_ix, df_item in enumerate(df):
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

        if introns_off:
            # compute threshold
            if isinstance(shrink_threshold, int):
                subdf["shrink_threshold"] = [shrink_threshold] * len(subdf)
            elif isinstance(shrink_threshold, float):
                subdf["shrink_threshold"] = [shrink_threshold] * len(subdf)
                subdf = subdf.groupby(
                    "Chromosome", group_keys=False, observed=True
                ).apply(
                    lambda x: compute_thresh(x, chrmd_df_grouped)
                    if not x.empty
                    else None
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

        genesmd_df["ycoord"] = genesmd_df["ycoord"] * feat_dict["v_space"]
        chrmd_df["pr_line"] = chrmd_df["pr_line"] * feat_dict["v_space"]

        if engine == "plt" or engine == "matplotlib":
            plot_exons_plt(
                subdf=subdf,
                vcf=vcf,
                tot_ngenes_l=tot_ngenes_l,
                feat_dict=feat_dict,
                genesmd_df=genesmd_df,
                chrmd_df=chrmd_df,
                chrmd_df_grouped=chrmd_df_grouped,
                ts_data=ts_data,
                max_shown=max_shown,
                id_col=id_col,
                transcript_str=transcript_str,
                showinfo=showinfo,
                legend=legend,
                id_ann=id_ann,
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
                vcf=vcf,
                tot_ngenes_l=tot_ngenes_l,
                feat_dict=feat_dict,
                genesmd_df=genesmd_df,
                chrmd_df=chrmd_df,
                chrmd_df_grouped=chrmd_df_grouped,
                ts_data=ts_data,
                max_shown=max_shown,
                id_col=id_col,
                transcript_str=transcript_str,
                showinfo=showinfo,
                legend=legend,
                id_ann=id_ann,
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
    except SystemExit as e:
        print("An error occured:", e)
