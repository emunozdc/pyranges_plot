import numpy as np
from intervaltree import IntervalTree
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pyranges as pr
from matplotlib.patches import Rectangle
import sys
import plotly.colors as pc
import matplotlib.cm as cm
from .core import get_engine, get_warnings, cumdelting
from .matplotlib_base._core import plt_popup_warning


############ COMPUTE INTRONS OFF THRESHOLD
def compute_thresh(df, chrmd_df):
    """Get threshold from limits"""

    chrom = df["Chromosome"].iloc[0]
    limit_range = chrmd_df.loc[chrom]["max"] - chrmd_df.loc[chrom]["min"]
    df["shrink_threshold"] = [int(df["shrink_threshold"].iloc[0] * limit_range)] * len(
        df
    )

    return df


############ SUBSET
def make_subset(df, id_col, max_ngenes):
    """xxx"""

    # create a column indexing all the genes in the df
    genesix_l = [i for i in enumerate(df[id_col].drop_duplicates())]
    genesix_d = {key: value for value, key in genesix_l}
    df["gene_index"] = df[id_col].map(
        genesix_d
    )  # create a column indexing all the genes in the df
    tot_ngenes = max(genesix_l)[0]

    # select maximum number of genes
    if max(df.gene_index) + 1 <= max_ngenes:
        subdf = df
    else:
        subdf = df[df.gene_index < max_ngenes]

    # remove the gene_index column from the original df
    df.drop("gene_index", axis=1, inplace=True)

    return subdf, tot_ngenes


############ GENESMD_DF


###packed
def _genesmd_packed(genesmd_df):
    """xxx"""

    # Initialize IntervalTree and used y-coordinates list
    trees = [IntervalTree()]

    def find_tree(row):
        for tree in trees:
            if not tree.overlaps(row["Start"], row["End"]):
                return tree
        trees.append(IntervalTree())
        return trees[-1]

    # Assign y-coordinates
    for idx, row in genesmd_df.iterrows():
        tree = find_tree(row)
        tree.addi(row["Start"], row["End"], idx)
        genesmd_df.at[idx, "ycoord"] = trees.index(tree)

    return genesmd_df


###colors for genes
def is_pltcolormap(colormap_string):
    """Checks whether the string given is a valid plt colormap name."""

    try:
        colormap = cm.get_cmap(colormap_string)
        if colormap is not None and colormap._isinit:
            return True
        else:
            return False

    except ValueError:
        return False


def is_plycolormap(colormap_string):
    """Checks whether the string given is a valid plotly color object name."""

    if hasattr(pc.sequential, colormap_string):
        return True
    elif hasattr(pc.diverging, colormap_string):
        return True
    elif hasattr(pc.cyclical, colormap_string):
        return True
    elif hasattr(pc.qualitative, colormap_string):
        return True


def get_plycolormap(colormap_string):
    """Provides the plotly color object corresponding to the string given."""

    if hasattr(pc.sequential, colormap_string):
        return getattr(pc.sequential, colormap_string)
    elif hasattr(pc.diverging, colormap_string):
        return getattr(pc.diverging, colormap_string)
    elif hasattr(pc.cyclical, colormap_string):
        return getattr(pc.cyclical, colormap_string)
    elif hasattr(pc.qualitative, colormap_string):
        return getattr(pc.qualitative, colormap_string)


def _genesmd_assigncolor(genesmd_df, colormap):
    """xxx"""

    color_tags = genesmd_df.color_tag.drop_duplicates()
    n_color_tags = len(color_tags)

    # 0-string to colormap object if possible
    if isinstance(colormap, str):
        try:
            if is_pltcolormap(colormap):
                colormap = plt.get_cmap(colormap)
            elif is_plycolormap(colormap):
                colormap = get_plycolormap(colormap)
            else:
                sys.exit(1)
        except SystemExit as e:
            print("The provided string does not match any plt or plotly colormap.", e)

    # 1-plt colormap to list
    if isinstance(colormap, mcolors.ListedColormap):
        colormap = list(colormap.colors)  # colors of plt object

    # 2-list to dict
    if isinstance(colormap, list):
        # adjust number of colors
        if n_color_tags > len(colormap):
            engine = get_engine()
            warnings = get_warnings()
            if engine in ["plt", "matplotlib"] and warnings:
                plt_popup_warning(
                    "The genes are colored by iterating over the given color list."
                )
            elif engine in ["ply", "plotly"] and warnings:
                genesmd_df["_iterwarning!"] = [1] * len(genesmd_df)
        else:
            colormap = colormap[:n_color_tags]
        # make plotly rgb colors compatible with plt
        if colormap[0][:3] == "rgb":
            numb_list = [
                rgb[rgb.find("(") + 1 : rgb.find(")")].split(",") for rgb in colormap
            ]
            colormap = [
                "#{:02x}{:02x}{:02x}".format(int(r), int(g), int(b))
                for r, g, b in numb_list
            ]
        # create dict of colors
        colormap = {
            str(color_tags[i]): colormap[i % len(colormap)] for i in range(n_color_tags)
        }

    # 3- Use dict to assign color to gene
    if isinstance(colormap, dict):
        genesmd_df["color_tag"] = genesmd_df["color_tag"].astype(str)
        genesmd_df["color"] = genesmd_df["color_tag"].map(colormap)
        if genesmd_df["color"].isna().any():
            engine = get_engine()
            warnings = get_warnings()
            if engine in ["plt", "matplotlib"] and warnings:
                plt_popup_warning(
                    "Some genes do not have a color assigned so they are colored in black."
                )
            elif engine in ["ply", "plotly"] and warnings:
                genesmd_df["_blackwarning!"] = [1] * len(genesmd_df)
            genesmd_df["color"].fillna("black", inplace=True)  # black for not specified

    return genesmd_df


def get_genes_metadata(df, id_col, color_col, packed, colormap):
    """xxx"""

    # Start df with chromosome and the column defining color
    if color_col == "Chromosome":
        genesmd_df = df.groupby(id_col, group_keys=False).agg(
            {"Chromosome": "first", "Start": "min", "End": "max"}
        )
    else:
        genesmd_df = df.groupby(id_col, group_keys=False).agg(
            {"Chromosome": "first", "Start": "min", "End": "max", color_col: "first"}
        )
    genesmd_df["chrix"] = genesmd_df["Chromosome"].copy()
    genesmd_df.rename(columns={color_col: "color_tag"}, inplace=True)
    genesmd_df["gene_ix_xchrom"] = genesmd_df.groupby(
        "chrix", group_keys=False
    ).cumcount()

    # Assign y-coordinate to genes
    if packed:
        genesmd_df["ycoord"] = -1
        genesmd_df = genesmd_df.groupby("chrix", group_keys=False).apply(
            _genesmd_packed
        )  # add packed ycoord column

    else:
        genesmd_df["ycoord"] = genesmd_df.loc[:, "gene_ix_xchrom"]

    # Assign color to each gene
    genesmd_df = _genesmd_assigncolor(genesmd_df, colormap)  # adds a column 'color'

    # Add legend info
    def create_legend_rect(color):
        return Rectangle((0, 0), 1, 1, color=color)

    # column with rectangle objects with same color as gene
    genesmd_df["legend_item"] = genesmd_df["color"].apply(create_legend_rect)

    return genesmd_df


############ CHRMD_DF


##limits
def _chrmd_limits(chrmd_df, limits):
    """xxx"""

    # 1- create min_max column containing (plot min, plot max)

    # no limits no info
    if limits is None:
        chrmd_df["min_max"] = [(np.nan, np.nan)] * len(chrmd_df)

    # one tuple for all chromosomes
    elif type(limits) is tuple:
        chrmd_df["min_max"] = [limits] * len(chrmd_df)

    # pyranges object
    elif type(limits) is pr.pyranges_main.PyRanges:
        # create dict to map limits
        limits_df = limits.df
        limits_chrmd_df = limits_df.groupby("Chromosome", group_keys=False).agg(
            {"Start": "min", "End": "max"}
        )
        limits_chrmd_dict = limits_chrmd_df.to_dict(orient="index")

        # function to get matching values from limits_chrmd_df
        def make_min_max(row):
            chromosome = str(row.name)
            limits = limits_chrmd_dict.get(chromosome)
            if limits:
                return (
                    limits["Start"],
                    limits["End"],
                )  # chromosome in both sets of data
            else:
                return (np.nan, np.nan)  # chromosome does not match

        # create limits column in plotting data
        chrmd_df["min_max"] = chrmd_df.apply(make_min_max, axis=1)

    # dictionary as limits
    else:
        chrmd_df["min_max"] = [
            limits.get(index) for index in chrmd_df.index
        ]  # fills with None the chromosomes not specified


def get_chromosome_metadata(df, id_col, limits, genesmd_df, ts_data=None):
    """xxx"""

    # Start df
    chrmd_df = df.groupby("Chromosome", group_keys=False).agg(
        {"Start": "min", "End": "max", id_col: "nunique"}
    )
    chrmd_df.dropna(inplace=True)  # remove chr not present in subset (NaN)
    chrmd_df.rename(
        columns={id_col: "n_genes", "Start": "min", "End": "max"}, inplace=True
    )

    if ts_data:
        chrmd_df["maxdelta"] = df.groupby("Chromosome", group_keys=False).agg(
            {"cumdelta": "max"}
        )
    # Add limits
    _chrmd_limits(chrmd_df, limits)  # unknown limits are nan

    def fill_min_max(row):
        minmax_t = row["min_max"]
        # deal with empty rows
        if minmax_t is None:
            minmax_t = (np.nan, np.nan)

        # check both items and put default if necessary
        minmax_l = list(minmax_t)

        # add default to lower limit
        if minmax_l[0] is None or np.isnan(minmax_l[0]):
            minmax_l[0] = row["min"]
        # add default to higher limit
        if minmax_l[1] is None or np.isnan(minmax_l[1]):
            minmax_l[1] = row["max"]
        # consider introns off for higher limit
        else:
            if len(row) == 5:
                new_upper_lim = cumdelting([minmax_l[1]], ts_data, row.name)
                minmax_l[1] = new_upper_lim[0]

        # put plot coordinates in min_max
        row["min_max"] = minmax_l
        return row

    chrmd_df = chrmd_df.apply(fill_min_max, axis=1)

    # Store plot y height
    chrmd_df["y_height"] = genesmd_df.groupby("chrix", group_keys=False).ycoord.max()
    chrmd_df["y_height"] += 1  # count from 1

    return chrmd_df
