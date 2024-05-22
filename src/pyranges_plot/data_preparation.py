import numpy as np
from intervaltree import IntervalTree
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pyranges as pr
from matplotlib.patches import Rectangle
import sys
import plotly.colors as pc
from .core import get_engine, get_warnings, cumdelting
from .matplotlib_base._core import plt_popup_warning


############ COMPUTE INTRONS OFF THRESHOLD
def compute_thresh(df, chrmd_df_grouped):
    """Get shrink threshold from limits"""

    chrom = df["Chromosome"].iloc[0]
    chrmd = chrmd_df_grouped.loc[chrom]
    limit_range = chrmd["max"] - chrmd["min"]
    df["shrink_threshold"] = [int(df["shrink_threshold"].iloc[0] * limit_range)] * len(
        df
    )

    return df


############ SUBSET
def make_subset(df, id_col, max_shown):
    """Reduce the number of genes to work with."""

    # create a column indexing all the genes in the df
    genesix_l = [i for i in enumerate(df[id_col].drop_duplicates())]
    genesix_d = {key: value for value, key in genesix_l}
    df["gene_index"] = df[id_col].map(
        genesix_d
    )  # create a column indexing all the genes in the df
    tot_ngenes = max(genesix_l)[0]

    # select maximum number of genes
    if max(df.gene_index) + 1 <= max_shown:
        subdf = df
    else:
        subdf = df[df.gene_index < max_shown]

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


def _update_y(genesmd_df):
    """xxx"""

    min_pr_ix = genesmd_df["pr_ix"].min()

    # Consider pr dividing lines spot
    genesmd_df["update_y_prix"] = genesmd_df.groupby("pr_ix").ngroup(ascending=False)

    # Consider the height of the previous pr to update y coords
    y_prev_df = (
        genesmd_df.groupby("pr_ix")["ycoord"]
        .max()
        .shift(-1, fill_value=-1)
        .apply(lambda x: x + 1)
        .loc[::-1]
        .cumsum()[::-1]
    )
    y_prev_df.name = "update_y_prev"
    genesmd_df = genesmd_df.join(y_prev_df, on="pr_ix")
    genesmd_df["ycoord"] += genesmd_df["update_y_prix"] + genesmd_df["update_y_prev"]

    return genesmd_df


###colors for genes
def is_pltcolormap(colormap_string):
    """Checks whether the string given is a valid plt colormap name."""
    try:
        colormap = plt.colormaps[colormap_string]
        if colormap is not None and isinstance(colormap, mcolors.Colormap):
            return True
        else:
            return False

    except KeyError:
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
    """Match color with gene"""

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
            str(color_tags.iloc[i]): colormap[i % len(colormap)]
            for i in range(n_color_tags)
        }

    # 3- Use dict to assign color to gene
    if isinstance(colormap, dict):
        genesmd_df["color_tag"] = genesmd_df["color_tag"].astype(str)
        genesmd_df["color"] = genesmd_df["color_tag"].map(colormap)

        # add black genes warning if needed
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


def get_genes_metadata(df, id_col, color_col, packed, colormap, v_space):
    """Create genes metadata df."""

    # Start df with chromosome and the column defining color
    if color_col == "Chromosome":
        genesmd_df = df.groupby([id_col, "pr_ix"], group_keys=False, observed=True).agg(
            {"Chromosome": "first", "Start": "min", "End": "max"}
        )
        genesmd_df["chromosome"] = genesmd_df["Chromosome"]
        color_col = "chromosome"
    else:
        genesmd_df = (
            df.groupby([id_col, "pr_ix"], group_keys=False, observed=True)
            .agg(
                {
                    "Chromosome": "first",
                    "Start": "min",
                    "End": "max",
                    color_col: "first",
                }
            )
            .reset_index(level=1)
        )

    # Sort bu pr_ix
    genesmd_df.sort_values(by="pr_ix", inplace=True)

    genesmd_df["chrix"] = genesmd_df.groupby(
        "Chromosome", group_keys=False, observed=True
    ).ngroup()
    genesmd_df.rename(columns={color_col: "color_tag"}, inplace=True)
    genesmd_df["gene_ix_xchrom"] = genesmd_df.groupby(
        ["chrix", "pr_ix"], group_keys=False, observed=True
    ).cumcount()

    # Assign y-coordinate to genes
    if packed:
        genesmd_df["ycoord"] = -1
        genesmd_df = genesmd_df.groupby(
            ["chrix", "pr_ix"], group_keys=False, observed=True
        ).apply(_genesmd_packed)  # add packed ycoord column
        genesmd_df = genesmd_df.groupby("Chromosome").apply(_update_y)
        genesmd_df.reset_index(level="Chromosome", drop=True, inplace=True)

    else:
        # one gene in each height
        genesmd_df.set_index("pr_ix", inplace=True, append=True)
        genesmd_df["ycoord"] = (
            genesmd_df.sort_values(by="pr_ix", ascending=False)
            .groupby("Chromosome", group_keys=False, observed=True)
            .cumcount()
        )
        genesmd_df = genesmd_df.assign(
            upd_yc=genesmd_df.groupby("Chromosome", group_keys=False).apply(
                lambda x: x.groupby("pr_ix").ngroup(ascending=False)
            )
        )
        genesmd_df["ycoord"] += genesmd_df["upd_yc"]
        ### now create col to update according to prev pr height
        genesmd_df.reset_index("pr_ix", inplace=True)

    # Assign color to each gene
    genesmd_df = _genesmd_assigncolor(genesmd_df, colormap)  # adds a column 'color'

    # Add legend info
    def create_legend_rect(color):
        return Rectangle((0, 0), 1, 1, color=color)

    # column with rectangle objects with same color as gene
    genesmd_df["legend_item"] = genesmd_df["color"].apply(create_legend_rect)

    # Update vertical space
    ##genesmd_df["ycoord_1base"] = genesmd_df["ycoord"]*v_space
    ##genesmd_df["ycoord"] = genesmd_df["ycoord"] * v_space

    return genesmd_df


############ CHRMD_DF


##limits
def _chrmd_limits(chrmd_df, limits):
    """Compute 'min_max' column for chromosome metadata"""

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
        limits_chrmd_df = limits_df.groupby(
            "Chromosome", group_keys=False, observed=True
        ).agg({"Start": "min", "End": "max"})
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


def _fill_min_max(row, ts_data):
    """Complete min_max column for chromosome metadata if needed."""

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
            new_upper_lim = cumdelting([minmax_l[1]], ts_data, row.name[0], row.name[1])
            minmax_l[1] = new_upper_lim[0]

    # put plot coordinates in min_max
    row["min_max"] = minmax_l
    return row


def get_chromosome_metadata(
    df, id_col, limits, genesmd_df, packed, v_space, ts_data=None
):
    """Create chromosome metadata df."""

    # Start df
    chrmd_df = (
        df.groupby(["Chromosome", "pr_ix"], observed=True).agg(
            {"Start": "min", "End": "max", id_col: "nunique"}
        )
        # .reset_index(level="pr_ix")
    )

    chrmd_df.rename(
        columns={id_col: "n_genes", "Start": "min", "End": "max"}, inplace=True
    )
    # Adjust limits in case +1 pr
    if len(df["pr_ix"].drop_duplicates()) > 1:
        chrmd_df["min"] = chrmd_df.groupby(
            "Chromosome", group_keys=False, observed=True
        )["min"].transform("min")
        chrmd_df["max"] = chrmd_df.groupby(
            "Chromosome", group_keys=False, observed=True
        )["max"].transform("max")

    # Add limits
    _chrmd_limits(chrmd_df, limits)  # unknown limits are nan
    chrmd_df = chrmd_df.apply(lambda x: _fill_min_max(x, ts_data), axis=1)

    chrmd_df_grouped = (
        chrmd_df.reset_index(level="pr_ix")
        .groupby("Chromosome", group_keys=False, observed=True)
        .agg({"min": "first", "max": "first", "min_max": "first", "pr_ix": "size"})
    )
    chrmd_df_grouped.rename(columns={"pr_ix": "n_pr_ix"}, inplace=True)

    # Store plot y height
    if packed:
        chrmd_df_grouped = chrmd_df_grouped.join(
            genesmd_df.groupby(["Chromosome"], group_keys=False, observed=True)[
                "ycoord"
            ].max()
            * v_space
        )
        chrmd_df_grouped.rename(columns={"ycoord": "y_height"}, inplace=True)
        chrmd_df_grouped["y_height"] += 1 * v_space  # count from 1

    else:
        y_height_df = chrmd_df.groupby("Chromosome", observed=True)["n_genes"].sum()
        chrmd_df_grouped = chrmd_df_grouped.join(y_height_df, on="Chromosome")
        chrmd_df_grouped.rename(columns={"n_genes": "y_height"}, inplace=True)
        chrmd_df_grouped["y_height"] += chrmd_df_grouped["n_pr_ix"]
        chrmd_df_grouped["y_height"] -= 1
        chrmd_df_grouped["y_height"] = chrmd_df_grouped["y_height"] * v_space

    # Obtain the positions of lines separating pr objects
    chrmd_df = chrmd_df.join(
        genesmd_df.groupby(["Chromosome", "pr_ix"], group_keys=False, observed=True)[
            "ycoord"
        ].max()
    )
    chrmd_df.rename(columns={"ycoord": "pr_line"}, inplace=True)

    chrmd_df["pr_line"] = chrmd_df.groupby("Chromosome")["pr_line"].shift(
        -1, fill_value=-1
    )
    chrmd_df["pr_line"] += 1
    ##chrmd_df["pr_line"] = chrmd_df["pr_line"]*v_space

    # Set chrom_ix to get the right association to the plot index
    chrmd_df_grouped["chrom_ix"] = chrmd_df_grouped.groupby(
        "Chromosome", group_keys=False, observed=True
    ).ngroup()

    return chrmd_df, chrmd_df_grouped
