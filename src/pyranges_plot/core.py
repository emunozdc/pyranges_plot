import pandas as pd
from .plot_features import plot_features_dict, plot_features_dict_in_use


# CORE FUNCTIONS
engine = None


def set_engine(name):
    """
    Defines the engine for the plots

    Parameters
    ----------
    name: str
        Indicates if Matplotlib ('plt', 'matplotlib') or Plotly ('ply', 'plotly') should be used.

    Examples
    --------
    >>> import pyranges_plot as prp

    >>> prp.set_engine('plt')

    """

    global engine
    engine = name


def get_engine():
    """Shows the current defined engine."""

    return engine


id_col = None


def set_idcol(name):
    """
    Defines the ID column for the data.

    Parameters
    ----------
    name: str

         Indicates the name of the ID column to be used when dealing with data.

    Examples
    --------
    >>> import pyranges_plot as prp

    >>> prp.set_idcol('gene_id')

    """

    global id_col
    id_col = name


def get_idcol():
    """Shows the current defined ID column (id_col)."""

    return id_col


warnings = True


def set_warnings(option):
    """
    Defines if the warnings should be shown.

    Parameters
    ----------
    option: bool

         True for showing the warnings and false to turn them off.

    Examples
    --------
    >>> import pyranges_plot as prp

    >>> prp.set_warnings(False)

    """

    global warnings
    warnings = option


def get_warnings():
    """Returns the current warnings state."""

    return warnings


# Related to default features


def set_default(varname, value):
    """
    Define some features of the plot layout.

    Parameters
    ----------
    varname: str

        Name of the variable to change.

    value:

        New value of the variable to be assigned.

    Examples
    --------
    >>> import pyranges_plot as prp

    >>> prp.set_default('plot_background', 'magenta')

    >>> prp.set_default('title_size', 20)

    """

    plot_features_dict_in_use[varname] = (
        value,
        plot_features_dict[varname][1],
        "*",
    )  # (value, description, modified tag)


def get_default(varname="all"):
    """
    Obtain the deafault value for a plot layout variable/s and its description.

    Parameters
    ----------
    varname: {str, list}, default 'all'

        Name of the variable/s to get the value and description.

    """

    # list of variables
    if isinstance(varname, list):
        vars_list = []
        for var in varname:
            vars_list.append(plot_features_dict_in_use[var][0])
        return vars_list

    # all variables
    elif varname == "all":
        return plot_features_dict_in_use

    # one variable
    else:
        try:
            if varname in plot_features_dict_in_use:
                return plot_features_dict_in_use[varname][0]
            else:
                raise Exception(
                    f"The variable you provided is not customizable. The customizable variables are: {list(plot_features_dict.keys())}"
                )
        except SystemExit as e:
            print("An error occured:", e)


def get_original_default():
    """Returns the dictionary with the original plot features."""

    return plot_features_dict


def reset_default(varname="all"):
    """
    Reset the deafault value for one, some or all plot layout variables to their original vlaue.

    Parameters
    ----------
    varname: {str, list}, default 'all'

        Name of the variable/ to reset the value.

    Examples
    --------
    >>> import pyranges_plot as prp

    >>> prp.reset_default()

    >>> prp.reset_default('all')

    >>> prp.reset_default('tag_background')

    >>> prp.reset_default(['title_dict_plt', 'tag_background'])

    >>> prp.reset_default('title_dict_ply')
    """

    plot_features_dict_in_use = get_default()
    plot_features_dict = get_original_default()

    # list of variables
    if isinstance(varname, list):
        for var in varname:
            plot_features_dict_in_use[var] = plot_features_dict[var]

    # all variables
    elif varname == "all":
        for var in plot_features_dict_in_use.keys():
            plot_features_dict_in_use[var] = plot_features_dict[var]

    # one variable
    else:
        try:
            if varname in plot_features_dict_in_use.keys():
                plot_features_dict_in_use[varname] = plot_features_dict[varname]
            else:
                raise Exception(
                    f"The variable you provided is not customizable. The customizable variables are: {list(plot_features_dict.keys())}"
                )
        except SystemExit as e:
            print("An error occured:", e)


def divide_desc(desc, cutoff):
    """Divide long feature description in lines."""

    lines_l = []
    while len(desc) > cutoff:
        for i in range(59, -1, -1):
            if desc[i] == " ":
                lines_l.append(desc[:i])
                desc = desc[i + 1 :]
                break
    lines_l.append(desc)

    return lines_l


def print_default(return_keys=False):
    """xxx"""

    # store data
    plot_features_dict_in_use = get_default()

    # prepare data to print
    if not return_keys:
        feat_df = pd.DataFrame.from_dict(
            plot_features_dict_in_use,
            orient="index",
            columns=["Value", "Description", "Modified"],
        )

        # Calculate column sizes
        name_sz = max([len(val) for val in plot_features_dict_in_use])
        value_sz = max([len(str(val)) for val in feat_df["Value"]])
        if value_sz < 5:  # value has a minimum of 5
            value_sz = 5
        mod_sz = 7  # according to "Edited?" length
        desc_sz = 60

        # Function to format row
        def format_row(key, value):
            if len(value[1]) <= 60:
                return f"| {key:^{name_sz}} | {str(value[0]):^{value_sz}} | {value[2]:^{mod_sz}} | {value[1]:<{desc_sz}} |"

            else:
                lines_l = divide_desc(value[1], cutoff=desc_sz)
                fstr = f"| {key:^{name_sz}} | {str(value[0]):^{value_sz}} | {value[2]:^{mod_sz}} | {lines_l[0]:<{desc_sz}} |"
                empty = " "
                for i in range(1, len(lines_l)):
                    fstr += f"\n| {empty:^{name_sz}} | {empty:^{value_sz}} | {empty:^{mod_sz}} | {lines_l[i]:<{desc_sz}} |"

                return fstr

        # Create table header
        header = f"+{'-' * (name_sz+2)}+{'-' * (value_sz+2)}+{'-' * (mod_sz+2)}+{'-' * (desc_sz+2)}+\n"
        header += f"| {'Feature':^{name_sz}} | {'Value':^{value_sz}} | {'Edited?':^{mod_sz}} | {'Description':^{desc_sz}} |\n"
        header += f"+{'-' * (name_sz+2)}+{'-' * (value_sz+2)}+{'-' * (mod_sz+2)}+{'-' * (desc_sz+2)}+"

        # Divide features
        print(feat_df)
        extragen_feat_df = feat_df[
            feat_df.index.isin(
                [
                    "tag_background",
                    "plot_background",
                    "plot_border",
                    "title_size",
                    "title_color",
                ]
            )
        ].copy()
        intragen_feat_df = feat_df[
            feat_df.index.isin(
                [
                    "exon_width",
                    "arrow_line_width",
                    "arrow_color",
                    "arrow_size",
                    "arrow_intron_threshold",
                ]
            )
        ].copy()
        other_feat_df = feat_df[
            feat_df.index.isin(["shrink_threshold", "plotly_port"])
        ].copy()

        # Create table rows
        rows_eg = "\n".join(
            [format_row(key, value) for key, value in extragen_feat_df.iterrows()]
        )
        rows_ig = "\n".join(
            [format_row(key, value) for key, value in intragen_feat_df.iterrows()]
        )
        rows_o = "\n".join(
            [format_row(key, value) for key, value in other_feat_df.iterrows()]
        )

        # Print table
        print(header)
        print(rows_eg)
        print(
            f"+{'-' * (name_sz+2)}+{'-' * (value_sz+2)}+{'-' * (mod_sz+2)}+{'-' * (desc_sz+2)}+"
        )
        print(rows_ig)
        print(
            f"+{'-' * (name_sz + 2)}+{'-' * (value_sz + 2)}+{'-' * (mod_sz + 2)}+{'-' * (desc_sz + 2)}+"
        )
        print(rows_o)
        print(
            f"+{'-' * (name_sz + 2)}+{'-' * (value_sz + 2)}+{'-' * (mod_sz + 2)}+{'-' * (desc_sz + 2)}+"
        )

    if return_keys:
        return set(plot_features_dict_in_use.keys())


def cumdelting(num_l, ts_data, chrom):
    """Update a list of numbers according to cumdelta."""

    for i in range(len(num_l)):
        cdel = 0
        # get proper cumdelta
        for ix, row in ts_data[chrom].iterrows():
            if row["End"] <= num_l[i]:
                cdel = row["cumdelta"]
            else:
                break
        num_l[i] -= cdel

    return num_l
