import numpy as np
import pandas as pd
from .plot_features import plot_features_dict, plot_features_dict_in_use
import tkinter as tk

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


def coord2inches(fig, ax, X0, X1, Y0, Y1):
    """Provides the inches length from the points given. Plt friendly"""

    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    fig_width, fig_height = fig.get_size_inches()

    x_scale = fig_width / (x_max - x_min)
    y_scale = fig_height / (y_max - y_min)

    inch_len = float(np.sqrt(((X1 - X0) * x_scale) ** 2 + ((Y1 - Y0) * y_scale) ** 2))

    return inch_len


def coord2percent(fig, trace, X0, X1):
    """Provides the plot percentage length from the points given. Plotly friendly"""

    if trace == 1:
        trace = ""
    else:
        trace = str(trace)

    x_min, x_max = fig["layout"]["xaxis" + trace]["range"]
    x_rang = x_max - x_min

    percent_size = float(((X1 - X0) / x_rang))

    return percent_size


def inches2coord(fig, ax, x_inches):
    """Provides the coordinates distance from the inches given. Plt friendly"""

    x_min, x_max = ax.get_xlim()
    fig_width, fig_height = fig.get_size_inches()
    x_scale = fig_width / (x_max - x_min)

    cord_size = float(x_inches / x_scale)

    return cord_size


def percent2coord(fig, trace, x_percent):
    """Provides the coordinates distance from the plot percentage given. Plt friendly"""

    if trace == 1:
        trace = ""
    else:
        trace = str(trace)

    x_min, x_max = fig["layout"]["xaxis" + trace]["range"]
    x_rang = x_max - x_min

    percent_coord = float(x_percent * x_rang)

    return percent_coord


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


def print_default(return_keys=False):
    """xxx"""

    # prepare data to print
    plot_features_dict_in_use = get_default()
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
    mod_sz = 8  # according to "Modified" length
    desc_sz = max([len(val) for val in feat_df["Description"]])

    # Function to format row
    def format_row(key, value):
        name_sz = max([len(val) for val in plot_features_dict_in_use])
        value_sz = max([len(str(val)) for val in feat_df["Value"]])
        if value_sz < 5:  # value has a minimum of 5
            value_sz = 5
        mod_sz = 8  # according to "Modified" length
        desc_sz = max([len(val) for val in feat_df["Description"]])

        return f"| {key:^{name_sz}} | {str(value[0]):^{value_sz}} | {value[2]:^{mod_sz}} | {value[1]:<{desc_sz}} |"

    # Create table header
    header = f"+{'-' * (name_sz+2)}+{'-' * (value_sz+2)}+{'-' * (mod_sz+2)}+{'-' * (desc_sz+2)}+\n"
    header += f"| {'Feature':^{name_sz}} | {'Value':^{value_sz}} | {'Modified':^{mod_sz}} | {'Description':^{desc_sz}} |\n"
    header += f"+{'-' * (name_sz+2)}+{'-' * (value_sz+2)}+{'-' * (mod_sz+2)}+{'-' * (desc_sz+2)}+"

    # Create table rows
    rows = "\n".join([format_row(key, value) for key, value in feat_df.iterrows()])

    # Print table
    print(header)
    print(rows)
    print(
        f"+{'-' * (name_sz+2)}+{'-' * (value_sz+2)}+{'-' * (mod_sz+2)}+{'-' * (desc_sz+2)}+"
    )

    if return_keys:
        return set(plot_features_dict_in_use.keys())


def plt_popup_warning(txt, bkg="#1f1f1f", txtcol="white", botcol="#D6AA00"):
    """Create warning window for Matplotlib plots."""

    warn = tk.Tk()

    # Title and background
    warn.wm_title("Warning!")
    warn.configure(background=bkg)

    # Label for warning text
    label = tk.Label(warn, text=txt, font=("Sans", 15), fg=txtcol, bg=bkg)
    label.pack(side="top", anchor="center", pady=10)

    # Button
    bot = tk.Button(
        warn, text="Got it", command=lambda: warn.destroy(), fg="black", bg=botcol
    )
    bot.pack(pady=10)

    # Start main loop
    warn.mainloop()
