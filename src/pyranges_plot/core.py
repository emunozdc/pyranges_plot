import numpy as np
import matplotlib.pyplot as plt
import plotly.subplots as sp
import matplotlib.cm as cm
import plotly.colors as pc
from intervaltree import Interval, IntervalTree
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
    >>> pyranges_plot.set_engine('plt')
    
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
    >>> pyranges_plot.set_idcol('gene_id')
    
    """

    global id_col
    id_col = name


def get_idcol():
    """Shows the current defined ID column (id_col)."""

    return id_col





def packed_for_genesmd(genesmd_df):
    """xxx"""

    # Sort the dataframe by Start values
    genesmd_df = genesmd_df.sort_values(by='Start')

    # Initialize IntervalTree and used y-coordinates list
    trees = [IntervalTree()]

    def find_tree(row):
        for tree in trees:
            if not tree.overlaps(row['Start'], row['End']):
                return tree
        trees.append(IntervalTree())
        return trees[-1]

    # Assign y-coordinates
    for idx, row in genesmd_df.iterrows():
        tree = find_tree(row)
        tree.addi(row['Start'], row['End'], idx)
        genesmd_df.at[idx, 'ycoord'] = trees.index(tree)

    return genesmd_df


def coord2inches(fig, ax, X0, X1, Y0, Y1):
    """Provides the inches length from the points given. Plt friendly"""
    
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    fig_width, fig_height = fig.get_size_inches()
    
    x_scale = fig_width / (x_max - x_min)
    y_scale = fig_height / (y_max - y_min)
    
    inch_len = float(np.sqrt( ((X1 - X0) * x_scale)**2 + ((Y1 - Y0) * y_scale)**2 ))
    
    return inch_len

def coord2percent(fig, trace, X0,X1):
    """Provides the plot percentage length from the points given. Plotly friendly"""

    if trace == 1:
        trace = ""
    else:
        trace = str(trace)
        
    x_min, x_max = fig['layout']["xaxis"+trace]['range']
    x_rang = x_max-x_min
    
    percent_size = float(((X1-X0)/x_rang) )
    
    return percent_size
    

def inches2coord(fig, ax, x_inches):
    """Provides the coordinates distance from the inches given. Plt friendly"""
    
    x_min, x_max = ax.get_xlim()
    fig_width, fig_height = fig.get_size_inches()
    x_scale = fig_width / (x_max - x_min)
    
    cord_size = float(x_inches/x_scale)
    
    return cord_size 


def percent2coord(fig, trace, x_percent):
    """Provides the coordinates distance from the plot percentage given. Plt friendly"""

    if trace == 1:
        trace = ""
    else:
        trace = str(trace)
        
    x_min, x_max = fig['layout']["xaxis"+trace]['range']
    x_rang = x_max-x_min
    
    percent_coord = float(x_percent * x_rang)
    
    return percent_coord





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



def on_hover_factory(fig, annotation, object, geneinfo):
    def on_hover(event):
        """Check if the mouse is over the object and show annotation if so."""
        visible = annotation.get_visible()
        contains_object = object.contains(event)  # Check if the mouse is over the object

        if contains_object:
            annotation.set_text(geneinfo)  # Add the information
            annotation.xy = (event.xdata, event.ydata)
            annotation.set_visible(True)
            fig.canvas.draw()
        elif visible:
            annotation.set_visible(False)
            fig.canvas.draw()

    return on_hover



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
    >>> pyranges_plot.set_default('plot_background', 'magenta')
    
    >>> pyranges_plot.set_default('title_dict_ply.size', 20)
    
    >>> pyranges_plot.set_default('title_dict_plt.color', 'red')
    
    """

    if '.' in varname:
        dictname = varname.split('.')[0]
        keyname = varname.split('.')[1]
        plot_features_dict_in_use[dictname][0][keyname] = value

    else:
        plot_features_dict_in_use[varname] = (value, plot_features_dict[varname][1])



def get_default(varname = 'all'):
    """
    Obtain the deafault value for a plot layout variable/s and its description.
    
    Parameters
    ----------
    varname: {str, list}, default 'all'
        
        Name of the variable/s to get the value and description.
    
    Examples
    --------
    >>> pyranges_plot.get_default()
    
    >>> pyranges_plot.get_default('all')
    
    >>> pyranges_plot.get_default('plot_border')
    
    >>> pyranges_plot.get_default(['title_dict_plt', 'plot_background'])
    
    >>> pyranges_plot.get_default('title_dict_ply')[0][color]
    
    """

    # list of variables
    if type(varname) is list:
        vars_dict = {}
        for var in varname:
            vars_dict[var] = plot_features_dict_in_use[var]
        return vars_dict
        
    # all variables
    elif varname == "all":
        return plot_features_dict_in_use
    
    # one variable
    else:
        try:
            if varname in plot_features_dict_in_use:
                return plot_features_dict_in_use[varname]
            else:
                raise Exception(f"The variable you provided is not customizable. The customizable variables are: {list(plot_features_dict.keys())}")
        except SystemExit as e:
            print("An error occured:", e) 
    


def get_original_default():
   """Returns the dictionary with the original plot features."""
   
   return plot_features_dict

 
    
def reset_default(varname = "all"):
    """
    Reset the deafault value for one, some or all plot layout variables to their original vlaue.
    
    Parameters
    ----------
    varname: {str, list}, default 'all'
        
        Name of the variable/ to reset the value.
    
    Examples
    --------
    >>> pyranges_plot.reset_default()
    
    >>> pyranges_plot.reset_default('all')
    
    >>> pyranges_plot.reset_default('tag_background')
    
    >>> pyranges_plot.reset_default(['title_dict_plt', 'tag_background'])
    
    >>> pyranges_plot.reset_default('title_dict_ply')
    """
    
    plot_features_dict_in_use = get_default()
    plot_features_dict = get_original_default()
    
    # list of variables
    if type(varname) is list:
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
                raise Exception(f"The variable you provided is not customizable. The customizable variables are: {list(plot_features_dict.keys())}")
        except SystemExit as e:
            print("An error occured:", e)
