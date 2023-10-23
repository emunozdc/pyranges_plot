import numpy as np
import matplotlib.pyplot as plt
import plotly.subplots as sp
import matplotlib.cm as cm
import plotly.colors as pc
from intervaltree import Interval, IntervalTree

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
