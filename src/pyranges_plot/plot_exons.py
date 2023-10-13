import pyranges
import plotly.colors
from .core import get_engine
from .plot_exons_plt.plot_exons_plt import plot_exons_plt
from .plot_exons_ply.plot_exons_ply import plot_exons_ply

colormap = plotly.colors.sequential.thermal

def plot_exons(df, engine = None, max_ngenes = 25, id_column = 'gene_id', color_column = None, colormap = colormap, custom_coords = None):

    """
    Create genes plot from PyRanges object DataFrame
    
    Parameters
    ----------
    df: {pyranges.pyranges_main.PyRanges, pandas.DataFrame}
    	
    	Pyranges or derived dataframe with genes' data.
    	
    engine: str, default None
    
    	Library in which the plot sould be built, it accepts either Matplotlib ['matplotlib'/'plt'] or 
    	Plotly ['ply'/'plotly']
    	    
    max_ngenes: int, default 20
    	
    	Maximum number of genes plotted in the dataframe order.
    
    id_column: str, default 'gene_id'
        
        Name of the column containing gene ID.
    
    color_column: str, default None
    	
    	Name of the column used to color the genes.
    
    colormap: {matplotlib.colors.ListedColormap, list, str, dict}, default plotly.colors.sequential.thermal
    
    	Sequence of colors for the genes, it can be provided as a Matplotlib colormap, 
    	a Plotly color sequence (built as lists), a string naming the previously mentioned
    	color objects from Matplotlib and Plotly, or a dictionary with the following 
    	structure {color_column_value1: color1, color_column_value2: color2, ...}. When a specific
    	color_column value is not specified in the dictionary it will be colored in black.
    	
    custom_coords: {None, dict}, default None
    
    	Customization of coordinates for the chromosome plots. As default it is defined as 
    	the minimun and maximum exon coordinate plotted plus a 5% of the range on each side.
    	It also accepts a dictionary with the following strunctue {chr_name1: (min_coord, max coord), 
    	chr_name2: (min_coord, max_coord), ...}. Note that in the dictionary not all the chromosomes 
    	have to be present and some coordinates can be indicated as None leading to the use of the 
    	default value.
    	
    Examples
    --------
    
    >>> plot_exons(df, engine='plt', max_ngenes=25, colormap='Set3')
    
    >>> plot_exons(df, engine='matplotlib', color_column='Strand', colormap={'+': 'green', '-': 'red'})
    
    >>> plot_exons(df, engine='ply', custom_coords = {'1': (1000, 50000), '2': None, '3': (10000, None)})
    
    >>> plot_exons(df, engine='plotly', colormap=plt.get_cmap('Dark2'))
    	

    """
    
    # Get dataframe if the provided object is pyranges
    if type(df) == pyranges.pyranges_main.PyRanges:
        df = df.df
    
    
    # Deal with column id
    try:
        if id_column not in df.columns:
            raise Exception("Please specify the name of the ID column with plot_exons(..., id_column = 'your_id_column')")
    except SystemExit as e:
        print("An error occured:", e)
        
    
    # Deal with engine
    if engine is None:
        engine = get_engine()
    
    try:
    	if engine == 'plt' or engine == 'matplotlib':
    	    plot_exons_plt(df, max_ngenes=max_ngenes, id_column = id_column, color_column = color_column, colormap = colormap, custom_coords = custom_coords)
    	elif engine == 'ply' or engine == 'plotly':
    	    plot_exons_ply(df, max_ngenes=max_ngenes, id_column = id_column, color_column = color_column, colormap = colormap, custom_coords = custom_coords)
    	else:
            raise Exception("Please define engine with set_engine().")
    except SystemExit as e:
        print("An error occured:", e)




        
        
        
        
