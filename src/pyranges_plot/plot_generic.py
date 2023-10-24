import pyranges
import plotly.colors
from .core import get_engine, get_idcol
from .plot_exons_plt.plot_exons_plt import plot_exons_plt
from .plot_exons_ply.plot_exons_ply import plot_exons_ply

colormap = plotly.colors.sequential.thermal

def plot(df, engine = None, max_ngenes = 25, id_col = None, color_col = None, colormap = colormap, 
		custom_coords = None, showinfo = None, disposition = 'packed', to_file = None):

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
    
    id_col: str, default 'gene_id'
        
        Name of the column containing gene ID.
    
    color_col: str, default None
    	
    	Name of the column used to color the genes.
    
    colormap: {matplotlib.colors.ListedColormap, list, str, dict}, default plotly.colors.sequential.thermal
    
    	Sequence of colors for the genes, it can be provided as a Matplotlib colormap, 
    	a Plotly color sequence (built as lists), a string naming the previously mentioned
    	color objects from Matplotlib and Plotly, or a dictionary with the following 
    	structure {color_column_value1: color1, color_column_value2: color2, ...}. When a specific
    	color_col value is not specified in the dictionary it will be colored in black.
    	
    custom_coords: {None, dict}, default None
    
    	Customization of coordinates for the chromosome plots. As default it is defined as 
    	the minimun and maximum exon coordinate plotted plus a 5% of the range on each side.
    	It also accepts a dictionary with the following strunctue {chr_name1: (min_coord, max coord), 
    	chr_name2: (min_coord, max_coord), ...}. Note that in the dictionary not all the chromosomes 
    	have to be present and some coordinates can be indicated as None leading to the use of the 
    	default value.
    	
    showinfo: list, default None
    
    	Dataframe information to show when placing the mouse over a gene. This must be provided as a list 
    	of column names. By default it shows the ID of the gene followed by its start and end position.
    
    disposition: str, default 'packed'
    
    	Select wether the genes should be presented in full display (one row each) using the 'full' option,
    	or if they should be presented in a packed (in the same line if they do not overlap) using 'packed'.
    	
    to_file: str, default None
    
    	Name of the file to export specifying the desired extension. The supported extensions are '.png' and '.pdf'.
    	
    Examples
    --------
    
    >>> plot_generic(df, engine='plt', max_ngenes=25, colormap='Set3')
    
    >>> plot_generic(df, engine='matplotlib', color_col='Strand', colormap={'+': 'green', '-': 'red'})
    
    >>> plot_generic(df, engine='ply', custom_coords = {'1': (1000, 50000), '2': None, '3': (10000, None)})
    
    >>> plot_generic(df, engine='plotly', colormap=plt.get_cmap('Dark2'), showinfo = ['feature1', 'feature3'])
    
    >>> plot_generic(df, engine='plt', color_col='Strand', disposition='full', to_file='my_plot.pdf')
    	

    """
    
    # Get dataframe if the provided object is pyranges
    if type(df) is pyranges.pyranges_main.PyRanges:
        df = df.df
    
    
    # Deal with export
    if to_file:
        ext = to_file[-4:]
        try:
            if ext not in ['.pdf', '.png']:
                raise Exception("Please specify the desired format to export the file including either '.png' or '.pdf' as an extension.")
        except SystemExit as e:
            print("An error occured:", e)

    
    # Deal with id column
    if id_col is None:
        id_col = get_idcol()
    
    try:
        if id_col is None or id_col not in df.columns:
            raise Exception("Please define the name of the ID column using either set_idcol() function or plot_generic parameter as plot_generic(..., id_col = 'your_id_col')")
    except SystemExit as e:
        print("An error occured:", e)
        
    
    # Deal with engine
    if engine is None:
        engine = get_engine()
    
    try:
    	if engine == 'plt' or engine == 'matplotlib':
    	    plot_exons_plt(df, max_ngenes=max_ngenes, id_col = id_col, color_col = color_col, colormap = colormap, 
    	    		custom_coords = custom_coords, showinfo = showinfo, disposition = disposition, to_file = to_file)
    	elif engine == 'ply' or engine == 'plotly':
    	    plot_exons_ply(df, max_ngenes=max_ngenes, id_col = id_col, color_col = color_col, colormap = colormap, 
    	    		custom_coords = custom_coords, showinfo = showinfo, disposition = disposition, to_file = to_file)
    	else:
            raise Exception("Please define engine with set_engine().")
    except SystemExit as e:
        print("An error occured:", e)


