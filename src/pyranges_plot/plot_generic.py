import pyranges
import plotly.colors
from .core import get_engine, get_idcol
from .plot_exons_plt.plot_exons_plt import plot_exons_plt
from .plot_exons_plt.plot_transcript_plt import plot_transcript_plt
from .plot_exons_ply.plot_exons_ply import plot_exons_ply
from .plot_exons_ply.plot_transcript_ply import plot_transcript_ply
from dash import Dash, dcc, html, Input, Output
import dash_bootstrap_components as dbc


colormap = plotly.colors.qualitative.Alphabet


def plot(df, engine = None, max_ngenes = 25, id_col = None, transcript_str = False, color_col = None, colormap = colormap, 
		limits = None, showinfo = None, chr_string = None, packed = True, to_file = None, file_size = None,
		**kargs):

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
        
    transcript_str: bool, default False
    
        Display differentially transcript regions belonging and not belonging to CDS. The CDS/exon information
        must be stored in the 'Feature' column of the PyRanges object or the dataframe.
    
    color_col: str, default None
    	
    	Name of the column used to color the genes.
    
    colormap: {matplotlib.colors.ListedColormap, list, str, dict}, default plotly.colors.sequential.thermal
    
    	Sequence of colors for the genes, it can be provided as a Matplotlib colormap, 
    	a Plotly color sequence (built as lists), a string naming the previously mentioned
    	color objects from Matplotlib and Plotly, or a dictionary with the following 
    	structure {color_column_value1: color1, color_column_value2: color2, ...}. When a specific
    	color_col value is not specified in the dictionary it will be colored in black.
    	
    limits: {None, dict, tuple, pyranges.pyranges_main.PyRanges}, default None
    
    	Customization of coordinates for the chromosome plots. 
    	- None: minimun and maximum exon coordinate plotted plus a 5% of the range on each side.
    	- dict: {chr_name1: (min_coord, max coord), chr_name2: (min_coord, max_coord), ...}. Not 
    	all the plotted chromosomes need to be specified in the dictionary and some coordinates 
    	can be indicated as None, both cases lead to the use of the default value. 
    	- tuple: the coordinate limits of all chromosomes will be defined as indicated.
    	- pyranges.pyranges_main.PyRanges: for each matching chromosome between the plotted data 
    	and the limits data, the limits will be defined by the minimum and maximum coordinates 
    	in the pyranges object defined as limits. If some plotted chromosomes are not present they 
    	will be left as default.
    	    	
    showinfo: list, default None
    
    	Dataframe information to show when placing the mouse over a gene. This must be provided as a list 
    	of column names. By default it shows the ID of the gene followed by its start and end position.
    	
    chr_string: string, default 
    
    	String naming the chromosome used as titles for the plots. An f string is used, where the variable 
    	chrom corresponds to the PyRanges object value for chromosome. Must contain chrom variable.
    
    packed: bool, default True
    
    	Disposition of the genes in the plot. Use True for a packed disposition (genes in the same line if
        they do not overlap) and False for unpacked (one row per gene).
    	
    to_file: str, default None
    
    	Name of the file to export specifying the desired extension. The supported extensions are '.png' and '.pdf'.
    	
    file_size: {list, tuple}, default None
    
    	Size of the plot to export defined by a sequence object like: (height, width). The default values 
    	make the height according to the number of genes and the width as 20 in Matplotlib and 1600 in Plotly.
    	
    **kargs: 
    
    	Customizable plot features can be defined using kargs. Use print_default() function to check the variables'
    	nomenclature, description and default values.

    
    	
    Examples
    --------
    
    >>> plot_generic(df, engine='plt', max_ngenes=25, colormap='Set3')
    
    >>> plot_generic(df, engine='matplotlib', color_col='Strand', colormap={'+': 'green', '-': 'red'})
    
    >>> plot_generic(df, engine='ply', limits = {'1': (1000, 50000), '2': None, '3': (10000, None)})
    
    >>> plot_generic(df, engine='plotly', colormap=plt.get_cmap('Dark2'), showinfo = ['feature1', 'feature3'])
    
    >>> plot_generic(df, engine='plt', color_col='Strand', packed='False', to_file='my_plot.pdf', background)
    	

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
        
    
    # Deal with transcript structure
    if transcript_str:
        try:
            if 'Feature' not in df.columns:
                raise Exception("The transcript structure information must be stored in 'Feature' column of the data.")
        except SystemExit as e:
            print("An error occured:", e)
    
    
    # Deal with engine
    if engine is None:
        engine = get_engine()
    
    try:  #call plot functions
    	kargs_l = [k for k in kargs]
    	if engine == 'plt' or engine == 'matplotlib':
    	    if not transcript_str:
    	        plot_exons_plt(df, max_ngenes=max_ngenes, id_col = id_col, color_col = color_col, 
    	                       colormap = colormap, limits = limits, showinfo = showinfo, chr_string = chr_string, packed = packed, 
    	                       to_file = to_file, file_size = file_size, **kargs)
    	    else:
    	        plot_transcript_plt(df, max_ngenes=max_ngenes, id_col = id_col, color_col = color_col, 
    	                            colormap = colormap, limits = limits, showinfo = showinfo, chr_string = chr_string, packed = packed, 
    	                            to_file = to_file, file_size = file_size, **kargs)                       
    	         
    	                       
    	elif engine == 'ply' or engine == 'plotly':
    	    if not transcript_str:
    	        #store figure
    	        fig = plot_exons_ply(df, max_ngenes=max_ngenes, id_col = id_col, color_col = color_col, 
    	                       colormap = colormap, limits = limits, showinfo = showinfo, chr_string = chr_string, packed = packed, 
    	                       to_file = to_file, file_size = file_size, **kargs)
    	        
    	        #initialize and run dash app
    	        app_instance = initialize_dash_app(fig, max_ngenes)
    	        app_instance.run_server()
           
    	    else:
    	        #store figure
    	        fig = plot_transcript_ply(df, max_ngenes=max_ngenes, id_col = id_col, color_col = color_col, 
    	                            colormap = colormap, limits = limits, showinfo = showinfo, chr_string = chr_string, packed = packed, 
    	                            to_file = to_file, file_size = file_size, **kargs)
    	        
    	        #initialize and run dash app
    	        app_instance = initialize_dash_app(fig, max_ngenes)
    	        app_instance.run_server()
    	                       
    	else:
            raise Exception("Please define engine with set_engine().")
    except SystemExit as e:
        print("An error occured:", e)




# Plotly - Function to initialize Dash app layout and callbacks
def initialize_dash_app(fig, max_ngenes):
    app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

    # create alert and graph components
    subdf_alert = dbc.Alert("The provided data contains more genes than the ones plotted.",
                            id="alert-subset", color="warning", dismissable=True, is_open=False)
    gr = dcc.Graph(id="genes-plot", figure=fig)

    # define layout
    app.layout = html.Div([
        dbc.Row([subdf_alert, gr], justify="around")
    ])

    # callback function
    @app.callback(
        Output("alert-subset", "is_open"), 
        Input("genes-plot", "figure"), 
        )
        
    def show_warning(grfig):
        n_genes = int(grfig['data'][0]['customdata'][0])
        print(n_genes)
        if n_genes > max_ngenes:
            return True

    return app
