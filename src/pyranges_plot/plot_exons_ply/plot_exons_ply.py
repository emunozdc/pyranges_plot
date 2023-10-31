import pyranges as pr
import plotly.express as px
import plotly.graph_objects as go
import plotly.subplots as sp
import plotly.colors
import plotly.io as pio
import matplotlib.pyplot as plt  ## priorizar no depender de las 2?
import matplotlib.colors as mcolors
import numpy as np
from ..core import coord2percent, percent2coord, is_pltcolormap, is_plycolormap, get_plycolormap, packed_for_genesmd, get_default
import dash
from dash import dcc, html, Output, Input



# plot parameters
exon_width = 0.4
transcript_utr_width = 0.2 * exon_width
colormap = plotly.colors.sequential.thermal
arrow_width = 1
arrow_color = "grey"
arrow_size_max = 0.02
arrow_size_min = 0.005
intron_threshold = 0.03




# PLOT_EXONS FUNCTIONS 

def plot_exons_ply(df, max_ngenes = 25, id_col = 'gene_id', transcript_str = False, color_col = None, colormap = colormap, 
		limits = None, showinfo = None, packed = True, to_file = None, file_size = None):
    """
    Create genes plot from PyRanges object DataFrame
    
    Parameters
    ----------
    df: pandas.DataFrame 
    	
    	Dataframe with genes.
    	    
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
    	structure {color_column_value1: color1, color_column_value2: color2, ...}
    	
    limits: {None, dict, tuple}, default None
    
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

    packed: bool, default True
    
    	Disposition of the genes in the plot. Use True for a packed disposition (genes in the same line if
        they do not overlap) and False for unpacked (one row per gene).
    	
    to_file: str, default None
    
    	Name of the file to export specifying the desired extension. The supported extensions are '.png' and '.pdf'.
    	
    file_size: {list, tuple}, default None
    
    	Size of the plot to export defined by a sequence object like: (height, width). The default values 
    	make the height according to the number of genes and the width as 1600.
    	
    Examples
    --------
    
    >>> plot_exons_ply(df, max_ngenes=25, colormap='Set3')
    
    >>> plot_exons_ply(df, color_col='Strand', colormap={'+': 'green', '-': 'red'})
    
    >>> plot_exons_ply(df, limits = {'1': (1000, 50000), '2': None, '3': (10000, None)})
    	

    """
    
    
    # Get default plot features
    tag_background = get_default('tag_background')[0]
    plot_background = get_default('plot_background')[0]
    plot_border = get_default('plot_border')[0]
    title_dict_ply = get_default('title_dict_ply')[0]

    
    # Make DataFrame subset if needed
    # create a cloumn indexing all the genes in the df
    genesix_l = [i for i in enumerate(df[id_col].drop_duplicates())]
    genesix_d = {}
    for ix, gene in genesix_l:
        genesix_d[gene] = ix
    df["gene_index"] = df[id_col].map(genesix_d) 
    
    # select maximun number of genes
    if max(df.gene_index)+1 <= max_ngenes:
        subdf = df
    else:
        subdf = df[df.gene_index < max_ngenes]
    
    # remove the gene_index column from the original df
    df.drop('gene_index', axis=1, inplace=True)
    
    
    # Create chromosome metadata DataFrame
    chrmd_df = subdf.groupby("Chromosome").agg({'Start': 'min', 'End': 'max', id_col: 'nunique'})
    chrmd_df.dropna(inplace=True) # remove chr not present in subset
    chrmd_df.rename(columns={id_col: 'n_genes',
                            'Start': 'min',
                            'End': 'max'}, 
                   inplace=True)
    nchrs = len(chrmd_df)
    
    # consider custom coordinates limits
    # 1- create min_max column containing (plot min, plot max)
    
    #no limits no info
    if limits is None:
    	chrmd_df['min_max'] = [(np.nan, np.nan)] * len(chrmd_df)
    
    #tuple for all chromosomes
    elif type(limits) is tuple:
        chrmd_df['min_max'] = [limits] * len(chrmd_df)
    
    #pyranges object
    elif type(limits) is pr.pyranges_main.PyRanges:
        #create dict to map limits
        limits_df = limits.df
        limits_chrmd_df = limits_df.groupby("Chromosome").agg({'Start': 'min', 'End': 'max'})
        limits_chrmd_dict = limits_chrmd_df.to_dict(orient='index')
        
        #function to get matching values from limits_chrmd_df
        def make_min_max(row):
            chromosome = str(row.name)
            limits = limits_chrmd_dict.get(chromosome)
            if limits:
                return (limits['Start'], limits['End'])  #chromosome in both sets of data
            else:
                return (np.nan, np.nan)  #chromosome does not match
        
        #create limits column in plotting data
        chrmd_df['min_max'] = chrmd_df.apply(make_min_max, axis=1)
    
    #dictionary as limits
    else:
    	chrmd_df['min_max'] = [limits_chrmd_df.get(index) for index in chrmd_df.index] # fills with None the chromosomes not specified
    
    # 2- not specified values are (np.nan, np.nan), get default from data
    def fill_min_max(row):
    	minmax_t = row['min_max']
    	#deal with empty rows
    	if minmax_t is None: 
    	    minmax_t = (np.nan, np.nan)
    
    	#check both items and put default if necessary
    	minmax_l = list(minmax_t) 
    	if minmax_l[0] is None or np.isnan(minmax_l[0]):
            minmax_l[0] = row['min']
    	if minmax_l[1] is None or np.isnan(minmax_l[1]):
            minmax_l[1] = row['max']
    
    	#put plot coordinates in min_max
    	row['min_max'] = minmax_l
    	return row

    chrmd_df = chrmd_df.apply(fill_min_max, axis=1)
    
    
    
    # Create genes metadata DataFrame
    if color_col is None:
        color_col = id_col

    # start df with chromosome and the column defining color
    genesmd_df = subdf.groupby(id_col).agg({'Chromosome': 'first', 'Start': 'min', 'End': 'max', color_col: 'first'})
    genesmd_df.dropna(inplace=True) # remove genes not present in subset (NaN)
    genesmd_df.rename(columns={color_col: 'color_tag'}, inplace=True)
    genesmd_df['gene_ix_xchrom'] = genesmd_df.groupby('Chromosome').cumcount()
    
    #Assign y-coordinate
    if packed:
        genesmd_df['ycoord'] = -1
        genesmd_df = genesmd_df.groupby(genesmd_df['Chromosome']).apply(packed_for_genesmd) # add packed ycoord column using intervaltree
        genesmd_df = genesmd_df.reset_index(level='Chromosome', drop=True)
    else:
        genesmd_df['ycoord'] = genesmd_df.loc[:, 'gene_ix_xchrom']
    
    #store plot y height
    chrmd_df['y_height'] = genesmd_df.groupby('Chromosome').ycoord.max() 
    chrmd_df['y_height'] += 1# count from 1

    #Assign colors to genes
    color_tags = genesmd_df.color_tag.drop_duplicates()
    n_color_tags = len(color_tags)
    
    # 0-string to colormap object if possible
    if type(colormap) == str:
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
        colormap = list(colormap.colors) #colors of plt oject
        #adjust number of colors
        if n_color_tags < len(colormap):
            colormap = colormap[:n_color_tags] 
        #make plt colors compatible with plotly
        colormap = [f'rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})' for r, g, b in colormap] 
    
    # 2-list to dict
    if type(colormap) == list: 
        #adjust number of colors
        if n_color_tags < len(colormap):
            colormap = colormap[:n_color_tags]
        #create dict of colors
        colormap = {color_tags[i]: colormap[i%len(colormap)] for i in range(n_color_tags)} #iterate over colors if necessary
    
    # 3- Use dict to assign color to gene
    if type(colormap) == dict: 
        genesmd_df['color'] = genesmd_df['color_tag'].map(colormap)  ## NOTE: when specifying color by dict, careful with color_col
        genesmd_df['color'].fillna('black', inplace=True) # not specified in dict will be colored as black
    
    
    # Create figure and chromosome plots
    titles = ["Chromosome %s" % chrom for chrom in chrmd_df.index]
    fig = sp.make_subplots(rows=nchrs, cols=1, row_heights=chrmd_df.y_height.to_list(), subplot_titles=titles)
    #print('\n\n\n' + str(sum(chrmd_df.y_height)+ nchrs) + '\n\n\n')
    
    # one subplot per chromosome
    for i in range(nchrs):
        chrom = chrmd_df.index[i]
        fig.add_trace(go.Scatter(x=[], y=[]), row=i+1, col=1)
        
        # set title format
        fig.layout.annotations[i].update(font=title_dict_ply) 
	
        # set x axis limits
        x_min, x_max = chrmd_df.iloc[i]['min_max']
        x_rang = x_max - x_min
        fig.update_xaxes(range=[x_min - 0.05*x_rang, x_max + 0.05*x_rang], tickformat='d', showgrid=False, zeroline=False, row=i+1, col=1) # add 5% to limit coordinates range
        #fig.update_xaxes(title_text="Chromosome %s" % chrom, row=i + 1, col=1)
        
        # set y axis limits
        y_min = 0
        y_max = chrmd_df.iloc[i].y_height
        y_ticks_val = []
        y_ticks_name = []
        if not packed:
            y_ticks_val = [i + 0.5 for i in range(int(y_max))]
            y_ticks_name = genesmd_df.groupby(genesmd_df['Chromosome']).groups[chrom]
        fig.update_yaxes(range=[y_min, y_max], tickvals=y_ticks_val, ticktext=y_ticks_name, showgrid=False, row=i+1, col=1)


    # Plot genes
    subdf.groupby(id_col).apply(lambda subdf: _gby_plot_exons(subdf, fig, chrmd_df, genesmd_df, id_col, transcript_str, showinfo, tag_background))
    
    
    # Adjust plot display
    if max_ngenes > 25:
        fig.update_layout(title_text="<span style='color:red;'>Warning! The plot integity might be compromised when displaying too many genes.</span>") #warning for too many genes
    fig.update_layout(plot_bgcolor=plot_background, font_color=plot_border, showlegend=False)  #label text color == plot border color
    fig.update_xaxes(showline=True, linewidth=1, linecolor=plot_border, mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor=plot_border, mirror=True)
    
    
    # Provide output
    if to_file == None:
        fig.show()
    else:
        if not file_size:
            fig.update_layout(width=1600, height=800)
        else:
            fig.update_layout(width=file_size[0], height=file_size[1])
        
        pio.write_image(fig, to_file)


    

def _gby_plot_exons(df, fig, chrmd_df, genesmd_df, id_col, transcript_str, showinfo, tag_background):

    """Plot elements corresponding to the df rows of one gene."""

    # Gene parameters
    genename = df[id_col].iloc[0]
    genelabel = genename ### select label info
    gene_ix = genesmd_df.loc[genename]["ycoord"] + 0.5
    exon_color = genesmd_df.loc[genename].color
    chrom = genesmd_df.loc[genename].Chromosome
    chrom_ix = chrmd_df.index.get_loc(chrom)
    n_exons = len(df)
    if 'Strand' in df.columns:
        strand = df["Strand"].unique()[0]
    else:
        strand = ''
    
    # Get the gene information to print on hover
    if strand:
        geneinfo = f"[{strand}] ({min(df.Start)}, {max(df.End)})<br>ID: {genename}" #default with strand
    else:
        geneinfo = f"({min(df.Start)}, {max(df.End)})<br>ID: {genename}" #default without strand
    showinfo_data = []
    if showinfo:
        for i in range(len(showinfo)):
           col = showinfo[i]
           showinfo_data.append(df[col].iloc[0]) # first by default, change for exons.
           geneinfo += f"<br>{showinfo[i]}: {showinfo_data[i]}"


    # Plot transcript structure
    if transcript_str:
        # transcript has CDS and exon
        if df.Feature.str.contains('CDS').any() and df.Feature.str.contains('exon').any(): 
            #get coordinates for utr and cds
            tr_start, cds_start = df.groupby('Feature').Start.apply(min)[['exon', 'CDS']]
            tr_end, cds_end = df.groupby('Feature').End.apply(max)[['exon', 'CDS']]
        
            #create start utr
            x0, x1 = tr_start, cds_start
            y0, y1 = gene_ix - transcript_utr_width/2, gene_ix + transcript_utr_width/2
            fig.add_trace(
                go.Scatter(x=[x0, x1, x1, x0, x0], y=[y0, y0, y1, y1, y0],
                           fill = "toself", fillcolor = exon_color, mode = 'lines', 
                           line = dict(color=exon_color), text = geneinfo),
                           row=chrom_ix+1, 
                           col=1)
                           
            #create end utr
            x0, x1 = cds_end, tr_end
            y0, y1 = gene_ix - transcript_utr_width/2, gene_ix + transcript_utr_width/2
            fig.add_trace(
                go.Scatter(x=[x0, x1, x1, x0, x0], y=[y0, y0, y1, y1, y0],
                           fill = "toself", fillcolor = exon_color, mode = 'lines', 
                           line = dict(color=exon_color), text = geneinfo),
                           row=chrom_ix+1, 
                           col=1)
                           
            #remove non-CDS from data
            df = df.groupby('Feature').get_group('CDS')
            
        # transcript only has CDS
        elif df.Feature.str.contains('CDS').any() and not df.Feature.str.contains('exon').any():
            print()
        
        # trancript only has exon    
        elif not df.Feature.str.contains('CDS').any() and df.Feature.str.contains('exon').any():
            #plot just as utr
            df.apply(_apply_gene, args=(fig, strand, genename, gene_ix, exon_color, chrom, chrom_ix, n_exons, genelabel, geneinfo, transcript_utr_width), axis=1)
            
        # transcript has neither, skip it
        else:
            return


    # Plot LINE binding the exons
    x0, x1 = min(df.Start), max(df.End)
    y0, y1 = gene_ix - exon_width/100, gene_ix + exon_width/100
    exon_line = go.Scatter(
        x=[x0, x1, x1, x0, x0], 
        y=[y0, y0, y1, y1, y0],
        fill="toself",
        fillcolor=exon_color,
        mode='lines',
        line=dict(color=exon_color, width=0.5),
        text=geneinfo)
    fig.add_trace(exon_line, row=chrom_ix+1, col=1)
    
    # Plot the gene rows
    # not transcript
    if not transcript_str:
        df.apply(_apply_gene, args=(fig, strand, genename, gene_ix, exon_color, chrom, chrom_ix, n_exons, genelabel, geneinfo, exon_width), axis=1)
    # trancript does not only have exon
    if transcript_str and df.Feature.str.contains('CDS').any() and not df.Feature.str.contains('exon').any():
        df.apply(_apply_gene, args=(fig, strand, genename, gene_ix, exon_color, chrom, chrom_ix, n_exons, genelabel, geneinfo, exon_width), axis=1)
    
    # Plot DIRECTION ARROW in INTRONS if strand is known
    sorted_exons = df[['Start', 'End']].sort_values(by = 'Start')
    
    if strand:
        # evaluate  each intron
        for i in range(len(sorted_exons)-1):
            start =  sorted_exons['End'].iloc[i] 
            stop = sorted_exons['Start'].iloc[i+1]
            intron_size = coord2percent(fig, chrom_ix+1, start, stop)
            incl = percent2coord(fig, chrom_ix+1, 0.003) #how long in the plot (OX)
            
            # create and plot lines
            if intron_size > intron_threshold:
                ##diagonal_line = OX arrow extension(intron middle point +- incl), OY arrow extension (intron middle point + half of exon width)
                top_plus = [(start+stop)/2 + incl, (start+stop)/2 - incl], [gene_ix, gene_ix + exon_width/2-0.01]
                bot_plus = [(start+stop)/2 - incl, (start+stop)/2 + incl], [gene_ix - exon_width/2+0.01, gene_ix]
                top_minus = [(start+stop)/2 + incl, (start+stop)/2 - incl], [gene_ix - exon_width/2+0.01, gene_ix]
                bot_minus = [(start+stop)/2 - incl, (start+stop)/2 + incl], [gene_ix, gene_ix + exon_width/2-0.01]
                
                if strand == '+':
                    arrow_bot = go.Scatter(x=bot_plus[0], y=bot_plus[1], mode='lines',
                                           line=go.scatter.Line(color=arrow_color, width=1), showlegend=False,
                                           name=genename, hoverinfo='skip')
                    arrow_top = go.Scatter(x=top_plus[0], y=top_plus[1], mode='lines',
                                           line=go.scatter.Line(color=arrow_color, width=1), showlegend=False,
                                           name=genename, hoverinfo='skip')
                    fig.add_trace(arrow_bot, row=chrom_ix+1, col=1)
                    fig.add_trace(arrow_top, row=chrom_ix+1, col=1)

                elif strand == '-':
                    arrow_bot = go.Scatter(x=bot_minus[0], y=bot_minus[1], mode='lines',
                                           line=go.scatter.Line(color=arrow_color, width=1), showlegend=False,
                                           name=genename, hoverinfo='skip')
                    arrow_top = go.Scatter(x=top_minus[0], y=top_minus[1], mode='lines',
                                           line=go.scatter.Line(color=arrow_color, width=1), showlegend=False,
                                           name=genename, hoverinfo='skip')
                    fig.add_trace(arrow_bot, row=chrom_ix+1, col=1)
                    fig.add_trace(arrow_top, row=chrom_ix+1, col=1)
    




def _apply_gene(row, fig, strand, genename, gene_ix, exon_color, chrom, chrom_ix, 
                n_exons, genelabel, geneinfo, exon_width):

    """Plot elements corresponding to one row of one gene."""
    
    # Exon start and stop
    start = int(row["Start"])
    stop = int(row["End"])
    # convert to coordinates for rectangle
    x0, x1 = start, stop
    y0, y1=gene_ix-exon_width/2, gene_ix+exon_width/2 ##gene middle point -+ half of exon size
    
    # Plot EXON as rectangle
    fig.add_trace(
        go.Scatter(x=[x0, x1, x1, x0, x0], y=[y0, y0, y1, y1, y0],
            fill = "toself", fillcolor = exon_color, mode = 'lines', 
            line = dict(color=exon_color), text = geneinfo, hoverinfo = "text"),
            row=chrom_ix+1, 
            col = 1
            )
    
    # Plot DIRECTION ARROW in EXON
    # decide about placing a direction arrow
    arrow_size = coord2percent(fig, chrom_ix+1, 0.05*start, 0.05*stop)
    #too small to plot
    if arrow_size <= arrow_size_min: 
        incl = 0
    #plot the arrow
    else:
        incl = percent2coord(fig, chrom_ix+1, 0.003) #how long in the plot (OX)
    
    # create and plot lines
    if incl:
        ##diagonal_line = OX arrow extension(exon middle point +- incl), OY arrow extension (exon middle point + half of exon width)
        top_plus = [(start+stop)/2 + incl, (start+stop)/2 - incl], [gene_ix, gene_ix + exon_width/2-0.01]
        bot_plus = [(start+stop)/2 - incl, (start+stop)/2 + incl], [gene_ix - exon_width/2+0.01, gene_ix]
        top_minus = [(start+stop)/2 + incl, (start+stop)/2 - incl], [gene_ix - exon_width/2+0.01, gene_ix]
        bot_minus = [(start+stop)/2 - incl, (start+stop)/2 + incl], [gene_ix, gene_ix + exon_width/2-0.01]

        if strand == '+':
            arrow_bot = go.Scatter(x=bot_plus[0], y=bot_plus[1], mode='lines',
                                   line=go.scatter.Line(color=arrow_color, width=arrow_width),
                                   showlegend=False, name=genename, hoverinfo='skip')
            arrow_top = go.Scatter(x=top_plus[0], y=top_plus[1], mode='lines',
                                   line=go.scatter.Line(color=arrow_color, width=arrow_width),
                                   showlegend=False, name=genename, hoverinfo='skip')

        elif strand == '-':
            arrow_bot = go.Scatter(x=bot_minus[0], y=bot_minus[1], mode='lines',
                                   line=go.scatter.Line(color=arrow_color, width=arrow_width),
                                   showlegend=False, name=genename, hoverinfo='skip')
            arrow_top = go.Scatter(x=top_minus[0], y=top_minus[1], mode='lines',
                                   line=go.scatter.Line(color=arrow_color, width=arrow_width),
                                   showlegend=False, name=genename, hoverinfo='skip')
            
        fig.add_trace(arrow_bot, row=chrom_ix+1, col=1)
        fig.add_trace(arrow_top, row=chrom_ix+1, col=1)















    
    

