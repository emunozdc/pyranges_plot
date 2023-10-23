tag_background = "grey"
plot_background = "white"
plot_border = "black"
#{'family': title_font, 'color': title_color, 'size': title_size}
title_dict_ply = {'family': 'Arial', 'color': 'goldenrod', 'size': 18} 
title_dict_plt = {'family': 'sans-serif', 'color': 'goldenrod', 'size': 13}


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
        globals()[dictname][keyname] = value

    else:
        globals()[varname] = value



def get_default(varname):
    """
    Obtain the deafault value for a plot layout variable.
    
    Parameters
    ----------
    varname: str
        
        Name of the variable to get the value.
    
    Examples
    --------
    >>> pyranges_plot.get_default('plot_border')
    
    >>> pyranges_plot.get_default('title_dict_plt')
    
    >>> pyranges_plot.get_default('title_dict_ply')[color]
    
    """

    return globals()[varname]
