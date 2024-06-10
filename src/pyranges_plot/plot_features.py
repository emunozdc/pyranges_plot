import copy


plot_features_dict = {
    "colormap": (
        "Alphabet",
        "Sequence of colors for the genes, it can be provided as a Matplotlib colormap,a Plotly color sequence (built as lists), a string naming the previously mentioned color objects from Matplotlib and Plotly, or a dictionary with the following structure {color_column_value1: color1, color_column_value2: color2, ...}. When a specific color_col value is not specified in the dictionary it will be colored in black.",
        " ",
    ),
    "tag_bkg": (
        "grey",
        "Background color of the tooltip annotation for the gene in Matplotlib.",
        " ",
    ),
    "fig_bkg": ("white", "Bakground color of the whole figure.", " "),
    "plot_bkg": ("white", "Background color for the chromosomes plots.", " "),
    "plot_border": ("black", "Color of the line defining the chromosome plots.", " "),
    "title_size": (18, "Size of the plots' titles.", " "),
    "title_color": ("black", "Color of the plots' titles.", " "),
    "grid_color": ("lightgrey", "Color of x coordinates grid lines.", " "),
    "exon_border": (None, "Color of the interval's rectangle border.", " "),
    "exon_width": (0.4, "Height of the exon rectangle in the plot.", " "),
    "text_size": (10, "Fontsize of the text annotation beside the intervals.", " "),
    "text_pad": (
        1,
        "Space where the id annotation is placed beside the interval. When text_pad is float, it represents the percentage of the plot space, while an int pad represents number of positions or base pairs.",
        " ",
    ),
    "v_space": (0.7, "Vertical distance between exons in different y heights.", " "),
    "shrink_threshold": (
        0.01,
        "Minimum length of an intron in order for it to be shrinked while using the introns_off feature. When threshold is float, it represents the percentage of the plot space, while an int threshold represents number of positions or base pairs.",
        " ",
    ),
    "plotly_port": (8050, "Port to run plotly app.", " "),
    "arrow_line_width": (
        1,
        "Line width of the arrow lines (for stranded PyRanges).",
        " ",
    ),
    "arrow_color": ("grey", "Direction arrow color (for stranded PyRanges).", " "),
    "arrow_size_min": (
        0.002,
        "Minimum size of the arrow to plot direction in exons if necessary. Provided as a float corresponding to the plot fraction or percentage.",
        " ",
    ),
    "arrow_size": (
        0.006,
        "Fraction or percentage of the plot occupied by a direction arrow.",
        " ",
    ),
    "arrow_intron_threshold": (
        0.04,
        "Minimum size of the intron to plot a direction arrow in it. Provided as a float corresponding to the plot fraction or percentage.",
        " ",
    ),
    "shrinked_bkg": (
        "lightyellow",
        "Color of the shrinked region background.",
        " ",
    ),
    "shrinked_alpha": (
        0.7,
        "Opacity of the shrinked region background color.",
        " ",
    ),
}

# Normal (light theme)
plot_features_dict_in_use = copy.deepcopy(plot_features_dict)
plot_features_dict_vals = {}
for key, val in plot_features_dict.items():
    plot_features_dict_vals[key] = val[0]

# Dark theme
theme_dark = {
    "colormap": "G10",
    "fig_bkg": "#1f1f1f",
    "plot_border": "white",
    "title_color": "goldenrod",
    "plot_bkg": "grey",
    "grid_color": "darkgrey",
    "arrow_color": "lightgrey",
    "shrinked_bkg": "lightblue",
    "shrinked_alpha": 0.4,
}

# Store themes
builtin_themes = {
    "light": plot_features_dict_vals,
    "dark": theme_dark,
}
