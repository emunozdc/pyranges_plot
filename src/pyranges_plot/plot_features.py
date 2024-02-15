import copy

plot_features_dict = {
    "tag_background": (
        "grey",
        "Background color of the tooltip annotation for the gene in Matplotlib.",
        " ",
    ),
    "plot_background": ("white", "Background color for the chromosomes plots.", " "),
    "plot_border": ("black", "Color of the line defining the chromosome plots.", " "),
    "title_size": (18, "Size of the plots' titles.", " "),
    "title_color": ("black", "Color of the plots' titles.", " "),
    "exon_width": (0.4, "Height of the exon rectangle in the plot.", " "),
    "shrink_threshold": (
        0.05,
        "Minimum lenght of an intron in order for it to be shrinked while using the introns_off feature. When threshold is float, it represents the percentage of the plot space, while an int threshold represents number of positions or base pairs.",
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
        "Minimum size of the arrow to plot direction in exons if necessary. Provided as a float correspondig to the plot fraction or percentage.",
        " ",
    ),
    "arrow_size": (
        "0.006",
        "Fraction or percentage of the plot occupied by a direction arrow.",
        " ",
    ),
    "arrow_intron_threshold": (
        0.04,
        "Minimum size of the intron to plot a direction arrow in it. Provided as a float correspondig to the plot fraction or percentage.",
        " ",
    ),
}


plot_features_dict_in_use = copy.deepcopy(plot_features_dict)
