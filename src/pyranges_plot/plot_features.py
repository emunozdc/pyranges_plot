import copy


plot_features_dict = {
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

# Normal (light mode)
plot_features_dict_in_use = copy.deepcopy(plot_features_dict)

# Dark mode
dark_plot_features_dict_in_use = copy.deepcopy(plot_features_dict)
dark_plot_features_dict_in_use["fig_bkg"] = (
    "#1f1f1f",
    "Background color of the whole figure.",
    " ",
)
dark_plot_features_dict_in_use["plot_border"] = (
    "white",
    "Color of the line defining the chromosome plots.",
    " ",
)
dark_plot_features_dict_in_use["title_color"] = (
    "goldenrod",
    "Color of the plots' titles.",
    " ",
)
dark_plot_features_dict_in_use["plot_bkg"] = (
    "grey",
    "Background color for the chromosomes plots.",
    " ",
)
dark_plot_features_dict_in_use["grid_color"] = (
    "darkgrey",
    "Color of x coordinates grid lines.",
    " ",
)
dark_plot_features_dict_in_use["arrow_color"] = (
    "lightgrey",
    "Direction arrow color (for stranded PyRanges).",
    " ",
)
