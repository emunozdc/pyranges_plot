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
    "title_color": ("goldenrod", "Color of the plots' titles.", " "),
    "exon_width": (0.4, "Height of the exon rectangle in the plot.", " "),
}


plot_features_dict_in_use = copy.deepcopy(plot_features_dict)
