import pyranges as pr, pyranges_plot as prp

# Start r object
r = pr.PyRanges(
    {
        "Chromosome": [1, 1, 1],
        "Start": [3, 11, 16],
        "End": [8, 15, 19],
        "transcript_id": ["r", "r", "r"],
    }
)

# r.window object
r_window = r.window(2)
r_window["transcript_id"] = ["r.window(2)"] * len(r_window)  # rename for plot

# tile_genome object
tile_g = pr.tile_genome(
    pr.PyRanges({"Chromosome": [1], "Start": [0], "End": [max(r["End"])]}), 2
)
tile_g["transcript_id"] = ["pr.tile_genome(2)"] * len(tile_g)  # rename for plot

# r.tile object
r_tile = r.tile(2)
r_tile["transcript_id"] = ["r.tile(2)"] * len(r_tile)  # rename for plot

# plot
prp.plot(
    [r, r_window, tile_g, r_tile],
    engine="ply",
    id_col="transcript_id",
    exon_border="black",
    # packed=False,
)
