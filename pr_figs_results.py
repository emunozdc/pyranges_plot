import pyranges as pr, pyranges_plot as prp


### prp simple
gr = pr.PyRanges(
    {
        "Chromosome": ["1"] * 10 + ["2"] * 10,
        "Strand": ["+", "+", "+", "+", "-", "-", "-", "-", "+", "+"]
        + ["+", "+", "+", "+", "-", "-", "-", "-", "+", "+"],
        "Start": [90, 61, 104, 228, 9, 142, 52, 149, 218, 151]
        + [6, 27, 37, 47, 1, 7, 42, 37, 60, 80],
        "End": [92, 64, 113, 229, 12, 147, 57, 155, 224, 153]
        + [8, 32, 40, 50, 5, 10, 46, 40, 70, 90],
        "transcript_id": ["t1", "t1", "t1", "t1", "t2", "t2", "t2", "t2", "t3", "t3"]
        + ["t4", "t4", "t4", "t4", "t5", "t5", "t5", "t5", "t6", "t6"],
    }
)

# "default" plot
prp.set_engine("ply")
prp.plot(gr)

# id, color col and cmap
prp.plot(
    gr,
    id_col="transcript_id",
    color_col="Strand",
    colormap={"+": "lightgreen", "-": "lightblue"},
)

# id, packed, title and default
prp.plot(
    gr,
    id_col="transcript_id",
    packed=False,
    title_chr="Chr{chrom}",
    title_size=30,
    exon_width=0.7,
)

### prp complex

gr = pr.example_data.ensembl_gtf

prp.plot(
    gr,
    id_col="gene_id",
)
