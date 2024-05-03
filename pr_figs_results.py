import pyranges as pr, pyranges_plot as prp


# Figure 1
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

gr_1 = gr[gr["Chromosome"] == "1"]

# "default" plot (Figure 1.1)
prp.set_engine("plt")
prp.plot(gr)

# id, color col and cmap (Figure 1.2)
prp.plot(
    gr_1,
    id_col="transcript_id",
    color_col="Strand",
    colormap={"+": "lightgreen", "-": "lightblue"},
    # to_file="fig1_2.png",
)


# Figure 2
gr = pr.example_data.ncbi_gff
grp = gr[gr.Feature.isin(["CDS", "exon"])]
grp = grp[
    grp.Parent.isin(["rna-DGYR_LOCUS12552-2", "rna-DGYR_LOCUS12552"])
]  # , "rna-DGYR_LOCUS13738", "rna-DGYR_LOCUS13739", "rna-DGYR_LOCUS13730"])]

grpp = grp[["Chromosome", "Feature", "Start", "End", "Strand", "Parent"]]

# show transcript_str (Figure 2.1)
prp.plot(
    grpp,
    id_col="Parent",
    transcript_str=True,
    id_ann=False,
)

# show introns off (Figure 2.2)
prp.plot(
    gr_1,
    id_col="transcript_id",
    introns_off=True,
)
