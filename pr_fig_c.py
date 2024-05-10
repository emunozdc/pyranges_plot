import pyranges as pr, pyranges_plot as prp

a = pr.PyRanges(
    {
        "Chromosome": [1] * 7,
        "Start": [2, 12, 17, 22, 27, 32, 33],
        "End": [5, 14, 20, 26, 29, 37, 36],
        "Strand": ["+"] * 2 + ["-"] * 3 + ["+"] * 2,
    }
)
.
b = pr.PyRanges(
    {
        "Chromosome": [1] * 5,
        "Start": [6, 11, 18, 24, 34],
        "End": [8, 13, 19, 28, 36],
        "Strand": ["+"] * 3 + ["-"] * 1 + ["+"] * 1,
    }
)

a_ov_b = a.overlap(b)
a_ov_b_slack = a.overlap(b, slack=2)
a_ov_b_nostrand = a.overlap(b, strand_behavior="ignore")
a_ov_b_opstrand = a.overlap(b, strand_behavior="opposite")

a_inters_b = a.intersect(b)
a_setinters_b = a.set_intersect(b)

a_subt_b = a.subtract_ranges(b)  ## is this a.subtract(the slack one)


prp.plot(
    [
        a,
        b,
        a_ov_b,
        a_ov_b_slack,
        a_ov_b_nostrand,
        a_ov_b_opstrand,
        a_inters_b,
        a_setinters_b,
        a_subt_b,
    ],
    engine="plt",
    title_chr=" ",
    id_ann=False,
    # limits=(-4, None),
    # to_file="fig3_3.png",
    file_size=(10, 8),
    colormap=["black"],
)
