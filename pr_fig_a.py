import pyranges as pr, pyranges_plot as prp

g = pr.PyRanges(
    {
        "Chromosome": [1, 1, 1],
        "Start": [4, 13, 19],
        "End": [10, 15, 23],
        "transcript_id": ["g", "g", "g"],
    }
)

g_ext = g.extend(1)
g_ext["transcript_id"] = ["g.extend(1)"] * len(g_ext)
g_ext_id = g.extend(1, transcript_id="transcript_id")
g_ext_id["transcript_id"] = ["g.extend(1, id)"] * len(g_ext_id)

b = g.boundaries("transcript_id")
b["transcript_id"] = ["g.boundries"] * len(b)
b_subseq = b.subsequence(0, 10)
b_subseq["transcript_id"] = ["g.boundries.subseq"] * len(b_subseq)

g_subseq_id = g.subsequence(0, 10, "transcript_id")
g_subseq_id["transcript_id"] = ["g.subsequence(0,10,id)"] * len(g_subseq_id)
g_subseq = g.subsequence(0, 3)
g_subseq["transcript_id"] = ["g.subsequence(0,3)"] * len(g_subseq)

g_spl_subseq_for = g.spliced_subsequence(1, 8, "transcript_id")
g_spl_subseq_for["transcript_id"] = ["g.spliced_subsequence(1,8,id)"] * len(
    g_spl_subseq_for
)
g_spl_subseq_rev = g.spliced_subsequence(-5, -1, "transcript_id")
g_spl_subseq_rev["transcript_id"] = ["g.spliced_subsequence(-5,-1,id)"] * len(
    g_spl_subseq_rev
)


prp.plot(
    [
        g,
        pr.concat([g_ext, g_ext_id]),
        pr.concat([b, b_subseq]),
        pr.concat([g_subseq, g_subseq_id]),
        pr.concat([g_spl_subseq_for, g_spl_subseq_rev]),
    ],
    engine="ply",
    id_col="transcript_id",
    # packed=False,
)
