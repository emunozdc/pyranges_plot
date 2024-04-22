import pyranges as pr, pyranges_plot as prp

p = pr.PyRanges(
    {
        "Chromosome": [1, 1, 1],
        "Start": [4, 13, 19],
        "End": [10, 15, 23],
        "transcript_id": ["a", "a", "a"],
    }
)

prp.plot(
    [p, p.extend(1), p.extend(1, transcript_id="transcript_id")],
    engine="plt",
    id_col="transcript_id",
)
