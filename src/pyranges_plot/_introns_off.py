import pyranges as pr
import pandas as pd


# def get_introns(self, id_col: str) -> "pr.PyRanges":
def get_introns(p, id_col) -> "pr.PyRanges":
    """Calculate introns from a PyRanges object.

     Parameters
    ----------
    id_col : str
        Name of column containing the ID infromation.

    Examples
    --------
    >>> gr = pr.PyRanges({"Chromosome": [1, 1, 1],
    ...                   "Start": [0, 20, 40],
    ...                   "End": [10, 35, 50],
    ...                   "transcript_id": ['t1', 't1' ,'t1']})
    >>> gr
      index  |      Chromosome    Start      End  transcript_id
      int64  |           int64    int64    int64  object
    -------  ---  ------------  -------  -------  ---------------
          0  |               1        0       10  t1
          1  |               1       20       35  t1
          2  |               1       40       50  t1
    PyRanges with 3 rows, 4 columns, and 1 index columns.
    Contains 1 chromosomes.

    >>> gr.get_introns("transcript_id")
      index  |      Chromosome      End    Start  transcript_id
      int64  |           int64    int64    int64  object
    -------  ---  ------------  -------  -------  ---------------
          1  |               1       20       10  t1
          2  |               1       40       35  t1
    PyRanges with 2 rows, 4 columns, and 1 index columns.
    Contains 1 chromosomes.
    """

    introns = (
        p.copy()
        .groupby(
            [id_col, "pr_ix"], group_keys=False, observed=True
        )  ## tambien por feature si esta
        .apply(pr.PyRanges.sort_by_position)
    )

    # intron start is exon end shifted
    if "pr_ix" in introns.columns:
        introns["End"] = introns.groupby(
            [id_col, "pr_ix"], group_keys=False, observed=True
        )["End"].shift()
    else:
        introns["End"] = introns.groupby(id_col, group_keys=False, observed=True)[
            "End"
        ].shift()

    # remove first exon rows
    introns.dropna(subset="End", inplace=True)
    introns.reset_index(inplace=True)

    # intron end is exon start
    introns.rename(columns={"Start": "End", "End": "Start"}, inplace=True)
    introns["Start"] = [int(i) for i in introns["Start"]]

    return introns


def introns_resize(df, ts_data, id_col):
    """Calculate intron resizes and provide info for plotting"""

    chrom = df["Chromosome"].iloc[0]
    p = df
    thresh = df["shrink_threshold"].iloc[0]

    # Calculate shrinkable intron ranges
    # get flexible introns
    exons = p.copy()
    introns = get_introns(p, id_col)
    to_shrink = pr.PyRanges()

    if not introns.empty:
        flex_introns = introns.subtract_ranges(exons, strand_behavior="ignore")

        # obtain shrinkable regions
        to_shrink = flex_introns.merge_overlaps(use_strand=False)  # unique ranges
        if not to_shrink.empty:
            to_shrink = to_shrink[
                to_shrink["End"] - to_shrink["Start"] > thresh
            ]  # filtered

    # nohing to shrink
    if to_shrink.empty:
        ts_data[chrom] = pd.DataFrame(
            columns=["Chromosome", "Start", "End", "Start_adj", "End_adj", "cumdelta"]
        )
        result = p
        result["Start_adj"] = result["Start"]
        result["End_adj"] = result["End"]
        result["delta"] = [0] * len(result)
        result["cumdelta"] = [0] * len(result)

        return result

    # get coordinate shift (delta) and cumulative coordinate shift (cumdelta)
    to_shrink["delta"] = (
        to_shrink["End"] - to_shrink["Start"]
    ) - thresh  # calculate coord shift considering margins
    assert to_shrink.sort_values("Start").equals(to_shrink), "PyRanges not sorted."
    to_shrink["cumdelta"] = to_shrink["delta"].cumsum()

    # store adjusted coord to plot shrinked intron regions
    to_shrink["Start_adj"] = to_shrink["Start"] - to_shrink.cumdelta.shift().fillna(0)
    to_shrink["End_adj"] = to_shrink["End"] - to_shrink.cumdelta

    # store to shrink data
    ts_data[chrom] = to_shrink

    # Calculate exons coordinate shift
    exons = pr.concat([exons, to_shrink])
    exons.sort_values("Start", inplace=True)
    exons["cumdelta"] = exons["cumdelta"].ffill()
    exons = exons.fillna({"cumdelta": 0})
    # match exons with its cumdelta
    result = exons.dropna(subset=[id_col])

    # Adjust coordinates
    result["Start_adj"] = result["Start"] - result["cumdelta"]
    result["End_adj"] = result["End"] - result["cumdelta"]

    # Provide result
    return result[list(p.columns) + ["Start_adj", "End_adj", "delta"]]


def recalc_axis(ts_data, tick_pos_d, ori_tick_pos_d):
    """Calculate shrinked axis values according to original coordinates."""

    for chrom in ts_data.keys():
        # add to-shrinked reagions limits to axis
        ori_tick_pos = []
        tick_pos = []

        if ts_data[chrom].empty:  # nothing to shrink
            tick_pos_d[chrom] = []
            ori_tick_pos_d[chrom] = []

        else:
            pos = [
                [a, b] for a, b in zip(ts_data[chrom]["Start"], ts_data[chrom]["End"])
            ]
            cdel = list(ts_data[chrom]["cumdelta"])

            # update tick positions for shrinked regions and keep original values as names
            for i in range(len(pos)):
                if i == 0:
                    tick_pos.append(pos[i][0])
                    tick_pos.append(pos[i][1] - cdel[i])
                else:
                    tick_pos.append(pos[i][0] - cdel[i - 1])
                    tick_pos.append(pos[i][1] - cdel[i])

                ori_tick_pos.append(pos[i][0])
                ori_tick_pos.append(pos[i][1])

            tick_pos_d[chrom] = tick_pos
            ori_tick_pos_d[chrom] = ori_tick_pos

    return tick_pos_d, ori_tick_pos_d
