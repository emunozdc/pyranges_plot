import pyranges as pr


def introns_shrink(df, thresh=0):
    p = pr.from_dict(df.to_dict())

    # Calculate shrinkable intron ranges
    # get flexible introns
    p = p.sort(by=['transcript_id', 'Start'])
    introns = p.features.introns(by='gene')
    exons = p[p.Feature == 'exon']
    flex_introns = introns.subtract(exons)

    # get coordinate shift (delta) and cumulative coordinate shift (cumdelta)
    to_shrink = flex_introns.merge(strand=False, count=True, count_col='count')  # unique
    to_shrink = to_shrink[to_shrink.End - to_shrink.Start > thresh]  # filtered
    to_shrink.delta = (to_shrink.df["End"] - to_shrink.df[
        "Start"]) - thresh  # calculate coord shift considering margins
    assert to_shrink.df.sort_values("Start").equals(to_shrink.df), "PyRanges not sorted."
    to_shrink.cumdelta = to_shrink.df["delta"].cumsum()

    # Match ranges where not shrinked intron lines must be placed
    # associate proper cumdelta to fixed intron lines
    inters = introns.intersect(exons, strandedness=False)
    inters = inters.df.merge(to_shrink.df, how='left', left_on=['Start'], right_on=['End'], suffixes=('', '_y')).fillna(
        {'cumdelta': 0, 'delta': 0})
    inters = inters.sort_values("Start")
    inters["cumdelta"] = inters["cumdelta"].replace(0, method='ffill')
    inters.drop(["Chromosome_y", "Start_y", "End_y", "count"], axis=1, inplace=True)
    inters = pr.from_dict(inters.to_dict())

    # relate fixed intron lines to each intron
    introns = introns.join(inters, how='left', suffix='_line').df
    introns.rename(columns={"cumdelta": "cumdelta_line"}, inplace=True)
    introns["delta"] = introns['delta'].replace(-1, 0)

    # adjust fixed intron coordinates
    introns["Start_line"] -= introns["cumdelta_line"]
    introns["End_line"] -= introns["cumdelta_line"]
    introns["i_lines"] = introns.apply(
        lambda row: [row["Start_line"], row["End_line"]] if row["Start_line"] != 0 and row["End_line"] != 0 else None,
        axis=1,
    )

    # store all fixed intron lines from one intron in a list to plot it later
    introns = introns.groupby(["Start", "End"], group_keys=False).agg({
        'Chromosome': 'first',
        'Strand': 'first',
        'transcript_id': 'first',
        'gene_id': 'first',
        'Feature': 'first',
        'i_lines': lambda lines: [interval for interval in lines if interval is not None],
        'delta': 'first',
    }).reset_index()

    introns = pr.from_dict(introns.to_dict())

    # Calculate coordinate shift for each intron
    # keep highest cumdelta for each intron
    intronsdf = introns.join(to_shrink, how='left').df
    intronsdf["cumdelta"] = intronsdf['cumdelta'].replace(-1, 0)

    intronsdf = intronsdf.groupby(["Start", "End"], group_keys=False).agg({
        'Chromosome': 'first',
        'Strand': 'first',
        'transcript_id': 'first',
        'gene_id': 'first',
        'Feature': 'first',
        'cumdelta': 'max',
        'delta': 'first',
        'i_lines': 'first',
    }).reset_index()

    # fill cumdelta empty values correctly
    intronsdf = intronsdf.sort_values(by="Start")
    intronsdf["cumdelta"] = intronsdf['cumdelta'].replace(-1, 0)
    intronsdf["cumdelta"] = intronsdf['cumdelta'].replace(0, method='ffill')
    intronsdf["num_exon"] = intronsdf.groupby("transcript_id", group_keys=False).cumcount() + 1

    # Match introns-exons and calculate coordinates shift for exons
    exons = exons.df
    exons["num_exon"] = exons.groupby("transcript_id", group_keys=False).cumcount()  # not needed
    result = exons.merge(intronsdf, how='left', left_on=['transcript_id', 'Start'], right_on=["transcript_id", "End"],
                         suffixes=('', '_y'))

    # fill cumdelta
    result = result.sort_values("Start")
    result = result.fillna({'cumdelta': 0, 'delta': 0, 'i_lines': 0})
    result["cumdelta"] = result['cumdelta'].replace(0, method='ffill').apply(int)
    # result["delta"] = result['delta'].replace(0, method='ffill').apply(int)
    # adjust coordinates
    # result["Start"] -= result["cumdelta"]
    # result["End"] -= result["cumdelta"]
    result["Start_adj"] = result["Start"] - result["cumdelta"]
    result["End_adj"] = result["End"] - result["cumdelta"]

    # return result
    return result[list(p.columns) + ["Start_adj", "End_adj", "cumdelta", "delta", "i_lines"]]



def introns_shrinkkk(df, thresh=4):

    p = pr.from_dict(df.to_dict())

    # Calculate shrinkable intron ranges
    # get flexible introns
    p = p.sort(by=['transcript_id', 'Start'])
    introns = p.features.introns(by='gene')
    exons = p[p.Feature == 'exon']
    flex_introns = introns.subtract(exons)

    # get coordinate shift (delta) and cumulative coordinate shift (cumdelta)
    to_shrink = flex_introns.merge(strand=False, count=True, count_col='count')  # unique
    to_shrink = to_shrink[to_shrink.End - to_shrink.Start > thresh]  # filtered
    to_shrink.delta = (to_shrink.df["End"] - to_shrink.df[
        "Start"]) - thresh  # calculate coord shift considering margins
    assert to_shrink.df.sort_values("Start").equals(to_shrink.df), "PyRanges not sorted."
    to_shrink.cumdelta = to_shrink.df["delta"].cumsum()

    # Calculate coordinate shift for each intron
    # keep highest cumdelta for each intron
    intronsdf = introns.join(to_shrink, how='left')
    intronsdf = intronsdf.df.sort_values("cumdelta", ascending=False).drop_duplicates(["Start", "End"])

    # fill cumdelta empty values correctly
    intronsdf = intronsdf.sort_values(by="Start")
    intronsdf["cumdelta"] = intronsdf['cumdelta'].replace(-1, 0)
    intronsdf["cumdelta"] = intronsdf['cumdelta'].replace(0, method='ffill')
    intronsdf["num_exon"] = intronsdf.groupby("transcript_id", group_keys=False).cumcount() + 1

    # Match introns-exons and calculate coordinates shift for exons
    exons = exons.df
    exons["num_exon"] = exons.groupby("transcript_id", group_keys=False).cumcount()  # not needed
    result = exons.merge(intronsdf, how='left', left_on=['transcript_id', 'Start'], right_on=["transcript_id", "End"],
                         suffixes=('', '_y'))
    # fill cumdelta
    result = result.sort_values("Start")
    result = result.fillna({'cumdelta': 0})
    result["cumdelta"] = result['cumdelta'].replace(0, method='ffill').apply(int)
    # adjust coordinates
    # result["Start"] -= result["cumdelta"]
    # result["End"] -= result["cumdelta"]
    result["Start_adj"] = result["Start"] - result["cumdelta"]
    result["End_adj"] = result["End"] - result["cumdelta"]

    # return result
    return result[list(p.columns) + ["Start_adj", "End_adj", "cumdelta"]]
