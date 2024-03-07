def _gby_plot_vcf(
    vcf,
    vcf_axes,
    chrmd_df,
):
    """Add vcf data to plots."""

    # Select axis
    chrom = vcf["Chromosome"].iloc[0]
    chrom_ix = chrmd_df.index.get_loc(chrom)
    ax = vcf_axes[chrom_ix]

    # adjust y limits?
    ax.set_ylim(0.5, 1.5)

    # Plot each variant
    for i, row in vcf.iterrows():
        x = int(row["Position"])
        y = 1  # provisional
        ax.scatter(x, y)
