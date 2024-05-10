import plotly.subplots as sp
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from pyranges_plot.core import cumdelting


def calculate_ticks(subdf, num_ticks=10):
    """Calculate tick values for a given data range."""

    # Calculate range and initial tick interval
    data_min = subdf["oriStart"].min()
    data_max = subdf["oriEnd"].max()
    data_range = data_max - data_min
    initial_interval = data_range / (num_ticks - 1)

    # Calculate a 'nice' interval between ticks
    # The exponent of the range
    exponent = np.floor(np.log10(initial_interval))
    # Fractional part of the range
    fractional_part = initial_interval / 10**exponent

    # Determine nice fractional part
    if fractional_part < 1.5:
        nice_fractional_part = 1
    elif fractional_part < 3:
        nice_fractional_part = 2
    elif fractional_part < 7:
        nice_fractional_part = 5
    else:
        nice_fractional_part = 10

    nice_interval = nice_fractional_part * 10**exponent

    # Calculate tick values
    tick_values = np.arange(
        start=data_min - (data_min % nice_interval), stop=data_max, step=nice_interval
    )
    # Adjust the last tick value if necessary
    if tick_values[-1] < data_max:
        tick_values = np.append(tick_values, tick_values[-1] + nice_interval)

    return tick_values


def create_fig(
    subdf,
    chrmd_df,
    chrmd_df_grouped,
    genesmd_df,
    ts_data,
    title_chr,
    title_dict_ply,
    grid_color,
    packed,
    plot_background,
    plot_border,
    tick_pos_d,
    ori_tick_pos_d,
    shrinked_bkg,
    shrinked_alpha,
    v_space,
):
    """Generate the figure and axes fitting the data."""

    # Unify titles and start figure
    titles = [title_chr.format(**locals()) for chrom in chrmd_df_grouped.index]
    titles = list(pd.Series(titles).drop_duplicates())
    fig = sp.make_subplots(
        rows=len(titles),
        cols=1,
        row_heights=chrmd_df_grouped["y_height"].to_list(),
        subplot_titles=titles,
    )

    # one subplot per chromosome
    # chrmd_df_grouped = chrmd_df.groupby(
    #     ["Chromosome"], group_keys=False, observed=True
    # ).agg({"n_genes": "sum", "min_max": "first", "y_height": "first"})

    for i in range(len(titles)):
        chrom = chrmd_df_grouped.index[i]
        fig.add_trace(go.Scatter(x=[], y=[]), row=i + 1, col=1)

        # set title format
        fig.layout.annotations[i].update(font=title_dict_ply)

        # set x axis limits
        x_min, x_max = chrmd_df_grouped.loc[chrom]["min_max"]
        x_rang = x_max - x_min
        fig.update_xaxes(
            range=[x_min - 0.05 * x_rang, x_max + 0.05 * x_rang],
            tickformat="d",
            showgrid=True,
            gridcolor=grid_color,
            griddash="dot",
            zeroline=False,
            row=i + 1,
            col=1,
        )  # add 5% to limit coordinates range

        # consider introns off
        if tick_pos_d:
            # get previous default ticks
            chrom_subdf = subdf[subdf["Chromosome"] == chrom]
            original_ticks = list(calculate_ticks(chrom_subdf))

            jump = original_ticks[1] - original_ticks[0]

            # find previous ticks that should be conserved
            to_add_val = []
            # there is data to shrink
            if ori_tick_pos_d[chrom]:
                to_add_val += [
                    tick
                    for tick in original_ticks
                    if tick < min(ori_tick_pos_d[chrom])
                    or tick > max(ori_tick_pos_d[chrom])
                ]
                for ii in range(1, len(ori_tick_pos_d[chrom]) - 1, 2):
                    not_shr0 = ori_tick_pos_d[chrom][ii]
                    not_shr1 = ori_tick_pos_d[chrom][ii + 1]
                    to_add_val += [
                        i for i in original_ticks if not_shr0 < i <= not_shr1
                    ]

            # nothing to shrink
            else:
                to_add_val += original_ticks

            # compute new coordinates of conserved previous ticks
            to_add = to_add_val.copy()
            to_add = cumdelting(to_add, ts_data, chrom)

            # set new ticks
            x_ticks_val = sorted(to_add)
            # do not add ticks beyond adjusted limits
            x_ticks_val = [
                num for num in x_ticks_val if num <= chrmd_df_grouped.loc[chrom]["max"]
            ]
            x_ticks_name = sorted(to_add_val)[: len(x_ticks_val)]

            # set new ticks
            fig.update_xaxes(
                tickvals=x_ticks_val,
                ticktext=x_ticks_name,
                row=i + 1,
                col=1,
            )

        # set y axis limits
        y_min = 0
        y_max = chrmd_df_grouped.loc[chrom]["y_height"]
        y_ticks_val = []
        y_ticks_name = []
        if not packed:
            y_ticks_val = [(i + 0.5) * v_space for i in range(int(y_max / v_space))]
            y_ticks_name_d = (
                genesmd_df[genesmd_df["Chromosome"] == chrom]
                .groupby("pr_ix", group_keys=False, observed=True)
                .groups
            )
            y_ticks_name_d = dict(sorted(y_ticks_name_d.items(), reverse=True))
            y_ticks_name = [list(id)[::-1] + [""] for id in y_ticks_name_d.values()]
            y_ticks_name = [item for sublist in y_ticks_name for item in sublist][:-1]
        fig.update_yaxes(
            range=[y_min - 0.5 * v_space, y_max + 0.5 * v_space],
            fixedrange=True,
            tickvals=y_ticks_val,
            ticktext=y_ticks_name,
            showgrid=False,
            row=i + 1,
            col=1,
        )

        # Add shrink rectangles
        if ts_data:
            rects_df = ts_data[chrom]
            rects_df["cumdelta_end"] = rects_df["cumdelta"]
            rects_df["cumdelta_start"] = rects_df["cumdelta"].shift(
                periods=1, fill_value=0
            )
            rects_df["Start"] -= rects_df["cumdelta_start"]
            rects_df["End"] -= rects_df["cumdelta_end"]

            for a, b, c, d in zip(
                rects_df["Start"],
                rects_df["End"],
                rects_df["cumdelta_start"],
                rects_df["cumdelta_end"],
            ):
                x0, x1 = a, b
                y0, y1 = y_min - 1, y_max + 1
                fig.add_trace(
                    go.Scatter(
                        x=[x0, x1, x1, x0, x0],
                        y=[y0, y0, y1, y1, y0],
                        fill="toself",
                        fillcolor=shrinked_bkg,
                        mode="lines",
                        line={"color": "lightyellow"},
                        text=f"Shrinked region:\n[{x0+c} - {x1+d}]",
                        hoverinfo="text",
                        opacity=shrinked_alpha,
                        line_width=0,
                    ),
                    row=i + 1,
                    col=1,
                )

        # Draw lines separating pr objects
        if chrmd_df["pr_line"].drop_duplicates().max() != 0:
            pr_line_y_l = chrmd_df.loc[chrom]["pr_line"].tolist()
            if isinstance(pr_line_y_l, int):
                pr_line_y_l = [pr_line_y_l]
            # separate items with horizontal lines
            for pr_line_y in pr_line_y_l:
                if pr_line_y != 0:
                    fig.add_trace(
                        go.Scatter(
                            x=[x_min - 0.1 * x_rang, x_max + 0.1 * x_rang],
                            y=[pr_line_y + 0.5 * v_space, pr_line_y + 0.5 * v_space],
                            mode="lines",
                            line=dict(color=plot_border, width=1, dash="solid"),
                            hoverinfo="skip",
                        ),
                        row=i + 1,
                        col=1,
                    )

    return fig


def create_fig_with_vfc(
    subdf,
    vcf,
    chrmd_df,
    genesmd_df,
    ts_data,
    title_chr,
    title_dict_ply,
    packed,
    plot_background,
    tick_pos_d,
    ori_tick_pos_d,
    shrinked_bkg,
    shrinked_alpha,
):
    """xxx"""

    # Store titles
    titles = [title_chr.format(**locals()) for chrom in chrmd_df.index]

    # Update height ratios
    hr = [
        num
        for pair in zip([1] * len(chrmd_df), list(chrmd_df.y_height))
        for num in pair
    ]

    # Start figure
    fig = sp.make_subplots(
        rows=len(chrmd_df) * 2,
        cols=1,
        row_heights=hr,
    )

    # Calculate and set the position for each subplot, and add titles
    for i in range(len(chrmd_df) * 2):
        # store vcf position
        if i % 2 == 0:
            vcfpos = fig["layout"][f"yaxis{i + 1}" if i else "yaxis"]["domain"]
        # move exon plot upwards to be next to its vcf plot
        else:
            exonpos = fig["layout"][f"yaxis{i + 1}" if i else "yaxis"]["domain"]
            incr = vcfpos[0] - exonpos[1]
            exonadj = [pos + incr for pos in exonpos]
            fig.update_layout(
                **{f"yaxis{i + 1}" if i else "yaxis": dict(domain=exonadj)}
            )

        # configure title
        title = titles[i // 2]
        if isinstance(title, str):
            titles[i // 2] = dict(
                xref="paper",
                yref="paper",
                x=0.5,
                y=vcfpos[1] + 0.05,
                xanchor="center",
                yanchor="middle",
                text=title,
                showarrow=False,
                font=dict(size=16),
            )
    fig.update_layout(annotations=titles)

    # One subplot per chromosome and vcf
    for i in range(len(chrmd_df) * 2):
        chrom = chrmd_df.index[i // 2]
        fig.add_trace(go.Scatter(x=[], y=[]), row=i + 1, col=1)

        # set title format
        fig.layout.annotations[i // 2].update(font=title_dict_ply)

        # set x axis limits
        x_min, x_max = chrmd_df.iloc[i // 2]["min_max"]
        x_rang = x_max - x_min
        fig.update_xaxes(
            range=[x_min - 0.05 * x_rang, x_max + 0.05 * x_rang],
            tickformat="d",
            showgrid=True,
            gridcolor="lightgrey",
            griddash="dot",
            zeroline=False,
            row=i + 1,
            col=1,
        )  # add 5% to limit coordinates range

        # consider introns off
        if tick_pos_d:
            # get previous default ticks
            chrom_subdf = subdf[subdf["Chromosome"] == chrom]
            original_ticks = list(calculate_ticks(chrom_subdf))

            jump = original_ticks[1] - original_ticks[0]

            # find previous ticks that should be conserved
            to_add_val = []
            # there is data to shrink
            if ori_tick_pos_d[chrom]:
                to_add_val += [
                    tick
                    for tick in original_ticks
                    if tick < min(ori_tick_pos_d[chrom])
                    or tick > max(ori_tick_pos_d[chrom])
                ]
                for ii in range(1, len(ori_tick_pos_d[chrom]) - 1, 2):
                    not_shr0 = ori_tick_pos_d[chrom][ii]
                    not_shr1 = ori_tick_pos_d[chrom][ii + 1]
                    to_add_val += [
                        i for i in original_ticks if not_shr0 < i <= not_shr1
                    ]

            # nothing to shrink
            else:
                to_add_val += original_ticks

            # compute new coordinates of conserved previous ticks
            to_add = to_add_val.copy()
            to_add = cumdelting(to_add, ts_data, chrom)

            # set new ticks
            x_ticks_val = sorted(to_add)
            # do not add ticks beyond adjusted limits
            x_ticks_val = [
                num for num in x_ticks_val if num <= chrmd_df.loc[chrom]["min_max"][1]
            ]

            x_ticks_name = sorted(to_add_val)[: len(x_ticks_val)]

            # set new ticks
            fig.update_xaxes(
                tickvals=x_ticks_val,
                ticktext=x_ticks_name,
                row=i + 1,
                col=1,
            )

        # hide x labels for vcf plots
        if i % 2 == 1:
            fig.update_xaxes(showticklabels=False, row=i, col=1)

        # set y axis limits
        y_min = 0
        y_max = chrmd_df.iloc[i // 2].y_height
        y_ticks_val = []
        y_ticks_name = []
        if not packed:
            y_ticks_val = [i + 0.5 for i in range(int(y_max))]
            y_ticks_name = genesmd_df.groupby(
                "chrix", group_keys=False, observed=True
            ).groups[chrom]
        fig.update_yaxes(
            range=[y_min, y_max],
            fixedrange=True,
            tickvals=y_ticks_val,
            ticktext=y_ticks_name,
            showgrid=False,
            row=i + 1,
            col=1,
        )

        # Add shrink rectangles
        if ts_data:
            rects_df = ts_data[chrom].copy()
            rects_df["cumdelta_end"] = rects_df["cumdelta"]
            rects_df["cumdelta_start"] = rects_df["cumdelta"].shift(
                periods=1, fill_value=0
            )
            rects_df["Start"] -= rects_df["cumdelta_start"]
            rects_df["End"] -= rects_df["cumdelta_end"]

            for a, b, c, d in zip(
                rects_df["Start"],
                rects_df["End"],
                rects_df["cumdelta_start"],
                rects_df["cumdelta_end"],
            ):
                x0, x1 = a, b
                y0, y1 = y_min - 1, y_max + 1
                fig.add_trace(
                    go.Scatter(
                        x=[x0, x1, x1, x0, x0],
                        y=[y0, y0, y1, y1, y0],
                        fill="toself",
                        fillcolor=shrinked_bkg,
                        mode="lines",
                        line={"color": "lightyellow"},
                        text=f"Shrinked region:\n[{x0 + c} - {x1 + d}]",
                        hoverinfo="text",
                        opacity=shrinked_alpha,
                        line_width=0,
                    ),
                    row=i + 1,
                    col=1,
                )

    return fig
