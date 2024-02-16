import plotly.subplots as sp
import plotly.graph_objects as go
import numpy as np
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
    genesmd_df,
    ts_data,
    title_chr,
    title_dict_ply,
    packed,
    plot_background,
    tick_pos_d,
    ori_tick_pos_d,
):
    """Generate the figure and axes fitting the data."""

    titles = [title_chr.format(**locals()) for chrom in chrmd_df.index]
    fig = sp.make_subplots(
        rows=len(chrmd_df),
        cols=1,
        row_heights=chrmd_df.y_height.to_list(),
        subplot_titles=titles,
    )

    # one subplot per chromosome
    for i in range(len(chrmd_df)):
        chrom = chrmd_df.index[i]
        fig.add_trace(go.Scatter(x=[], y=[]), row=i + 1, col=1)

        # set title format
        fig.layout.annotations[i].update(font=title_dict_ply)

        # set x axis limits
        x_min, x_max = chrmd_df.iloc[i]["min_max"]
        x_rang = x_max - x_min
        fig.update_xaxes(
            range=[x_min - 0.05 * x_rang, x_max + 0.05 * x_rang],
            tickformat="d",
            showgrid=True,
            gridcolor="grey",
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
            # x_ticks_val = sorted(tick_pos_d[chrom] + to_add)
            x_ticks_val = sorted(to_add)
            # do not add ticks beyond adjusted limits
            x_ticks_val = [
                num for num in x_ticks_val if num <= chrmd_df.loc[chrom]["min_max"][1]
            ]
            # x_ticks_name = sorted(ori_tick_pos_d[chrom] + to_add_val)[
            #     : len(x_ticks_val)
            # ]
            x_ticks_name = to_add_val[: len(x_ticks_val)]

            # set new ticks
            fig.update_xaxes(
                tickvals=x_ticks_val,
                ticktext=x_ticks_name,
                row=i + 1,
                col=1,
            )

        # set y axis limits
        y_min = 0
        y_max = chrmd_df.iloc[i].y_height
        y_ticks_val = []
        y_ticks_name = []
        if not packed:
            y_ticks_val = [i + 0.5 for i in range(int(y_max))]
            y_ticks_name = genesmd_df.groupby("chrix", group_keys=False).groups[chrom]
        fig.update_yaxes(
            range=[y_min, y_max],
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
                        fillcolor="lightyellow",
                        mode="lines",
                        line={"color": "lightyellow"},
                        text=f"Shrinked region:\n[{x0+c} - {x1+d}]",
                        hoverinfo="text",
                        opacity=0.7,
                        line_width=0,
                    ),
                    row=i + 1,
                    col=1,
                )

    return fig
