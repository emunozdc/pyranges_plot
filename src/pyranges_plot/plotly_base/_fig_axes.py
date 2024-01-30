import plotly.subplots as sp
import plotly.graph_objects as go


def create_fig(
    chrmd_df, genesmd_df, ts_data, chr_string, title_dict_ply, packed, plot_background
):
    """Generate the figure and axes fitting the data."""

    titles = [chr_string.format(**locals()) for chrom in chrmd_df.index]
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
            print(ts_data)
            rects_df = ts_data[str(i + 1)]  #### THIS HAS TO BE CHROMOSOME VALUE not i+1
            rects_df["cumdelta_end"] = rects_df["cumdelta"]
            rects_df["cumdelta_start"] = rects_df["cumdelta"].shift(
                periods=1, fill_value=0
            )
            rects_df["Start"] -= rects_df["cumdelta_start"]
            rects_df["End"] -= rects_df["cumdelta_end"]

            for a, b in zip(rects_df["Start"], rects_df["End"]):
                x0, x1 = a, b
                y0, y1 = y_min - 1, y_max + 1
                fig.add_trace(
                    go.Scatter(
                        x=[x0, x1, x1, x0, x0],
                        y=[y0, y0, y1, y1, y0],
                        fill="tonexty",
                        fillcolor=plot_background,
                        fillpattern_shape="/",
                        fillpattern_solidity=0.3,
                        line={"color": "whitesmoke"},
                        hoverinfo=None,
                        # opacity=0.1,
                        line_width=0,
                    ),
                    row=i + 1,
                    col=1,
                )

    return fig
