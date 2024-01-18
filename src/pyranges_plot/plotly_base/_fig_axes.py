import plotly.subplots as sp
import plotly.graph_objects as go


def create_fig(chrmd_df, genesmd_df, chr_string, title_dict_ply, packed):
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

    return fig
