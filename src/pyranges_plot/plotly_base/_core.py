
def coord2percent(fig, trace, X0, X1):
    """Provides the plot percentage length from the points given. Plotly friendly"""

    if trace == 1:
        trace = ""
    else:
        trace = str(trace)

    x_min, x_max = fig["layout"]["xaxis" + trace]["range"]
    x_rang = x_max - x_min

    percent_size = float(((X1 - X0) / x_rang))

    return percent_size



def percent2coord(fig, trace, x_percent):
    """Provides the coordinates distance from the plot percentage given. Plt friendly"""

    if trace == 1:
        trace = ""
    else:
        trace = str(trace)

    x_min, x_max = fig["layout"]["xaxis" + trace]["range"]
    x_rang = x_max - x_min

    percent_coord = float(x_percent * x_rang)

    return percent_coord
