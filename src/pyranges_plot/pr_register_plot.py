from pyranges import PyRanges
from .core import set_engine
from .plot_main import plot


def register_plot(engine=None):
    """Register the plot function as a method to PyRanges."""

    if engine is not None:
        set_engine(engine)

    # Attach the wrapper as a method to PyRanges
    PyRanges.plot = plot
