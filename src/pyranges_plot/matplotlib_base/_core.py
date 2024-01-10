import numpy as np
import tkinter as tk


def coord2inches(fig, ax, X0, X1, Y0, Y1):
    """Provides the inches length from the points given. Plt friendly"""

    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    fig_width, fig_height = fig.get_size_inches()

    x_scale = fig_width / (x_max - x_min)
    y_scale = fig_height / (y_max - y_min)

    inch_len = float(np.sqrt(((X1 - X0) * x_scale) ** 2 + ((Y1 - Y0) * y_scale) ** 2))

    return inch_len


def inches2coord(fig, ax, x_inches):
    """Provides the coordinates distance from the inches given. Plt friendly"""

    x_min, x_max = ax.get_xlim()
    fig_width, fig_height = fig.get_size_inches()
    x_scale = fig_width / (x_max - x_min)

    cord_size = float(x_inches / x_scale)

    return cord_size


def plt_popup_warning(txt, bkg="#1f1f1f", txtcol="white", botcol="#D6AA00"):
    """Create warning window for Matplotlib plots."""

    warn = tk.Tk()

    # Title and background
    warn.wm_title("Warning!")
    warn.configure(background=bkg)

    # Label for warning text
    label = tk.Label(warn, text=txt, font=("Sans", 15), fg=txtcol, bg=bkg)
    label.pack(side="top", anchor="center", pady=10)

    # Button
    bot = tk.Button(
        warn, text="Got it", command=lambda: warn.destroy(), fg="black", bg=botcol
    )
    bot.pack(pady=10)

    # Start main loop
    warn.mainloop()
