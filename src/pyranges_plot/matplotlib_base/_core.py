import tkinter as tk


def coord2percent(ax, X0, X1):
    """Provides the plot percentage length from the points given. Matplotlib friendly"""

    x_min, x_max = ax.get_xlim()
    x_rang = x_max - x_min
    percent_size = float(((X1 - X0) / x_rang))

    return percent_size


def percent2coord(ax, x_percent):
    """Provides the coordinates distance from the plot percentage given. Matplotlib friendly"""

    x_min, x_max = ax.get_xlim()
    x_rang = x_max - x_min
    percent_coord = float(x_percent * x_rang)

    return percent_coord


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
