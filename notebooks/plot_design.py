import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

for f in fm.findSystemFonts():
    if "Aptos" in f:
        fm.fontManager.addfont(f)


# Set scientific publication style
plt.rcParams.update({
    # Font
    'font.size': 14,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'font.family': ['Aptos'],
    'mathtext.fontset': 'custom',
    'mathtext.it': 'Helvetica:italic',
    'mathtext.bf': 'Helvetica:bold',

    # Axes and lines
    'lines.linewidth': 1.5,
    'lines.markersize': 4,
    'axes.linewidth': 1.5,

    # Ticks
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'xtick.major.width': 1,
    'ytick.major.width': 1,
    'xtick.minor.width': 0.5,
    'ytick.minor.width': 0.5,
    'xtick.major.size': 4,
    'ytick.major.size': 4,
    'xtick.minor.size': 1.5,
    'ytick.minor.size': 1.5,

    # Legend
    'legend.frameon': False,
    'legend.handlelength': 2,
    'legend.handletextpad': 0.5,

    # Savefig
    'figure.dpi': 100,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02,
})