import sys, os

# add the parent directory to PYTHONPATH 
# KEEP THIS IN FOR NOW UNTIL WE MAKE A TOML FILE SO WE CAN JUST INSTALL THE PACKAGE
sys.path.append(os.path.abspath(".."))
from approx.approximate import * 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Set scientific publication style for Energy Storage plot
plt.rcParams.update({
    # Font
    'font.size': 14,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'font.family': ['Helvetica'],
    'mathtext.fontset': 'custom',
    'mathtext.it': 'Helvetica:italic',
    'mathtext.bf': 'Helvetica:bold',

    # Axes and lines
    'lines.linewidth': 2,
    'lines.markersize': 6,
    'axes.linewidth': 2,

    # Ticks
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'xtick.major.width': 2,
    'ytick.major.width': 2,
    'xtick.minor.width': 1,
    'ytick.minor.width': 1,
    'xtick.major.size': 6,
    'ytick.major.size': 6,
    'xtick.minor.size': 3,
    'ytick.minor.size': 3,

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


def generate_alkali_values(alkali_amount, interval, max_range):
    """Generate array of alkali metal values to test."""
    return np.arange(interval, min(alkali_amount, max_range), interval)


def get_colormap_colors(num_items, colormap='viridis', cmap_range=(0.1, 0.9)):
    """Generate colors from a colormap within a specified range."""
    cmap = plt.cm.get_cmap(colormap)
    
    if num_items == 1:
        return [cmap((cmap_range[0] + cmap_range[1]) / 2)]
    else:
        return [cmap(cmap_range[0] + (cmap_range[1] - cmap_range[0]) * i / (num_items - 1)) 
                for i in range(num_items)]


def setup_figure(figsize=(6, 5)):
    """Create figure and axes with standard settings."""
    return plt.subplots(figsize=figsize)


def finalize_plot(ax, xlabel, ylabel, invert_x=False):
    """Apply common plot formatting and display."""
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if invert_x:
        ax.invert_xaxis()
    ax.legend(loc='best')
    plt.tight_layout()
    plt.show()


def calculate_charge_balances(analyzer, parsed_formula, alkali_metal, values):
    """Calculate charge balance for each alkali metal value."""
    charge_balances = []
    
    for val in values:
        try:
            modified_formula = {**parsed_formula, alkali_metal: val}
            result = analyzer.charge_balance(modified_formula, return_format='object')
            charge_balances.append(result.final_charge if result else np.nan)
        except Exception:
            charge_balances.append(np.nan)
    
    return charge_balances


def extract_transition_metal_data(result, alkali_metal):
    """Extract transition metal data from charge balance result, excluding alkali metals and oxygen."""
    alkali_metals = {'Li', 'Na', 'K', 'Rb', 'Cs', 'O'}
    transition_metals = {}
    
    for element, data in result['elements'].items():
        if element not in alkali_metals and element != alkali_metal:
            transition_metals[element] = data
    
    return transition_metals


def calculate_weighted_average_oxidation_state(data):
    """Calculate weighted average oxidation state from element data."""
    if 'states' in data:
        total_quantity = sum(state['quantity'] for state in data['states'])
        if total_quantity == 0:
            return np.nan
        weighted_sum = sum(state['oxidation_state'] * state['quantity'] 
                         for state in data['states'])
        return weighted_sum / total_quantity
    else:
        return data['oxidation_state']


def add_panel_label(ax, label, x=0.05, y=0.95, fontsize=16, fontweight='bold'):
    """
    Add a panel label (a, b, c, etc.) to the top-left inside corner of an axis.
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axis to add the label to
    label : str
        The label text (e.g., '(a)', '(b)')
    x : float
        X position in axis coordinates (0-1)
    y : float
        Y position in axis coordinates (0-1)
    fontsize : int
        Font size for the label
    fontweight : str
        Font weight ('bold', 'normal', etc.)
    """
    ax.text(x, y, label, transform=ax.transAxes, 
            fontsize=fontsize, fontweight=fontweight, va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                     edgecolor='none', alpha=0.8))


def create_multi_panel_figure(nrows=1, ncols=2, figsize=None, 
                               width_ratios=None, height_ratios=None,
                               hspace=0.3, wspace=0.3):
    """
    Create a multi-panel figure for side-by-side plots.
    
    Parameters:
    -----------
    nrows : int
        Number of rows
    ncols : int
        Number of columns
    figsize : tuple, optional
        Figure size (width, height). If None, calculated automatically
    width_ratios : list, optional
        Relative widths of columns
    height_ratios : list, optional
        Relative heights of rows
    hspace : float
        Height spacing between subplots
    wspace : float
        Width spacing between subplots
    
    Returns:
    --------
    fig : matplotlib.figure.Figure
    axes : array of matplotlib.axes.Axes
    """
    if figsize is None:
        figsize = (6 * ncols, 5 * nrows)
    
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(nrows, ncols, figure=fig, 
                  width_ratios=width_ratios, height_ratios=height_ratios,
                  hspace=hspace, wspace=wspace)
    
    axes = np.array([[fig.add_subplot(gs[i, j]) for j in range(ncols)] 
                     for i in range(nrows)])
    
    if nrows == 1 and ncols == 1:
        return fig, axes[0, 0]
    elif nrows == 1:
        return fig, axes[0]
    elif ncols == 1:
        return fig, axes[:, 0]
    else:
        return fig, axes


def plot_charge_balance_simple_ax(ax, analyzer, formula='NaFeO2', alkali_metal='Na', 
                                   interval=0.01, max_range=2.0):
    """Plot charge balance on provided axis (for multi-panel figures)."""
    parsed_formula = analyzer.parse_formula(formula)
    alkali_amount = parsed_formula.get(alkali_metal, 0)
    values = generate_alkali_values(alkali_amount, interval, max_range)
    charge_balances = calculate_charge_balances(analyzer, parsed_formula, alkali_metal, values)
    
    ax.plot(values, charge_balances, 'o', markersize=4, color='#e55c30', 
            label=f'{alkali_metal} Charge Balance')
    ax.axhline(y=0, color='grey', linestyle='--', label='Charge Neutral Line', alpha=0.8)
    
    ax.set_xlabel(f'Amount of {alkali_metal}')
    ax.set_ylabel('Charge Balance')
    ax.legend(loc='best')


def plot_charge_balance_substitution_ax(ax, analyzer, base_formula, metal_A, metal_B, 
                                        fractions, alkali_metal='Na', interval=0.015, 
                                        max_range=2.0, cmap_name='viridis'):
    """Plot charge balance substitution on provided axis (for multi-panel figures)."""
    parsed_base = analyzer.parse_formula(base_formula)
    if not parsed_base:
        raise ValueError(f"Could not parse base formula '{base_formula}'")
    if metal_A not in parsed_base:
        raise ValueError(f"{metal_A} not found in base formula.")
    
    base_amount_A = parsed_base[metal_A]
    fractions = np.array(fractions)
    colors = get_colormap_colors(len(fractions), cmap_name)
    
    for idx, frac in enumerate(fractions):
        new_amount_A = base_amount_A * (1 - frac)
        new_amount_B = base_amount_A * frac
        
        formula_dict = parsed_base.copy()
        formula_dict[metal_A] = new_amount_A
        formula_dict[metal_B] = new_amount_B
        
        label = f"{metal_A}{new_amount_A:.2f}{metal_B}{new_amount_B:.2f}"
        
        values = generate_alkali_values(formula_dict.get(alkali_metal, 0), interval, max_range)
        charge_balances = calculate_charge_balances(analyzer, formula_dict, alkali_metal, values)
        
        ax.plot(values, charge_balances, 'o', markersize=4, color=colors[idx], label=label)
    
    ax.axhline(0, color='gray', linestyle='--')
    ax.set_xlabel(f'Amount of {alkali_metal}')
    ax.set_ylabel('Charge Balance')
    ax.legend(loc='best')


def plot_avg_oxidation_states_ax(ax, analyzer, formula='NaMnO2', alkali_metal='Na', 
                                  interval=0.015, max_range=2.0, colormap='viridis', 
                                  cmap_range=(0.1, 0.9)):
    """Plot average oxidation states on provided axis (for multi-panel figures)."""
    parsed_formula = analyzer.parse_formula(formula)
    alkali_amount = parsed_formula.get(alkali_metal, 0)
    values = generate_alkali_values(alkali_amount, interval, max_range)
    
    transition_metals_data = {}
    
    for val in values:
        try:
            modified_formula = {**parsed_formula, alkali_metal: val}
            result = analyzer.charge_balance(modified_formula, return_format='dict')
            
            if result:
                tm_data = extract_transition_metal_data(result, alkali_metal)
                
                for element, data in tm_data.items():
                    if element not in transition_metals_data:
                        transition_metals_data[element] = []
                    
                    avg_ox_state = calculate_weighted_average_oxidation_state(data)
                    transition_metals_data[element].append(avg_ox_state)
        except Exception:
            for element in transition_metals_data:
                if len(transition_metals_data[element]) < len(values):
                    transition_metals_data[element].append(np.nan)
    
    colors = get_colormap_colors(len(transition_metals_data), colormap, cmap_range)
    
    for i, (metal, ox_states) in enumerate(transition_metals_data.items()):
        ax.plot(values[:len(ox_states)], ox_states, 'o-', markersize=4, 
                color=colors[i], label=metal, linewidth=2, alpha=0.6)
    
    ax.set_xlabel(f'Fractional amount of {alkali_metal}')
    ax.set_ylabel('Average Oxidation State')
    ax.invert_xaxis()
    ax.legend(loc='best')


def plot_individual_oxidation_states_ax(ax, analyzer, formula='NaMnO2', alkali_metal='Na', 
                                        interval=0.01, max_range=2.0, colormap='viridis', 
                                        cmap_range=(0.1, 0.9)):
    """Plot individual oxidation states on provided axis (for multi-panel figures)."""
    parsed_formula = analyzer.parse_formula(formula)
    alkali_amount = parsed_formula.get(alkali_metal, 0)
    values = generate_alkali_values(alkali_amount, interval, max_range)
    
    oxidation_state_data = {}
    
    for val in values:
        try:
            modified_formula = {**parsed_formula, alkali_metal: val}
            result = analyzer.charge_balance(modified_formula, return_format='dict')
            
            if result:
                tm_data = extract_transition_metal_data(result, alkali_metal)
                
                for element, data in tm_data.items():
                    if 'states' in data:
                        for state in data['states']:
                            key = (element, state['oxidation_state'])
                            if key not in oxidation_state_data:
                                oxidation_state_data[key] = []
                            oxidation_state_data[key].append(state['quantity'])
                    else:
                        key = (element, data['oxidation_state'])
                        if key not in oxidation_state_data:
                            oxidation_state_data[key] = []
                        oxidation_state_data[key].append(data['quantity'])
                
                current_length = len([v for v in values if v <= val])
                for key in oxidation_state_data:
                    while len(oxidation_state_data[key]) < current_length:
                        oxidation_state_data[key].append(0)
                        
        except Exception as e:
            print(f"Error at {alkali_metal}={val}: {e}")
    
    max_length = len(values)
    for key in oxidation_state_data:
        while len(oxidation_state_data[key]) < max_length:
            oxidation_state_data[key].append(0)
    
    colors = get_colormap_colors(len(oxidation_state_data), colormap, cmap_range)
    
    for i, ((element, ox_state), quantities) in enumerate(sorted(oxidation_state_data.items())):
        label = f'{element}$^{{{ox_state:+d}}}$'
        ax.plot(values[:len(quantities)], quantities, 'o-', markersize=4,
                color=colors[i], label=label, linewidth=2, alpha=0.6)
    
    ax.set_ylim(-0.03, 1.003)
    ax.set_xlabel(f'Fractional amount of {alkali_metal}')
    ax.set_ylabel('Fractional amount for each oxidation state')
    ax.invert_xaxis()
    ax.legend(loc='best')


# Original standalone plotting functions remain unchanged
def plot_charge_balance_simple(analyzer, formula='NaFeO2', alkali_metal='Na', 
                               interval=0.01, max_range=2.0):
    """Plot charge balance vs alkali metal content for a single formula."""
    fig, ax = setup_figure()
    plot_charge_balance_simple_ax(ax, analyzer, formula, alkali_metal, interval, max_range)
    plt.tight_layout()
    plt.show()


def plot_charge_balance_substitution(analyzer, base_formula, metal_A, metal_B, 
                                     fractions, alkali_metal='Na', interval=0.015, 
                                     max_range=2.0, cmap_name='viridis'):
    """Plot charge balance for different metal substitution fractions."""
    fig, ax = setup_figure()
    plot_charge_balance_substitution_ax(ax, analyzer, base_formula, metal_A, metal_B,
                                        fractions, alkali_metal, interval, max_range, cmap_name)
    plt.tight_layout()
    plt.show()


def plot_avg_oxidation_states(analyzer, formula='NaMnO2', alkali_metal='Na', 
                              interval=0.015, max_range=2.0, colormap='viridis', 
                              cmap_range=(0.1, 0.9)):
    """Plot average oxidation states of transition metals vs alkali metal content."""
    fig, ax = setup_figure()
    plot_avg_oxidation_states_ax(ax, analyzer, formula, alkali_metal, 
                                  interval, max_range, colormap, cmap_range)
    plt.tight_layout()
    plt.show()


def plot_individual_oxidation_states(analyzer, formula='NaMnO2', alkali_metal='Na', 
                                     interval=0.01, max_range=2.0, colormap='viridis', 
                                     cmap_range=(0.1, 0.9)):
    """Plot individual oxidation state quantities vs alkali metal content."""
    fig, ax = setup_figure()
    plot_individual_oxidation_states_ax(ax, analyzer, formula, alkali_metal,
                                        interval, max_range, colormap, cmap_range)
    plt.tight_layout()
    plt.show()