from collections import defaultdict
import warnings
from pymatgen.core import Element, Species

def weighted_average(element_list, value_fn):
    """
    element_list: [(element, ox, qty)]
    value_fn: function(element, ox) → dict
    """
    if len(element_list) == 0:
        return {}

    totals = defaultdict(float)
    total_qty = sum(q for _, _, q in element_list)

    for el, ox, qty in element_list:
        values = value_fn(el, ox)  # must be dict
        for k, v in values.items():
            totals[k] += v * qty

    # Convert to averages
    return {k: totals[k] / total_qty for k in totals}

def lookup_ionic_radius(element, ox, default=0.0):
    """
    Safe ionic radius lookup:
    1. Try Element(element).ionic_radii[ox]  -> no warnings
    2. Try Species(element, ox).ionic_radius -> suppress warnings
    3. Return default if missing

    Always returns a float in picometers (pm).
    """
    try:
        el = Element(element)
        r = el.ionic_radii.get(ox, None)
        
        if r is not None:
            # r may be FloatWithUnit or float
            if hasattr(r, "to"):
                return float(r.to("pm"))
            return float(r)
    except Exception:
        pass

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = Species(element, ox).ionic_radius
            
        if r is not None:
            # again: r may be FloatWithUnit or float
            if hasattr(r, "to"):
                return float(r.to("pm")) #type: ignore[attr-defined]
            return float(r)
    except Exception:
        pass

    return float(default)

