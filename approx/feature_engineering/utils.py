from collections import defaultdict
from functools import lru_cache
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

def _radius_to_pm(value):
    if value is None:
        return None

    if hasattr(value, "to"):
        return float(value.to("pm"))

    numeric = float(value)

    # Plain floats from pymatgen may come back in angstrom.
    # Typical ionic radii in pm are ~50-200, while angstrom values are usually < 10.
    if numeric < 10:
        return numeric * 100.0

    return numeric

@lru_cache(maxsize=512)
def _lookup_ionic_radius_cached(element, ox):
    """
    Safe ionic radius lookup:
    1. Try Element(element).ionic_radii[ox]  -> no warnings
    2. Try Species(element, ox).ionic_radius -> suppress warnings

    Always returns a float in picometers (pm).
    """
    try:
        el = Element(element)
        r = el.ionic_radii.get(ox, None)
        
        if r is not None:
            normalized = _radius_to_pm(r)
            if normalized is not None:
                return normalized
    except Exception:
        pass

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = Species(element, ox).ionic_radius
            
        if r is not None:
            normalized = _radius_to_pm(r)
            if normalized is not None:
                return normalized
    except Exception:
        pass

    return None

def lookup_ionic_radius(element, ox, default=0.0):
    value = _lookup_ionic_radius_cached(element, ox)
    if value is None:
        return float(default)
    return float(value)

    

