from collections import defaultdict
from math import sqrt

from .base import FeatureModule
from .ionic_radius_features import lookup_ionic_radius
from .registry import register_feature
from .stats_expander import StatsExpander


#TODO: Fix pauling EN

@register_feature
class ElectronegativityModule(FeatureModule):
    IP_SCALE = 12.0
    EA_SCALE = 8.0
    PAULING_Q_SLOPE = 0.03

    def __init__(self, approx, ptable):
        super().__init__(approx, ptable)
        self._zeff_cache = {}
        self._valence_electron_cache = {}

    @staticmethod
    def _orbital_occupancy(config, period):
        return {
            "s": config.conf.get((period, "s"), 0),
            "p": config.conf.get((period, "p"), 0),
            "d": config.conf.get((period - 1, "d"), 0),
            "f": config.conf.get((period - 1, "f"), 0),
        }

    @staticmethod
    def _fill_anion_valence(occupancy, electrons_to_add):
        capacities = {"s": 2, "p": 6, "d": 10, "f": 14}

        for orbital in ("s", "p", "d", "f"):
            current = occupancy[orbital]
            room = max(capacities[orbital] - current, 0)
            added = min(electrons_to_add, room)
            occupancy[orbital] += added
            electrons_to_add -= added

            if electrons_to_add == 0:
                break

        return occupancy

    def _Z_eff(self, element):
        if element in self._zeff_cache:
            return self._zeff_cache[element]
        row = self.get_element_row(element)
        Z = int(row["atomic_number"])
        value = Z - sqrt(Z)
        self._zeff_cache[element] = value
        return value

    def _valence_electrons(self, element, ox):
        cache_key = (element, int(ox))
        if cache_key in self._valence_electron_cache:
            return self._valence_electron_cache[cache_key]

        row = self.get_element_row(element)
        config = self.get_electronic_configuration(element)
        period = int(row["period"])

        if ox >= 0:
            ion = config.ionize(ox)
            occupancy = self._orbital_occupancy(ion, period)
        else:
            occupancy = self._orbital_occupancy(config, period)
            occupancy = self._fill_anion_valence(occupancy, abs(int(ox)))

        value = (
            occupancy["s"] +
            occupancy["p"] +
            occupancy["d"] +
            occupancy["f"]
        )
        self._valence_electron_cache[cache_key] = value
        return value

    def compute_en(self, element, ox):
        row = self.get_element_row(element)

        r = lookup_ionic_radius(element, ox, default=100.0)
        Z_eff = self._Z_eff(element)
        n_val = self._valence_electrons(element, ox)

        pauling_0 = float(row["en_pauling"] or 0.0)
        pauling = pauling_0 + self.PAULING_Q_SLOPE * ox

        allred_rochow = (Z_eff / (r ** 2)) * 0.359 + 0.744
        gordy = Z_eff / r
        metallic_bond = n_val / r

        return {
            "pauling": pauling,
            "allred_rochow": allred_rochow,
            "gordy": gordy,
            "mb": metallic_bond,
        }

    def _expand_group(self, elems, prefix):
        values = defaultdict(list)
        weights = []

        for el, ox, qty in elems:
            en = self.compute_en(el, ox)
            for k, v in en.items():
                values[k].append(v)
            weights.append(qty)

        features = {}
        for k, vals in values.items():
            features |= StatsExpander.expand(
                vals,
                weights,
                prefix=f"{prefix}{k}_en_"
            )

        return features

    def get_features(self, formula):
        all_elems, var_elems = self.get_parsed_elements(formula)

        features = {}
        features |= self._expand_group(all_elems, "all_")

        if var_elems:
            features |= self._expand_group(var_elems, "var_")
        else:
            features |= self.zero_fill_from_all(features)

        return features
