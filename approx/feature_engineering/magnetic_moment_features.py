from math import sqrt

from .base import FeatureModule
from .registry import register_feature
from .stats_expander import StatsExpander


@register_feature
class MagneticMomentModule(FeatureModule):
    """
    Computes oxidation-state-dependent magnetic moment features.

    Spin-only magnetic moment:
        mu = sqrt(n_unpaired * (n_unpaired + 2))
    """

    def __init__(self, approx, ptable):
        super().__init__(approx, ptable)
        self._unpaired_cache = {}

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

    def count_unpaired_electrons(self, element: str, ox: int) -> int:
        cache_key = (element, int(ox))
        if cache_key in self._unpaired_cache:
            return self._unpaired_cache[cache_key]

        row = self.get_element_row(element)
        config = self.get_electronic_configuration(element)
        period = int(row["period"])

        if ox >= 0:
            ion = config.ionize(ox)
            occupancy = self._orbital_occupancy(ion, period)
        else:
            occupancy = self._orbital_occupancy(config, period)
            occupancy = self._fill_anion_valence(occupancy, abs(int(ox)))

        electrons = []

        electrons += self.fill_orbitals(occupancy["d"], 5)
        electrons += self.fill_orbitals(occupancy["p"], 3)
        electrons += self.fill_orbitals(occupancy["s"], 1)
        electrons += self.fill_orbitals(occupancy["f"], 7)

        result = sum(1 for e in electrons if e == 1)
        self._unpaired_cache[cache_key] = result
        return result

    @staticmethod
    def fill_orbitals(nelect: int, norb: int):
        occ = [0] * norb

        for i in range(min(nelect, norb)):
            occ[i] = 1

        rem = nelect - norb
        if rem > 0:
            for i in range(min(rem, norb)):
                occ[i] += 1

        return occ

    def _expand_group(self, elems, prefix):
        features = {}

        n_vals, weights = [], []
        mu_vals = []

        for el, ox, qty in elems:
            n = self.count_unpaired_electrons(el, ox)
            mu = sqrt(n * (n + 2))

            n_vals.append(n)
            mu_vals.append(mu)
            weights.append(qty)

        features |= StatsExpander.expand(
            n_vals, weights, prefix=f"{prefix}unpaired_electrons_"
        )

        features |= StatsExpander.expand(
            mu_vals, weights, prefix=f"{prefix}spin_moment_"
        )

        return features

    def get_features(self, formula: str):
        all_elems, var_elems = self.get_parsed_elements(formula)

        features = {}
        features |= self._expand_group(all_elems, "all_")

        if var_elems:
            features |= self._expand_group(var_elems, "var_")
        else:
            features |= self.zero_fill_from_all(features)

        return features
