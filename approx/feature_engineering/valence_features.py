from .base import FeatureModule
from .registry import register_feature
from .stats_expander import StatsExpander


@register_feature
class ValenceFeatureModule(FeatureModule):
    def __init__(self, approx, ptable):
        super().__init__(approx, ptable)
        self._valence_cache = {}

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

    def get_valence(self, element, ox):
        cache_key = (element, int(ox))
        if cache_key in self._valence_cache:
            return self._valence_cache[cache_key]

        row = self.get_element_row(element)
        config = self.get_electronic_configuration(element)
        period = int(row["period"])

        if ox >= 0:
            ion = config.ionize(ox)
            occupancy = self._orbital_occupancy(ion, period)
        else:
            occupancy = self._orbital_occupancy(config, period)
            occupancy = self._fill_anion_valence(occupancy, abs(int(ox)))

        result = occupancy | {"total": sum(occupancy.values())}
        self._valence_cache[cache_key] = result
        return result

    def get_unfilled(self, val):
        maxc = {"s": 2, "p": 6, "d": 10, "f": 14}
        return {orb: maxc[orb] - val[orb] for orb in maxc} | {
            "total": sum(maxc[o] - val[o] for o in maxc)
        }

    def _expand_group(self, elems, prefix):
        out = {}
        weights = [qty for _, _, qty in elems]
        filled_values = []
        unfilled_values = []

        for el, ox, _ in elems:
            filled = self.get_valence(el, ox)
            filled_values.append(filled)
            unfilled_values.append(self.get_unfilled(filled))

        for channel in ["s", "p", "d", "f", "total"]:
            out |= StatsExpander.expand(
                [vals[channel] for vals in filled_values],
                weights,
                prefix=f"{prefix}valence_{channel}_"
            )

            out |= StatsExpander.expand(
                [vals[channel] for vals in unfilled_values],
                weights,
                prefix=f"{prefix}unfilled_valence_{channel}_"
            )

        return out

    def get_features(self, formula):
        all_elems, var_elems = self.get_parsed_elements(formula)

        features = {}
        features |= self._expand_group(all_elems, "all_")

        if var_elems:
            features |= self._expand_group(var_elems, "var_")
        else:
            features |= self.zero_fill_from_all(features)

        return features
