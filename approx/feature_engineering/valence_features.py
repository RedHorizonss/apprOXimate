from .base import FeatureModule
from .registry import register_feature
from .stats_expander import StatsExpander

from mendeleev.econf import ElectronicConfiguration


@register_feature
class ValenceFeatureModule(FeatureModule):

    def get_valence(self, element, ox):
        row = self.ptable.query("symbol == @element").iloc[0]
        config = ElectronicConfiguration(row["electronic_configuration"])
        period = int(row["period"])

        ion = config.ionize(ox)

        s = ion.conf.get((period, "s"), 0)
        p = ion.conf.get((period, "p"), 0)
        d = ion.conf.get((period - 1, "d"), 0)
        f = ion.conf.get((period - 1, "f"), 0)

        return {"s": s, "p": p, "d": d, "f": f, "total": s + p + d + f}

    def get_unfilled(self, val):
        maxc = {"s": 2, "p": 6, "d": 10, "f": 14}
        return {orb: maxc[orb] - val[orb] for orb in maxc} | {
            "total": sum(maxc[o] - val[o] for o in maxc)
        }

    def _expand_group(self, elems, prefix):
        """
        elems: list of (element, ox, qty)
        prefix: 'all_' or 'var_'
        """
        out = {}

        for channel in ["s", "p", "d", "f", "total"]:
            # filled
            values, weights = [], []
            for el, ox, qty in elems:
                v = self.get_valence(el, ox)[channel]
                values.append(v)
                weights.append(qty)

            out |= StatsExpander.expand(
                values,
                weights,
                prefix=f"{prefix}valence_{channel}_"
            )

            # unfilled
            values, weights = [], []
            for el, ox, qty in elems:
                u = self.get_unfilled(self.get_valence(el, ox))[channel]
                values.append(u)
                weights.append(qty)

            out |= StatsExpander.expand(
                values,
                weights,
                prefix=f"{prefix}unfilled_valence_{channel}_"
            )

        return out

    def get_features(self, formula):
        result = self.approx.charge_balance(formula, return_format="dict")
        all_elems, var_elems = self.parse_result(result)

        features = {}

        # ALL elements
        features |= self._expand_group(all_elems, "all_")

        # VAR elements
        if var_elems:
            features |= self._expand_group(var_elems, "var_")
        else:
            # if no variable elements → zero-fill var stats
            for k, v in self._expand_group(all_elems, "all_").items():
                features[k.replace("all_", "var_")] = 0.0

        return features
