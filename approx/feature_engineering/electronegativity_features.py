from math import sqrt
from collections import defaultdict

from .base import FeatureModule
from .registry import register_feature
from .stats_expander import StatsExpander
from .ionic_radius_features import lookup_ionic_radius
from mendeleev.econf import ElectronicConfiguration


@register_feature
class ElectronegativityModule(FeatureModule):

    IP_SCALE = 12.0
    EA_SCALE = 8.0
    PAULING_Q_SLOPE = 0.03

    def _Z_eff(self, element):
        row = self.ptable.query("symbol == @element").iloc[0]
        Z = int(row["atomic_number"])
        return Z - sqrt(Z)

    def _valence_electrons(self, element, ox):
        row = self.ptable.query("symbol == @element").iloc[0]
        config = ElectronicConfiguration(row["electronic_configuration"])
        period = int(row["period"])
        ion = config.ionize(ox)

        return (
            ion.conf.get((period, "s"), 0) +
            ion.conf.get((period, "p"), 0) +
            ion.conf.get((period - 1, "d"), 0) +
            ion.conf.get((period - 1, "f"), 0)
        )

    def compute_en(self, element, ox):
        row = self.ptable.query("symbol == @element").iloc[0]

        r = lookup_ionic_radius(element, ox, default=100.0)
        Z_eff = self._Z_eff(element)
        n_val = self._valence_electrons(element, ox)

        # Pauling
        χP0 = float(row["en_pauling"] or 0.0)
        χP = χP0 + self.PAULING_Q_SLOPE * ox

        # Mulliken (DOES NOT WORK CURRENTLY)
        # EA0 = float(row["electron_affinity"] or 0.0)
        # IP = 5.0 + (Z_eff / r) * self.IP_SCALE
        # EA = EA0 - self.EA_SCALE * ox / r
        # χM = (IP + EA) / 2.0

        # Allred–Rochow
        χAR = (Z_eff / (r ** 2)) * 0.359 + 0.744

        # Gordy
        χG = Z_eff / r

        # Metallic bond EN
        χMB = n_val / r

        return {
            "pauling": χP,
            # "mulliken": χM,
            "allred_rochow": χAR,
            "gordy": χG,
            "mb": χMB,
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
        result = self.approx.charge_balance(formula, return_format="dict")
        all_elems, var_elems = self.parse_result(result)

        features = {}

        # ALL
        features |= self._expand_group(all_elems, "all_")

        # VAR
        if var_elems:
            features |= self._expand_group(var_elems, "var_")
        else:
            zero = self._expand_group(all_elems, "all_")
            for k in zero:
                features[k.replace("all_", "var_")] = 0.0

        return features
