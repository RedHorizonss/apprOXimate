from math import sqrt
from .base import FeatureModule
from .registry import register_feature
from .stats_expander import StatsExpander
from mendeleev.econf import ElectronicConfiguration


@register_feature
class MagneticMomentModule(FeatureModule):
    """
    Computes oxidation-state-dependent magnetic moment features.

    Spin-only magnetic moment:
        μ = sqrt(n_unpaired * (n_unpaired + 2))
    """

    def count_unpaired_electrons(self, element: str, ox: int) -> int:
        row = self.ptable.query("symbol == @element").iloc[0]
        config = ElectronicConfiguration(row["electronic_configuration"])
        period = int(row["period"])

        ion = config.ionize(ox)

        electrons = []

        # d orbitals
        d_e = ion.conf.get((period - 1, "d"), 0)
        electrons += self.fill_orbitals(d_e, 5)

        # p orbitals
        p_e = ion.conf.get((period, "p"), 0)
        electrons += self.fill_orbitals(p_e, 3)

        # s orbitals
        s_e = ion.conf.get((period, "s"), 0)
        electrons += self.fill_orbitals(s_e, 1)

        # f orbitals
        f_e = ion.conf.get((period - 1, "f"), 0)
        electrons += self.fill_orbitals(f_e, 7)

        return sum(1 for e in electrons if e == 1)

    @staticmethod
    def fill_orbitals(nelect: int, norb: int):
        occ = [0] * norb

        # Hund's rule: single-fill first
        for i in range(min(nelect, norb)):
            occ[i] = 1

        # Pair remaining electrons
        rem = nelect - norb
        if rem > 0:
            for i in range(min(rem, norb)):
                occ[i] += 1

        return occ

    def _expand_group(self, elems, prefix):
        """
        elems: list of (element, ox, qty)
        prefix: 'all_' or 'var_'
        """
        features = {}

        # --- unpaired electrons ---
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
        result = self.approx.charge_balance(formula, return_format="dict")
        all_elems, var_elems = self.parse_result(result)

        features = {}

        # ALL elements
        features |= self._expand_group(all_elems, "all_")

        # VARIABLE elements
        if var_elems:
            features |= self._expand_group(var_elems, "var_")
        else:
            # zero-fill var stats if none exist
            zeroed = self._expand_group(all_elems, "all_")
            for k in zeroed:
                features[k.replace("all_", "var_")] = 0.0

        return features
