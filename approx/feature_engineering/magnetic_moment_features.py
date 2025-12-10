from collections import defaultdict
from .base import FeatureModule
from math import sqrt
from .registry import register_feature

#Might be worth making the config, period and ion variables into a utility function?
from mendeleev.econf import ElectronicConfiguration

@register_feature
class MagneticMomentModule(FeatureModule):
    """
    Computes magnetic-moment-related features based on oxidation-state–adjusted
    electronic configurations.

    -------------------------------------------
    THEORY / EQUATIONS USED
    -------------------------------------------

    (1) Determine number of unpaired electrons:
        n_unpaired = number of singly-occupied orbitals
                      in d, s, p, f shells after ionization

    (2) Spin quantum number:
        S = n_unpaired / 2

    (3) Spin-only magnetic moment (μ_B):
        μ = sqrt(n_unpaired * (n_unpaired + 2))

        This is the standard approximation used for transition metals,
        assuming quenching of orbital angular momentum (valid for oxides).

    These quantities are oxidation-state dependent because the electronic
    configuration changes after ionization.
    """

    # --------------------------------------------
    # Extract unpaired-electron count
    # --------------------------------------------
    def count_unpaired_electrons(self, element: str, ox: int) -> int:
        """
        Count unpaired electrons after oxidation-state ionization.
        Uses your ElectronicConfiguration class.

        Returns:
            n_unpaired (int)
        """
        cfg_str = self.ptable.query("symbol == @element")["electronic_configuration"].values[0]
        period = self.ptable.query("symbol == @element")["period"].values[0]

        config = ElectronicConfiguration(cfg_str)
        ion = config.ionize(ox)

        # Collect electrons from orbitals relevant to magnetism
        electrons = []

        # d orbitals (period - 1)
        d_e = ion.conf.get((period - 1, "d"), 0)
        electrons += self.fill_orbitals(d_e, 5)

        # p orbitals (period)
        p_e = ion.conf.get((period, "p"), 0)
        electrons += self.fill_orbitals(p_e, 3)

        # s orbitals (period)
        s_e = ion.conf.get((period, "s"), 0)
        electrons += self.fill_orbitals(s_e, 1)

        # f orbitals rarely relevant here, but included for completeness
        f_e = ion.conf.get((period - 1, "f"), 0)
        electrons += self.fill_orbitals(f_e, 7)

        # Count unpaired electrons
        n_unpaired = sum(1 for e in electrons if e == 1)
        return n_unpaired
    
    @staticmethod
    def fill_orbitals(nelect: int, norb: int):
        """
        Fills orbitals following Hund's rule.

        Args:
            nelect: total electrons in the subshell
            norb: number of orbitals (1=s, 3=p, 5=d, 7=f)

        Returns:
            A list of orbital occupancies:
            [1,1,1,0,0] etc.
        """
        occ = [0] * norb

        # First pass: singly fill
        for i in range(min(nelect, norb)):
            occ[i] = 1

        # Second pass: pair up remaining electrons
        rem = nelect - norb
        if rem > 0:
            for i in range(min(rem, norb)):
                occ[i] += 1

        return occ

    def get_features(self, formula: str, mode: str = "both"):
        result = self.approx.charge_balance(formula, return_format="dict")
        all_elems, var_elems = self.parse_result(result)

        features = {}

        # Weighted helper
        def weighted_vals(elems, prefix):
            total_qty = sum(q for _, _, q in elems)
            if total_qty == 0:
                return {
                    f"{prefix}unpaired_electrons": 0.0,
                    f"{prefix}spin_moment": 0.0,
                }

            total_unpaired = 0.0
            total_moment = 0.0

            for el, ox, qty in elems:
                n_unpaired = self.count_unpaired_electrons(el, ox)
                # Spin-only magnetic moment
                mu = sqrt(n_unpaired * (n_unpaired + 2))

                total_unpaired += n_unpaired * qty
                total_moment += mu * qty

            return {
                f"{prefix}unpaired_electrons": total_unpaired / total_qty,
                f"{prefix}spin_moment": total_moment / total_qty,
            }

        # --- ALL elements ---
        if mode in ("all", "both"):
            features.update(weighted_vals(all_elems, "all_"))

        # --- VARIABLE oxidation-state elements ---
        if mode in ("var", "both"):
            if var_elems:
                features.update(weighted_vals(var_elems, "var_"))
            else:
                features.update({
                    "var_unpaired_electrons": 0.0,
                    "var_spin_moment": 0.0,
                })

        return features
