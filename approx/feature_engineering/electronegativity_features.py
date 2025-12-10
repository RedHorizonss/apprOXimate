from collections import defaultdict
from math import sqrt
from .registry import register_feature

from .base import FeatureModule
from .ionic_radius_features import lookup_ionic_radius
from .valence_features import ValenceFeatureModule

@register_feature
class ElectronegativityModule(FeatureModule):
    """
    Oxidation-state–sensitive electronegativity features.

    Computed:
        - Pauling EN (oxidation-corrected)
        - Mulliken EN (approx IP + EA)/2
        - Allred–Rochow EN
        - Gordy EN
        - Metallic-bond EN

    Your mendeleev table DOES NOT contain:
        - ionization energies
        - Allred–Rochow
        - Gordy
        - Mulliken

    So we compute these using ionic radius + Z_eff approximations.
    """

    # Factors for Mulliken EN approximation
    IP_SCALE = 12.0
    EA_SCALE = 8.0

    PAULING_Q_SLOPE = 0.03  # empirical charge correction

    # --------------------------
    # Helper: compute Z_eff
    # --------------------------
    def _Z_eff(self, element):
        row = self.ptable.query("symbol == @element")
        Z = int(row["atomic_number"].iloc[0])
        S = sqrt(Z)
        return Z - S


    # --------------------------
    # Main per-element EN computation
    # --------------------------
    def compute_en(self, element, ox):

        row = self.ptable.query("symbol == @element")

        # ionic radius (pm)
        r = lookup_ionic_radius(element, ox, default=100.0)

        # effective nuclear charge
        Z_eff = self._Z_eff(element)

        # valence electrons
        vmod: ValenceFeatureModule = ValenceFeatureModule(self.approx, self.ptable)
        val = vmod.get_valence(element, ox) # type: ignore
        n_valence = val["total"]

        # ---------------------------
        # Pauling EN (corrected)
        # ---------------------------
        χP0 = float(row["en_pauling"].values[0] or 0.0)
        χP = χP0 + self.PAULING_Q_SLOPE * ox

        # ---------------------------
        # Mulliken EN
        # ---------------------------
        EA0 = float(row["electron_affinity"].values[0] or 0.0)

        # approximate IP using Z_eff and ionic radius
        IP_approx = 5.0 + (Z_eff / r) * self.IP_SCALE
        EA_corr = EA0 - self.EA_SCALE * ox / r

        χM = (IP_approx + EA_corr) / 2.0

        # ---------------------------
        # Allred–Rochow EN
        # χ_AR = (Z_eff / r^2)*0.359 + 0.744
        # ---------------------------
        χAR = (Z_eff / (r ** 2)) * 0.359 + 0.744

        # ---------------------------
        # Gordy EN
        # χ_G = Z_eff / r
        # ---------------------------
        χG = Z_eff / r

        # ---------------------------
        # Metallic EN
        # χ_MB = n_valence / r
        # ---------------------------
        χMB = n_valence / r

        return {
            "pauling_en": χP,
            "mulliken_en": χM,
            "allred_rochow_en": χAR,
            "gordy_en": χG,
            "mb_en": χMB,
        }

    # --------------------------
    # Weighted averaging
    # --------------------------
    def weighted(self, elems, prefix):
        total = sum(q for _, _, q in elems)
        if total == 0:
            return {f"{prefix}{k}": 0.0 for k in
                    ["pauling_en", "mulliken_en", "allred_rochow_en", "gordy_en", "mb_en"]}

        sums = defaultdict(float)

        for el, ox, qty in elems:
            vals = self.compute_en(el, ox)
            for k, v in vals.items():
                sums[f"{prefix}{k}"] += v * qty

        return {k: v / total for k, v in sums.items()}

    # --------------------------
    # Public API
    # --------------------------
    def get_features(self, formula, mode="both"):
        result = self.approx.charge_balance(formula, return_format="dict")
        all_elems, var_elems = self.parse_result(result)

        out = {}

        if mode in ("all", "both"):
            out |= self.weighted(all_elems, "all_")

        if mode in ("var", "both"):
            if var_elems:
                out |= self.weighted(var_elems, "var_")
            else:
                out |= {
                    "var_pauling_en": 0.0,
                    "var_mulliken_en": 0.0,
                    "var_allred_rochow_en": 0.0,
                    "var_gordy_en": 0.0,
                    "var_mb_en": 0.0,
                }

        return out

