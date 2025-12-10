from pymatgen.core import Species
from .base import FeatureModule
from .registry import register_feature
from .utils import lookup_ionic_radius

@register_feature
class TransitionMetalPotentialModule(FeatureModule):
    """
    Computes transition-metal ionic potential, cation/anion potentials,
    cationic potential, and radial potential.

    These features are **only meaningful for layered NaMO2, LiMO2, etc.**
    MODULE IS OPTIONAL.

    User must specify cation="Na" and anion="O" (or their own choice).
    """

    def __init__(self, approx, ptable, cation="Na", anion="O"):
        super().__init__(approx, ptable)
        self.cation = cation
        self.anion = anion

    # --------- TM ionic potential ---------
    def TM_weighted_ionic_pot(self, elems):
        total = 0.0
        for el, ox, qty in elems:
            if el not in (self.cation, self.anion):
                total += lookup_ionic_radius(el, ox) * qty
        return total

    # --------- Cation / Anion ionic potentials ---------
    def cat_an_ionic_pots(self, elems):
        cat_pot = an_pot = 0.0

        for el, ox, qty in elems:
            if el == self.cation:
                cat_pot += lookup_ionic_radius(el, ox) * qty
            elif el == self.anion:
                an_pot += lookup_ionic_radius(el, ox) * qty

        return cat_pot, an_pot

    # --------- TM weighted ionic radii ---------
    def TM_weighted_radius(self, elems):
        weighted = 0.0
        total_qty = 0.0

        for el, ox, qty in elems:
            if el not in (self.cation, self.anion):
                r = Species(el, ox).ionic_radius
                if r:
                    weighted += r * qty
                    total_qty += qty

        return weighted / total_qty if total_qty else 0.0

    # --------- radial potential ---------
    def radial_potential(self, elems):
        TM_r = self.TM_weighted_radius(elems)
        cat_p, an_p = self.cat_an_ionic_pots(elems)

        if an_p == 0:
            return 0.0

        return (TM_r * cat_p) / an_p

    def get_features(self, formula):
        result = self.approx.charge_balance(formula, return_format="dict")
        all_elems, _ = self.parse_result(result)

        TM_pot = self.TM_weighted_ionic_pot(all_elems)
        cat_pot, an_pot = self.cat_an_ionic_pots(all_elems)
        catonic_potential = (TM_pot * cat_pot) / an_pot if an_pot else 0.0
        radial_pot = self.radial_potential(all_elems)

        return {
            "m_weighted_ionic_potential": TM_pot,
            "cation_ionic_potential": cat_pot,
            "anion_ionic_potential": an_pot,
            "cation_anion_relative_potential": catonic_potential,
            "m_radial_potential": radial_pot,
        }
