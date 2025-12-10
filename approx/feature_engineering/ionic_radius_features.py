# approx/feature_engineering/ionic_radius.py
from pymatgen.core.periodic_table import Element
from .base import FeatureModule
from .utils import weighted_average, lookup_ionic_radius
from .registry import register_feature


@register_feature
class IonicRadiusModule(FeatureModule):
    def __init__(self, approx, ptable):
        super().__init__(approx, ptable)


    def get_features(self, formula):
        result = self.approx.charge_balance(formula, return_format="dict")
        all_elems, var_elems = self.parse_result(result)

        all_rad = weighted_average(
            all_elems,
            lambda el, ox: {"r": lookup_ionic_radius(el, ox)}
        )
        all_out = {"all_ionic_radius": all_rad.get("r", 0.0)}

        if var_elems:
            var_rad = weighted_average(
                var_elems,
                lambda el, ox: {"r": lookup_ionic_radius(el, ox)}
            )
            var_out = {"var_ionic_radius": var_rad.get("r", 0.0)}
        else:
            var_out = {"var_ionic_radius": 0.0}

        return all_out | var_out
