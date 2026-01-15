# approx/feature_engineering/ionic_radius.py

from .base import FeatureModule
from .registry import register_feature
from .stats_expander import StatsExpander
from .utils import lookup_ionic_radius


@register_feature
class IonicRadiusModule(FeatureModule):

    def _expand_group(self, elems, prefix):
        values = []
        weights = []

        for el, ox, qty in elems:
            r = lookup_ionic_radius(el, ox, default=100.0)
            values.append(r)
            weights.append(qty)

        return StatsExpander.expand(
            values,
            weights,
            prefix=f"{prefix}ionic_radius_"
        )

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
            # zero-fill
            zero = self._expand_group(all_elems, "all_")
            for k in zero:
                features[k.replace("all_", "var_")] = 0.0

        return features
