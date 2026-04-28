from .base import FeatureModule
from .element_properties import DEFAULT_ELEMENT_PROPERTIES, get_property_value
from .registry import register_feature
from .stats_expander import StatsExpander


@register_feature
class ElementPropertyModule(FeatureModule):
    """
    Generic composition feature module for elemental property statistics.
    """

    def __init__(self, approx, ptable, property_names=None):
        super().__init__(approx, ptable)
        self.property_names = property_names or DEFAULT_ELEMENT_PROPERTIES

    def _expand_group(self, elems, prefix):
        features = {}

        for property_name in self.property_names:
            values = []
            weights = []

            for el, ox, qty in elems:
                row = self.get_element_row(el)
                values.append(get_property_value(property_name, el, ox, row))
                weights.append(qty)

            features |= StatsExpander.expand(
                values,
                weights,
                prefix=f"{prefix}{property_name}_",
            )

        return features

    def get_features(self, formula):
        all_elems, var_elems = self.get_parsed_elements(formula)

        features = {}
        features |= self._expand_group(all_elems, "all_")

        if var_elems:
            features |= self._expand_group(var_elems, "var_")
        else:
            features |= self.zero_fill_from_all(features)

        return features
