from .base import FeatureModule
from .utils import weighted_average
from .registry import register_feature

#Might be worth making the config, period and ion variables into a utility function?
from mendeleev.econf import ElectronicConfiguration

@register_feature
class ValenceFeatureModule(FeatureModule):

    def get_valence(self, element, ox):
        config = ElectronicConfiguration(
            self.ptable.query("symbol == @element")["electronic_configuration"].values[0]
        )
        period = self.ptable.query("symbol == @element")["period"].values[0]
        # Remove the number of electrons for the given element
        ion = config.ionize(ox)

        s = ion.conf.get((period, "s"), 0)
        p = ion.conf.get((period, "p"), 0)
        d = ion.conf.get((period - 1, "d"), 0)
        f = ion.conf.get((period - 1, "f"), 0)

        return {"s": s, "p": p, "d": d, "f": f, "total": s + p + d + f}

    def get_unfilled(self, val):
        maxc = {"s": 2, "p": 6, "d": 10, "f": 14}
        # Return a dictionary with the unfilled orbitals (valence - max configuration)
        return {
            orb: maxc[orb] - val[orb] for orb in maxc
        } | {
            "total": sum(maxc[o] - val[o] for o in maxc)
        }

    def get_features(self, formula):
        result = self.approx.charge_balance(formula, return_format="dict")
        all_elems, var_elems = self.parse_result(result)

        # ALL
        all_val = weighted_average(all_elems, lambda e, ox: self.get_valence(e, ox))
        all_unf = weighted_average(all_elems, lambda e, ox: self.get_unfilled(self.get_valence(e, ox)))

        all_val = {f"all_valence_{k}": v for k, v in all_val.items()}
        all_unf = {f"all_unfilled_valence_{k}": v for k, v in all_unf.items()}

        # VAR
        if len(var_elems) > 0:
            var_val = weighted_average(var_elems, lambda e, ox: self.get_valence(e, ox))
            var_unf = weighted_average(var_elems, lambda e, ox: self.get_unfilled(self.get_valence(e, ox)))

            var_val = {f"var_valence_{k}": v for k, v in var_val.items()}
            var_unf = {f"var_unfilled_valence_{k}": v for k, v in var_unf.items()}
        else:
            # No variable oxidation states → return zeros
            var_val = {f"var_valence_{k}": 0 for k in all_val}
            var_unf = {f"var_unfilled_valence_{k}": 0 for k in all_unf}
            
        # Return a dictionary of all the variable calculated
        return all_val | all_unf | var_val | var_unf
