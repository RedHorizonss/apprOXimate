class FeatureModule:
    def __init__(self, approx, ptable):
        self.approx = approx
        self.ptable = ptable

    def parse_result(self, result_dict):
        all_elems = []
        var_elems = []

        for el, data in result_dict["elements"].items():
            if "states" in data:
                for st in data["states"]:
                    ox, qty, fixed = st["oxidation_state"], st["quantity"], st["is_fixed"]
                    all_elems.append((el, ox, qty))
                    if not fixed:
                        var_elems.append((el, ox, qty))
            else:
                ox, qty, fixed = data["oxidation_state"], data["quantity"], data["is_fixed"]
                all_elems.append((el, ox, qty))
                if not fixed:
                    var_elems.append((el, ox, qty))

        return all_elems, var_elems

    def filter_mode(self, features, mode):
        if mode == "both":
            return features
        if mode == "all":
            return {k: v for k, v in features.items() if k.startswith("all_")}
        if mode == "var":
            return {k: v for k, v in features.items() if k.startswith("var_")}
        raise ValueError("mode must be: all, var, both")

    def get_features(self, formula):
        raise NotImplementedError
