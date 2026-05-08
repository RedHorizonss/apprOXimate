from .registry import FEATURE_REGISTRY

class MaterialFeatureExtractor:
    """
    Main feature extraction engine.
    Loads modules from registry, runs each module, and merges results.
    """

    def __init__(
        self,
        approx,
        ptable,
        mode="both",
        modules="all",
        enable_tm_potential=False,
        tm_cation="Na",
        tm_anion="O",
        ionic_radius_unit="pm",
        element_property_names=None,
    ):
        """
        Args:
            approx: ApprOXimate instance
            ptable: periodic table dataframe
            mode: "all", "var", "both"
            modules: list of module names OR "all"
            enable_tm_potential: enable TM-specific feature module?
            tm_cation: TM module cation symbol
            tm_anion: TM module anion symbol
            ionic_radius_unit: "pm" or "angstrom" for ionic-radius-based features
            element_property_names: optional list of elemental property names. Use
                None for defaults, or an empty list to skip ElementPropertyModule.
        """

        self.approx = approx
        self.ptable = ptable
        self.mode = mode
        self.ionic_radius_unit = ionic_radius_unit
        self.element_property_names = element_property_names
        self.modules = []

        # Load modules dynamically
        for name, ModuleClass in FEATURE_REGISTRY.items():

            # Skip TM module unless enabled
            if name == "TransitionMetalPotentialModule" and not enable_tm_potential:
                continue

            if name == "ElementPropertyModule" and element_property_names == []:
                continue

            # If user specified a custom list of modules
            if modules != "all" and name not in modules:
                continue

            # TM module needs extra arguments
            if name == "TransitionMetalPotentialModule":
                module = ModuleClass(
                    approx,
                    ptable,
                    cation=tm_cation,
                    anion=tm_anion,
                    ionic_radius_unit=ionic_radius_unit,
                )
            elif name == "ElementPropertyModule":
                module = ModuleClass(
                    approx,
                    ptable,
                    property_names=element_property_names,
                )
            elif name in {"IonicRadiusModule", "ElectronegativityModule"}:
                module = ModuleClass(
                    approx,
                    ptable,
                    ionic_radius_unit=ionic_radius_unit,
                )
            else:
                module = ModuleClass(approx, ptable)

            self.modules.append(module)

    def get_features(self, formula):
        """
        Run all feature modules on a single formula and merge results.
        Applies "all", "var", or "both" filtering to feature names.
        """
        features = {"formula": formula}

        for module in self.modules:
            mod_feats = module.get_features(formula)

            # Filter feature names based on mode
            mod_feats = module.filter_mode(mod_feats, self.mode)

            features.update(mod_feats)

        return features

    def featurize_many(self, formulas):
        """
        Extract features for multiple formulas → list of dicts
        """
        rows = []
        for f in formulas:
            try:
                rows.append(self.get_features(f))
            except Exception as e:
                print(f"Error in {f}: {e}")
                rows.append({"formula": f})
        return rows
