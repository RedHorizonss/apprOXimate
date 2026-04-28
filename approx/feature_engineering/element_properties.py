from dataclasses import dataclass
from math import isnan


METALLOID_SYMBOLS = {"B", "Si", "Ge", "As", "Sb", "Te", "Po"}
NONMETAL_SYMBOLS = {"H","C", "N", "O", "P", "S", "Se", "F", "Cl", "Br", "I","At","He","Ne","Ar","Kr","Xe","Rn","Og"}

FAMILY_CODES = {
    "alkali_metal": 1.0,
    "alkaline_earth_metal": 2.0,
    "transition_metal": 3.0,
    "post_transition_metal": 4.0,
    "metalloid": 5.0,
    "nonmetal": 6.0,
    "halogen": 7.0,
    "noble_gas": 8.0,
    "lanthanoid": 9.0,
    "actinoid": 10.0,
    "unknown": 0.0,
}


@dataclass(frozen=True)
class PropertySpec:
    column_candidates: tuple[str, ...] = ()
    default: float = 0.0
    helper: str | None = None


PROPERTY_SPECS = {
    "Atomic_Number": PropertySpec(column_candidates=("atomic_number",)),
    "Atomic_Weight": PropertySpec(
        column_candidates=("atomic_weight", "atomic_mass", "mass")
    ),
    "Period": PropertySpec(column_candidates=("period",)),
    "group": PropertySpec(column_candidates=("group_id", "group")),
    "families": PropertySpec(helper="families"),
    "Metal": PropertySpec(helper="metal"),
    "Nonmetal": PropertySpec(helper="nonmetal"),
    "Metalloid": PropertySpec(helper="metalloid"),
    "Mendeleev_Number": PropertySpec(
        column_candidates=("mendeleev_number", "mendeleev_no")
    ),
    "Atomic_Radius": PropertySpec(
        column_candidates=("atomic_radius", "atomic_radius_rahm", "atomic_radius_cordero")
    ),
    "Covalent_Radius": PropertySpec(
        column_candidates=(
            "covalent_radius",
            "covalent_radius_pyykko",
            "covalent_radius_cordero",
        )
    ),
    "Density_(g/mL)": PropertySpec(column_candidates=("density",)),
    "specific_heat_(J/g_K)_": PropertySpec(
        column_candidates=("specific_heat", "specific_heat_capacity")
    ),
    "heat_of_fusion_(kJ/mol)_": PropertySpec(
        column_candidates=("heat_of_fusion", "fusion_heat")
    ),
    "heat_of_vaporization_(kJ/mol)_": PropertySpec(
        column_candidates=("heat_of_vaporization", "evaporation_heat")
    ),
    "thermal_conductivity_(W/(m_K))_": PropertySpec(
        column_candidates=("thermal_conductivity",)
    ),
}

DEFAULT_ELEMENT_PROPERTIES = list(PROPERTY_SPECS)


def _coerce_float(value, default):
    if value is None:
        return float(default)

    try:
        if isnan(value):
            return float(default)
    except TypeError:
        pass

    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def _first_present(row, column_candidates):
    for column in column_candidates:
        if column in row.index:
            value = row[column]
            if value is not None:
                try:
                    if isnan(value):
                        continue
                except TypeError:
                    return value
                return value
    return None


def _group_id(row):
    raw_group = _first_present(row, ("group_id", "group"))
    return int(_coerce_float(raw_group, 0.0))


def _series_text(row):
    for column in ("series", "series_name", "classification", "category", "block"):
        if column in row.index and row[column] is not None:
            return str(row[column]).strip().lower()
    return ""


def _is_metalloid(symbol, row):
    series = _series_text(row)
    return symbol in METALLOID_SYMBOLS or "metalloid" in series


def _is_nonmetal(symbol, row):
    series = _series_text(row)
    return symbol in NONMETAL_SYMBOLS or "nonmetal" in series or "noble gas" in series


def _is_metal(symbol, row):
    if _is_metalloid(symbol, row) or _is_nonmetal(symbol, row):
        return False

    series = _series_text(row)
    if "metal" in series:
        return True

    group_id = _group_id(row)
    atomic_number = int(_coerce_float(_first_present(row, ("atomic_number",)), 0.0))

    if group_id in {1, 2} and symbol != "H":
        return True
    if 3 <= group_id <= 12:
        return True
    if group_id in {13, 14, 15, 16} and atomic_number > 0 and symbol not in NONMETAL_SYMBOLS:
        return True
    if 57 <= atomic_number <= 80 or 89 <= atomic_number <= 112:
        return True

    return False


def _family_code(symbol, row):
    if _is_metalloid(symbol, row):
        return FAMILY_CODES["metalloid"]
    if _is_nonmetal(symbol, row):
        group_id = _group_id(row)
        if group_id == 17:
            return FAMILY_CODES["halogen"]
        if group_id == 18:
            return FAMILY_CODES["noble_gas"]
        return FAMILY_CODES["nonmetal"]

    atomic_number = int(_coerce_float(_first_present(row, ("atomic_number",)), 0.0))
    group_id = _group_id(row)

    if group_id == 1 and symbol != "H":
        return FAMILY_CODES["alkali_metal"]
    if group_id == 2:
        return FAMILY_CODES["alkaline_earth_metal"]
    if 3 <= group_id <= 12:
        return FAMILY_CODES["transition_metal"]
    if 57 <= atomic_number <= 71:
        return FAMILY_CODES["lanthanoid"]
    if 89 <= atomic_number <= 103:
        return FAMILY_CODES["actinoid"]
    if _is_metal(symbol, row):
        return FAMILY_CODES["post_transition_metal"]

    return FAMILY_CODES["unknown"]


HELPERS = {
    "families": _family_code,
    "metal": lambda symbol, row: 1.0 if _is_metal(symbol, row) else 0.0,
    "nonmetal": lambda symbol, row: 1.0 if _is_nonmetal(symbol, row) else 0.0,
    "metalloid": lambda symbol, row: 1.0 if _is_metalloid(symbol, row) else 0.0,
}


def get_property_value(property_name, symbol, ox, row):
    del ox

    if property_name not in PROPERTY_SPECS:
        raise KeyError(f"Unknown property: {property_name}")

    spec = PROPERTY_SPECS[property_name]

    if spec.helper is not None:
        helper = HELPERS[spec.helper]
        return _coerce_float(helper(symbol, row), spec.default)

    value = _first_present(row, spec.column_candidates)
    return _coerce_float(value, spec.default)
