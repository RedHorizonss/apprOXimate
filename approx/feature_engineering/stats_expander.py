import numpy as np
from collections import Counter

class StatsExpander:
    """
    Expands element-level properties into
    sum / avg / std / min / max / range / mode
    """

    @staticmethod
    def expand(values, weights=None, prefix=""):
        values = np.array(values, dtype=float)

        if weights is None:
            weights = np.ones_like(values)
        else:
            weights = np.array(weights, dtype=float)

        total_weight = weights.sum()

        avg = np.average(values, weights=weights)
        var = np.average((values - avg) ** 2, weights=weights)
        std = np.sqrt(var)

        mode_val = Counter(values).most_common(1)[0][0]

        return {
            f"{prefix}sum": np.sum(values * weights),
            f"{prefix}avg": avg,
            f"{prefix}dev": std,
            f"{prefix}min": np.min(values),
            f"{prefix}max": np.max(values),
            f"{prefix}range": np.max(values) - np.min(values),
            f"{prefix}mode": mode_val,
        }
