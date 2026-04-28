import numpy as np
from collections import defaultdict

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
        var = np.average(np.square(values - avg), weights=weights)
        std = np.sqrt(var)

        weighted_counts = defaultdict(float)
        for value, weight in zip(values, weights):
            weighted_counts[float(value)] += float(weight)
        mode_val = max(weighted_counts.items(), key=lambda item: item[1])[0]

        return {
            f"{prefix}sum": np.sum(values * weights),
            f"{prefix}avg": avg,
            f"{prefix}dev": std,
            f"{prefix}min": np.min(values),
            f"{prefix}max": np.max(values),
            f"{prefix}range": np.max(values) - np.min(values),
            f"{prefix}mode": mode_val,
        }
