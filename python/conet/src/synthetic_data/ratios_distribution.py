import random
import numpy as np


def _check_parameters(weights, means, variances, var_0):
    lengths = {len(weights), len(means), len(variances)}
    if len(lengths) > 1:
        raise RuntimeError("Received different lengths of vectors!")
    if len(weights) < 1:
        raise RuntimeError("Mixture must have positive size!")
    if len(list(filter(lambda x: x <= 0, means + weights + variances + [var_0]))) > 0:
        raise RuntimeError("Means, weights and variances must be positive!")


class RatiosDistribution:
    def __init__(self, weights, means, variances, var_0):
        _check_parameters(weights, means, variances, var_0)
        self.weights = weights
        self.means = means
        self.variances = variances
        self.var_0 = var_0

    def sample_breakpoint_ratio(self):
        index = random.choices(range(0, len(self.weights)), k=1, weights=self.weights)[0]
        res = np.random.normal(self.means[index], self.variances[index] ** 0.5, 1)
        if res > 0:
            return res[0]
        else:
            return self.sample_breakpoint_ratio()

    def sample_non_breakpoint_ratio(self):
        return abs(np.random.normal(0, self.var_0 ** 0.5, 1)[0])