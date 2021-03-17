import random
import numpy as np
from scipy.stats import norm
import math


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
        if res[0] > 0:
            return res[0]
        else:
            return self.sample_breakpoint_ratio()

    def sample_non_breakpoint_ratio(self):
        return abs(np.random.normal(0, self.var_0 ** 0.5, 1)[0])

    def apply_density_at_points(self, points, breakpoint):
        if not breakpoint:
            return math.sqrt(self.var_0) * norm.pdf(points)
        else:
            result = 2 * self.weights[0] * math.sqrt(self.variances[0]) * norm.pdf(points, self.means[0])
            for i in range(1, len(self.means)):
                result += self.weights[i] * math.sqrt(self.variances[i]) * norm.pdf(points, self.means[i])
            return result / self.__estimate_norming_constant()

    def __estimate_norming_constant(self):
        result = 0
        for i in range(0, 100000):
            index = random.choices(range(0, len(self.weights)), k=1, weights=self.weights)[0]
            sample = np.random.normal(self.means[index], self.variances[index] ** 0.5, 1)[0]
            if sample > 0:
                result += 1
        return result / 100000
