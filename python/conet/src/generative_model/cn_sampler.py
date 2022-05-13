from dataclasses import dataclass
from typing import Dict

import numpy as np

from src.generative_model.utils import sample_conditionally


@dataclass
class CNSampler:
    CN_to_weight: Dict[int, float]
    CN_to_sd: Dict[int, float]
    neutral_cn: int
    noise_prob: float

    def sample_non_neutral_CN(self):
        return sample_conditionally(
            sampler=lambda: np.random.choice(list(self.CN_to_weight.keys()), 1, p=list(self.CN_to_weight.values()))[0],
            condition=lambda cn: cn != self.neutral_cn)

    def sample_corrected_count(self, count: int):
        if np.random.uniform() <= self.noise_prob:
            cn = np.random.choice(list(self.CN_to_weight.keys()), 1, p=list(self.CN_to_weight.values()))[0]
            return self.__sample_truncated_normal(cn, self.CN_to_sd[int(cn)])

        return self.__sample_truncated_normal(count, self.CN_to_sd[count])

    @staticmethod
    def __sample_truncated_normal(mean: float, sd: float) -> float:
        return sample_conditionally(sampler=lambda: np.random.normal(mean, sd, 1), condition=lambda x: x > 0)
