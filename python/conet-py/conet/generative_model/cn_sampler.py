from dataclasses import dataclass
from typing import Dict

import numpy as np

from conet.generative_model.utils import sample_conditionally


@dataclass
class CNSampler:
    """
    This class provides copy number sampling functionalities for per-bin generative model.
    """
    CN_to_weight: Dict[int, float]  # CNs will be sampled using weights provided by the user
    CN_to_sd: Dict[int, float]  # Each CN is assigned a sd - CC for this CN will be sampled by adding gaussian noise
    # with this deviation
    neutral_cn: int
    noise_prob: float

    def sample_non_neutral_CN(self) -> int:
        return sample_conditionally(
            sampler=lambda: np.random.choice(list(self.CN_to_weight.keys()), 1, p=list(self.CN_to_weight.values()))[0],
            condition=lambda cn: cn != self.neutral_cn)

    def sample_corrected_count(self, count: int) -> float:
        if np.random.uniform() <= self.noise_prob:
            cn = np.random.choice(list(self.CN_to_weight.keys()), 1, p=list(self.CN_to_weight.values()))[0]
            return self.__sample_truncated_normal(cn, self.CN_to_sd[int(cn)])

        return self.__sample_truncated_normal(count, self.CN_to_sd[count])

    @classmethod
    def create_default_sampler(cls) -> 'CNSampler':
        """
        Create sampler corresponding to the model used in the article
        """
        return cls({0: 0.020240121, 1: 0.203724532, 3: 0.050340118, 4: 0.038828672, 2: 0.686866557},
                   {0: 0.449, 1: 0.116, 2: 0.187, 3: 0.114, 4: 0.279, 5: 0.0957, 6: 0.4833, 7: 0.2760, 8: 6.15780,
                    9: 4.72105270}, 2, 0.1)

    @staticmethod
    def __sample_truncated_normal(mean: float, sd: float) -> float:
        return sample_conditionally(sampler=lambda: np.random.normal(mean, sd, 1), condition=lambda x: x > 0)
