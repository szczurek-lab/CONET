import numpy as np
import random

class CNSampler:
    def __init__(self, CN_to_weight, CN_to_sd):
        self.CN_to_weight = CN_to_weight
        sum_ = sum(list(CN_to_weight.values()))
        for key in self.CN_to_weight.keys():
            self.CN_to_weight[key] = self.CN_to_weight[key] / sum_
        self.CN_to_sd = CN_to_sd
        self.NEUTRAL_CN = 2 
        self.NOISE_PROB = 0.1
        
    def sample_CN(self):
        cn = np.random.choice(list(self.CN_to_weight.keys()), 1, p = list(self.CN_to_weight.values()))[0]
        while cn == self.NEUTRAL_CN:
            cn = np.random.choice(list(self.CN_to_weight.keys()), 1, p = list(self.CN_to_weight.values()))[0]
        return cn
    
    def sample_corrected_count(self, count): 
        if np.random.uniform() <= self.NOISE_PROB:
            cn = np.random.choice(list(self.CN_to_weight.keys()), 1, p = list(self.CN_to_weight.values()))[0]
            cc = cn + np.random.normal(0, self.CN_to_sd[int(count)], 1)
            while cc < 0:
                cc = cn + np.random.normal(0, self.CN_to_sd[int(count)], 1)
            return cc
            
        cc = count + np.random.normal(0, self.CN_to_sd[int(count)], 1)
        while cc < 0:
            cc = count + np.random.normal(0, self.CN_to_sd[int(count)], 1)
        return cc
