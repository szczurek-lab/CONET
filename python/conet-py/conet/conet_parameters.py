import dataclasses
from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class CONETParameters:
    data_dir: str = "./"
    param_inf_iters: int = 100000
    pt_inf_iters: int = 100000
    counts_penalty_s1: float = 0.0
    counts_penalty_s2: float = 0.0
    event_length_penalty_k0: float = 1.0
    tree_structure_prior_k1: float = 1.0
    use_event_lengths_in_attachment: bool = True
    seed: int = 12312
    mixture_size: int = 4
    num_replicas: int = 5
    threads_likelihood: int = 4
    verbose: bool = True
    neutral_cn: float = 2.0
    output_dir: str = "./"

    def to_arg_value_pairs(self) -> List[Tuple[str, str]]:
        return [(f"--{key}", f"{value}") for key, value in dataclasses.asdict(self).items()]
