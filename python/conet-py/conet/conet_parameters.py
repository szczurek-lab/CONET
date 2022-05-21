from dataclasses import dataclass


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
    parameter_resampling_frequency: int = 10
    moves_between_swaps: int = 10
    burn_in: int = 10000
    verbose: bool = True
    neutral_cn: float = 2.0
    output_dir: str = "./"

    def to_string(self):
        return [
            str(self.data_dir),
            str(self.param_inf_iters),
            str(self.pt_inf_iters),
            str(self.counts_penalty_s1),
            str(self.counts_penalty_s2),
            str(self.event_length_penalty_k0),
            str(self.tree_structure_prior_k1),
            str((lambda x: 1 if x else 0)(self.use_event_lengths_in_attachment)),
            str(self.seed),
            str(self.mixture_size),
            str(self.num_replicas),
            str(self.threads_likelihood),
            str(self.parameter_resampling_frequency),
            str(self.moves_between_swaps),
            str(self.burn_in),
            str((lambda x: 1 if x else 0)(self.verbose)),
            str(self.neutral_cn),
            str(self.output_dir)]
