

class CONETParameters:
    #
    # @data_dir - path to output directory
    # @param_inf_iters - number of iterations for parameters inference
    # @pt_inf_iters - number of iterations for tree search with parallel tempering

    #
    # @mixture_size - maximal size of mixture which defines difference distribution
    # @num_replicas - number of replicas in parallel tempering
    # @threads_likelihood - number of threads which shall be used in likelihood calculation during parameter inference
    #
    def __init__(self,
                 data_dir="./",
                 param_inf_iters=100000,
                 pt_inf_iters=100000,
                 counts_penalty_s1=0.0,
                 counts_penalty_s2=0.0,
                 event_length_penalty_k0=1.0,
                 tree_structure_prior_k1=1.0,
                 use_event_lengths_in_attachment=True,
                 seed=12312,
                 mixture_size=4,
                 num_replicas=5,
                 threads_likelihood=4,
                 parameter_resampling_frequency=10,
                 moves_between_swaps=10,
                 burn_in=10000,
                 verbose=True,
                 neutral_cn=2.0,
                 output_dir="./"
                 ):
        self.data_dir = data_dir
        self.param_inf_iters = param_inf_iters
        self.pt_inf_iters = pt_inf_iters
        self.counts_penalty_s1 = counts_penalty_s1
        self.counts_penalty_s2 = counts_penalty_s2
        self.event_length_penalty_k0 = event_length_penalty_k0
        self.tree_structure_prior_k1 = tree_structure_prior_k1
        self.use_event_lengths_in_attachment = (lambda x: 1 if x else 0)(use_event_lengths_in_attachment)
        self.seed = seed
        self.mixture_size = mixture_size
        self.num_replicas = num_replicas
        self.threads_likelihood = threads_likelihood
        self.parameter_resampling_frequency = parameter_resampling_frequency
        self.moves_between_swaps = moves_between_swaps
        self.burn_in = burn_in
        self.verbose = (lambda x: 1 if x else 0)(verbose)
        self.neutral_cn = neutral_cn
        self.output_dir = output_dir

    def to_string(self):
        return [
            str(self.data_dir),
            str(self.param_inf_iters),
            str(self.pt_inf_iters),
            str(self.counts_penalty_s1),
            str(self.counts_penalty_s2),
            str(self.event_length_penalty_k0),
            str(self.tree_structure_prior_k1),
            str(self.use_event_lengths_in_attachment),
            str(self.seed),
            str(self.mixture_size),
            str(self.num_replicas),
            str(self.threads_likelihood),
            str(self.parameter_resampling_frequency),
            str(self.moves_between_swaps),
            str(self.burn_in),
            str(self.verbose),
            str(self.neutral_cn),
            str(self.output_dir)]
