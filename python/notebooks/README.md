# CoNET

## Contents
This directory contains two notebooks that present examples of CONET usage for:
### Synthetic data
  Contains notebook for synthetic data generation, inference and result scoring.
### Biological_data
  Applies CONET to SA501X3F xenograft breast cancer data.
## Usage details
### CONET user-defined parameters
 CONET depends on CONETParameters object of class 
| Parameter name | Description | Default value |
| ---- | -------- | --- |
| **data_dir** | Path to directory containing input files. Inference results will be saved there.  |"./" |
| **param_inf_iters** | Number of MCMC iterations for inference of difference distribution. | 100000 |
| **pt_inf_iters** | Number of MCMC iterations for the search of MAP event tree. | 100000 |
| **counts_penalty_c** | Constant controlling impact of penalty for latge discrepancies between inferred and real count matrices. | 0.0 |
| **event_length_penalty_c** | Constant controlling impact of penalty for long inferred events. | 1.0 |
| **data_size_prior_c** | Constant controlling impact of data size prior. | 1.0 |
| **use_event_lengths_in_attachment** | If True cell attachment probability will depend on average event length in the history, otherwise it will be uniform.| True |
| **seed** | Seed for C++ RNG | 12312 |
| **mixture_size** | Initial number of components in difference distribution for breakpoint loci. This value may be decreased in the course of inference but will never be increased.| 4 |
| **num_replicas** | Number of tempered chain replicas in MAP event tree search. | 5 |
| **threads_likelihood** | Number of threads which will be used for the most demanding likelihood calculations. | 4 |
| **parameter_resampling_frequency** | Number of tree MCMC moves for each parameter MCMC move. | 10 |
| **moves_between_swaps** | Number of MCMC moves done by each replica before swap move is attempted. | 10 |
| **burn_in** | Number of MCMC iterations which should be skipped before statistics gathering. | 10000 |
| **verbose** | True if CONET should print messages during inference. | True |
