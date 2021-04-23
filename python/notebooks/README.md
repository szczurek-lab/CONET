# CoNET



## Running the script (synthetic data)
```bash
set @bin_dir variable to path to directory which contains CoNET executable
```

## Running the script (biological data)
```bash
set @bin_dir variable to path to directory which contains CoNET executable
Use provided R script to plot inferred tree and count matrix:
  -
```
## Usage details
### CONET user-defined parameters

| Parameter name | Description | Default value |
| ---- | -------- | --- |
| **sdfdsf** | fadfs  | "" |
| **sdfsdf** | sdfdsf | - |



 data_dir="./",
                 param_inf_iters=100000,
                 pt_inf_iters=100000,
                 counts_penalty_c=0.0,
                 event_length_penalty_c=1.0,
                 data_size_prior_c=1.0,
                 use_event_lengths_in_attachment=True,
                 seed=12312,
                 mixture_size=4,
                 num_replicas=5,
                 threads_likelihood=4,
                 parameter_resampling_frequency=10,
                 moves_between_swaps=10,
                 burn_in=10000,
                 verbose=True
