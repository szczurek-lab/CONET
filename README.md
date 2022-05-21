# CONET

## About
### CONET: Copy number event tree model of evolutionary tumor history for single-cell data
#### Authors:
Magda Markowska, Tomasz Cąkała, Błażej Miasojedow, Bogac Aybey, Dilafruz Juraeva, Johanna Mazur, Edith Ross, Eike Staub and Ewa Szczurek
#### Questions about C++ and python implementation
Tomasz Cąkała, tc360950@gmail.com
#### Questions about R implementation and general questions
Magda Markowska, magda.markowska@gmail.com
#### Manuscript's corresponding author:
Ewa Szczurek, szczurek@mimuw.edu.pl
 
 
## Requirements
### Necessary
* C++ compiler that supports C++14 standard
* Python 3.6 or higher
* GNU Make
### Additional
* R 4.0 or higher (for output plotting)

# Usage 
We advise new users to examine notebooks found in python/notebooks directory.

# Contents 

This project consists of 3 main components:
* CONET cpp sources which define CONET executable 
* CONET executable should not be used directly but with aid of conet-py python package (*python/conet-py*)
* R script for advanced plots of inference results 

# Installation

## In container 
Use image defined in *CONET.Dockerfile*. It installs conet-py and compiles cpp CONET into executable 
*~/conet-py/CONET*. If you want to install CONET locally it's easy to mimic steps executed in the image. 


# Input Data

### Corrected counts matrix 
Basic input data should be provided in the form of corrected counts matrix. With subsequent bins in rows and cells in columns.
The matrix should contain 5 additional columns (placed at positions 1,2,3,4,5 in the matrix):

#### chr
Bin's chromosome number - should always be an integer (please change X to 23 and Y to 24).
#### start
Bin's start locus
#### end
Bin's end locus
#### width 
Bin's width
#### candidate_brkp
binary breakpoint indicator:
1 -- if the start locus of the bin is a candidate breakpoint
0 -- otherwise

```bash
Example of input matrix for SA501X3F xenograft breast cancer data is contained in CONET/python/notebooks/biological_data/data/SA501X3F_filtered_corrected_counts.csv
and for TN2 breast cancer data -- in CONET/R/TN2_corrected_counts_with_indices_50cells.csv
```

## CONET should be used with the aid of provided Python scripts. 

Examples and details are provided in three notebooks:
#### Biological_data
Applies CONET to SA501X3F xenograft breast cancer data (DLP sequencing)
* python/notebooks/biological_data/biological_data.ipynb
#### Synthetic data
Contains notebook for synthetic data generation, inference and result scoring.
* python/notebooks/per_bin_generative_model/generative_model.ipynb
* python/notebooks/per_breakpoint_generative_model/generative_model.ipynb


## Usage details
### CONET user-defined parameters
CONET depends on a number of user-defined parameters which are represented by objects of class CONETParameters. 

| Parameter name                      | Description                                                                                                                                                      | Default value |
|-------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|
| **data_dir**                        | Path to directory containing input file.                                                                                                                         | "./"          |
| **output_dir**                      | Path to output directory. Inference results will be saved there.                                                                                                 | "./output"    |
| **param_inf_iters**                 | Number of MCMC iterations for joint tree and model parameters inference.                                                                                         | 100000        |
| **pt_inf_iters**                    | Number of MCMC iterations for tree inference.                                                                                                                    | 100000        |
| **counts_penalty_s1**               | Constant controlling impact of penalty for large discrepancies between inferred and real count matrices.                                                         | 0.0           |
| **counts_penalty_s2**               | Constant controlling impact of penalty for inferring clusters with changed copy number equal to basal ploidy.                                                    | 0.0           |
| **event_length_penalty_k0**         | Constant controlling impact of penalty for long inferred events.                                                                                                 | 1.0           |
| **tree_structure_prior_k1**         | Constant controlling impact of data size part of tree structure prior.                                                                                           | 1.0           |
| **use_event_lengths_in_attachment** | If True cell attachment probability will depend on average event length in the history, otherwise it will be uniform.                                            | True          |
| **seed**                            | Seed for C++ RNG                                                                                                                                                 | 12312         |
| **mixture_size**                    | Initial number of components in difference distribution for breakpoint loci. This value may be decreased in the course of inference but will never be increased. | 4             |
| **num_replicas**                    | Number of tempered chain replicas in MAP event tree search.                                                                                                      | 5             |
| **threads_likelihood**              | Number of threads which will be used for the most demanding likelihood calculations.                                                                             | 4             |
| **parameter_resampling_frequency**  | Number of tree MCMC moves for each parameter MCMC move.                                                                                                          | 10            |
| **moves_between_swaps**             | Number of MCMC moves done by each replica before swap move is attempted.                                                                                         | 10            |
| **burn_in**                         | Number of MCMC iterations which should be skipped before statistics gathering.                                                                                   | 10000         |
| **verbose**                         | True if CONET should print messages during inference.                                                                                                            | True          |

### Guide to parameter settings

For more details please refer to Additional File 1: S7 A recommended procedure for setting CONET regularization parameters.

| Parameter name              | Recommendation                                                                                                                                  | Initial value |
|-----------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------|---------------|
| **param_inf_iters**         | Start with initial value, save and plot likelihood to check convergence. Depends on input size - number of cells and candidate breakpoint loci. | 250000        |
| **pt_inf_iters**            | Start with initial value, save and plot likelihood to check convergence. Depends on input size - number of cells and candidate breakpoint loci. | 500000        |
| **event_length_penalty_k0** | Start with initial value and increase if you want to penalize trees inferring long events.                                                      | 1.0           |
| **tree_structure_prior_k1** | Start with initial value and try increasing/decreasing if quality measures are not satisfactory.                                                | 0.0           |
| **counts_penalty_s1**       | Start with initial value and try increasing/decreasing if quality measures are not satisfactory.                                                | 100000.0      |
| **counts_penalty_s2**       | Start with initial value and try increasing/decreasing if quality measures are not satisfactory.                                                | 100000.0      |
| **seed**                    | Try using different seed to make sure you do not stuck in local optima.                                                                         | 12312         |

We recommend that all other parametr values are left at default.


## Output of CONET
### inferred_tree
Structure of inferred CONET in a format readable by readTree.R script. Can be saved as Newick using this script.
#### inferred_attachment
Cells attachment to CONET nodes. Readable by readTree.R script.
#### inferred_breakpoints
Binary matrix with inferred breakpoints per genomic locus and cell (coressponding to input file).
#### inferred_distribution
Model parameters inffered by CONET. 
On first line mean and variance of no_breakpoint distribution (R+ truncated normal).
On the following lines: weigth; mean; variance of the components of breakpoint distribution (R+ truncated mixed normal). One line per each component.
#### edge_confidence
Inferred CONET edges with number of iterations (after burn in phase) they appeared in.

## Final CN calling and output visualization 
Final CN matrix can be inferred with provided R script. It also allows the user to visualize results of CONET model.
### Input
CONET input, output plus all dependencies provided at CONET/R directory.
### Usage
All details in CONET/R/readTree.R
### Biological example
CONET/R/TN2example Illustration of the results from CONET applied to 100 cells from TN2 breast cancer data (ACT sequencing).
### Output
#### tree_with_genes.pdf
Plot of CONET with chr, start and end breakpoint loci, cancer genes and number of attached cells. 
#### inferred_tree_Newick
CONET in the Newick format
#### inferred_CN_matrix.csv
Final inferred CN matrix in the same format as input corrected counts matrix
#### CONET_quality_measures.csv
Quality measures calculated for the inferred CONET and CN matrix (described in manuscript, Additonal File 1: Section S6.1 )
#### heatmap_CCs.tiff
Heatmap of inferred CNs (genomic loci on X axis, cells on Y axis)
#### heatmap_CNs.tiff
Heatmap of corrected counts (genomic loci on X axis, cells on Y axis), plotted in the same cells order as CN heatmap.
