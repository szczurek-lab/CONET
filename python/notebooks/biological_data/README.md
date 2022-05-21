# About 
This notebook applies CONET to SA501X3F xenograft breast cancer data. It is intended as a guideline on how to 
use CONET for other datasets. 
# Running the notebook 
## In container 

Run image defined in *CONET.Jupyter.Dockerfile* (in the root directory)
and expose its port 8889 to 8889 on the host.

Access jupyter notebook through url displayed on the stdout of the container. 
The notebook can be run there without any additional steps.

## Locally
1. Compile CONET 
2. Install conet-py 
3. set @bin_dir variable in the notebook to path to directory which contains CoNET executable
4. Run the notebook
## Input Data

The only required input is corrected counts matrix which should be read into memory in the form of 
pandas.DataFrame.

### Corrected Counts matrix
Basic input data should be provided in the form of corrected counts matrix. With subsequent bins in rows and cells in columns.
The matrix should contain 5 additional columns (placed at positions 1,2,3,4,5 in the matrix):

### Chromosome
Bin's chromosome number - should always be an integer.
### Start
Bin's start locus
### End
Bin's end locus
### Width 
Bin's width
### Breakpoint candidate
Should be equal to 1 if bin should be considered as a breakpoint candidate, to 0 otherwise.
(Potential events on the CONET tree are formed from pairs of bins with this column set to 1)

**Example of input matrix for SA501X3F xenograft breast cancer data is contained in data/SA501X3F_filtered_corrected_counts.csv**

## Plotting
Set RESULTS_DIR variable to path to directory containing inference result files.

Set READS_DIR variable to path to directory containing corrected counts matrix.

The script will output:
* inferred tree plot (tree.html)
* inferred counts heatmap (heatmap_CNs.tiff)
* real corrected counts heatmap (heatmap_CCs.tiff)
