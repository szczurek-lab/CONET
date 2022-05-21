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


## Outputs
We end the computation by dumping results with command:
``` 
InferenceResult.dump_results_to_dir
```
Results consist of three files:
* inferred_tree 
* inferred_attachment
* inferred_counts 

*inferred_counts* is a n x m matrix where n is equal to number of cells and m is equal to number of loci in input corrected counts matrix. 
It stores inferred Copy Number profiles. 

*inferred_tree* contains lines in the format:

```
   (<chr>_<bin_start>,<chr>_<bin_start2>)-(<chr2>_<bin_start3>,<chr2>_<bin_start4>)
```
and should be interpreted as:
Node with event spanning segment *[bin_start, bin_start2)* in chromosome *chr* is a parent of node with event spanning segment 
*[bin_start3, bin_start4)* in chromosome *chr2*.

*inferred_attachment* contains lines in the format:
    
```
            <cell_name>;<chr>_<bin_start>;<chr>_<bin_start2>
```
And should be interpreted as:
Cell with name <cell_name> (cell names are taken from corrected counts matrix headers) is attached to node with event 
spanning *[bin_start, bin_start2)* in chromosome *chr*.

## Plotting
Set RESULTS_DIR variable to path to directory containing inference result files.

Set READS_DIR variable to path to directory containing corrected counts matrix.

The script will output:
* inferred tree plot (tree.html)
* inferred counts heatmap (heatmap_CNs.tiff)
* real corrected counts heatmap (heatmap_CCs.tiff)

