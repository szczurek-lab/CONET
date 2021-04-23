
## Running the script (biological data)
```bash
set @bin_dir variable to path to directory which contains CoNET executable
Use provided R script to plot inferred tree and count matrix:
  -
```
## Input Data

Input data should be passed to object from DataConverter class. 

### Corrected Counts matrix
Basic input data should be provided in the form of corrected counts matrix. With subsequent bins in rows and cells in columns.
The matrix should contain 4 additional columns (placed at positions 1,2,3,4 in the matrix):

### Chromosome
Bin's chromosome number - should always be an integer.
### Start
Bin's start loci
### End
Bin's end loci
### Width 
Bin's width
```bash
Example of input matrix for SA501X3F xenograft breast cancer data is contained in data/SA501X3F_filtered_corrected_counts.csv
```
### Breakpoint candidate loci
Apart from the corrected counts matrix CONET expects indices of breakpoint candidate loci (only pairs of loci from this set will constitute potential event set).
Each candidate loci should be given in the form of corresponding bin's index in the corrected counts matrix (indices should start at 0).

### neutral_cn (defaults to 2.0)
Neutral copy number. 

### default_bin_width
Number which will be used for artifical bin's corresponding to ends of chromosomes.

### event_length_normalizer
Real constant which will be used to normalize event lengths (in the notebook this constant is set to total length of all chromosomes). 

### add_chr_ends_to_indices (True / False)
It is recommended to add artifical chromosome ends to the set of candidate loci with neutral CN. This can either be done manually by the user by insertion of additonal bins to corrected counts matrix (add_chr_ends_to_indices=False) or automatically (add_chr_ends_to_indices=True).

Option add_chr_ends_to_indices from method DataConverter.create_CoNET_input_files should be set to the same value as add_chr_ends_to_indices.


For the user's convenience DataConverter.create_CoNET_input_files allows passing of @chromosomes list which contains numbers of chromosomes for which the inference will be conducted. Note that the plotting script always outputs inferred heatmap for the whole genome, hence for chromosomes not contained in the list every bin will have neutral copy number. 
## Plotting
Set RESULTS_DIR variable to path to directory containing inference result files.

Set READS_DIR variable to path to directory containing corrected counts matrix.

The script will output:
* inferred tree plot (tree.html)
* inferred counts heatmap (heatmap_CNs.tiff)
* real corrected counts heatmap (heatmap_CCs.tiff)
