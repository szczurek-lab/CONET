
## Running the script (biological data)
```bash
set @bin_dir variable to path to directory which contains CoNET executable
Use provided R script to plot inferred tree and count matrix:
  -
```

## Input Data

Input data should be provided in the form of corrected counts matrix. With subsequent bins in rows and cells in columns.
The matrix should contain 4 additional columns (placed at positions 1,2,3,4 in the matrix):

### Chromosome
Bin's chromosome number - should always be an integer.
### Start
Bin's start loci
### End
Bin's end loci
### Width 
Bin's width

Example of input matrix for SA501X3F xenograft breast cancer data is contained in data/SA501X3F_filtered_corrected_counts.csv
