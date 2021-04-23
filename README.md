# CONET

## About
 CONET is a copy number event tree inference and copy number calling bayesian procedure from whole-genome DNA sequencing data.
 
 CONET should be used with the aid of provided Python scripts. 
 Examples and details are provided in two notebooks:
 * python/notebooks/biological_data/biological_data.ipynb
 * python/notebooks/synthetic_data/synthetic_data.ipynb
    
 
## Requirements
* C++ compiler that supports C++14 standard
* Python 3.0 or higher
* GNU Make
* R 4.0 or higher (for output plotting)
## Installation

### C++
```bash
git clone https://github.com/tc360950/CONET.git # Clone the repository
cd CONET/src
make                                          # Build the executables
```

### Python wrapper
```bash
cd python/conet
pip install .
```
