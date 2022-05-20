# CONET

## Contents
This directory contains two notebooks that present examples of CONET usage for:
### Synthetic data
  Notebook for synthetic data generation and inference. 
### Biological_data
  Applies CONET to SA501X3F xenograft breast cancer data.
## Usage
For users' convenience we have supplied *CONET.Jupyter.Dockerfile* (may be found in the repository root directory) which 
defines image with running jupyter notebook (port 8889) with conet-py and CONET installed. 

Run the image and expose its port 8889 to 8889 on the host. 
Access jupyter notebook through URL given in the stdout of the container.
Notebooks can be run there without any additional steps.
