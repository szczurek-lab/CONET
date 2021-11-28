#Set this variable to directory containing CONET executable
bin_dir = './'

import sys
sys.path.append('../..')
import pandas as pd
import conet 
import conet.src.data_converter.data_converter as dc
import conet.src.conet as c
import conet.src.conet_parameters as cp
import conet.src.inference_result as ir


# Use DataConverter class to convert corrected counts matrix into CONET specific input files
data_converter = dc.DataConverter("data/SA501X3F_filtered_corrected_counts_2.csv", 
                                  delimiter= ',', 
                                  default_bin_length = 150000, 
                                  event_length_normalizer = 3095677412,
                                  add_chromosome_ends = True,
                                  neutral_cn = 2.0)
                                  
                                  
                                  
# DataConverter expects list of potential breakpoint candidates, here we load precalculated set of candidates
breakpoint_candidates_indices = pd.read_csv('data/indices.csv', header=None, sep = ' ')[1].tolist()
breakpoint_candidates_indices = list(map(lambda x : x - 1, breakpoint_candidates_indices))




# Converts corrected counts matrix to CONET input files. @chromosomes parameter can be set to restrict inference to 
# a subset of chromosomes
data_converter.create_CoNET_input_files( bin_dir, chromosomes=[17,18, 20, 23], add_chr_ends_to_indices=True)




# this may take up to 10 minutes
conet = c.CONET(bin_dir + "CONET")
params = cp.CONETParameters(data_size_prior_c = 0.5, data_dir = bin_dir, counts_penalty_c=200000, 
                            param_inf_iters=30000, seed = 21567, mixture_size=2, pt_inf_iters=200000, neutral_cn =2.0, output_dir = "./output/")
conet.infer_tree(params)




