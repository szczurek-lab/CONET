#Set this variable to directory containing CONET executable
bin_dir = './'

import sys
sys.path.append('../..')
import pandas as pd
import conet
import conet.src.data_converter.data_converter as dc
import conet.src.conet as c
import conet.src.conet_parameters as cp
from  conet.src.data_converter.data_converter_2 import CorrectedCounts, DataConverter2

# Use DataConverter class to convert corrected counts matrix into CONET specific input files
#data_converter = dc.DataConverter("data/SA501X3F_filtered_corrected_counts.csv",
 #                                 delimiter= ',',
  #                                default_bin_length = 150000,
   #                               event_length_normalizer = 3095677412,
    #                              add_chromosome_ends = True,
     #                             neutral_cn = 2.0)

#data_converter.create_CoNET_input_files(bin_dir, chromosomes=[17,18, 20, 23], add_chr_ends_to_indices=True)

cc = CorrectedCounts("data/SA501X3F_filtered_corrected_counts.csv",
                                  delimiter= ',')
cc.add_chromosome_ends(2, 150000)
data_converter2 = DataConverter2(event_length_normalizer = 3095677412)
data_converter2.create_CoNET_input_files(bin_dir, chromosomes=[17,18, 20, 23], corrected_counts=cc)




