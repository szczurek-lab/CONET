# Set this variable to directory containing CONET executable
bin_dir = './'

from conet import CorrectedCounts, DataConverter, CONET, CONETParameters

# Use DataConverter class to convert corrected counts matrix into CONET specific input files

cc = CorrectedCounts("data/SA501X3F_filtered_corrected_counts.csv", delimiter=',')
cc.add_chromosome_ends(neutral_cn=2, end_bin_length=150000)
data_converter2 = DataConverter(event_length_normalizer=3095677412)
data_converter2.create_CoNET_input_files(bin_dir, chromosomes=[17, 18, 20, 23], corrected_counts=cc)


# this may take up to 10 minutes
conet = CONET(bin_dir + "CONET")
params = CONETParameters(tree_structure_prior_k1 = 0.5, data_dir = bin_dir, counts_penalty_s1=200000, counts_penalty_s2=200000,
                            param_inf_iters=30000, seed = 21567, mixture_size=2, pt_inf_iters=200000, neutral_cn =2.0, output_dir = "./output/")
conet.infer_tree(params)