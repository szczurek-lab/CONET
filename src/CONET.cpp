#include <utility>
#include <chrono>
#include <cstdlib>
#include <tuple>
#include <chrono>
#include <fstream>

#include "src/tree/event_tree.h"
#include "src/input_data/csv_reader.h"
#include "src/tree_sampler_coordinator.h"
#include "src/utils/random.h"
#include "src/parallel_tempering_coordinator.h"
#include "src/tree/tree_formatter.h"
#include "src/conet_result.h"

#include <boost/program_options.hpp>

using namespace std;

namespace po = boost::program_options;

int main(int argc, char **argv) {
	po::options_description description("MyTool Usage");

	description.add_options()
		("data_dir", po::value<string>()->required(),  "Path to directory containing input files.")
		("output_dir",  po::value<string>()->required(), "Path to output directory. Inference results will be saved there.")
		("param_inf_iters",  po::value<int>()->default_value(100000), "Number of MCMC iterations for joint tree and model parameters inference.")
		("pt_inf_iters",  po::value<int>()->default_value(100000), "Number of MCMC iterations for tree inference.")
		("counts_penalty_s1",  po::value<double>()->default_value(0.0), "Constant controlling impact of penalty for large discrepancies between inferred and real count matrices.")
		("counts_penalty_s2",  po::value<double>()->default_value(0.0), "Constant controlling impact of penalty for inferring clusters with changed copy number equal to basal ploidy.")
		("event_length_penalty_k0",  po::value<double>()->default_value(1.0), "Constant controlling impact of penalty for long inferred events.")
		("tree_structure_prior_k1",  po::value<double>()->default_value(1.0), "Constant controlling impact of data size part of tree structure prior.")
		("use_event_lengths_in_attachment",  po::value<bool>()->default_value(true), "If True cell attachment probability will depend on average event length in the history, otherwise it will be uniform.")
		("seed",  po::value<int>()->default_value(12312), "Seed for C++ RNG")
		("mixture_size",  po::value<size_t>()->default_value(4), "Initial number of components in difference distribution for breakpoint loci.")
		("num_replicas",  po::value<size_t>()->default_value(5), "Number of tempered chain replicas in MAP event tree search.")
		("threads_likelihood",  po::value<size_t>()->default_value(4), "Number of threads which will be used for the most demanding likelihood calculations.")
		("verbose",  po::value<bool>()->default_value(true), "True if CONET should print messages during inference.")
		("neutral_cn",  po::value<double>()->default_value(2.0), "Neutral copy number");
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(description).run(), vm);
	po::notify(vm);
	
	auto data_dir = vm["data_dir"].as<string>();
	auto output_dir = vm["output_dir"].as<string>();
    if(data_dir.back() != '/') data_dir.push_back('/');
    if(output_dir.back() != '/') output_dir.push_back('/');

	auto param_inf_iters = vm["param_inf_iters"].as<int>();
	auto pt_inf_iters = vm["pt_inf_iters"].as<int>();

	COUNTS_SCORE_CONSTANT_0 = vm["counts_penalty_s1"].as<double>();
	COUNTS_SCORE_CONSTANT_1 = vm["counts_penalty_s2"].as<double>();
	EVENTS_LENGTH_PENALTY = vm["event_length_penalty_k0"].as<double>();
	DATA_SIZE_PRIOR_CONSTANT = vm["tree_structure_prior_k1"].as<double>();
	USE_EVENT_LENGTHS_IN_ATTACHMENT = vm["use_event_lengths_in_attachment"].as<bool>();
	SEED = vm["seed"].as<int>();
	MIXTURE_SIZE = vm["mixture_size"].as<size_t>();
	NUM_REPLICAS = vm["num_replicas"].as<size_t>();
	THREADS_LIKELIHOOD = vm["threads_likelihood"].as<size_t>();
    VERBOSE = vm["verbose"].as<bool>();
    NEUTRAL_CN = vm["neutral_cn"].as<double>();

	Random<double> random(SEED);
    CONETInputData<double> provider = create_from_file(string(data_dir).append("ratios"), string(data_dir).append("counts"), string(data_dir).append("counts_squared"), ';');
    
    log("Input files have been loaded successfully");
    ParallelTemperingCoordinator<double> PT(provider, random);
	CONETInferenceResult<double> result = PT.simulate(param_inf_iters, pt_inf_iters);
	log("Tree inference has finished");

	std::ofstream tree_file{ string(output_dir).append("inferred_tree") };
	tree_file << TreeFormatter::to_string_representation(result.tree);

	std::ofstream attachment_file{ string(output_dir).append("inferred_attachment") };
	attachment_file << result.attachment;
    return 0;
}
