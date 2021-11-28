#include <utility>
#include <chrono>
#include <cstdlib>
#include <tuple>
#include <chrono>

#include "src/tree/pointer_tree.h"
#include "src/cell_provider/csv_reader/csv_reader.h"
#include "src/tree_sampler_coordinator.h"
#include "src/utils/random.h"
#include "src/likelihood/implementations/utils/csv_reader.h"
#include "src/parallel_tempering_coordinator.h"

using namespace std;

using NodeHandle = PointerTree::NodeHandle;

void save_hmm_matrix(std::string path, PointerTree &tree, std::vector<BreakpointPair> attachment, size_t no_of_loci) {
	std::ofstream file{ path };
    std::vector<NodeHandle> nodes = tree.getNodes();
	std::map<BreakpointPair, NodeHandle> breakpoint_to_node;
    std::map<NodeHandle, std::set<size_t>> ancestors;

	for (auto node : nodes) {
		breakpoint_to_node[tree.getNodeBreakpoints(node)] = node;
		std::unordered_set<size_t> breakpoints;
		tree.gatherAncestorsBreakpoints(node, breakpoints);
		ancestors[node] = std::set<size_t>(breakpoints.begin(), breakpoints.end());
	}

    for (size_t c = 0; c < attachment.size(); c++) {
        auto brkps = ancestors[breakpoint_to_node[attachment[c]]];
        for (size_t l = 0; l < no_of_loci; l++) {
			std::string value = brkps.find(l) == brkps.end() ? "0.0" : "1.0";
			std::string separator = l < no_of_loci - 1 ? ";" : "\n";
			file << value << separator;
        }
    }
}

void save_attachment(std::string data_path, std::string path, std::map<size_t, std::string> loci_names, std::vector<BreakpointPair> attachment) {

	std::ifstream in(data_path.append("cell_names"));
	std::string str;
	std::vector<std::string> cells;
	while (std::getline(in, str))
	{
	    if(str.size() > 0)
		cells.push_back(str);
	}
	
	std::ofstream file{ path };
	for (size_t j = 0; j < attachment.size(); j++)
	{
		if (attachment[j].first == 0 && attachment[j].second == 0) {
				file << cells[j] << ";" << j << ";" << attachment[j].first << ";" << attachment[j].second << "\n";
		} else {
			file << cells[j] << ";" << j << ";" << loci_names[attachment[j].first] << ";" << loci_names[attachment[j].second] << "\n";
		}
	}
}

void save_edge_confidence(std::string path, PointerTree &tree, std::map<std::pair<BreakpointPair, BreakpointPair>, size_t> &edge_count, std::map<size_t, std::string> loci_to_name) {
	std::ofstream file{ path };
    auto edges = tree.gatherEdges();

	for (auto edge : edges) {
		file << tree.get_node_label(edge.first, loci_to_name) << ";" << tree.get_node_label(edge.second, loci_to_name) << ";" << edge_count[edge] << " \n";
    }
}

void save_parameters(string path, NormalMixtureLikelihood<double> likelihood) {
    std::ofstream file{ path };
    file << likelihood.noBrkpNormal.mean << ";" << likelihood.noBrkpNormal.sd << "\n";
    for (size_t i = 0; i < likelihood.mixture.weights.size(); i++) {
        file << likelihood.mixture.weights[i] << ";" << -likelihood.mixture.gaussians[i].mean << ";" << likelihood.mixture.gaussians[i].sd << "\n";
    }
}

/**
Input parameters
   1. data dir - path to the directory where results should be saved
   2. number of iterations for parameter estimation
   3. number of iterations for parallel tempering
   4. use regularization penalty 
   THREADS_NUM;
   2. events length penalty - 
   3. attachment probability version -
   4. regularization penalty constant - 
   5.

   extern int PARAMETER_RESAMPLING_FREQUENCY;
extern int NUMBER_OF_MOVES_BETWEEN_SWAPS;
**/


const int PARAMETERS_COUNT = 18;

std::tuple<std::string, size_t, size_t, std::string> read_parameters(char **argv) {
	std::string data_dir{ argv[1] };
	size_t iterations_parameters = (size_t) std::stoi(argv[2]);
	size_t iterations_pt = (size_t)std::stoi(argv[3]);
	
	
	COUNTS_SCORE_CONSTANT = std::stod(argv[4]);
	EVENTS_LENGTH_PENALTY = std::stod(argv[5]);
	DATA_SIZE_PRIOR_CONSTANT = std::stod(argv[6]);
	USE_EVENT_LENGTHS_IN_ATTACHMENT = std::stoi(argv[7]);
	
	SEED = std::stol(argv[8]);
	MIXTURE_SIZE = (size_t) std::stoi(argv[9]);

	THREADS_NUM = (size_t) std::stoi(argv[10]);
	THREADS_LIKELIHOOD = (size_t)std::stoi(argv[11]);
	PARAMETER_RESAMPLING_FREQUENCY = (size_t)std::stoi(argv[12]);
	NUMBER_OF_MOVES_BETWEEN_SWAPS = (size_t)std::stoi(argv[13]);

	BURNIN = (size_t)std::stoi(argv[14]);
    VERBOSE = std::stoi(argv[15]);
    NEUTRAL_CN = std::stod(argv[16]);
    std::string output_dir{ argv[17] };
    
    if(data_dir.back() != '/')
    	data_dir.push_back('/');
    if(output_dir.back() != '/')
    	output_dir.push_back('/');
    return make_tuple(data_dir, iterations_parameters, iterations_pt, output_dir);
}

int main(int argc, char **argv)
{
	if (argc != PARAMETERS_COUNT) {
		logErr("Invalid number of parameters");
		logErr("Received: ", argc, " expected: ", PARAMETERS_COUNT); 
		return -1;
	}
	
	auto parameters = read_parameters(argv);
	auto data_dir = get<0>(parameters);
	auto output_dir = get<3>(parameters);
	std::ofstream tree_file{ string(output_dir).append("inferred_tree") };

	Random<double> random(SEED);
    VectorCellProvider<double> provider = readFile(string(data_dir).append("ratios"), string(data_dir).append("counts"), string(data_dir).append("counts_squared"), ';', COUNTS_SCORE_CONSTANT != 0.0);
    
    log("Input files have been loaded succesfully");
    ParallelTemperingCoordinator<double> PT(provider, random);
	auto result = PT.simulate(get<1>(parameters), get<2>(parameters));
	log("Tree inference has finished");

	PointerTree inferred_tree = std::get<0>(std::get<0>(result));
	tree_file << inferred_tree.toString(provider.get_loci_to_name());

    save_hmm_matrix(string(output_dir).append("inferred_breakpoints"), get<0>(get<0>(result)), get<1>(get<0>(result)), provider.getLociCount());
    save_edge_confidence(string(output_dir).append("edge_confidence"), std::get<0>(get<0>(result)), std::get<3>(result), provider.get_loci_to_name());
    save_attachment(data_dir, string(output_dir).append("inferred_attachment"), provider.get_loci_to_name(), std::get<1>(get<0>(result)));
    save_parameters(string(output_dir).append("inferred_distribution"), get<4>(result));
    return 0;
}
