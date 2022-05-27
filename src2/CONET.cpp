#include <utility>
#include <chrono>
#include <cstdlib>
#include <tuple>
#include <chrono>

#include "src/tree/event_tree.h"
#include "src/cell_provider/csv_reader/csv_reader.h"
#include "src/tree_sampler_coordinator.h"
#include "src/utils/random.h"
#include "src/parallel_tempering_coordinator.h"
#include "src/tree/tree_formatter.h"
#include "src/conet_result.h"

using namespace std;

using NodeHandle = EventTree::NodeHandle;

void save_attachment(std::string data_path, std::string path, std::vector<Event> attachment) {

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
				file << cells[j] << ";" << j << ";" << attachment[j].first << ";" << attachment[j].second << "\n";
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


const int PARAMETERS_COUNT = 19;

std::tuple<std::string, size_t, size_t, std::string> read_parameters(char **argv) {
	std::string data_dir{ argv[1] };
	size_t iterations_parameters = (size_t) std::stoi(argv[2]);
	size_t iterations_pt = (size_t)std::stoi(argv[3]);
	
	
	COUNTS_SCORE_CONSTANT_0 = std::stod(argv[4]);
	COUNTS_SCORE_CONSTANT_1 = std::stod(argv[5]);
	EVENTS_LENGTH_PENALTY = std::stod(argv[6]);
	DATA_SIZE_PRIOR_CONSTANT = std::stod(argv[7]);
	USE_EVENT_LENGTHS_IN_ATTACHMENT = std::stoi(argv[8]);
	
	SEED = std::stol(argv[9]);
	MIXTURE_SIZE = (size_t) std::stoi(argv[10]);

	THREADS_NUM = (size_t) std::stoi(argv[11]);
	THREADS_LIKELIHOOD = (size_t)std::stoi(argv[12]);
	PARAMETER_RESAMPLING_FREQUENCY = (size_t)std::stoi(argv[13]);
	NUMBER_OF_MOVES_BETWEEN_SWAPS = (size_t)std::stoi(argv[14]);

	BURNIN = (size_t)std::stoi(argv[15]);
    VERBOSE = std::stoi(argv[16]);
    NEUTRAL_CN = std::stod(argv[17]);
    std::string output_dir{ argv[18] };
    
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
    VectorCellProvider<double> provider = create_from_file(string(data_dir).append("ratios"), string(data_dir).append("counts"), string(data_dir).append("counts_squared"), ';');
    
    log("Input files have been loaded succesfully");
    ParallelTemperingCoordinator<double> PT(provider, random);
	CONETInferenceResult<double> result = PT.simulate(get<1>(parameters), get<2>(parameters));
	log("Tree inference has finished");

	tree_file << TreeFormatter::to_string_representation(result.tree);

    save_attachment(data_dir, string(output_dir).append("inferred_attachment"), result.attachment);
    return 0;
}
