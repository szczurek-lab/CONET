#include <string>
#include<fstream>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include <exception>
#include <stdexcept>

#include "csv_reader.h" 


std::vector<std::vector<std::string>> getData(std::string path, char delimiter)
{
	std::ifstream file(path);

	std::vector<std::vector<std::string>> dataList;

	std::string line = "";
	while (getline(file, line))
	{
		std::vector<std::string> vec;
		std::istringstream f(line);
		std::string s;
		while (getline(f, s, delimiter)) {
			vec.push_back(s);
		}
		dataList.push_back(vec);
	}
	file.close();

	return dataList;
}

std::vector<std::vector<double>> stringDataToDouble(std::vector<std::vector<std::string>> data) {
	std::vector<std::vector<double>> result(data.size());
	for (size_t i = 0; i < result.size(); i++) {
		result[i].resize(data[0].size());
        try {
		std::transform(data[i].begin(), data[i].end(), result[i].begin(),
				[](std::string s) -> double {  
                            return std::stold(s);
                     });
        } catch(std::invalid_argument& e) {
                            std::cout << e.what();
                            throw e;
                    }
	};
	return result;
}

std::vector<size_t> getChromosomeMarkers(std::vector<double> data) {
	std::vector<size_t> result;
	double previous = data[0];
	for (size_t counter = 0; counter < data.size(); counter++) {
		if (data[counter] != previous) {
			result.push_back(counter);
			previous = data[counter];
		}
	}
	result.push_back(data.size());
	return result;
}

VectorCellProvider<double> dataToCellProvider(std::vector<std::vector<double>> data, std::vector<std::string> names, std::vector<std::string> chromosomes) {
	VectorCellProvider<double> provider(data[0].size());
	provider.setChromosomeMarkers(getChromosomeMarkers(data[0]));
    provider.setBetweenLengths(data[1]);
	for (size_t i = 2; i < data.size(); i++) {
		for (auto &el: data[i])
		{
			el = -std::abs(el);
		}
		provider.postCell(data[i]);
	}
	
	provider.set_loci_to_name_map(names, chromosomes);
	return provider;
}

VectorCellProvider<double> getCellProvider(std::string path, char delimiter) {
	auto string_data = getData(path, delimiter);
	auto names = string_data[0];
	auto chromosomes = string_data[1];
	string_data.erase(string_data.begin());
	auto data = stringDataToDouble(string_data);
	return dataToCellProvider(data, names, chromosomes);
}

void read_counts_penalty_files(VectorCellProvider<double> &cells, std::string sum_counts_path, std::string squared_counts_path, char delimiter) {
	auto counts_sum = getData(sum_counts_path, delimiter);
	auto squared_counts = getData(squared_counts_path, delimiter);
	auto regions_sizes = counts_sum[0];
	counts_sum.erase(counts_sum.begin());
	squared_counts.erase(squared_counts.begin());
	std::vector<std::vector<std::string>> regions_sizes_2{ regions_sizes };
	auto region_sizes_double = stringDataToDouble(regions_sizes_2);
	auto counts_sum_double = stringDataToDouble(counts_sum);
	auto squared_counts_double = stringDataToDouble(squared_counts);
	cells.post_counts_dispersion_data(region_sizes_double[0], counts_sum_double, squared_counts_double);
}

VectorCellProvider<double> readFile(std::string path, std::string sum_counts_path, std::string squared_counts_path, char delimiter, bool use_counts_penalty_files) {
	auto cell_provider = getCellProvider(path, delimiter);
	if (use_counts_penalty_files) {
		read_counts_penalty_files(cell_provider, sum_counts_path, squared_counts_path, delimiter);
	}
	return cell_provider;
}
