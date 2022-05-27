#include <string>
#include<fstream>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>

#include "csv_reader.h" 

/**
 * @brief Split files by delimiter 
 * 
 * @param path Path to file
 * @param delimiter 
 * @return std::vector<std::vector<std::string>> Each row corresponds to one line, columns are lines split by the delimiter 
 */

std::vector<std::vector<std::string>> split_file_by_delimiter(std::string path, char delimiter)
{
	std::ifstream file(path);
	std::vector<std::vector<std::string>> data;
	std::vector<std::string> line_placeholder;
	std::string buffer, line;

	while (getline(file, line))
	{
		std::istringstream f(line);
		while (getline(f, buffer, delimiter)) {
			line_placeholder.push_back(buffer);
		}
		data.push_back(line_placeholder);
		line_placeholder.clear();
	}
	return data;
}

std::vector<std::vector<double>> string_matrix_to_double(std::vector<std::vector<std::string>> data) {
	std::vector<std::vector<double>> double_data(data.size());
	for (size_t i = 0; i < double_data.size(); i++) {
		double_data[i].resize(data[0].size());
		std::transform(data[i].begin(), data[i].end(), double_data[i].begin(), [](std::string s) -> double {  return std::stold(s); });
	};
	return double_data;
}

/**
 * @brief Convert vector of chromosomes to chromosome markers
 * 
 * Create vector of positive numbers, where each element from the vector @i is such that 
 * @chromosomes[@i] != @chromosomes[@i - 1]. 
 */
std::vector<size_t> get_chromosome_markers(std::vector<double> chromosomes) {
	std::vector<size_t> result;
	double previous = chromosomes[0];
	for (size_t counter = 0; counter < chromosomes.size(); counter++) {
		if (chromosomes[counter] != previous) {
			result.push_back(counter);
			previous = chromosomes[counter];
		}
	}
	result.push_back(chromosomes.size());
	return result;
}

VectorCellProvider<double> diff_matrix_to_cell_provider(std::vector<std::vector<double>> data) {
	VectorCellProvider<double> provider(data[0].size(), get_chromosome_markers(data[0]), data[1]);
	for (size_t i = 2; i < data.size(); i++) {
		std::for_each(data[i].begin(), data[i].end(), [](double &r){r = -std::abs(r);});  
		provider.post_cell(data[i]);
	}
	return provider;
}

void read_counts_penalty_files(VectorCellProvider<double> &provider, std::string summed_counts_path, std::string squared_counts_path, char delimiter) {
	auto summed_counts = string_matrix_to_double(split_file_by_delimiter(summed_counts_path, delimiter));
	auto squared_counts = string_matrix_to_double(split_file_by_delimiter(squared_counts_path, delimiter));
	auto regions_sizes = summed_counts[0];
	summed_counts.erase(summed_counts.begin());
	squared_counts.erase(squared_counts.begin());
	provider.post_counts_dispersion_data(regions_sizes, summed_counts, squared_counts);
}

VectorCellProvider<double> create_from_file(std::string path, std::string summed_counts_path, std::string squared_counts_path, char delimiter) {
	VectorCellProvider<double> provider = diff_matrix_to_cell_provider(string_matrix_to_double(split_file_by_delimiter(path, delimiter)));
	read_counts_penalty_files(provider, summed_counts_path, squared_counts_path, delimiter);
	return provider;
}
