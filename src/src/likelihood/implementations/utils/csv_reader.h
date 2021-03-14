#ifndef PARAMETERS_CSV_READER_H
#define PARAMETERS_CSV_READER_H

#include <utility>
#include <vector>

#include <string>
#include<fstream>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "../utils/gaussian_mixture.h"
#include "../utils/gaussian.h"
#include "../../../utils/random.h"

namespace Gauss {
	template <class Real_t> std::vector<std::vector<Real_t>> stringDataToDouble(std::vector<std::vector<std::string>> data) {
		std::vector<std::vector<Real_t>> result(data.size());
		for (size_t i = 0; i < result.size(); i++) {
			result[i].resize(data[0].size());
			std::transform(data[i].begin(), data[i].end(), result[i].begin(),
				[](std::string s) -> Real_t { return std::stold(s); });
		};
		return result;
	}

	std::vector<std::vector<std::string>> getData(std::string path, char delimiter)
	{
		std::ifstream file(path);

		std::vector<std::vector<std::string> > dataList;

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

	template<class Real_t> std::pair<GaussianMixture<Real_t>, Gaussian<Real_t>> readCSV(std::string path, char delimiter, Random<Real_t> &random) {
		auto data = stringDataToDouble<Real_t>(getData(path, delimiter));
		Gauss::GaussianMixture<Real_t> mixt(data[0], data[1], data[2], random);
		Gauss::Gaussian<Real_t> gauss(data[3][0], data[3][1], random);
		return std::make_pair(mixt, gauss);
	}
}

#endif // !PARAMETERS_CSV_READER_H
