#ifndef CSV_READER_H
#define CSV_READER_H
#include "../vector_cell_provider.h"



VectorCellProvider<double> readFile(std::string path, std::string sum_counts_path, std::string squared_counts_path, char delimiter, bool use_counts_penalty_files);




#endif // !CSV_READER_H
