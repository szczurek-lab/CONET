#ifndef CSV_READER_H
#define CSV_READER_H
#include "input_data.h"

CONETInputData<double> create_from_file(std::string path,
                                        std::string summed_counts_path,
                                        std::string squared_counts_path,
                                        char delimiter);

#endif // !CSV_READER_H
