#ifndef MATRIX_H
#define MATRIX_H

#include <vector> 
#include <algorithm> 

namespace Matrix {
    
    template <class Type> std::vector<std::vector<Type>> create_2d_matrix(size_t rows, size_t columns, Type value) {
        std::vector<std::vector<Type>> result; 
        result.resize(rows);
        for (auto &row : result) {
            row.resize(columns);
            std::fill(row.begin(), row.end(), value);
        }
        return result;
    }


}


#endif 