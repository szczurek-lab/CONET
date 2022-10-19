#ifndef MATRIX_H
#define MATRIX_H

#include <algorithm>
#include <vector>

namespace Matrix {

template <class Type>
std::vector<std::vector<Type>> create_2d_matrix(size_t rows, size_t columns,
                                                Type value) {
  std::vector<std::vector<Type>> result;
  result.resize(rows);
  for (auto &row : result) {
    row.resize(columns);
    std::fill(row.begin(), row.end(), value);
  }
  return result;
}

template <class Type>
std::vector<Type> create_1d_matrix(size_t size, Type value) {
  std::vector<Type> result;
  result.resize(size);
  std::fill(result.begin(), result.end(), value);
  return result;
}

template <class Real_t>
void normalize_1d_matrix_elements(std::vector<Real_t> &matrix) {
  Real_t sum = std::accumulate(matrix.begin(), matrix.end(), (Real_t)0.0);
  for (auto &entry : matrix) {
    entry = entry / sum;
  }
}
} // namespace Matrix

#endif