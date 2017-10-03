#ifndef REDUCTION_H
#define REDUCTION_H

#include <vector>

std::vector<double> Sweeping (const std::vector<double> &A,
                              const std::vector<double> &C,
                              const std::vector<double> &F);

std::vector<double> UoverU (size_t k, size_t n, const std::vector<double> &A, const std::vector<double> &C, const std::vector<double> &F);

std::vector<double> UUoverU(size_t k, size_t l, size_t n, const std::vector<double> &A, const std::vector<double> &C, const std::vector<double> &F);

std::vector<std::vector<double>> straight (std::vector<std::vector<double>> &P, const std::vector<double> &A, const std::vector<double> &C);

std::vector<std::vector<double>> reverse  (std::vector<std::vector<double>> &P, const std::vector<double> &A, const std::vector<double> &C);

#endif //REDUCTION_H
