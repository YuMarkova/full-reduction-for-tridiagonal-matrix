#ifndef VECTOR_OPERATION_H
#define VECTOR_OPERATION_H

#include <vector>

std::vector<double> sum (const std::vector<double> &vec1,
                         const std::vector<double> &vec2);

std::vector<double> sum (double x, const std::vector<double> &vec);

std::vector<double> mult (double x, const std::vector<double> &vec);

#endif //VECTOR_OPERATION_H
