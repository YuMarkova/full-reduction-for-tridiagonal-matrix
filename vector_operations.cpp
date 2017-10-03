#include <vector>
#include <omp.h>
#include <iostream>
#include "vector_operations.h"

#define MIN_PAR_SIZE 1000

inline void threads_configure() {
    int nThreads = std::max(2, omp_get_max_threads());
    omp_set_num_threads(nThreads);
}


std::vector<double> sum (const std::vector<double> &vec1,
                         const std::vector<double> &vec2)
{
    std::vector<double> vecSum = vec1;

    if (vec1.size() != vec2.size())
        std::cerr << "Vectors are not equalent" << std::endl;
    else {
        threads_configure();
        #pragma omp parallel for if (vec1.size() > MIN_PAR_SIZE)
        for (size_t i = 0; i < vec1.size(); ++i)
            vecSum[i] += vec2[i];
    }

    return vecSum;
}

std::vector<double> sum (double x, const std::vector<double> &vec) {
    std::vector<double> vecSum = vec;

    threads_configure();
    #pragma omp parallel for if (vec.size() > MIN_PAR_SIZE)
    for (size_t i = 0; i < vec.size(); ++i)
        vecSum[i] += x;

    return vecSum;
}

std::vector<double> mult (double x, const std::vector<double> &vec) {
    std::vector<double> vecSum = vec;

    threads_configure();
    #pragma omp parallel for if (vec.size() > MIN_PAR_SIZE)
    for (size_t i = 0; i < vec.size(); ++i)
        vecSum[i] *= x;

    return vecSum;
}

