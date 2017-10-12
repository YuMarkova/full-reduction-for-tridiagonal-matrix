#include <iostream>
#include <omp.h>
#include "input_output.h"
#include "reduction.h"

std::vector<double> A;
std::vector<double> C;
std::vector<std::vector<double>> P;


void Init (size_t n, size_t m, double value) {
    P.clear();
    for (size_t i = 0; i < n; ++i) {
        std::vector<double> tmp;
        for (size_t j = 0; j < m; ++j)
            tmp.push_back(value);
        P.push_back(tmp);
    }
}
void Init (size_t m) {
    A.clear();
    C.clear();
    for (size_t i = 0; i < m; ++i) {
        A.push_back(-1);
        C.push_back(4);
    }
}


int main()
{
    std::ifstream in("input.txt");
    std::ofstream out("test.txt");
    size_t bMatrixSize, bSize;
    double value = 16;

    std::cin >> bMatrixSize >> bSize;

    Init(bSize);
    Init(bMatrixSize, bSize, value);

    float start1  = omp_get_wtime();
    straight(P, A, C);
//    float finish1 = omp_get_wtime();
//    std::cout << 2 << std::endl;
//    float start2  = omp_get_wtime();
    reverse(P, A, C);
    float finish2 = omp_get_wtime();

    std::cout << finish2 - start1 << std::endl;
//    std::cout << finish2 - start2 << std::endl;
    out << P << '\n';

    in.close();
    out.close();
    return 0;
}

