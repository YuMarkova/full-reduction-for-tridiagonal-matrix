#include <iostream>
#include <omp.h>
#include "input_output.h"
#include "reduction.h"

std::vector<double> A;
std::vector<double> C;
std::vector<std::vector<double>> P;


void Init (size_t n){
    A.clear();
    C.clear();
    P.clear();
    for (size_t i = 0; i < n; ++i) {
        A.push_back(-1);
        C.push_back(4);
    }
    for (size_t i = 0; i < n; ++i){
        std::vector<double> tmp;
        for (size_t j = 0; j < n; ++j)
            tmp.push_back(16);
        P.push_back(tmp);
    }
}


int main()
{
    std::ifstream in("input.txt");
    std::ofstream out("test.txt");
    size_t n;
    std::cin >> n;

    for (int i = 0; i < 200; i++){

        Init(n);
        //float start1  = omp_get_wtime();
        straight(P, A, C);
        //float finish1 = omp_get_wtime();

        //float start2  = omp_get_wtime();
        reverse(P, A, C);
        //float finish2 = omp_get_wtime();

        out << P << '\n';
    }

    in.close();
    out.close();
    return 0;
}

