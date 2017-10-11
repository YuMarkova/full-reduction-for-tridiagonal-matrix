#include <vector>
#include <algorithm>
#include <cmath>
#include <omp.h>
#include <iostream>
#include "reduction.h"
#include "vector_operations.h"

#define MIN_PAR_SIZE 1

inline void threads_configure() {
    int nThreads = std::max(2, omp_get_max_threads());
    omp_set_num_threads(nThreads);
}


std::vector<double> Sweeping (const std::vector<double> &A,
                              const std::vector<double> &C,
                              const std::vector<double> &F)
{
    std::vector<double> alp, bet, X;

    alp.push_back(-A[0] / C[0]);

    bet.push_back(F[0] / C[0]);
    size_t N = A.size();

    for (size_t i = 1; i < N - 1; ++i) {
        alp.push_back(-A[i] /
                     (A[i] * alp[i-1] + C[i]));

        bet.push_back((F[i] - A[i] * bet[i-1]) /
                      (A[i] * alp[i-1] + C[i]));
    }

    X.push_back((F[N-1] - A[N-1] * bet[N-2]) /
                (C[N-1] + A[N-1] * alp[N-2]));

    for (size_t i = 0; i < N - 1; ++i)
        X.push_back(X[i] * alp[N-i-2] + bet[N-i-2]);

    std::reverse(&X[0], &X[N]);
    return X;
}

std::vector<double> UoverU(const size_t k, const size_t n, const std::vector<double> &A, const std::vector<double> &C, const std::vector<double> &F)
{
    std::vector<double> R;
    for (size_t i = 0; i < F.size(); ++i)
        R.push_back(0);

    for (size_t s = 0; s < n; s++) {
        if (((k + 1) * (s + 1)) % (n + 1) == 0) continue;
        double a = 2 * pow(-1, s) *
                sin((k + 1) * M_PI * (s + 1) / (n + 1)) *
                sin(M_PI * (s + 1) / (n + 1)) / (n + 1);
        double b = - 2 * cos(M_PI * (s + 1) / (n + 1));

        R = sum(R, Sweeping(A, sum(b, C), mult(a, F)));
    }

    return R;
}

std::vector<double> UUoverU(size_t k, size_t l, size_t n, const std::vector<double> &A, const std::vector<double> &C, const std::vector<double> &F) {
    std::vector<double> R;
    for (size_t i = 0; i < F.size(); ++i)
        R.push_back(0);

    for(size_t s = 0; s < n; ++s){
        if (!(((k + 1) * (s + 1)) % (n + 1)) || !(((l + 1) * (s + 1)) % (n + 1))) continue;
        double a = 2 * pow(-1, s) *
                   sin((k + 1) * M_PI * (s + 1) / (n + 1)) *
                   sin((l + 1) * M_PI * (s + 1) / (n + 1)) / (n + 1);
        double b = - 2 * cos(M_PI * (s + 1) / (n + 1));

        R = sum(R, Sweeping(A, sum(b, C), mult(a, F)));
    }
    return R;
}

std::vector<std::vector<double>> straight (std::vector<std::vector<double>> &P, const std::vector<double> &A, const std::vector<double> &C)
{
    int n = P[0].size();
    int m = log(n) / log(2);

    for (int k = 0; k < m; ++k) {
        int twoPowK = pow(2, k);
        int step = pow(2, k + 1);

        threads_configure();
        #pragma omp parallel for if (n > MIN_PAR_SIZE)
        for (int i = twoPowK - 1; i < n; i += step) {

            int l = i - twoPowK;

            int r = std::min(n, i + twoPowK);

            if (l >= 0) {
                std::vector<double> V;
                V = UoverU(r - i - 1, r - l - 1, A, C, P[i]);
                for (size_t j = 0; j < P[l].size(); ++j)
                    P[l][j] += V[j];
            }
            if (r < n) {
                std::vector<double> V;
                V = UoverU(i - l - 1, r - l - 1, A, C, P[i]);
                for (size_t j = 0; j < P[r].size(); ++j)
                    P[r][j] += V[j];

            }
        }
    }

    return P;
}

std::vector<std::vector<double>> reverse (std::vector<std::vector<double>> &P, const std::vector<double> &A, const std::vector<double> &C) {
    int n = P[0].size();
    int m = log(n) / log(2);

    for (int k = m; k >= 0; --k) {
        int twoPowK = pow(2, k);
        int step = pow(2, k + 1);

        threads_configure();
        #pragma omp parallel for if (n > MIN_PAR_SIZE)
        for (int i = twoPowK - 1; i < n; i += step) {
            int l = i - twoPowK;
            int r = std::min(i + twoPowK, n);

            P[i] = UUoverU(r - i - 1, i - l - 1, r - l - 1, A,  C, P[i]);

            if (l >= 0)
                P[i] = sum(P[i], UoverU(r - i - 1, r - l - 1, A, C, P[l]));
            if (r < n)
                P[i] = sum(P[i], UoverU(i - l - 1, r - l - 1, A, C, P[r]));
        }
    }

    return P;
}

