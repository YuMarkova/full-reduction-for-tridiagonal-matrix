#include <fstream>
#include <iomanip>
#include "input_output.h"

std::ostream &operator<<(std::ostream &out, const std::vector<std::vector<double>> &obj) {

    for (size_t i = 0; i < obj.size(); ++i){
        for (size_t j = 0; j < obj[0].size(); j++)
            out << std::fixed << std::setprecision(3) << '\t' << obj[i][j];
        out << std::endl;
    }
    return out;
}

std::istream &operator>>(std::istream &in, std::vector<std::vector<double>> &obj) {
    size_t N;
    in >> N;

    for (size_t i = 0; i < N; ++i) {
        std::vector<double> tmp;
        double var;
        for (size_t j = 0; j < N; ++j) {
            in >> var;
            tmp.push_back(var);
        }
        obj.push_back(tmp);
    }
    return in;
}
