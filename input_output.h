#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include <fstream>
#include <vector>

std::ostream &operator<<(std::ostream &out, const std::vector <std::vector<double>> &obj);

std::istream &operator>>(std::istream &in, std::vector <std::vector<double>> &obj);

#endif // #ifndef INPUT_OUTPUT_H
