#ifndef T_INTERIORSIDE_H
#define T_INTERIORSIDE_H
#include <fstream>
#include <vector>

namespace interiorSide {
    void T_IO(int bcGrp, std::ifstream &FileFlux);

    void correctT(double &varM, double varP);

    void correctGradT(std::vector<double> &gradM, const std::vector<double> &gradP);
}
#endif // T_INTERIORSIDE_H
