#ifndef P_INTERIORSIDE_H
#define P_INTERIORSIDE_H
#include <fstream>
#include <vector>
namespace interiorSide {
    void p_IO(int bcGrp, std::ifstream &FileFlux);

    void correctP(double &pM, double pP);

    void correctGradP(std::vector<double> &gradM, const std::vector<double> &gradP);
}

#endif // P_INTERIORSIDE_H
