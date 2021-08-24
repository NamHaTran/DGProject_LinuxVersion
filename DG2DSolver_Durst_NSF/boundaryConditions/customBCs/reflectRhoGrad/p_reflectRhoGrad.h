#ifndef P_REFLECTRHOGRAD_H
#define P_REFLECTRHOGRAD_H
#include <fstream>
#include <vector>
namespace reflectRhoGrad {
    void p_IO(int bcGrp, std::ifstream &FileFlux);

    void correctP(double &pM, double pP);

    void correctGradP(std::vector<double> &gradM, const std::vector<double> &gradP);
}

#endif // P_REFLECTRHOGRAD_H
