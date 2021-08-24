#ifndef U_INTERIORSIDE_H
#define U_INTERIORSIDE_H
#include <fstream>
#include <vector>

namespace interiorSide {
    void u_IO(int bcGrp, std::ifstream &FileFlux);

    void correctU(std::vector<double> &varM, std::vector<double> varP);

    void correctGradU(std::vector<double> &gradM, const std::vector<double> &gradP);
}

#endif // U_INTERIORSIDE_H
