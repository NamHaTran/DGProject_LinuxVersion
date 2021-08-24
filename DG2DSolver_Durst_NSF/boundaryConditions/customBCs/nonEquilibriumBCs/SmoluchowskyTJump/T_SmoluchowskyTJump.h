#ifndef T_SMOLUCHOWSKYTJUMP_H
#define T_SMOLUCHOWSKYTJUMP_H
#include <fstream>
#include <vector>

namespace SmoluchowskyTJump {
    void T_IO(int bcGrp, std::ifstream &FileFlux);

    //void correctT(int edge, int edgeGrp, double &varM, double varP);

    //void correctGradT(std::vector<double> &gradM, const std::vector<double> &gradP);

    void calcTJump_DGTypeExplicit(int edge, int edgeGrp);

    void calcTJump_FDMTypeImplicit(int edge, int edgeGrp);
}
#endif // T_SMOLUCHOWSKYTJUMP_H
