#ifndef U_MAXWELLSLIP_H
#define U_MAXWELLSLIP_H
#include <fstream>
#include <vector>

namespace MaxwellSlip {
    void u_IO(int bcGrp, std::ifstream &FileFlux);

    //void correctU(int edge, int edgeGrp, std::vector<double> &varM, std::vector<double> varP);

    //void correctGradU(std::vector<double> &gradM, const std::vector<double> &gradP);

    void calcUSlip_DGTypeExplicit(int edge, int edgeGrp);

    void calcUSlip_FDMTypeImplicit(int edge, int edgeGrp);
}

#endif // U_MAXWELLSLIP_H
