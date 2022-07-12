#ifndef P_ZERORHOGRADUNCORRECTP_H
#define P_ZERORHOGRADUNCORRECTP_H
#include <fstream>
#include <vector>
namespace zeroRhoGradUncorectP {
    void p_IO(int bcGrp, std::ifstream &FileFlux);

    void correctP(double &pM, double pP);

    void correctGradP(std::vector<double> &gradM, const std::vector<double> &gradP);
}

#endif // P_ZERORHOGRADUNCORRECTP_H
