#ifndef RHO_REFLECTRHOGRAD_H
#define RHO_REFLECTRHOGRAD_H
#include <vector>
namespace reflectRhoGrad
{
    void correctRho(double &varM, double varP);

    void correctGradRho(int edgeGrp, std::vector<double> &gradM, const std::vector<double> &gradP, const std::vector<double> &n);
}
#endif // RHO_REFLECTRHOGRAD_H
