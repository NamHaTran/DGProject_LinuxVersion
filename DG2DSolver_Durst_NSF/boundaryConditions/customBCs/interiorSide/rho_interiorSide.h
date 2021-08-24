#ifndef RHO_INTERIORSIDE_H
#define RHO_INTERIORSIDE_H
#include <vector>
namespace interiorSide
{
    void correctRho(double &varM, double varP);

    void correctGradRho(std::vector<double> &gradM, const std::vector<double> &gradP, const std::vector<double> &n, bool isStrongMeth);
}
#endif // RHO_INTERIORSIDE_H
