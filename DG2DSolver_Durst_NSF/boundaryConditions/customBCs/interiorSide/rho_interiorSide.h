#ifndef RHO_INTERIORSIDE_H
#define RHO_INTERIORSIDE_H
#include <vector>
namespace interiorSide
{
    /**
     * @brief Function corrects density following interiorSide condition (using zeroGradient_scalar function).
     * @param varM: minus side density
     * @param varP: plus side density
     */
    void correctRho(double &varM, double varP);

    /**
     * @brief Function corrects gradient of density following interiorSide condition.
     *
     * Note that current function removes normal gradient component.
     *
     * @param gradM: minus side of grad(rho).
     * @param gradP: plus side of grad(rho).
     * @param n: normal unit vector.
     * @param isStrongMeth: flag of application method.
     */
    void correctGradRho(std::vector<double> &gradM, const std::vector<double> &gradP, const std::vector<double> &n, bool isStrongMeth);
}
#endif // RHO_INTERIORSIDE_H
