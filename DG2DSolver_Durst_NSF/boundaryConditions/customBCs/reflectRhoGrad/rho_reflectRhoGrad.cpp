#include "rho_reflectRhoGrad.h"
#include "./boundaryConditions/zeroGradient.h"
#include <vector>
#include "./boundaryConditions/bcVariables.h"

namespace reflectRhoGrad
{
    /* Cac dieu kien cua rho se overwrite cac dieu kien cua p */
    void correctRho(double &varM, double varP)
    {
        varM = zeroGradient_scalar(varP);
    }

    void correctGradRho(int edgeGrp, std::vector<double> &gradM, const std::vector<double> &gradP, const std::vector<double> &n)
    {

        double tangGrad(0.0), nx(n[0]), ny(n[1]);

        if (BCVars::NewmannAppMethGradGeneralBCStrong[edgeGrp-1])
        {
            /* For reference
            dUXMinus[i] = dUXPlus[i] - 2 * (dUXPlus[i] * nx + dUYPlus[i] * ny)*nx;
            dUYMinus[i] = dUYPlus[i] - 2 * (dUXPlus[i] * nx + dUYPlus[i] * ny)*ny;*/

            gradM[0] = gradP[0] - 2 * (gradP[0] * nx + gradP[1] * ny)*nx;
            gradM[1] = gradP[1] - 2 * (gradP[0] * nx + gradP[1] * ny)*ny;
        }
        else
        {
            /* For reference
            tangGrad = -dUXPlus[i]*ny+dUYPlus[i]*nx;
            dUXMinus[i]=-tangGrad*ny;
            dUYMinus[i]=tangGrad*nx;*/

            tangGrad = -gradP[0]*ny+gradP[1]*nx;
            gradM[0]=-tangGrad*ny;
            gradM[1]=tangGrad*nx;
        }
    }
}
