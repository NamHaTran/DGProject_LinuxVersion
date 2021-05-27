#include "./extNSFEqns/FranzDurst/boundaryConditions.h"
#include <tuple>
#include <vector>
#include "DGAuxUltilitiesLib.h"
#include "dynamicVarDeclaration.h"
#include "./boundaryConditions/bcVariables.h"
#include "DGMath.h"

namespace bcForExtNSF_Durst {
double tempdTx(0.0), tempdTy(0.0);

    void dropNormDiffVel(int edge, int edgeGrp, double &dRhoXM, double &dRhoYM, double dRhoXP, double dRhoYP, double rhoBC, double dTBCx, double dTBCy, const std::vector<double> &n)
    {
        int UType(bcValues::UBcType[edgeGrp - 1]), //Lay UType de xac dinh kieu dk bien
            loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));

        double TBC(SurfaceBCFields::TBc[loc]), nx(n[0]), ny(n[1]), dRhoBCx(0.0), dRhoBCy(0.0);

        //Chi drop normal diffusive velocity khi bien dang set la wall
        if (UType==BCVars::velocityBCId::noSlip
                || UType==BCVars::velocityBCId::movingWall
                || UType==BCVars::velocityBCId::slipWall)
        {
            double muBC(math::CalcVisCoef(TBC));
            //Calculate coefficients
            double a1=nx*TBC*material::R,
                    b1=ny*TBC*material::R,
                    c1=material::R*rhoBC*(nx*dTBCx + ny*dTBCy),
                    a2=-(muBC/rhoBC)*nx/rhoBC,
                    b2=-(muBC/rhoBC)*ny/rhoBC,
                    c2=-(muBC/rhoBC)*(nx*dTBCx + ny*dTBCy)/(2*TBC);

            dRhoBCy=(c2*a1-c1*a2)/(b1*a2-b2*1);
            dRhoBCx=(-b1*dRhoBCy-c1)/a1;

            //Correct div(rho)
            dRhoXM = 2*dRhoBCx-dRhoXP;
            dRhoYM = 2*dRhoBCy-dRhoYP;
        }
    }
}
