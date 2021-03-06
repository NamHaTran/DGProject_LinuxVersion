#include "DGBCsLib.h"
#include <vector>
#include "DGMath.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <tuple>
#include "DGMessagesLib.h"
#include <math.h>
#include "DGProcLib.h"
#include <algorithm>
#include "BCSupportFncs.h"
#include "bcVariables.h"
#include "./parallelFunctions/generalParallelFuncs.h"
#include "./parallelFunctions/parallelVariables.h"

//Matched BC type
#include "matched.h"
#include "symmetry.h"

//Debug
#include <iostream>
#include "debuggingFuncs.h"

/* NOTES:
- All boundary functions must return numerical fluxes (rho, rhou, rhov, rhoE fluxes) at surface Gauss points!
*/

std::vector<double> NSFEqBCsImplement(int element, int edge, int nG)
{
    std::vector<double> Fluxes(4, 0.0);
    int edgeGrp(auxUlti::getGrpOfEdge(edge)),
            UType(bcValues::UBcType[edgeGrp - 1]),
            TType(bcValues::TBcType[edgeGrp - 1]),
            pType(bcValues::pBcType[edgeGrp - 1]);

    if (UType == BCVars::generalBCId::matched && TType == BCVars::generalBCId::matched && pType == BCVars::generalBCId::matched)
    {
        matchedBC_NSFEqs(Fluxes,element,edge,nG);
    }
    else if (UType == BCVars::generalBCId::symmetry && TType == BCVars::generalBCId::symmetry && pType == BCVars::generalBCId::symmetry)
    {
        symmetryBC_NSFEqs(Fluxes,element,edge,edgeGrp,nG);
    }
    else
    {
        NSFEqBCsForNormalBoundary(Fluxes,element,edge,edgeGrp,nG);
    }

    return Fluxes;
}

std::vector<std::vector<double>> auxEqBCsImplement(int element, int edge, int nG)
{
    //Bien phu duoc dat la mu*divU --> can tra ve gia tri mu*U sau khi chay ham nay
    std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
    int edgeGrp(auxUlti::getGrpOfEdge(edge)),
            UType(bcValues::UBcType[edgeGrp - 1]),
            TType(bcValues::TBcType[edgeGrp - 1]),
            pType(bcValues::pBcType[edgeGrp - 1]);

    if (UType == BCVars::generalBCId::matched && TType == BCVars::generalBCId::matched && pType == BCVars::generalBCId::matched)
    {
        if (flowProperties::viscous)
        {
            matchedBC_auxEqs(Fluxes,edge,nG);
        }
    }
    else if (UType == BCVars::generalBCId::symmetry && TType == BCVars::generalBCId::symmetry && pType == BCVars::generalBCId::symmetry)
    {
        symmetryBC_auxEqs(Fluxes,element,edge,edgeGrp,nG);
    }
    else
    {
        auxEqBCsForNormalBoundary(Fluxes,element,edge,edgeGrp,nG);
    }
    return Fluxes;
}

void auxEqBCsForNormalBoundary(std::vector<std::vector<double>> &Fluxes, int element, int edge, int edgeGrp, int nG)
{
    std::vector<double>
        UPlus(4, 0.0), UMean(4, 0.0),
        UMinus(4, 0.0),
        priVarsPlus(5, 0.0), //priVarsPlus = [rho u v p T]
        priVarsMean(5, 0.0),
        priVarsMinus(5, 0.0),
        n{auxUlti::getNormVectorComp(element, edge, 1), auxUlti::getNormVectorComp(element, edge, 2)};
    //double a(0.0), b(0.0);
    //std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

    //Compute plus values (TPlus cung duoc tinh khi excute ham nay)
    BCSupportFncs::auxilaryBCs::getUPlus(edge,nG,UPlus);

    //Compute mean values
    BCSupportFncs::auxilaryBCs::calcUMean(element,UMean);

    //Decode plus pri vars from UPlus (khong can tinh lai TPlus)
    BCSupportFncs::decompseU(priVarsPlus,UPlus,0,0,0,false);

    //Decode mean pri vars from UMean (khong can tinh lai TMean)
    BCSupportFncs::decompseU(priVarsMean,UMean,0,0,0,false);

    //Cap nhat TPlus va pPlus vao priVarsPlus
    priVarsPlus[4] = surfaceFields::T[edge][nG];
    priVarsPlus[3] = math::CalcP(priVarsPlus[4], priVarsPlus[0]);

    //Cap nhat TMean va pMean vao priVarsMean
    priVarsMean[4] = math::calcMeanPriVar(element,1); //==> ham nay tinh TMean cua cac diem Gauss volume cua cell
    priVarsMean[3] = math::CalcP(priVarsMean[4], priVarsMean[0]);

    //Check inflow/outflow
    bool inflow(BCSupportFncs::checkInflow(priVarsMean[1], priVarsMean[2],n[0],n[1]));

    //Correct priVarsMinus
    BCSupportFncs::correctPriVars(edge, edgeGrp, priVarsMinus, priVarsPlus, priVarsMean, n, inflow);

    //Save TM
    surfaceFields::T[edge][nG+mathVar::nGauss+1]=priVarsMinus[4];

    //Reconstruct advective U
    BCSupportFncs::reconstructConvectiveU(UMinus,priVarsMinus);

    //Nhan mu truoc khi return
    double muP(math::CalcVisCoef(priVarsPlus[4])), muM(math::CalcVisCoef(priVarsMinus[4]));
    for (int i = 0; i < 4; i++)
    {
        Fluxes[i][0] = UPlus[i]*muP;
        Fluxes[i][1] = UMinus[i]*muM;
    }

    //Save U+ and U- at boundary to arrays
    auxUlti::saveUAtBCToSurfaceFields(edge,nG,UPlus,UMinus);
    /* UMinus va Plus o day trong tr.hop mass diffusion ON la Uconvective.
     * Utotal se duoc correct tu Uconvective o phan apply BC cua pt NSF*/
}

void NSFEqBCsForNormalBoundary(std::vector<double> &Fluxes, int element, int edge, int edgeGrp, int nG)
{
    std::vector<double> UMinus(4, 0.0),
        UPlus(4, 0.0),
        UMean(4, 0.0),

        dUXPlus(4, 0.0),
        dUYPlus(4, 0.0),
        dpriXPlus(4, 0.0),
        dpriYPlus(4, 0.0),

        dUXMean(4, 0.0),
        dUYMean(4, 0.0),
        dpriXMean(4, 0.0),
        dpriYMean(4, 0.0),

        dUXMinus(4, 0.0),
        dUYMinus(4, 0.0),
        dpriXMinus(4, 0.0),
        dpriYMinus(4, 0.0),

        norm(2, 0.0);
    double a(0.0), b(0.0), TPlus(surfaceFields::T[edge][nG]),TMinus(surfaceFields::T[edge][nG+mathVar::nGauss+1]), TMean(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
    std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
    norm[0] = nx;
    norm[1] = ny;

    //Tinh U+/-
    //Khi co mass diffusion, U+[3] va U-[3] chi la convective energy, vi vay can correct U+/-
    for (int i = 0; i < 4; i++)
    {
        std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
    }
    //Tinh UMean
    BCSupportFncs::auxilaryBCs::calcUMean(element,UMean);

    //Tinh TMean
    TMean = math::calcMeanPriVar(element,1);

    if (flowProperties::viscous)
    {
        //Tinh dUX/dUY plus
        BCSupportFncs::NSFEqBCs::calcdUMean(edge,element,dUXMean,dUYMean);
        BCSupportFncs::NSFEqBCs::getdUPlus(edge,nG,dUXPlus,dUYPlus);

        //Decompose dUX/Y plus/mean
        BCSupportFncs::decompsedU(dpriXPlus,UPlus,dUXPlus,TPlus);
        BCSupportFncs::decompsedU(dpriYPlus,UPlus,dUYPlus,TPlus);

        BCSupportFncs::decompsedU(dpriXMean,UMean,dUXMean,TMean);
        BCSupportFncs::decompsedU(dpriYMean,UMean,dUYMean,TMean);

        //Check inflow
        bool inflow(BCSupportFncs::checkInflow(UMean[1]/UMean[0], UMean[2]/UMean[0],norm[0],norm[1]));

        //Correct grad(priVars)
        BCSupportFncs::correctPriVarsGrad(edgeGrp,dpriXMinus,dpriYMinus,dpriXPlus,dpriYPlus,UMinus,TMinus,norm,inflow);

        //Reconstruct grad(U)
        BCSupportFncs::reconstructdU(dUXMinus,dpriXMinus,UMinus,TMinus);
        BCSupportFncs::reconstructdU(dUYMinus,dpriYMinus,UMinus,TMinus);
        //std::cout<<"dPriX "<<dpriXMinus[0]<<" "<<dpriXMinus[1]<<" "<<dpriXMinus[2]<<" "<<dpriXMinus[3]<<" "<<dpriXMinus[4]<<"\n";
        //std::cout<<"dPriY "<<dpriYMinus[0]<<" "<<dpriYMinus[1]<<" "<<dpriYMinus[2]<<" "<<dpriYMinus[3]<<" "<<dpriYMinus[4]<<"\n";

        //Neu mass diffusion on, correct U+/-[3] thanh total energy
        if (flowProperties::massDiffusion)
        {
            UPlus[3]=math::pointValue(element,a,b,4,2);
            UMinus[3]=BCSupportFncs::NSFEqBCs::calcTotalEnergy(UMinus,TMinus,dpriXMinus[0],dpriYMinus[0]);
        }
    }

    BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge, Fluxes, TPlus, TMinus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
}

/**
 * @brief Function applies boundary conditions for Rho (only when mass diffusion = ON).
 *
 * If edge has type "matched":
 * - Calculate Gauss points coordinates at side '-' then calculate /f$\rho_-/f$ as normal.
 * For another type, apply zeroGradient BC.
 *
 * @param element: element Id.
 * @param edge: edge Id.
 * @param nG: Gauss point Id.
 * @return
 */
std::tuple<double, double> rhoBCsImplement(int edge, int nG)
{
    int edgeGrp(auxUlti::getGrpOfEdge(edge));
    int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]), pType(bcValues::pBcType[edgeGrp - 1]);
    double rhoP(0.0), rhoM(0.0);

    //Tinh rhoP
    rhoP = surfaceFields::rho[edge][nG];

    //Tinh rhoM
    if (UType == BCVars::generalBCId::matched && TType == BCVars::generalBCId::matched && pType == BCVars::generalBCId::matched)
    {
        //int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        rhoM = surfaceFields::rho[edge][mathVar::nGauss+nG+1];
    }
    else
    {
        rhoM = rhoP;
    }
    return std::make_tuple(rhoP,rhoM);
}
