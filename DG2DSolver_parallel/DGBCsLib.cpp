#include "DGBCsLib.h"
#include <vector>
#include "DGMath.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <tuple>
#include "DGMessagesLib.h"
#include <math.h>

//Debug
#include <iostream>

/*
 * REFERENCE: cac dieu kien bien trong file nay tham khao theo tai lieu "A Guide to the Implementation of
 * Boundary Conditions in Compact High-Order Methods for Compressible Aerodynamics" - Mengaldo et al
 *
 * NOTES:
- All boundary functions must return numerical fluxes (rho, rhou, rhov, rhoE fluxes) at surface Gauss points!
- Method index: 1: weak weak Riemann

Boundary conditions compatibility
        Boundary conditions compatibility
        |U					|T					|p					|
        +-------------------+-------------------+-------------------+
        |1. inFlow			|1. inFlow			|1. inFlow			|
        |	Value u v w		|	Value T			|	Value p			|
        +-------------------+-------------------+-------------------+
        |2. noSlip			|2. WallIsothermal	|2. zeroGradient	|
        |					|	Value T			|					|
        +-------------------+-------------------+-------------------+
        |2. noSlip			|3. WallAdiabatic	|2. zeroGradient	|
        +-------------------+-------------------+-------------------+
        |7.	symmetry		|7. symmetry		|7. symmetry		|
        +-------------------+-------------------+-------------------+
        |4. outFlow			|4. outFlow			|4. outFlow			|
        |	Value u v w		|	Value T			|	Value p			|
        +-------------------+-------------------+-------------------+
        |5.	slip    		|6. temperatureJump	|2. zeroGradient    |
        |   v_wall u v w    |   T_wall T        |                   |
        +-------------------+-------------------+-------------------+
        U:
        + 3:
        movingWall
        v_wall        u v w
*/

std::vector<std::vector<double>> NSFEqBCsImplement(int element, int edge, int nG)
{
    /*Fluxes array has the following form:
    - column 0: advective fluxes
    - column 1: diffusive fluxes*/
    std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
    int edgeGrp(auxUlti::getGrpOfEdge(edge));
    int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]), pType(bcValues::pBcType[edgeGrp - 1]);
    if (UType == 1 && TType == 1 && pType == 1)
    {
        Fluxes = NSFEqBCs::patch::inFlow(element, edge, edgeGrp, nG);
    }
    else if (UType == 4 && TType == 4 && pType == 4)
    {
        Fluxes = NSFEqBCs::patch::outFlow(element, edge, edgeGrp, nG);
    }
    else if ((UType == 2 || UType == 3) && TType == 2 && pType == 2)
    {
        Fluxes = NSFEqBCs::wall::wallIsoThermal(element, edge, edgeGrp, nG);
    }
    else if ((UType == 2 || UType == 3) && TType == 3 && pType == 2)
    {
        Fluxes = NSFEqBCs::wall::wallAdiabatic(element, edge, edgeGrp, nG);
    }
    else if (UType == 7 && TType == 7 && pType == 7)
    {
        Fluxes = NSFEqBCs::Symmetry(element, edge, nG);
    }
    else if (UType == 10 && TType == 10 && pType == 10)
    {
        Fluxes = NSFEqBCs::matched(element, edge, nG);
    }
    else if (UType == 5 && TType == 6 && pType == 2)
    {
        Fluxes = NSFEqBCs::wall::wall_MaxwellSmoluchowski(element, edge, nG);
    }
    else
    {
        std::string errorStr = message::BcCompatibleError(edgeGrp);
        message::writeLog(systemVar::pwd, systemVar::caseName, errorStr);
    }
    return Fluxes;
}

std::vector<std::vector<double>> auxEqBCsImplement(int element, int edge, int nG)
{
    //Bien phu duoc dat la mu*divU --> can tra ve gia tri mu*U sau khi chay ham nay
    std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
    int edgeGrp(auxUlti::getGrpOfEdge(edge));
    int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]), pType(bcValues::pBcType[edgeGrp - 1]);

    if (UType == 1 && TType == 1 && pType == 1)
    {
        Fluxes = auxilaryBCs::patch::inFlow(element, edge, edgeGrp, nG);
    }
    else if (UType == 4 && TType == 4 && pType == 4)
    {
        Fluxes = auxilaryBCs::patch::outFlow(element, edge, edgeGrp, nG);
    }
    else if ((UType == 2 || UType == 3) && TType == 2 && pType == 2)
    {
        Fluxes = auxilaryBCs::wall::wallIsoThermal(element, edge, edgeGrp, nG);
    }
    else if ((UType == 2 || UType == 3) && TType == 3 && pType == 2)
    {
        Fluxes = auxilaryBCs::wall::wallAdiabatic(element, edge, edgeGrp, nG);
    }
    else if (UType == 7 && TType == 7 && pType == 7)
    {
        Fluxes = auxilaryBCs::Symmetry(element, edge, nG);
    }
    else if (UType == 10 && TType == 10 && pType == 10)
    {
        Fluxes = auxilaryBCs::matched(element, edge, nG);
    }
    else if (UType == 5 && TType == 6 && pType == 2)
    {
        Fluxes = auxilaryBCs::wall::wall_MaxwellSmoluchowski(element, edge, edgeGrp, nG);
    }
    else
    {
        std::string errorStr = message::BcCompatibleError(edgeGrp);
        message::writeLog(systemVar::pwd, systemVar::caseName, errorStr);
    }

    return Fluxes;
}

//Implement bondary condition of Rho (use when massDiffusion is on)
//Method weakRiemann is used
std::tuple<double, double> rhoBCsImplement(int element, int edge, int nG)
{
    int edgeGrp(auxUlti::getGrpOfEdge(edge));
    int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]), pType(bcValues::pBcType[edgeGrp - 1]);
    double rhoP(0.0), rhoM(0.0), a(0.0), b(0.0);
    std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

    //Tinh rhoP
    rhoP = math::pointValue(element, a, b, 1, 2);

    //Tinh rhoM
    if (UType == 1 && TType == 1 && pType == 1)
    {
        rhoM = (bcValues::pBCFixed[edgeGrp - 1] / (material::R*bcValues::TBCFixed[edgeGrp - 1]));
    }
    else if ((UType == 4 && TType == 4 && pType == 4)
             || ((UType == 2 || UType == 3) && (TType == 2 || TType == 3) && pType == 2)
             || (UType == 7 && TType == 7 && pType == 7)
             || (UType == 5 && TType == 6 && pType == 2) //slip & temperature jump conditions
             )
    {
        rhoM = rhoP;
    }
    else if (UType == 10 && TType == 10 && pType == 10)
    {
        /* Neu la bien matched, tinh lai gia tri a, b cua phia ben -
         * sau do tinh rhoM nhu binh thuong
        */
        int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::tie(a,b)=auxUlti::functionsOfParallelComputing::getGaussPointCoorsOfNeighborCell(loc,nG);
        math::basisFc(a, b, 3);
        for (int iorder = 1; iorder <= mathVar::orderElem; iorder++)
        {
            rhoM += parallelBuffer::rho[loc][iorder]*mathVar::B[iorder]*parallelBuffer::theta2[loc]*parallelBuffer::theta1[loc];
        }
        rhoM += parallelBuffer::rho[loc][0];
    }
    else
    {
        std::string errorStr = message::BcCompatibleError(edgeGrp);
        message::writeLog(systemVar::pwd, systemVar::caseName, errorStr);
    }
    return std::make_tuple(rhoP,rhoM);
}

namespace BCSupportFncs
{
    bool checkInflow(double u, double v, double nx, double ny)
    {
        bool inflow(true);
        double normUMag(0.0);
        std::vector < double > U(2, 0.0);
        std::vector < double > normVector(2, 0.0);

        U[0] = u;
        U[1] = v;
        normVector[0] = nx;
        normVector[1] = ny;

        normUMag = math::vectorDotProduct(U, normVector);

        if ((normUMag >= 0) || fabs(normUMag) <=0.00001)
        {
            inflow = false;
        }

        return inflow;
    }

    namespace auxilaryBCs {
        void calcUPlus(int element, int edge, int nG, std::vector<double> &UPlus)
        {
            double a(0.0), b(0.0), TP(0.0);
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            UPlus=math::pointUVars(element, a, b);

            if (flowProperties::massDiffusion)
            {
                /* Neu mass diffusion on: tinh T bang ham CalcTFromConsvVar_massDiff,
                 * Modify lai gia tri rhoE = rhoE_convective
                */
                double rhoXP(0.0), rhoYP(0.0);
                if (systemVar::auxVariables==1)
                {
                    rhoXP=math::pointAuxValue(element,a,b,1,1);
                    rhoYP=math::pointAuxValue(element,a,b,1,2);
                }
                else if (systemVar::auxVariables==2)
                {
                    rhoXP=math::BR2Fncs::pointAuxValue_sur(edge,element,a,b,1,1);
                    rhoYP=math::BR2Fncs::pointAuxValue_sur(edge,element,a,b,1,2);
                }

                TP=math::CalcTFromConsvVar_massDiff(UPlus[0],UPlus[1],UPlus[2],UPlus[3],rhoXP,rhoYP,surfaceFields::T[edge][nG]);

                UPlus[3]=UPlus[0]*material::Cv*TP+0.5*(UPlus[1]*UPlus[1]+UPlus[2]*UPlus[2])/UPlus[0];
            }
            else {
                TP=math::CalcTFromConsvVar(UPlus[0],UPlus[1],UPlus[2],UPlus[3]);
            }
            surfaceFields::T[edge][nG]=TP;
        }
    }

    namespace NSFEqBCs {
    void calcdUPlus(int edge, int element, double a, double b, std::vector<double> &dUXPlus, std::vector<double> &dUYPlus)
    {
        //Compute dU+
        /*
        for (int i = 0; i < 4; i++)
        {
            dUXPlus[i] = math::pointAuxValue(element, a, b, i + 1, 1);
            dUYPlus[i] = math::pointAuxValue(element, a, b, i + 1, 2);
        }
        */
        dUXPlus=math::pointSVars(edge,element,a,b,1,2);
        dUYPlus=math::pointSVars(edge,element,a,b,2,2);
    }

    std::tuple<double, double, double, double> calcTotalVelocity(int BCType, double rhoP, double rhoM, double uP, double uM, double vP, double vM, double mudRhoXP, double mudRhoXM, double mudRhoYP, double mudRhoYM)
    {
        double umP(0.0), vmP(0.0), umM(0.0), vmM(0.0);
        umP=math::massDiffusionFncs::calcTotalVelocity(rhoP,uP,mudRhoXP);
        vmP=math::massDiffusionFncs::calcTotalVelocity(rhoP,vP,mudRhoYP);

        if (BCType==1 && auxUlti::checkTimeVaryingBCAvailable()) //type wall with Maxwell-Smoluchowski BC
        {
            umM=-uM;
            vmM=-vM;
        }
        else if (BCType==1 && !auxUlti::checkTimeVaryingBCAvailable())
        {
            umM=-umP;
            vmM=-vmP;
        }
        else
        {
            umM=math::massDiffusionFncs::calcTotalVelocity(rhoM,uM,mudRhoXM);
            vmM=math::massDiffusionFncs::calcTotalVelocity(rhoM,vM,mudRhoYM);
        }

        return std::make_tuple(umP,umM,vmP,vmM);
    }

    std::vector<std::vector<double>> NSFEqFluxes(int edge, int BCType, double TPlus, double TMinus, std::vector<double> &UPlus, std::vector<double> &UMinus, std::vector<double> &dUXPlus, std::vector<double> &dUXMinus, std::vector<double> &dUYPlus, std::vector<double> &dUYMinus, std::vector<double> &normVector)
    {
        /*Fluxes array has the following form:
        - column 0: advective fluxes
        - column 1: diffusive fluxes*/
        std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));

        /*StressHeat matrix has form:
        [tauXx		tauXy		Qx]
        [tauYx		tauYy		Qy]
        */
        std::vector<std::vector<double>> StressHeatP(2, std::vector<double>(3, 0.0));
        std::vector<std::vector<double>> StressHeatM(2, std::vector<double>(3, 0.0));

        double rhoPlus(UPlus[0]), rhouPlus(UPlus[1]), rhovPlus(UPlus[2]), rhoEPlus(UPlus[3]),
            rhoMinus(UMinus[0]), rhouMinus(UMinus[1]), rhovMinus(UMinus[2]), rhoEMinus(UMinus[3]),
            nx(normVector[0]), ny(normVector[1]),
                umPlus(0.0), umMinus(0.0), vmPlus(0.0), vmMinus(0.0);

        double
            uPlus(rhouPlus / rhoPlus),
            uMinus(rhouMinus / rhoMinus),

            vPlus(rhovPlus / rhoPlus),
            vMinus(rhovMinus / rhoMinus),

            totalEPlus(rhoEPlus / rhoPlus),
            totalEMinus(rhoEMinus / rhoMinus),

            pPlus(0.0),
            pMinus(0.0);

        double
            termX1P(0.0), termX1M(0.0),  //(rho*u)					or 0
            termX2P(0.0), termX2M(0.0),  //(rho*u^2 + p)			or tauxx
            termX3P(0.0), termX3M(0.0),  //(rho*u*v)				or tauxy
            termX4P(0.0), termX4M(0.0),  //(rho*totalE + p)*u		or tauxx*u + tauxy*v + Qx

            termY1P(0.0), termY1M(0.0),  //(rho*v)					or 0
            termY2P(0.0), termY2M(0.0),  //(rho*u*v)				or tauxy
            termY3P(0.0), termY3M(0.0),  //(rho*v^2 + p)			or tauyy
            termY4P(0.0), termY4M(0.0);  //(rho*totalE + p)*v		or tauxy*u + tauyy*v + Qy

        double C(0.0), Beta(0.0),
            uMagP(0.0),
            uMagM(0.0),
            aP(0.0),
            aM(0.0);

        /*INVISCID TERMS*/
        /*- Calculate total velocity components u_m*/
        if (flowProperties::massDiffusion)
        {
            std::tie(umPlus,umMinus,vmPlus,vmMinus)=BCSupportFncs::NSFEqBCs::calcTotalVelocity(BCType,rhoPlus,rhoMinus,uPlus,uMinus,vPlus,vMinus,dUXPlus[0],dUXMinus[0],dUYPlus[0],dUYMinus[0]);
        }
        else {
            umPlus=uPlus;
            umMinus=uMinus;
            vmPlus=vPlus;
            vmMinus=vMinus;
        }

        /*calculate velocity magnitude*/
        uMagP = sqrt(pow(umPlus, 2) + pow(vmPlus, 2));
        uMagM = sqrt(pow(umMinus, 2) + pow(vmMinus, 2));

        /*calculate T and P*/
        pPlus = math::CalcP(TPlus, rhoPlus);
        pMinus = math::CalcP(TMinus, rhoMinus);

        /*calculate speed of sound*/
        aP = math::CalcSpeedOfSound(TPlus);
        aM = math::CalcSpeedOfSound(TMinus);

        if (auxUlti::getBCType(edge) == 0)
        {
            C = LxFConst[edge];
        }
        else
        {
            C = math::numericalFluxes::constantC(uMagP, uMagM, aP, aM);
        }

        /*Calculate inviscid terms*/
        std::tie(termX1P, termX2P, termX3P, termX4P) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoPlus, uPlus, umPlus, vPlus, vmPlus, totalEPlus, pPlus, 1);
        std::tie(termY1P, termY2P, termY3P, termY4P) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoPlus, uPlus, umPlus, vPlus, vmPlus, totalEPlus, pPlus, 2);

        std::tie(termX1M, termX2M, termX3M, termX4M) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoMinus, uMinus, umMinus, vMinus, vmMinus, totalEMinus, pMinus, 1);
        std::tie(termY1M, termY2M, termY3M, termY4M) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoMinus, uMinus, umMinus, vMinus, vmMinus, totalEMinus, pMinus, 2);

        /*Calculate fluxes*/
        Fluxes[0][0] = math::numericalFluxes::advectiveFlux(termX1P, termX1M, rhoPlus, rhoMinus, C, nx, true) + math::numericalFluxes::advectiveFlux(termY1P, termY1M, rhoPlus, rhoMinus, C, ny, true);
        Fluxes[1][0] = math::numericalFluxes::advectiveFlux(termX2P, termX2M, rhouPlus, rhouMinus, C, nx, false) + math::numericalFluxes::advectiveFlux(termY2P, termY2M, rhouPlus, rhouMinus, C, ny, false);
        Fluxes[2][0] = math::numericalFluxes::advectiveFlux(termX3P, termX3M, rhovPlus, rhovMinus, C, nx, false) + math::numericalFluxes::advectiveFlux(termY3P, termY3M, rhovPlus, rhovMinus, C, ny, false);
        Fluxes[3][0] = math::numericalFluxes::advectiveFlux(termX4P, termX4M, rhoEPlus, rhoEMinus, C, nx, false) + math::numericalFluxes::advectiveFlux(termY4P, termY4M, rhoEPlus, rhoEMinus, C, ny, false);

        /*VISCOUS TERMS*/
        /*calculate viscous terms*/
        StressHeatP = math::viscousTerms::calcStressTensorAndHeatFlux(UPlus, dUXPlus, dUYPlus, TPlus);
        StressHeatM = math::viscousTerms::calcStressTensorAndHeatFlux(UMinus, dUXMinus, dUYMinus, TMinus);

        if (BCType==1)
        {
            std::tie(termX1P, termX2P, termX3P, termX4P) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uPlus, vPlus, 0.0, 1);
            std::tie(termY1P, termY2P, termY3P, termY4P) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uPlus, vPlus, 0.0, 2);

            std::tie(termX1M, termX2M, termX3M, termX4M) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uMinus, vMinus, 0.0, 1);
            std::tie(termY1M, termY2M, termY3M, termY4M) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uMinus, vMinus, 0.0, 2);
        }
        else
        {
            std::tie(termX1P, termX2P, termX3P, termX4P) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uPlus, vPlus, dUXPlus[0]/rhoPlus, 1);
            std::tie(termY1P, termY2P, termY3P, termY4P) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uPlus, vPlus, dUYPlus[0]/rhoPlus, 2);

            std::tie(termX1M, termX2M, termX3M, termX4M) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uMinus, vMinus, dUXMinus[0]/rhoMinus, 1);
            std::tie(termY1M, termY2M, termY3M, termY4M) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uMinus, vMinus, dUYMinus[0]/rhoMinus, 2);
        }

        /*Calculate fluxes*/
        Fluxes[0][1] = math::numericalFluxes::diffusiveFlux(termX1M, termX1P, rhoPlus, rhoMinus, Beta, nx) + math::numericalFluxes::diffusiveFlux(termY1M, termY1P, rhoPlus, rhoMinus, Beta, ny);
        Fluxes[1][1] = math::numericalFluxes::diffusiveFlux(termX2M, termX2P, rhouPlus, rhouMinus, Beta, nx) + math::numericalFluxes::diffusiveFlux(termY2M, termY2P, rhouPlus, rhouMinus, Beta, ny);
        Fluxes[2][1] = math::numericalFluxes::diffusiveFlux(termX3M, termX3P, rhovPlus, rhovMinus, Beta, nx) + math::numericalFluxes::diffusiveFlux(termY3M, termY3P, rhovPlus, rhovMinus, Beta, ny);
        Fluxes[3][1] = math::numericalFluxes::diffusiveFlux(termX4M, termX4P, rhoEPlus, rhoEMinus, Beta, nx) + math::numericalFluxes::diffusiveFlux(termY4M, termY4P, rhoEPlus, rhoEMinus, Beta, ny);

        return Fluxes;
    }
    }

    namespace parallel {
    double calcUMinus_total(int edge, int nG, std::vector<double> &UMinus)
    {
        double a(0.0), b(0.0), dRhoX(0.0), dRhoY(0.0), TM(0.0), rhoETotal(0.0);
        int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::tie(a,b)=auxUlti::functionsOfParallelComputing::getGaussPointCoorsOfNeighborCell(loc,nG);
        math::basisFc(a, b, 3);
        for (int iorder = 1; iorder <= mathVar::orderElem; iorder++)
        {
            UMinus[0] += parallelBuffer::rho[loc][iorder]*mathVar::B[iorder]*parallelBuffer::theta2[loc]*parallelBuffer::theta1[loc];
            UMinus[1] += parallelBuffer::rhou[loc][iorder]*mathVar::B[iorder]*parallelBuffer::theta2[loc];
            UMinus[2] += parallelBuffer::rhov[loc][iorder]*mathVar::B[iorder]*parallelBuffer::theta2[loc];
            UMinus[3] += parallelBuffer::rhoE[loc][iorder]*mathVar::B[iorder]*parallelBuffer::theta2[loc];
            dRhoX += parallelBuffer::drhoX[loc][iorder]*mathVar::B[iorder];
            dRhoY += parallelBuffer::drhoY[loc][iorder]*mathVar::B[iorder];
        }
        UMinus[0] += parallelBuffer::rho[loc][0];
        UMinus[1] += parallelBuffer::rhou[loc][0];
        UMinus[2] += parallelBuffer::rhov[loc][0];
        UMinus[3] += parallelBuffer::rhoE[loc][0];
        dRhoX += parallelBuffer::drhoX[loc][0];
        dRhoY += parallelBuffer::drhoY[loc][0];

        rhoETotal=UMinus[3]; //gia tri nay la rhoE total, tra ve de su dung cho phan giai pt NSF

        TM=math::CalcTFromConsvVar_massDiff_implicit(UMinus[0],UMinus[1],UMinus[2],UMinus[3],dRhoX,dRhoY);
        UMinus[3]=UMinus[0]*material::Cv*TM+0.5*(UMinus[1]*UMinus[1]+UMinus[2]*UMinus[2])/UMinus[0];
        surfaceFields::T[edge][nG+mathVar::nGauss+1]=TM;

        return rhoETotal;
    }

    void calcUMinus_convective(int edge, int nG, std::vector<double> &UMinus)
    {
        double a(0.0), b(0.0), TM(0.0);
        int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::tie(a,b)=auxUlti::functionsOfParallelComputing::getGaussPointCoorsOfNeighborCell(loc,nG);
        math::basisFc(a, b, 3);
        for (int iorder = 1; iorder <= mathVar::orderElem; iorder++)
        {
            UMinus[0] += parallelBuffer::rho[loc][iorder]*mathVar::B[iorder]*parallelBuffer::theta2[loc]*parallelBuffer::theta1[loc];
            UMinus[1] += parallelBuffer::rhou[loc][iorder]*mathVar::B[iorder]*parallelBuffer::theta2[loc];
            UMinus[2] += parallelBuffer::rhov[loc][iorder]*mathVar::B[iorder]*parallelBuffer::theta2[loc];
            UMinus[3] += parallelBuffer::rhoE[loc][iorder]*mathVar::B[iorder]*parallelBuffer::theta2[loc];
        }
        UMinus[0] += parallelBuffer::rho[loc][0];
        UMinus[1] += parallelBuffer::rhou[loc][0];
        UMinus[2] += parallelBuffer::rhov[loc][0];
        UMinus[3] += parallelBuffer::rhoE[loc][0];

        TM=math::CalcTFromConsvVar(UMinus[0],UMinus[1],UMinus[2],UMinus[3]);
        surfaceFields::T[edge][nG+mathVar::nGauss+1]=TM;
    }

    void calcdUMinus(int edge, int nG, std::vector<double> &dUMinus, int dir)
    {
        double a(0.0), b(0.0);
        int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::tie(a,b)=auxUlti::functionsOfParallelComputing::getGaussPointCoorsOfNeighborCell(loc,nG);
        math::basisFc(a, b, 3);
        if (systemVar::auxVariables==1)
        {
            if (dir==1)
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    dUMinus[0] += parallelBuffer::drhoX[loc][iorder]* mathVar::B[iorder];
                    dUMinus[1] += parallelBuffer::drhouX[loc][iorder]* mathVar::B[iorder];
                    dUMinus[2] += parallelBuffer::drhovX[loc][iorder]* mathVar::B[iorder];
                    dUMinus[3] += parallelBuffer::drhoEX[loc][iorder]* mathVar::B[iorder];
                }
            }
            else {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    dUMinus[0] += parallelBuffer::drhoY[loc][iorder]* mathVar::B[iorder];
                    dUMinus[1] += parallelBuffer::drhouY[loc][iorder]* mathVar::B[iorder];
                    dUMinus[2] += parallelBuffer::drhovY[loc][iorder]* mathVar::B[iorder];
                    dUMinus[3] += parallelBuffer::drhoEY[loc][iorder]* mathVar::B[iorder];
                }
            }
        }
    }
    }
}

namespace NSFEqBCs
{
    namespace wall
    {
        std::vector <std::vector<double>> wallIsoThermal(int element, int edge, int edgeGrp, int nG)
        {
            std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
            std::vector<double> UMinus(4, 0.0),
                UPlus(4, 0.0),
                dUXPlus(4, 0.0),
                dUYPlus(4, 0.0),
                dUXMinus(4, 0.0),
                dUYMinus(4, 0.0),
                norm(2, 0.0);
            double a(0.0), b(0.0), TPlus(surfaceFields::T[edge][nG]),TMinus(bcValues::TBCFixed[edgeGrp - 1]), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            norm[0] = nx;
            norm[1] = ny;

            //Compute U+
            for (int i = 0; i < 4; i++)
            {
                std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
            }

            //Modify UMinus[3] = UPlus[3] theo tai lieu cua Mengaldo
            UMinus[3] = UPlus[3];

            //Compute dU+
            BCSupportFncs::NSFEqBCs::calcdUPlus(edge,element,a,b,dUXPlus,dUYPlus);
            dUXMinus=dUXPlus;
            dUYMinus=dUYPlus;

            Fluxes = BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge, 1, TPlus, TMinus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
            return Fluxes;
        }

        std::vector <std::vector<double>> wallAdiabatic(int element, int edge, int edgeGrp, int nG)
        {
            std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
            std::vector<double> UMinus(4, 0.0),
                UPlus(4, 0.0),
                dUXPlus(4, 0.0), dUXMinus(4, 0.0),
                dUYPlus(4, 0.0), dUYMinus(4, 0.0),
                norm(2, 0.0),
                dRhoPlus(2, 0.0);
            double a(0.0), b(0.0), TPlus(surfaceFields::T[edge][nG]), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            norm[0] = nx;
            norm[1] = ny;

            //Compute U+
            for (int i = 0; i < 4; i++)
            {
                std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
            }

            //Compute dU+
            BCSupportFncs::NSFEqBCs::calcdUPlus(edge,element,a,b,dUXPlus,dUYPlus);

            //zero normal temperature gradient (sua lai theo cach tong quat)
            dUXMinus = dUXPlus;
            dUYMinus = dUYPlus;
            dUXMinus[3] = 0;
            dUYMinus[3] = 0;

            Fluxes = BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge,1,TPlus, TPlus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
            return Fluxes;
        }


        std::vector <std::vector<double>> wall_MaxwellSmoluchowski(int element, int edge, int nG)
        {
            int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
            std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
            std::vector<double> UMinus(4, 0.0),
                UPlus(4, 0.0),
                dUXPlus(4, 0.0),
                dUYPlus(4, 0.0),
                dUXMinus(4, 0.0),
                dUYMinus(4, 0.0),
                norm(2, 0.0);

            //Ham nay chi khac ham wallIsothermal o cho TMinus = SurfaceBCFields::TBc[loc][nG]
            double a(0.0), b(0.0), TPlus(0.0),TMinus(SurfaceBCFields::TBc[loc][nG]), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            norm[0] = nx;
            norm[1] = ny;

            //Compute U+
            for (int i = 0; i < 4; i++)
            {
                std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
            }
            TPlus=surfaceFields::T[edge][nG];

            //Compute dU+
            BCSupportFncs::NSFEqBCs::calcdUPlus(edge,element,a,b,dUXPlus,dUYPlus);
            dUXMinus=dUXPlus;
            dUYMinus=dUYPlus;

            Fluxes = BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge, 1, TPlus, TMinus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
            return Fluxes;
        }

        /*
        std::vector <std::vector<double>> wall_MaxwellSmoluchowski(int element, int edge, int nG)
        {
            int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
            std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
            std::vector<double> UMinus(4, 0.0),
                UPlus(4, 0.0),
                dUXPlus(4, 0.0),
                dUYPlus(4, 0.0),
                dUXMinus(4, 0.0),
                dUYMinus(4, 0.0),
                norm(2, 0.0);

            //Ham nay chi khac ham wallIsothermal o cho TMinus = SurfaceBCFields::TBc[loc][nG]
            double a(0.0), b(0.0), TPlus(SurfaceBCFields::TBc[loc][nG]),TMinus(SurfaceBCFields::TBc[loc][nG]), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            norm[0] = nx;
            norm[1] = ny;

            //Compute U+
            for (int i = 0; i < 4; i++)
            {
                std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
            }

            //Compute dU+
            BCSupportFncs::NSFEqBCs::calcdUPlus(edge,element,a,b,dUXPlus,dUYPlus);
            dUXMinus=dUXPlus;
            dUYMinus=dUYPlus;

            Fluxes = BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge, 1, TPlus, TMinus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
            return Fluxes;
        }*/
    }

    namespace patch
    {
        std::vector <std::vector<double>> inFlow(int element, int edge, int edgeGrp, int nG)
        {
            std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
            std::vector<double>
                UPlus(4, 0.0),
                UMinus(4, 0.0),
                dUXPlus(4, 0.0),
                dUYPlus(4, 0.0),
                dUXMinus(4, 0.0),
                dUYMinus(4, 0.0),
                norm(2, 0.0);
            double a(0.0), b(0.0), TPlus(0.0), TMinus(bcValues::TBCFixed[edgeGrp - 1]), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            norm[0] = nx;
            norm[1] = ny;

            //Compute U+
            for (int i = 0; i < 4; i++)
            {
                std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
            }
            TPlus=surfaceFields::T[edge][nG];

            //Compute dU+
            BCSupportFncs::NSFEqBCs::calcdUPlus(edge,element,a,b,dUXPlus,dUYPlus);
            dUXMinus=dUXPlus;
            dUYMinus=dUYPlus;

            Fluxes = BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge,2, TPlus, TMinus, UPlus, UMinus, dUXPlus, dUXPlus, dUYPlus, dUYPlus, norm);
            return Fluxes;
        }

        std::vector <std::vector<double>> outFlow(int element, int edge, int edgeGrp, int nG)
        {
            std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
            std::vector<double>
                UPlus(4, 0.0),
                UMinus(4, 0.0),
                dUXPlus(4, 0.0),
                dUYPlus(4, 0.0),
                dUXMinus(4, 0.0),
                dUYMinus(4, 0.0),
                norm(2, 0.0);
            double a(0.0), b(0.0), TPlus(0.0), TMinus(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            norm[0] = nx;
            norm[1] = ny;

            //Compute U+
            for (int i = 0; i < 4; i++)
            {
                std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
            }
            TPlus=surfaceFields::T[edge][nG];
            if(refValues::subsonic)
            {
                TMinus=bcValues::TBCFixed[edgeGrp-1];
            }
            else {
                TMinus=TPlus;
            }

            //Compute dU+
            BCSupportFncs::NSFEqBCs::calcdUPlus(edge,element,a,b,dUXPlus,dUYPlus);
            dUXMinus=dUXPlus;
            dUYMinus=dUYPlus;

            Fluxes = BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge,2, TPlus, TMinus, UPlus, UMinus, dUXPlus, dUXPlus, dUYPlus, dUYPlus, norm);
            return Fluxes;
        }
    }

    std::vector <std::vector<double>> Symmetry(int element, int edge, int nG)
    {
        std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
        std::vector<double>
            UPlus(4, 0.0),
            UMinus(4, 0.0),
            dUXPlus(4, 0.0),
            dUYPlus(4, 0.0),
            dUXMinus(4, 0.0),
            dUYMinus(4, 0.0),
            norm(2, 0.0);
        double a(0.0), b(0.0), TPlus(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
        norm[0] = nx;
        norm[1] = ny;

        //Compute U+
        for (int i = 0; i < 4; i++)
        {
            std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
        }
        TPlus=surfaceFields::T[edge][nG];

        //Compute dU+/-
        BCSupportFncs::NSFEqBCs::calcdUPlus(edge,element,a,b,dUXPlus,dUYPlus);
        for (int i = 0; i < 4; i++)
        {
            dUXMinus[i] = dUXPlus[i] - 2 * (dUXPlus[i] * nx + dUYPlus[i] * ny)*nx;
            dUYMinus[i] = dUYPlus[i] - 2 * (dUXPlus[i] * nx + dUYPlus[i] * ny)*ny;
        }

        Fluxes = BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge,3, TPlus, TPlus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
        return Fluxes;
    }

    std::vector <std::vector<double>> matched(int element, int edge, int nG)
    {
        std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
        std::vector<double>
            UPlus(4, 0.0),
            UMinus(4, 0.0),
            dUXPlus(4, 0.0),
            dUYPlus(4, 0.0),
            dUXMinus(4, 0.0),
            dUYMinus(4, 0.0),
            norm(2, 0.0);
        double a(0.0), b(0.0), TPlus(surfaceFields::T[edge][nG]), TMinus(surfaceFields::T[edge][nG+mathVar::nGauss+1]), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
        norm[0] = nx;
        norm[1] = ny;

        //Compute U+
        for (int i = 0; i < 4; i++)
        {
            std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
        }

        //Compute dU+/-
        BCSupportFncs::NSFEqBCs::calcdUPlus(edge,element,a,b,dUXPlus,dUYPlus);
        BCSupportFncs::parallel::calcdUMinus(edge,nG,dUXMinus,1);
        BCSupportFncs::parallel::calcdUMinus(edge,nG,dUYMinus,2);

        Fluxes = BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge,4, TPlus, TMinus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
        return Fluxes;
    }
}

namespace timeVaryingBCs
{
    void MaxwellSmoluchowski(int edge, int edgeGrp)
    {
        //Dieu kien bien Maxwell-Smoluchowski
        //Hien tai ham nay set cac gia tri tai cac diem Gauss tren 1 edge deu bang nhau
        int nG=0;

        int element(0), tempE(0), localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::tie(element,tempE)=auxUlti::getMasterServantOfEdge(edge);
        //Tinh cac gia tri can thiet
        double aConst(0.0), dux(0.0), duy(0.0), dvx(0.0), dvy(0.0), lambda(0.0),
                TC, uC, vC,
                rhoC(rho[element][0]),
                rhouC(rhou[element][0]),
                rhovC(rhov[element][0]),
                rhoEC(rhoE[element][0]),
                TVal(SurfaceBCFields::TBc[localEdgeId][nG]),
                uVal(SurfaceBCFields::uBc[localEdgeId][nG]),
                vVal(SurfaceBCFields::vBc[localEdgeId][nG]),
                delta(meshVar::distanceFromCentroidToBCEdge[localEdgeId]),
                muVal(0.0), nx, ny, dTx, dTy,dTn, dun, dvn;

        double TWall(bcValues::TBCFixed[edgeGrp - 1]), TJump(0.0),
        uWall(bcValues::uBCFixed[edgeGrp - 1]), vWall(bcValues::vBCFixed[edgeGrp - 1]), a, b;
        std::tie(a,b)=auxUlti::getGaussSurfCoor(edge,element,nG);
        double rhoBC(math::pointValue(element,a,b,1,2));

        nx = auxUlti::getNormVectorComp(element, edge, 1);
        ny = auxUlti::getNormVectorComp(element, edge, 2);

        //Gia tri trung binh (tai centroid)
        TC=math::CalcTFromConsvVar(rhoC,rhouC,rhovC,rhoEC);
        uC=rhouC/rhoC;
        vC=rhovC/rhoC;

        //Gia tri normal gradient
        dTn=(TVal-TC)/delta;
        dun=(uVal-uC)/delta;
        dvn=(vVal-vC)/delta;

        //Tinh cac thanh phan dao ham theo 2 phuong
        dTx=dTn*nx;
        dTy=dTn*ny;

        dux=dun*nx;
        duy=dun*ny;

        dvx=dvn*nx;
        dvy=dvn*ny;

        //1. Smoluchowski temperature jump-----------------------------------------------------------------------------
        muVal=math::CalcVisCoef(TVal);

        lambda=math::calcMeanFreePath(muVal,rhoBC,TVal);
        aConst=(2-bcValues::sigmaT)*2*material::gamma*lambda/(bcValues::sigmaT*material::Pr*(material::gamma+1));

        TJump=TWall-aConst*dTn;
        if (TJump<0)
        {
            TJump=TVal;
        }
        //Update to surfaceBCfields
        for (int i=0; i<mathVar::nGauss+1; i++)
        {
            SurfaceBCFields::TBc[localEdgeId][i]=TJump;
        }

        //2. Maxwell slip velocity--------------------------------------------------------
        //Tinh lai gia tri muVal va lambda tu TJump moi cap nhat
        muVal=math::CalcVisCoef(TJump);
        lambda=math::calcMeanFreePath(muVal,rhoBC,TJump);

        //Trasnformation tensor components
        double S11(1-nx*nx),
                S12(-nx*ny),
                S21(-nx*ny),
                S22(1-ny*ny);

        //- Term div_n(S.u)
        double div_n_u1(dux*nx+dvx*ny), div_n_u2(duy*nx+dvy*ny);
        std::vector<double> div_n_Sdotu={
            S11*div_n_u1 + S12*div_n_u2,
            S21*div_n_u1 + S22*div_n_u2
        };

        //- Term S.(n.PI_mc)
        double PImc_n1(nx*(dux-(2/3)*(dux+dvy))+ny*duy), PImc_n2(nx*dvx+ny*(dvy-(2/3)*(dux+dvy)));
        std::vector<double> SdotndotPImc={
            S11*PImc_n1 + S12*PImc_n2,
            S21*PImc_n1 + S22*PImc_n2
        };

        //- Term S.div(T)
        std::vector<double> SdotDivT={
            S11*dTx + S12*dTy,
            S21*dTx + S22*dTy
        };

        //Giai slip velocity
        double c=(2-bcValues::sigmaU)/bcValues::sigmaU;
        double
                uSlip=uWall-c*lambda*SdotndotPImc[0]/muVal-(3/4)*SdotDivT[0]/(rhoBC)-c*lambda*div_n_Sdotu[0],
                vSlip=vWall-c*lambda*SdotndotPImc[1]/muVal-(3/4)*SdotDivT[1]/(rhoBC)-c*lambda*div_n_Sdotu[1];

        //Update
        for (int i=0; i<mathVar::nGauss+1; i++)
        {
            SurfaceBCFields::uBc[localEdgeId][i]=uSlip;
            SurfaceBCFields::vBc[localEdgeId][i]=vSlip;
        }
    }

    void MaxwellSmoluchowski_explicit(int edge, int edgeGrp)
    {
        //Dieu kien bien Maxwell-Smoluchowski
        //Hien tai ham nay set cac gia tri tai cac diem Gauss tren 1 edge deu bang nhau
        int nG=0;

        int element(0), tempE(0), localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::tie(element,tempE)=auxUlti::getMasterServantOfEdge(edge);

        //Tinh cac gia tri can thiet-----------------------------------------------------
        //Lay cac gia tri da luu
        double aConst(0.0), dux(0.0), duy(0.0), dvx(0.0), dvy(0.0), lambda(0.0), dEx, dEy,
                TVal(SurfaceBCFields::TBc[localEdgeId][nG]),
                muVal(0.0), nx, ny, dTx, dTy,dTn;

        double TWall(bcValues::TBCFixed[edgeGrp - 1]), TJump(0.0),
        uWall(bcValues::uBCFixed[edgeGrp - 1]), vWall(bcValues::vBCFixed[edgeGrp - 1]), a, b;
        nx = auxUlti::getNormVectorComp(element, edge, 1);
        ny = auxUlti::getNormVectorComp(element, edge, 2);

        //Tinh gia tri tai normal projection cua center xuong BC
        a=meshVar::normProjectionOfCenterToBCEdge_standardSysCoor[localEdgeId][0];
        b=meshVar::normProjectionOfCenterToBCEdge_standardSysCoor[localEdgeId][1];
        double rhoBC(math::pointValue(element,a,b,1,2)),
                rhouBC(math::pointValue(element,a,b,2,2)),
                rhovBC(math::pointValue(element,a,b,3,2)),
                rhoEBC(math::pointValue(element,a,b,4,2)),
                dRhoBCx(math::pointAuxValue(element,a,b,1,1)),
                dRhoBCy(math::pointAuxValue(element,a,b,1,2)),
                dRhouBCx(math::pointAuxValue(element,a,b,2,1)),
                dRhouBCy(math::pointAuxValue(element,a,b,2,2)),
                dRhovBCx(math::pointAuxValue(element,a,b,3,1)),
                dRhovBCy(math::pointAuxValue(element,a,b,3,2)),
                dRhoEBCx(math::pointAuxValue(element,a,b,4,1)),
                dRhoEBCy(math::pointAuxValue(element,a,b,4,2));

        //Tinh cac thanh phan dao ham theo 2 phuong
        //dTx=math::calcDivTFromAuxVariable(dRhoEBCx,dRhoBCx,rhoEBC,rhoBC);
        //dTy=math::calcDivTFromAuxVariable(dRhoEBCy,dRhoBCy,rhoEBC,rhoBC);

        dux=math::calcRhouvEDeriv(dRhouBCx,dRhoBCx,rhouBC,rhoBC);
        duy=math::calcRhouvEDeriv(dRhouBCy,dRhoBCy,rhouBC,rhoBC);

        dvx=math::calcRhouvEDeriv(dRhovBCx,dRhoBCx,rhovBC,rhoBC);
        dvy=math::calcRhouvEDeriv(dRhovBCy,dRhoBCy,rhovBC,rhoBC);

        dEx=math::calcRhouvEDeriv(dRhoEBCx,dRhoBCx,rhoEBC,rhoBC);
        dEy=math::calcRhouvEDeriv(dRhoEBCy,dRhoBCy,rhoEBC,rhoBC);

        dTx=math::calcTDeriv(dEx,dux,dvx,rhouBC/rhoBC,rhovBC/rhoBC);
        dTy=math::calcTDeriv(dEy,duy,dvy,rhouBC/rhoBC,rhovBC/rhoBC);

        //1. Smoluchowski temperature jump-----------------------------------------------------------------------------
        dTn=dTx*nx+dTy*ny;
        muVal=math::CalcVisCoef(TVal);

        lambda=math::calcMeanFreePath(muVal,rhoBC,TVal)/muVal;
        aConst=(2-bcValues::sigmaT)*2*material::gamma*lambda/(bcValues::sigmaT*material::Pr*(material::gamma+1));

        TJump=TWall-aConst*dTn;
        if (TJump<0)
        {
            TJump=TVal;
        }
        //Update to surfaceBCfields
        for (int i=0; i<mathVar::nGauss+1; i++)
        {
            SurfaceBCFields::TBc[localEdgeId][i]=TJump;
        }

        //2. Maxwell slip velocity--------------------------------------------------------
        //Tinh lai gia tri muVal va lambda tu TJump moi cap nhat
        muVal=math::CalcVisCoef(TJump);
        lambda=math::calcMeanFreePath(muVal,rhoBC,TJump)/muVal;

        //Transformation tensor components
        double S11(1-nx*nx),
                S12(-nx*ny),
                S21(-nx*ny),
                S22(1-ny*ny);

        //- Term div_n(S.u)
        double div_n_u1(dux*nx+dvx*ny), div_n_u2(duy*nx+dvy*ny);
        std::vector<double> div_n_Sdotu={
            S11*div_n_u1 + S12*div_n_u2,
            S21*div_n_u1 + S22*div_n_u2
        };

        //- Term S.(n.PI_mc)
        double PImc_n1(nx*(dux-(2/3)*(dux+dvy))+ny*duy), PImc_n2(nx*dvx+ny*(dvy-(2/3)*(dux+dvy)));
        std::vector<double> SdotndotPImc={
            S11*PImc_n1 + S12*PImc_n2,
            S21*PImc_n1 + S22*PImc_n2
        };

        //- Term S.div(T)
        std::vector<double> SdotDivT={
            S11*dTx + S12*dTy,
            S21*dTx + S22*dTy
        };

        //Giai slip velocity
        double c=(2-bcValues::sigmaU)/bcValues::sigmaU;
        double
                uSlip=uWall-c*lambda*SdotndotPImc[0]/muVal-(3/4)*SdotDivT[0]/(rhoBC*TJump)-c*lambda*div_n_Sdotu[0],
                vSlip=vWall-c*lambda*SdotndotPImc[1]/muVal-(3/4)*SdotDivT[1]/(rhoBC*TJump)-c*lambda*div_n_Sdotu[1];

        //Update
        for (int i=0; i<mathVar::nGauss+1; i++)
        {
            SurfaceBCFields::uBc[localEdgeId][i]=uSlip;
            SurfaceBCFields::vBc[localEdgeId][i]=vSlip;
        }
    }

    void MaxwellSmoluchowski_implicit2ndOrder(int edge, int edgeGrp)
    {
        //Dieu kien bien Maxwell-Smoluchowski
        //Hien tai ham nay set cac gia tri tai cac diem Gauss tren 1 edge deu bang nhau
        int nG=0;

        int element(0), tempE(0), localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
        std::tie(element,tempE)=auxUlti::getMasterServantOfEdge(edge);
        //Tinh cac gia tri can thiet
        double lambda(0.0),
                TC, uC, vC,
                rhoC(rho[element][0]),
                rhouC(rhou[element][0]),
                rhovC(rhov[element][0]),
                rhoEC(rhoE[element][0]),
                TVal(SurfaceBCFields::TBc[localEdgeId][nG]),
                delta(meshVar::distanceFromCentroidToBCEdge[localEdgeId]),
                muVal(0.0), nx, ny, dTx, dTy,dTn;

        double TWall(bcValues::TBCFixed[edgeGrp - 1]), TJump(0.0),
        uWall(bcValues::uBCFixed[edgeGrp - 1]), vWall(bcValues::vBCFixed[edgeGrp - 1]), a, b;
        std::tie(a,b)=auxUlti::getGaussSurfCoor(edge,element,nG);
        double rhoBC(math::pointValue(element,a,b,1,2));

        nx = auxUlti::getNormVectorComp(element, edge, 1);
        ny = auxUlti::getNormVectorComp(element, edge, 2);

        //Gia tri trung binh (tai centroid)
        TC=math::CalcTFromConsvVar(rhoC,rhouC,rhovC,rhoEC);
        uC=rhouC/rhoC;
        vC=rhovC/rhoC;

        //1. Smoluchowski temperature jump-----------------------------------------------------------------------------
        double AT(material::As*pow(3.1416/(2*material::R),0.5)/rhoBC),
                BT(2*material::gamma*(2-bcValues::sigmaT)/(bcValues::sigmaT*material::Pr*(material::gamma+1)));
        double ATT(delta+AT*BT),
                BTT(material::Ts*delta-TWall*delta-AT*BT*TC),
                CTT(-TWall*material::Ts*delta), T1, T2;

        bool realSolution(false);

        std::tie(realSolution, T1, T2) = math::solvePolynomialsEq::polynominal2d(ATT,BTT,CTT);
        if (realSolution)
        {
            if (T1>0)
            {
                TJump=T1;
            }
            else
            {
                TJump=T2;
            }
        }
        else
        {
            TJump=TVal;
        }
        //Update to surfaceBCfields
        for (int i=0; i<mathVar::nGauss+1; i++)
        {
            SurfaceBCFields::TBc[localEdgeId][i]=TJump;
        }

        //Gia tri normal gradient
        dTn=(TJump-TC)/delta;

        //Tinh cac thanh phan dao ham theo 2 phuong
        dTx=dTn*nx;
        dTy=dTn*ny;

        muVal=math::CalcVisCoef(TJump);
        lambda=math::calcMeanFreePath(muVal,rhoBC,TJump);

        //2. Maxwell slip velocity--------------------------------------------------------
        uC=uC/delta;
        vC=vC/delta;

        //Transformation tensor components
        double S11(1-nx*nx),
                S12(-nx*ny),
                S21(-nx*ny),
                S22(1-ny*ny);

        double A1=(2-bcValues::sigmaU)*lambda/bcValues::sigmaU;
        //- Term (-3/4)*mu*S.div(T)/(rho*T)
        std::vector<double> A2={
            -3*muVal*(S11*dTx + S12*dTy)/(4*rhoBC*TJump),
            -3*muVal*(S21*dTx + S22*dTy)/(4*rhoBC*TJump)
        };

        //He so B
        double B1(-(uC*nx*nx+vC*nx*ny)),
                B2(-(uC*nx*ny+vC*ny*ny)),
                B3(nx*nx/delta),
                B4(nx*ny/delta),
                B5(ny*ny/delta);

        //He so C
        double C1(S11*B3+S12*B4),
                C2(S11*B4+S12*B5),
                C3(S11*B1+S12*B2),
                C4(S21*B3+S22*B4),
                C5(S21*B4+S22*B5),
                C6(S21*B1+S22*B2);

        //He so D
        double D1(nx*nx/(3*delta)+ny*ny/delta),
                D2(-2*nx*ny/(3*delta)),
                D3(-uC*ny*ny-uC*nx*nx/3+2*nx*ny*vC/3),
                D4(-2*nx*ny/(3*delta)),
                D5(ny*ny/(3*delta)+nx*nx/delta),
                D6(-vC*nx*nx-vC*ny*ny/3+2*nx*ny*uC/3);

        //He so E
        double E1(S11*D1+S12*D4),
                E2(S11*D2+S12*D5),
                E3(S11*D3+S12*D6),
                E4(S21*D1+S22*D4),
                E5(S21*D2+S22*D5),
                E6(S21*D3+S22*D6);

        //He so F
        double F1(1+A1*C1+A1*E1/muVal),
                F2(A1*C2+A1*E2/muVal),
                F3(uWall-(A1/muVal)*E3-A1*C3+A2[0]),
                F4(A1*C4+A1*E4/muVal),
                F5(1+A1*C5+A1*E5/muVal),
                F6(vWall-(A1/muVal)*E6-A1*C6+A2[1]);

        //Giai uSlip, vSlip
        double uSlip, vSlip;
        vSlip=(F6-F4*F3/F1)/(F5-F4*F2/F1);
        uSlip=(F3-F2*vSlip)/F1;

        //Update
        for (int i=0; i<mathVar::nGauss+1; i++)
        {
            SurfaceBCFields::uBc[localEdgeId][i]=uSlip;
            SurfaceBCFields::vBc[localEdgeId][i]=vSlip;
        }
    }
}

namespace auxilaryBCs
{
    namespace wall
    {
        std::vector <std::vector<double>> wallIsoThermal(int element, int edge, int edgeGrp, int nG)
        {
            //columns 0, 1 are plus, minus values
            std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
            std::vector<double>
                UPlus(4, 0.0),
                UMinus(4, 0.0);
            double a(0.0), b(0.0), muP, muM;
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

            //Compute plus values--------------------------------------------------------
            //Trong ham nay da co tinh T va luu vao  surfaceFields::T[edge][nG]
            BCSupportFncs::auxilaryBCs::calcUPlus(element,edge,nG,UPlus);

            //Compute minus values--------------------------------------------------------
            UMinus[0] = UPlus[0];

            if (bcValues::UBcType[edgeGrp - 1] == 2)
            {
                UMinus[1] = -UPlus[1];
                UMinus[2] = -UPlus[2];
                UMinus[3] = UMinus[0]*material::Cv*bcValues::TBCFixed[edgeGrp - 1]+0.5*(pow(UMinus[1],2)+pow(UMinus[2],2))/UMinus[0];
            }
            else if (bcValues::UBcType[edgeGrp - 1] == 3)
            {
                UMinus[1] = bcValues::uBCFixed[edgeGrp - 1]*UMinus[0];
                UMinus[2] = bcValues::vBCFixed[edgeGrp - 1]*UMinus[0];
                UMinus[3] = UMinus[0]*(material::Cv*bcValues::TBCFixed[edgeGrp - 1] + 0.5*(pow(bcValues::uBCFixed[edgeGrp - 1],2)+pow(bcValues::vBCFixed[edgeGrp - 1],2)));
            }
            //-----------------------------------------------------------------------------

            muP=math::CalcVisCoef(surfaceFields::T[edge][nG]);
            muM=math::CalcVisCoef(bcValues::TBCFixed[edgeGrp - 1]);
            for (int i = 0; i < 4; i++)
            {
                Fluxes[i][0] = UPlus[i]*muP;
                Fluxes[i][1] = UMinus[i]*muM;
            }

            //Save U+ and U- at boundary to arrays
            //Correct rhoEm+/-
            if (flowProperties::massDiffusion)
            {
                UPlus[3]=math::pointValue(element,a,b,4,2);
                //UMinus[3]=(UPlus[3]-UPlus[0]*material::Cv*surfaceFields::T[edge][nG])+UMinus[0]*material::Cv*bcValues::TBCFixed[edgeGrp - 1];
            }

            auxUlti::saveUAtBCToSurfaceFields(edge,nG,UPlus,UMinus);
            return Fluxes;
        }

        std::vector <std::vector<double>> wallAdiabatic(int element, int edge, int edgeGrp, int nG)
        {
            //columns 0, 1 are plus, minus values
            std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
            std::vector<double>
                UPlus(4, 0.0),
                UMinus(4, 0.0);
            double a(0.0), b(0.0), muP;
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

            //Compute plus values--------------------------------------------------------
            BCSupportFncs::auxilaryBCs::calcUPlus(element,edge,nG,UPlus);

            //Compute minus values--------------------------------------------------------
            UMinus[0] = UPlus[0];
            if (bcValues::UBcType[edgeGrp - 1] == 2)
            {
                UMinus[1] = -UPlus[1];
                UMinus[2] = -UPlus[2];
                UMinus[3] = UPlus[3];
            }
            else if (bcValues::UBcType[edgeGrp - 1] == 3)
            {
                UMinus[1] = bcValues::uBCFixed[edgeGrp - 1]*UMinus[0];
                UMinus[2] = bcValues::vBCFixed[edgeGrp - 1]*UMinus[0];
                UMinus[3] = UMinus[0]*0.5*(pow(bcValues::uBCFixed[edgeGrp - 1],2)+pow(bcValues::vBCFixed[edgeGrp - 1],2))+UPlus[3];
            }
            //-----------------------------------------------------------------------------

            //Nhan mu vao vector flux truoc khi return
            muP=math::CalcVisCoef(surfaceFields::T[edge][nG]);
            for (int i = 0; i < 4; i++)
            {
                Fluxes[i][0] = UPlus[i]*muP;
                Fluxes[i][1] = UMinus[i]*muP;
            }
            //Save U+ and U- at boundary to arrays
            //Recompute rhoEm
            if (flowProperties::massDiffusion)
            {
                UPlus[3]=math::pointValue(element,a,b,4,2);
                UMinus[3]=UPlus[3];
            }
            auxUlti::saveUAtBCToSurfaceFields(edge,nG,UPlus,UMinus);
            return Fluxes;
        }


        std::vector <std::vector<double>> wall_MaxwellSmoluchowski(int element, int edge, int edgeGrp, int nG)
        {
            /*
             * Dieu kien bien temperature jump:
             * - vi co thanh phan uslip tai be mat, rhou va rhov cua UMinus se khac 0:
             * rhou=rhoPlus*uSlip
             * rhov=rhoPlus*vSlip
             * rhoMinus = rhoPlus
            */

            //columns 0, 1 are plus, minus values
            std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
            std::vector<double>
                UPlus(4, 0.0),
                UMinus(4, 0.0);
            double a(0.0), b(0.0), muPlus, muMinus;
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            //Lay local id cua edge tren BC
            int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));

            //Compute plus values--------------------------------------------------------
            BCSupportFncs::auxilaryBCs::calcUPlus(element,edge,nG,UPlus);

            //Compute minus values--------------------------------------------------------
            UMinus[0] = UPlus[0];
            UMinus[1] = 2*SurfaceBCFields::uBc[loc][nG]*UMinus[0]-UPlus[1];
            UMinus[2] = 2*SurfaceBCFields::vBc[loc][nG]*UMinus[0]-UPlus[2];
            UMinus[3] = UMinus[0]*material::Cv*SurfaceBCFields::TBc[loc][nG]+0.5*(pow(UMinus[1],2)+pow(UMinus[2],2))/UMinus[0];
            //-----------------------------------------------------------------------------

            muPlus=math::CalcVisCoef(surfaceFields::T[edge][nG]);
            muMinus=math::CalcVisCoef(SurfaceBCFields::TBc[loc][nG]);

            //Nhan mu vao vector flux truoc khi return
            for (int i = 0; i < 4; i++)
            {
                Fluxes[i][0] = UPlus[i]*muPlus;
                Fluxes[i][1] = UMinus[i]*muMinus;
            }

            //Save U+ and U- at boundary to arrays
            //Recompute rhoEm
            if (flowProperties::massDiffusion)
            {
                UPlus[3]=math::pointValue(element,a,b,4,2);
            }
            //UMinus[3]=UPlus[3];
            auxUlti::saveUAtBCToSurfaceFields(edge,nG,UPlus,UMinus);
            return Fluxes;
        }

        /*
        std::vector <std::vector<double>> wall_MaxwellSmoluchowski(int element, int edge, int edgeGrp, int nG)
        {
            /
             * Dieu kien bien temperature jump:
             * - vi co thanh phan uslip tai be mat, rhou va rhov cua UMinus se khac 0:
             * rhou=rhoPlus*uSlip
             * rhov=rhoPlus*vSlip
             * rhoMinus = rhoPlus
            /

            //columns 0, 1 are plus, minus values
            std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
            std::vector<double>
                UPlus(4, 0.0),
                UMinus(4, 0.0);
            double a(0.0), b(0.0), muPlus, muMinus, rhoPlus;
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            //Lay local id cua edge tren BC
            int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));

            //Compute plus values--------------------------------------------------------
            rhoPlus=math::pointValue(element,a,b,1,2);

            //Compute minus values--------------------------------------------------------
            UMinus[0] = rhoPlus;
            UMinus[1] = SurfaceBCFields::uBc[loc][nG]*UMinus[0];
            UMinus[2] = SurfaceBCFields::vBc[loc][nG]*UMinus[0];
            UMinus[3] = UMinus[0]*material::Cv*SurfaceBCFields::TBc[loc][nG]+0.5*(pow(UMinus[1],2)+pow(UMinus[2],2))/UMinus[0];
            //-----------------------------------------------------------------------------
            UPlus=UMinus;
            surfaceFields::T[edge][nG]=SurfaceBCFields::TBc[loc][nG];

            muMinus=math::CalcVisCoef(SurfaceBCFields::TBc[loc][nG]);
            muPlus=muMinus;

            //Nhan mu vao vector flux truoc khi return
            for (int i = 0; i < 4; i++)
            {
                Fluxes[i][0] = UPlus[i]*muPlus;
                Fluxes[i][1] = UMinus[i]*muMinus;
            }

            //Save U+ and U- at boundary to arrays
            auxUlti::saveUAtBCToSurfaceFields(edge,nG,UPlus,UMinus);
            return Fluxes;
        }*/
    }

    namespace patch
    {
        std::vector <std::vector<double>> inFlow(int element, int edge, int edgeGrp, int nG)
        {
            std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
            std::vector<double>
                UPlus(4, 0.0),
                UMinus(4, 0.0);
            double a(0.0), b(0.0), muP, muM;
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

            //Compute plus values--------------------------------------------------------
            BCSupportFncs::auxilaryBCs::calcUPlus(element,edge,nG,UPlus);

            //Compute minus values--------------------------------------------------------
            //Apply weak Riemann infinite value
            UMinus[0] = (bcValues::pBCFixed[edgeGrp - 1] / (material::R*bcValues::TBCFixed[edgeGrp - 1]));
            UMinus[1] = UMinus[0] * bcValues::uBCFixed[edgeGrp - 1];
            UMinus[2] = UMinus[0] * bcValues::vBCFixed[edgeGrp - 1];
            UMinus[3] = UMinus[0] * (bcValues::TBCFixed[edgeGrp - 1] * material::Cv + 0.5*(pow(bcValues::uBCFixed[edgeGrp - 1], 2) + pow(bcValues::vBCFixed[edgeGrp - 1], 2)));
            //-----------------------------------------------------------------------------

            //Nhan mu vao vector flux truoc khi return
            muP=math::CalcVisCoef(surfaceFields::T[edge][nG]);
            muM=math::CalcVisCoef(bcValues::TBCFixed[edgeGrp - 1]);
            for (int i = 0; i < 4; i++)
            {
                Fluxes[i][0] = UPlus[i]*muP;
                Fluxes[i][1] = UMinus[i]*muM;
            }

            //Save U+ and U- at boundary to arrays
            //Recompute rhoEm
            if (flowProperties::massDiffusion)
            {
                UPlus[3]=math::pointValue(element,a,b,4,2);
            }
            auxUlti::saveUAtBCToSurfaceFields(edge,nG,UPlus,UMinus);
            return Fluxes;
        }

        std::vector <std::vector<double>> outFlow(int element, int edge, int edgeGrp, int nG)
        {
            /* Theo Mengaldo, khi tinh auxilary variable, q_auxBC = q_in.
             * Do do gan U- = U+ khi giai pt phu, sau do modify lai U- truoc khi
             * luu vao surfaceField
            */
            std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
            std::vector<double>
                UPlus(4, 0.0),
                UMinus(4, 0.0);
            double a(0.0), b(0.0), pInternal(0), muP, muM;
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

            //Compute plus values--------------------------------------------------------
            BCSupportFncs::auxilaryBCs::calcUPlus(element,edge,nG,UPlus);

            //Compute minus values--------------------------------------------------------
            UMinus = UPlus;

            //---------------------------------------------------------------------------
            //Nhan mu vao vector flux truoc khi return
            muP=math::CalcVisCoef(surfaceFields::T[edge][nG]);
            if (refValues::subsonic)
            {
                muM=bcValues::TBCFixed[edgeGrp - 1];
            }
            else
            {
                muM=muP;
            }

            for (int i = 0; i < 4; i++)
            {
                Fluxes[i][0] = UPlus[i]*muP;
                Fluxes[i][1] = UMinus[i]*muM;
            }

            //Modify lai U- truoc khi luu vao surfaceField
            //Apply PNR (2), R (1)
            double uPlus(UPlus[1] / UPlus[0]), vPlus(UPlus[2] / UPlus[0]);
            int implementation(2);
            switch (implementation)
            {
            case 1: //R
            {
                if (refValues::subsonic)
                {
                    UMinus[3] = bcValues::pBCFixed[edgeGrp - 1] / (material::gamma - 1) + 0.5*UPlus[0] * (pow(uPlus, 2) + pow(vPlus, 2));
                }
            }
            break;
            case 2: //PNR
            {
                pInternal = UPlus[0] * material::R*math::CalcTFromConsvVar(UPlus[0],UPlus[1],UPlus[2],UPlus[3]);
                if (refValues::subsonic)
                {
                    UMinus[3] = (2 * bcValues::pBCFixed[edgeGrp - 1] - pInternal) / (material::gamma - 1) + 0.5*UPlus[0] * (pow(uPlus, 2) + pow(vPlus, 2));
                }
            }
            break;
            default:
                break;
            }

            //Save U+ and U- at boundary to arrays
            //Recompute rhoEm
            if (flowProperties::massDiffusion)
            {
                UPlus[3]=math::pointValue(element,a,b,4,2);
                if (!refValues::subsonic)
                {
                    UMinus[3] = UPlus[3];
                }
            }
            auxUlti::saveUAtBCToSurfaceFields(edge,nG,UPlus,UMinus);
            return Fluxes;
        }
    }

    std::vector <std::vector<double>> Symmetry(int element, int edge, int nG)
    {
        std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
        std::vector<double>
            UPlus(4, 0.0),
            UMinus(4, 0.0),
            norm(2, 0.0);
        double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2)), muP;
        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
        norm[0] = nx;
        norm[1] = ny;

        //Compute plus values--------------------------------------------------------
        BCSupportFncs::auxilaryBCs::calcUPlus(element,edge,nG,UPlus);

        //Compute minus values--------------------------------------------------------
        UMinus[0] = UPlus[0];
        UMinus[1] = UPlus[1] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*nx;
        UMinus[2] = UPlus[2] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*ny;
        UMinus[3] = UPlus[3];
        //----------------------------------------------------------------------------

        //Nhan mu vao vector flux truoc khi return
        muP=math::CalcVisCoef(surfaceFields::T[edge][nG]);
        for (int i = 0; i < 4; i++)
        {
            Fluxes[i][0] = UPlus[i]*muP;
            Fluxes[i][1] = UMinus[i]*muP;
        }

        //Save U+ and U- at boundary to arrays
        //Recompute rhoEm
        if (flowProperties::massDiffusion)
        {
            UPlus[3]=math::pointValue(element,a,b,4,2);
            UMinus[3]=UPlus[3];
        }
        auxUlti::saveUAtBCToSurfaceFields(edge,nG,UPlus,UMinus);
        return Fluxes;
    }

    std::vector <std::vector<double>> matched(int element, int edge, int nG)
    {
        std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
        std::vector<double>
            UPlus(4, 0.0),
            UMinus(4, 0.0);
        double a(0.0), b(0.0), rhoETotalMinus(0.0), muP, muM;
        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

        //Compute plus values--------------------------------------------------------
        BCSupportFncs::auxilaryBCs::calcUPlus(element,edge,nG,UPlus);

        //Compute minus values--------------------------------------------------------
        if (flowProperties::massDiffusion)
        {
            rhoETotalMinus=BCSupportFncs::parallel::calcUMinus_total(edge,nG,UMinus);
        }
        else
        {
            BCSupportFncs::parallel::calcUMinus_convective(edge,nG,UMinus);
        }

        //Nhan mu vao vector flux truoc khi return
        muP=math::CalcVisCoef(surfaceFields::T[edge][nG]);
        muM=math::CalcVisCoef(surfaceFields::T[edge][nG+mathVar::nGauss+1]);
        for (int i = 0; i < 4; i++)
        {
            Fluxes[i][0] = UPlus[i]*muP;
            Fluxes[i][1] = UMinus[i]*muM;
        }

        //Save U+ and U- at boundary to arrays
        //Correct rhoEm
        if (flowProperties::massDiffusion)
        {
            UPlus[3]=math::pointValue(element,a,b,4,2);
            UMinus[3]=rhoETotalMinus;
        }
        auxUlti::saveUAtBCToSurfaceFields(edge,nG,UPlus,UMinus);
        return Fluxes;
    }
}
