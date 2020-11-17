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

std::vector<double> NSFEqBCsImplement(int element, int edge, int nG)
{
    /*Fluxes array has the following form:
    - column 0: advective fluxes
    - column 1: diffusive fluxes*/
    std::vector<double> Fluxes(4, 0.0);
    int edgeGrp(auxUlti::getGrpOfEdge(edge));
    int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]), pType(bcValues::pBcType[edgeGrp - 1]);
    if (UType == 1 && TType == 1 && pType == 1)
    {
        NSFEqBCs::patch::inFlow(Fluxes, element, edge, edgeGrp, nG);
    }
    else if (UType == 4 && TType == 4 && pType == 4)
    {
        NSFEqBCs::patch::outFlow(Fluxes, element, edge, edgeGrp, nG);
    }
    else if ((UType == 2 || UType == 3) && TType == 2 && pType == 2)
    {
        NSFEqBCs::wall::wallIsoThermal(Fluxes, element, edge, edgeGrp, nG);
    }
    else if ((UType == 2 || UType == 3) && TType == 3 && pType == 2)
    {
        NSFEqBCs::wall::wallAdiabatic(Fluxes, element, edge, nG);
    }
    else if (UType == 7 && TType == 7 && pType == 7)
    {
        NSFEqBCs::Symmetry(Fluxes, element, edge, nG);
    }
    else if (UType == 10 && TType == 10 && pType == 10)
    {
        NSFEqBCs::matched(Fluxes, element, edge, nG);
    }
    else if (UType == 5 && TType == 6 && pType == 2)
    {
        NSFEqBCs::wall::wall_NonEquilibrium(Fluxes, element, edge, nG);
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
        Fluxes = auxilaryBCs::wall::wall_NonEquilibrium(element, edge, nG);
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

        /*Ham apply dieu kien back flow tai vi tri bi back (reversed) flow*/
        void calcUReversedFlow(int edge, int nG, std::vector<double> &U)
        {
            int edgeGrp(auxUlti::getGrpOfEdge(edge));
            U[0] = (bcValues::pBCFixed[edgeGrp - 1] / (material::R*bcValues::TBCFixed[edgeGrp - 1]));
            U[1] = U[0] * bcValues::uBCFixed[edgeGrp - 1];
            U[2] = U[0] * bcValues::vBCFixed[edgeGrp - 1];
            U[3] = U[0] * (bcValues::TBCFixed[edgeGrp - 1] * material::Cv + 0.5*(pow(bcValues::uBCFixed[edgeGrp - 1], 2) + pow(bcValues::vBCFixed[edgeGrp - 1], 2)));
            surfaceFields::T[edge][nG]=bcValues::TBCFixed[edgeGrp - 1];
        }

        void calcUMean(int element, int edge, int nG, std::vector<double> &UMean)
        {
            double TP(0.0);
            auxUlti::getVectorUMeanOfCell(element, UMean);

            if (flowProperties::massDiffusion)
            {
                /* Neu mass diffusion on: tinh T bang ham CalcTFromConsvVar_massDiff,
                 * Modify lai gia tri rhoE = rhoE_convective
                */
                double rhoXP(BR1Vars::rhoX[element][0]), rhoYP(BR1Vars::rhoY[element][0]);

                TP=math::CalcTFromConsvVar_massDiff(UMean[0],UMean[1],UMean[2],UMean[3],rhoXP,rhoYP,surfaceFields::T[edge][nG]);

                UMean[3]=UMean[0]*material::Cv*TP+0.5*(UMean[1]*UMean[1]+UMean[2]*UMean[2])/UMean[0];
            }
            else {
                TP=math::CalcTFromConsvVar(UMean[0],UMean[1],UMean[2],UMean[3]);
            }
            //Khong luu T mean
            //surfaceFields::T[edge][nG]=TP;
        }
    }

    namespace NSFEqBCs {
    void findMaxLxFConstantOnEdge(int element, int edge)
    {
        /* Ham tinh he so C tren 1 edge trong truong hop Flux type la LxF,
         * cell id input vao ham la master cell cua edge
        */

        double TM, TP;
        std::vector<double> UM(4), UP(4);
        std::vector<double> CArray(mathVar::nGauss+1);
        for (int nG = 0; nG <= mathVar::nGauss; nG++)
        {
            TP=surfaceFields::T[edge][nG];
            TM=surfaceFields::T[edge][nG+mathVar::nGauss+1];

            for (int i = 0; i < 4; i++)
            {
                //std::tie(UMaster[i], USlave[i]) = math::internalSurfaceValue(iedge, masterCell, nG, i + 1, 2);
                std::tie(UP[i], UM[i]) = auxUlti::getUAtInterfaces(edge, element, nG, i + 1);
            }

            double uMagP(sqrt(pow(UP[1]/UP[0],2)+pow(UP[2]/UP[0],2))),
                    uMagM(sqrt(pow(UM[1]/UM[0],2)+pow(UM[2]/UM[0],2))),
                    aP(math::CalcSpeedOfSound(TP)),
                    aM(math::CalcSpeedOfSound(TM));

            CArray[nG]=(math::numericalFluxes::constantC(uMagP,uMagM,aP,aM));
        }
        LxFConst[edge]=*max_element(CArray.begin(), CArray.end());
        //-----------------------------------------------
    }

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

    void calcdUMean(int edge, int element, std::vector<double> &dUXPlus, std::vector<double> &dUYPlus)
    {
        //Compute dU+
        /*
        for (int i = 0; i < 4; i++)
        {
            dUXPlus[i] = math::pointAuxValue(element, a, b, i + 1, 1);
            dUYPlus[i] = math::pointAuxValue(element, a, b, i + 1, 2);
        }
        */
        dUXPlus=math::pointSVars(edge,element,0,0,1,2);
        dUYPlus=math::pointSVars(edge,element,0,0,2,2);
    }

    std::tuple<double, double> calcInnerCellLocalInviscidFlux(std::vector<double> &UG, std::vector<double> &FLocalX, std::vector<double> &FLocalY, std::vector<double> &UL, std::vector<double> &n, double drhoX, double drhoY, double T)
    {
        /* Ham tinh inviscous term tu U,
         * - Xoay U ve he toa do local truoc khi tinh
         * - Gia tri xuat ra la gia tri F_invs da xoay ve he local
         * - Chu y: neu flux type la Lax Friedrichs thi k can xoay
         *
         * Cac bien tra ve (thong qua pass by reference):
         * - FlocalX: inviscid flux theo phuong Ox local, hay la phuong normal
         * - FlocalY: inviscid flux theo phuong Oy local, hay la phuong tangential
         * - UL: conservative variable da duoc xoay ve he local
         * - um, vm: van toc mass diffusion global
        */

        double nx(n[0]), ny(n[1]);
        double rho(UG[0]), rhouG(UG[1]), rhovG(UG[2]), rhoE(UG[3]), rhouG_M, rhovG_M; //G = global, M = mass diffusion
        double rhouL, //L = local, rhouL tuong duong normal direction
                rhovL, //rhovL tuong duong tangential direction
                rhouL_M,
                rhovL_M;

        //Tinh um, vm trong truong hop mass diffusion
        if (flowProperties::massDiffusion)
        {
            rhouG_M=rho*math::massDiffusionFncs::calcTotalVelocity(rho,rhouG/rho,drhoX);
            rhovG_M=rho*math::massDiffusionFncs::calcTotalVelocity(rho,rhovG/rho,drhoY);
        }
        else
        {
            rhouG_M=rhouG;
            rhovG_M=rhovG;
        }

        //Xoay ve he local
        if (DGSchemes::fluxControl::LxF)
        {
            rhouL=rhouG;
            rhovL=rhovG;
            rhouL_M=rhouG_M;
            rhovL_M=rhovG_M;
        }
        else
        {
            std::tie(rhouL,rhovL)=math::rotateToLocalCoordinateSystem(rhouG,rhovG,nx,ny);
            std::tie(rhouL_M,rhovL_M)=math::rotateToLocalCoordinateSystem(rhouG_M,rhovG_M,nx,ny);
        }


        //Cap nhat vector UL
        UL[0]=rho;
        UL[1]=rhouL;    //----> khong tra ve bien rhouL_M
        UL[2]=rhovL;    //----> khong tra ve bien rhovL_M
        UL[3]=rhoE;

        //Tinh inviscid flux
        double
        uL = (rhouL / rho),
        vL = (rhovL / rho),
        uL_M = (rhouL_M / rho),
        vL_M = (rhovL_M / rho),
        totalE = (rhoE / rho),

        /*calculate P*/
        p = math::CalcP(T, rho);

        /*calculate inviscid terms*/
        std::tie(FLocalX[0], FLocalX[1], FLocalX[2], FLocalX[3]) = math::inviscidTerms::calcInvisTermsFromPriVars(rho, uL, uL_M, vL, vL_M, totalE, p, 1);
        std::tie(FLocalY[0], FLocalY[1], FLocalY[2], FLocalY[3]) = math::inviscidTerms::calcInvisTermsFromPriVars(rho, uL, uL_M, vL, vL_M, totalE, p, 2);

        return std::make_tuple(rhouG_M/rho, rhovG_M/rho);
    }

    void calcGhostCellLocalInviscidFlux(int BCType, std::vector<double> &UG, std::vector<double> &FLocalX, std::vector<double> &FLocalY, std::vector<double> &UL, std::vector<double> &n, double drhoX, double drhoY, double um_interior, double vm_interior, double T)
    {
        /* Ham tinh inviscous term tu U,
         * - Xoay U ve he toa do local truoc khi tinh
         * - Gia tri xuat ra la gia tri F_invs da xoay ve he local
         * - Chu y: neu flux type la Lax Friedrichs thi k can xoay
         *
         * Cac bien tra ve (thong qua pass by reference):
         * - FlocalX: inviscid flux theo phuong Ox local, hay la phuong normal
         * - FlocalY: inviscid flux theo phuong Oy local, hay la phuong tangential
         * - UL: conservative variable da duoc xoay ve he local
        */

        double nx(n[0]), ny(n[1]);
        double rho(UG[0]), rhouG(UG[1]), rhovG(UG[2]), rhoE(UG[3]), rhouG_M, rhovG_M; //G = global, M = mass diffusion
        double rhouL, //L = local, rhouL tuong duong normal direction
                rhovL, //rhovL tuong duong tangential direction
                rhouL_M,
                rhovL_M;

        /* Hieu chinh lai van toc tai bien truoc khi dua vao ham tinh flux
         * - Neu bien dang xet la type wall, van toc mass diffusion um, vm deu lay
         * bang -(van toc convective o phia +).
        */
        if (flowProperties::massDiffusion)
        {
            if (BCType==1 && auxUlti::checkTimeVaryingBCAvailable()) //type wall with Maxwell-Smoluchowski BC
            {
                //Van toc mass diffusion bang van toc truot
                rhouG_M=rhouG;
                rhovG_M=rhovG;
            }
            else if (BCType==1 && !auxUlti::checkTimeVaryingBCAvailable())
            {
                //Van toc mass diffusion dao nguoc cua van toc mass diffusion phia +
                rhouG_M=-rho*um_interior;
                rhovG_M=-rho*vm_interior;
            }
            else
            {
                //Tinh van toc mass diffusion binh thuong
                rhouG_M=rho*math::massDiffusionFncs::calcTotalVelocity(rho,rhouG/rho,drhoX);
                rhovG_M=rho*math::massDiffusionFncs::calcTotalVelocity(rho,rhovG/rho,drhoY);
            }
        }
        else
        {
            rhouG_M=rhouG;
            rhovG_M=rhovG;
        }

        //Xoay ve he local
        if (DGSchemes::fluxControl::LxF)
        {
            rhouL=rhouG;
            rhovL=rhovG;
            rhouL_M=rhouG_M;
            rhovL_M=rhovG_M;
        }
        else
        {
            std::tie(rhouL,rhovL)=math::rotateToLocalCoordinateSystem(rhouG,rhovG,nx,ny);
            std::tie(rhouL_M,rhovL_M)=math::rotateToLocalCoordinateSystem(rhouG_M,rhovG_M,nx,ny);
        }


        //Cap nhat vector UL
        UL[0]=rho;
        UL[1]=rhouL;    //----> khong tra ve bien rhouL_M
        UL[2]=rhovL;    //----> khong tra ve bien rhovL_M
        UL[3]=rhoE;

        //Tinh inviscid flux
        double
        uL = (rhouL / rho),
        vL = (rhovL / rho),
        uL_M = (rhouL_M / rho),
        vL_M = (rhovL_M / rho),
        totalE = (rhoE / rho),

        /*calculate P*/
        p = math::CalcP(T, rho);

        /*calculate inviscid terms*/
        std::tie(FLocalX[0], FLocalX[1], FLocalX[2], FLocalX[3]) = math::inviscidTerms::calcInvisTermsFromPriVars(rho, uL, uL_M, vL, vL_M, totalE, p, 1);
        std::tie(FLocalY[0], FLocalY[1], FLocalY[2], FLocalY[3]) = math::inviscidTerms::calcInvisTermsFromPriVars(rho, uL, uL_M, vL, vL_M, totalE, p, 2);
    }

    void NSFEqFluxes(int edgeId, std::vector<double> &Fluxes, int BCType, double TPlus, double TMinus, std::vector<double> &UPlus, std::vector<double> &UMinus, std::vector<double> &dUXPlus, std::vector<double> &dUXMinus, std::vector<double> &dUYPlus, std::vector<double> &dUYMinus, std::vector<double> &normVector)
    {
        /*StressHeat matrix has form:
        [tauXx		tauXy		Qx]
        [tauYx		tauYy		Qy]
        */
        std::vector<std::vector<double>> StressHeatP(2, std::vector<double>(3, 0.0));
        std::vector<std::vector<double>> StressHeatM(2, std::vector<double>(3, 0.0));

        //Local coordinates
        std::vector<double> ULPlus(4, 0.0), ULMinus(4, 0.0);

        std::vector<double>
                //Local fluxes
                invisFluxLocalXPlus(4,0.0),
                invisFluxLocalYPlus(4,0.0),
                visFluxLocalXPlus(4,0.0),
                visFluxLocalYPlus(4,0.0),

                invisFluxLocalXMinus(4,0.0),
                invisFluxLocalYMinus(4,0.0),
                visFluxLocalXMinus(4,0.0),
                visFluxLocalYMinus(4,0.0),

                //vector lay gia tri flux o he toa do global
                invisFluxGlobal(4,0.0),
                visFluxGlobal(4,0.0);

        double umP, vmP;

        /* Voi ham calcLocalInviscidFlux, neu flux type la LxF hoac central, cac bien
         * invisFluxLocalXMaster,.., ULMaster,.. deu luu gia tri tai he toa do global.
         * Neu flux type la Roe hoac HLL, cac bien nay luu gia tri tai he local.
        */
        std::tie(umP,vmP)=BCSupportFncs::NSFEqBCs::calcInnerCellLocalInviscidFlux(UPlus,invisFluxLocalXPlus,invisFluxLocalYPlus,ULPlus,normVector,dUXPlus[0],dUYPlus[0],TPlus);
        BCSupportFncs::NSFEqBCs::calcGhostCellLocalInviscidFlux(BCType,UMinus,invisFluxLocalXMinus,invisFluxLocalYMinus,ULMinus,normVector,dUXMinus[0],dUYMinus[0],umP,vmP,TMinus);

        /* Voi ham calcLocalViscousFlux, vi flux type luon la central, cac bien
         * invisFluxLocalXMaster,.., ULMaster,.. deu luu gia tri tai he toa do global.
        */
        if (flowProperties::viscous)
        {
            process::NSFEq::calcLocalViscousFlux(UPlus,dUXPlus,dUYPlus,visFluxLocalXPlus,visFluxLocalYPlus,normVector,TPlus);
            process::NSFEq::calcLocalViscousFlux(UMinus,dUXMinus,dUYMinus,visFluxLocalXMinus,visFluxLocalYMinus,normVector,TMinus);
        }

        //Tinh flux
        //Vi dang tinh cho cell master cua edge nen phia trong (phia +) la cua cell master
        math::numericalFluxes::calConvectiveFluxes(edgeId, invisFluxGlobal,
                                                   invisFluxLocalXPlus, invisFluxLocalXMinus,
                                                   invisFluxLocalYPlus, invisFluxLocalYMinus,
                                                   ULPlus, ULMinus,
                                                   TPlus, TMinus,
                                                   normVector);
        //Viscous flux tinh bang central flux
        if (flowProperties::viscous)
        {
            math::numericalFluxes::centralFlux(visFluxGlobal,
                                               visFluxLocalXPlus, visFluxLocalXMinus,
                                               visFluxLocalYPlus, visFluxLocalYMinus,
                                               normVector);
        }

        //Luu flux vao vector
        for (int i=0; i<4; i++)
        {
            Fluxes[i]=invisFluxGlobal[i]+visFluxGlobal[i];
        }
    }
	
	std::tuple<double, double> zeroNormGradient(double gradXP, double gradYP, const std::vector<double> &n)
    {
        /* Nguyen tac: normal grad minus = 0
         *             tangential grad minus = tangential grad plus
        */
        double tangGrad(gradXP*n[1]+gradYP*n[0]);
        double gradXM=-tangGrad*n[1];
        double gradYM=tangGrad*n[0];

        return std::make_tuple(gradXM,gradYM);
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
        void wallIsoThermal(std::vector<double> &Fluxes, int element, int edge, int edgeGrp, int nG)
        {
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

            //Compute dU+
            if (flowProperties::viscous)
                BCSupportFncs::NSFEqBCs::calcdUPlus(edge,element,a,b,dUXPlus,dUYPlus);
            dUXMinus=dUXPlus;
            dUYMinus=dUYPlus;

            BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge, Fluxes, 1, TPlus, TMinus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
        }

        void wallAdiabatic(std::vector<double> &Fluxes, int element, int edge, int nG)
        {
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
            if (flowProperties::viscous)
                BCSupportFncs::NSFEqBCs::calcdUPlus(edge,element,a,b,dUXPlus,dUYPlus);

            //zero normal temperature gradient (sua lai theo cach tong quat)
            dUXMinus = dUXPlus;
            dUYMinus = dUYPlus;
            dUXMinus[3] = 0;
            dUYMinus[3] = 0;

            BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge, Fluxes, 1,TPlus, TPlus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
        }


        void wall_NonEquilibrium(std::vector<double> &Fluxes, int element, int edge, int nG)
        {
            int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
            std::vector<double> UMinus(4, 0.0),
                UPlus(4, 0.0),
                dUXPlus(4, 0.0),
                dUYPlus(4, 0.0),
                dUXMinus(4, 0.0),
                dUYMinus(4, 0.0),
                norm(2, 0.0);

            //Ham nay chi khac ham wallIsothermal o cho TMinus = SurfaceBCFields::TBc[loc][nG]
            double a(0.0), b(0.0), TPlus(0.0),TMinus, nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            norm[0] = nx;
            norm[1] = ny;

            //Compute U+
            for (int i = 0; i < 4; i++)
            {
                std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
            }
            TPlus=surfaceFields::T[edge][nG];
            TMinus=SurfaceBCFields::TBc[loc];

            //Compute dU+
            if (flowProperties::viscous)
                BCSupportFncs::NSFEqBCs::calcdUPlus(edge,element,a,b,dUXPlus,dUYPlus);
            dUXMinus=dUXPlus;
            dUYMinus=dUYPlus;

            BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge, Fluxes, 1, TPlus, TMinus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
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
        void inFlow(std::vector<double> &Fluxes, int element, int edge, int edgeGrp, int nG)
        {
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
            if (flowProperties::viscous)
                BCSupportFncs::NSFEqBCs::calcdUPlus(edge,element,a,b,dUXPlus,dUYPlus);

            dUXMinus=dUXPlus;
            dUYMinus=dUYPlus;

            BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge, Fluxes, 2, TPlus, TMinus, UPlus, UMinus, dUXPlus, dUXPlus, dUYPlus, dUYPlus, norm);
        }

        void outFlow(std::vector<double> &Fluxes, int element, int edge, int edgeGrp, int nG)
        {
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
            if (flowProperties::viscous)
                BCSupportFncs::NSFEqBCs::calcdUPlus(edge,element,a,b,dUXPlus,dUYPlus);
            //BCSupportFncs::NSFEqBCs::calcdUMean(edge,element,dUXPlus,dUYPlus);
            //dUXMinus=dUXPlus;
            //dUYMinus=dUYPlus;
            for (int i = 0; i < 4; i++)
            {
                std::tie(dUXMinus[i],dUYMinus[i])=BCSupportFncs::NSFEqBCs::zeroNormGradient(dUXPlus[i],dUYPlus[i],norm);
            }

            BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge, Fluxes, 2, TPlus, TMinus, UPlus, UMinus, dUXPlus, dUXPlus, dUYPlus, dUYPlus, norm);
        }
    }

    void Symmetry(std::vector<double> &Fluxes, int element, int edge, int nG)
    {
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
        if (flowProperties::viscous)
            BCSupportFncs::NSFEqBCs::calcdUPlus(edge,element,a,b,dUXPlus,dUYPlus);
        //BCSupportFncs::NSFEqBCs::calcdUMean(edge,element,dUXPlus,dUYPlus);
        for (int i = 0; i < 4; i++)
        {
            dUXMinus[i] = dUXPlus[i] - 2 * (dUXPlus[i] * nx + dUYPlus[i] * ny)*nx;
            dUYMinus[i] = dUYPlus[i] - 2 * (dUXPlus[i] * nx + dUYPlus[i] * ny)*ny;
        }

        BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge, Fluxes, 3, TPlus, TPlus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
    }

    void matched(std::vector<double> &Fluxes, int element, int edge, int nG)
    {
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
        if (flowProperties::viscous)
        {
            BCSupportFncs::NSFEqBCs::calcdUPlus(edge,element,a,b,dUXPlus,dUYPlus);
            BCSupportFncs::parallel::calcdUMinus(edge,nG,dUXMinus,1);
            BCSupportFncs::parallel::calcdUMinus(edge,nG,dUYMinus,2);
        }

        BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge, Fluxes, 4, TPlus, TMinus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
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
            double a(0.0), b(0.0), muP, muM, TMinus;
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

            //Compute plus values--------------------------------------------------------
            //Trong ham nay da co tinh T va luu vao  surfaceFields::T[edge][nG]
            BCSupportFncs::auxilaryBCs::calcUPlus(element,edge,nG,UPlus);

            //Compute minus values--------------------------------------------------------
            UMinus[0] = UPlus[0];
            //Tinh TMinus sao cho TBC = 0.5(T+ + T-)
            TMinus=2*bcValues::TBCFixed[edgeGrp - 1]-surfaceFields::T[edge][nG];
            if (TMinus<=0)
                TMinus=0.1;

            if (bcValues::UBcType[edgeGrp - 1] == 2)
            {
                UMinus[1] = -UPlus[1];
                UMinus[2] = -UPlus[2];
                UMinus[3] = UMinus[0]*material::Cv*TMinus+0.5*(pow(UMinus[1],2)+pow(UMinus[2],2))/UMinus[0];
            }
            else if (bcValues::UBcType[edgeGrp - 1] == 3)
            {
                UMinus[1] = 2*bcValues::uBCFixed[edgeGrp - 1]*UMinus[0]-UPlus[1];
                UMinus[2] = 2*bcValues::vBCFixed[edgeGrp - 1]*UMinus[0]-UPlus[2];
                UMinus[3] = UMinus[0]*material::Cv*TMinus+0.5*(pow(UMinus[1],2)+pow(UMinus[2],2))/UMinus[0];
            }
            //-----------------------------------------------------------------------------

            muP=math::CalcVisCoef(surfaceFields::T[edge][nG]);
            muM=math::CalcVisCoef(TMinus);
            //Luu TMinus (vi isothermal wall nen TMinus = TBC)
            surfaceFields::T[edge][nG+mathVar::nGauss+1]=TMinus;

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
                UMinus[1] = 2*bcValues::uBCFixed[edgeGrp - 1]*UMinus[0]-UPlus[1];
                UMinus[2] = 2*bcValues::vBCFixed[edgeGrp - 1]*UMinus[0]-UPlus[2];
                UMinus[3] = UMinus[0]*material::Cv*surfaceFields::T[edge][nG]+0.5*(pow(UMinus[1],2)+pow(UMinus[2],2))/UMinus[0];
            }
            //-----------------------------------------------------------------------------

            //Nhan mu vao vector flux truoc khi return
            muP=math::CalcVisCoef(surfaceFields::T[edge][nG]);
            //Luu TMinus (vi adiabatic wall nen TMinus = TPlus)
            surfaceFields::T[edge][nG+mathVar::nGauss+1]=surfaceFields::T[edge][nG];

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


        std::vector <std::vector<double>> wall_NonEquilibrium(int element, int edge, int nG)
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
            double a(0.0), b(0.0), muPlus, muMinus, TMinus;
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            //Lay local id cua edge tren BC
            int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));

            //Compute plus values--------------------------------------------------------
            BCSupportFncs::auxilaryBCs::calcUPlus(element,edge,nG,UPlus);

            //Compute minus values--------------------------------------------------------
            //Modify sao cho TBC = 0.5(T+ + T-)
            TMinus=2*SurfaceBCFields::TBc[loc]-surfaceFields::T[edge][nG];
            if (TMinus<=0)
                TMinus=0.1;

            UMinus[0] = UPlus[0];
            UMinus[1] = 2*SurfaceBCFields::uBc[loc]*UMinus[0]-UPlus[1];
            UMinus[2] = 2*SurfaceBCFields::vBc[loc]*UMinus[0]-UPlus[2];
            //UMinus[3] = UMinus[0]*material::Cv*SurfaceBCFields::TBc[loc]+0.5*(pow(UMinus[1],2)+pow(UMinus[2],2))/UMinus[0];
            UMinus[3] = UMinus[0]*material::Cv*TMinus+0.5*(pow(UMinus[1],2)+pow(UMinus[2],2))/UMinus[0];
            //-----------------------------------------------------------------------------

            muPlus=math::CalcVisCoef(surfaceFields::T[edge][nG]);
            //muMinus=math::CalcVisCoef(SurfaceBCFields::TBc[loc]);
            muMinus=math::CalcVisCoef(TMinus);

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
            //Luu TMinus (TMinus = TBC)
            surfaceFields::T[edge][nG+mathVar::nGauss+1]=bcValues::TBCFixed[edgeGrp - 1];
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
                UPlus(4, 0.0), UMean(4, 0.0),
                UMinus(4, 0.0);
            double a(0.0), b(0.0), pInternal(0), muP, muM, nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

            //Compute plus values--------------------------------------------------------
            BCSupportFncs::auxilaryBCs::calcUPlus(element,edge,nG,UPlus);
            BCSupportFncs::auxilaryBCs::calcUMean(element,edge,nG,UMean);

            //Compute minus values--------------------------------------------------------
            //Modify u, v if reversed flow is detected
            if (BCSupportFncs::checkInflow(UPlus[1]/UPlus[0],UPlus[2]/UPlus[0], nx, ny) || BCSupportFncs::checkInflow(UMean[1]/UMean[0],UMean[2]/UMean[0], nx, ny)) //Tai outlet neu co inflow thi bi reversed flow
            {
                if (!warningFlag::reversedFlowOccur)
                {
                    warningFlag::reversedFlowOccur=true;
                }
                BCSupportFncs::auxilaryBCs::calcUReversedFlow(edge,nG,UPlus);
            }

            UMinus = UPlus;
            //---------------------------------------------------------------------------
            //Nhan mu vao vector flux truoc khi return
            muP=math::CalcVisCoef(surfaceFields::T[edge][nG]);
            if (refValues::subsonic)
            {
                muM=math::CalcVisCoef(bcValues::TBCFixed[edgeGrp - 1]);

                //Luu TMinus (TMinus = TBC)
                surfaceFields::T[edge][nG+mathVar::nGauss+1]=bcValues::TBCFixed[edgeGrp - 1];
            }
            else
            {
                muM=muP;

                //Luu TMinus (TMinus = TPlus)
                surfaceFields::T[edge][nG+mathVar::nGauss+1]=surfaceFields::T[edge][nG];
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
        //BCSupportFncs::auxilaryBCs::calcUMean(element,edge,nG,UPlus);

        //Compute minus values--------------------------------------------------------
        UMinus[0] = UPlus[0];
        UMinus[1] = UPlus[1] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*nx;
        UMinus[2] = UPlus[2] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*ny;
        UMinus[3] = UPlus[3];
        //----------------------------------------------------------------------------

        //Nhan mu vao vector flux truoc khi return
        muP=math::CalcVisCoef(surfaceFields::T[edge][nG]);
        //Luu TMinus (TMinus = TPlus)
        surfaceFields::T[edge][nG+mathVar::nGauss+1]=surfaceFields::T[edge][nG];

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
