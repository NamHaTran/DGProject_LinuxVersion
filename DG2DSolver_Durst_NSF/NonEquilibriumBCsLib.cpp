#include "NonEquilibriumBCsLib.h"
#include "DGMath.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <math.h>
#include <iostream>

void updateTimeVaryingBCs()
{
    /*
     * Ham update gia tri tren cac field cua surfaceBCFields, chi dung cho cac BC bien thien theo thoi gian.
     * Hien tai, ham su dung cho temperatureJump va slip conditions
    */
    if (auxUlti::checkTimeVaryingBCAvailable() && flowProperties::viscous)
    {
        if (systemVar::currentProc==0)
        {
            std::cout<<"Updating time varying BCs.\n";
        }

        int globleEdge(0);
        for (int ilocalEdge=0; ilocalEdge<meshVar::numBCEdges; ilocalEdge++)
        {
            globleEdge=auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(ilocalEdge);
            int edgeGrp(auxUlti::getGrpOfEdge(globleEdge));
            int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]);
            if (UType == 5 && TType == 6)
            {
                nonequilibriumBCs::MaxwellSmoluchowskiBC::implicit2ndOrderMethod(globleEdge,edgeGrp);
            }
        }
    }
}

namespace nonequilibriumBCs
{
    namespace MaxwellSmoluchowskiBC
    {
        void explicitMethod_FVMType(int edge, int edgeGrp)
        {
            //Dieu kien bien Maxwell-Smoluchowski

            int element(0), tempE(0), localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
            std::tie(element,tempE)=auxUlti::getMasterServantOfEdge(edge);
            //Tinh cac gia tri can thiet
            double aConst(0.0), dux(0.0), duy(0.0), dvx(0.0), dvy(0.0), lambda(0.0),
                    TC, uC, vC,
                    rhoC(rho[element][0]),
                    rhouC(rhou[element][0]),
                    rhovC(rhov[element][0]),
                    rhoEC(rhoE[element][0]),
                    TVal(SurfaceBCFields::TBc[localEdgeId]),
                    uVal(SurfaceBCFields::uBc[localEdgeId]),
                    vVal(SurfaceBCFields::vBc[localEdgeId]),
                    delta(meshVar::distanceFromCentroidToBCEdge[localEdgeId]),
                    muVal(0.0), nx, ny, dTx, dTy,dTn, dun, dvn;

            double TWall(bcValues::TBCFixed[edgeGrp - 1]), TJump(0.0),
            uWall(bcValues::uBCFixed[edgeGrp - 1]), vWall(bcValues::vBCFixed[edgeGrp - 1]), a, b;

            //Tinh gia tri tai normal projection cua center xuong BC
            a=meshVar::normProjectionOfCenterToBCEdge_standardSysCoor[localEdgeId][0];
            b=meshVar::normProjectionOfCenterToBCEdge_standardSysCoor[localEdgeId][1];
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
            //muVal=math::CalcVisCoef(TVal);

            //lambda=math::calcMeanFreePath(muVal,rhoBC,TVal);
            aConst=(2-bcValues::sigmaT)*2*material::gamma*lambda/(bcValues::sigmaT*material::Pr*(material::gamma+1));

            TJump=TWall-aConst*dTn;
            if (TJump<0)
            {
                //TJump=TVal;
            }
            //Update to surfaceBCfields
            //SurfaceBCFields::TBc[localEdgeId]=TJump;

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
            double PImc_n1(muVal*(nx*(dux-(2/3)*(dux+dvy))+ny*duy)), PImc_n2(muVal*(nx*dvx+ny*(dvy-(2/3)*(dux+dvy))));
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
                    uSlip=uWall-c*lambda*SdotndotPImc[0]/muVal-(3/4)*SdotDivT[0]*muVal/(rhoBC*TJump)-c*lambda*div_n_Sdotu[0],
                    vSlip=vWall-c*lambda*SdotndotPImc[1]/muVal-(3/4)*SdotDivT[1]*muVal/(rhoBC*TJump)-c*lambda*div_n_Sdotu[1];

            //Update
            SurfaceBCFields::uBc[localEdgeId]=uSlip;
            SurfaceBCFields::vBc[localEdgeId]=vSlip;
        }

        void explicitMethod_DGType(int edge, int edgeGrp)
        {
            //Dieu kien bien Maxwell-Smoluchowski
            int element(0), tempE(0), localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
            std::tie(element,tempE)=auxUlti::getMasterServantOfEdge(edge);

            //Tinh cac gia tri can thiet-----------------------------------------------------
            //Lay cac gia tri da luu
            double aConst(0.0), dux(0.0), duy(0.0), dvx(0.0), dvy(0.0), lambda(0.0), dEx, dEy,
                    TVal(SurfaceBCFields::TBc[localEdgeId]),
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
                //TJump=TVal;
            }
            //Update to surfaceBCfields
            SurfaceBCFields::TBc[localEdgeId]=TJump;

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
                    uSlip=uWall-c*lambda*SdotndotPImc[0]-(3/4)*SdotDivT[0]/(rhoBC*TJump)-c*lambda*div_n_Sdotu[0],
                    vSlip=vWall-c*lambda*SdotndotPImc[1]-(3/4)*SdotDivT[1]/(rhoBC*TJump)-c*lambda*div_n_Sdotu[1];

            //Update
            SurfaceBCFields::uBc[localEdgeId]=uSlip;
            SurfaceBCFields::vBc[localEdgeId]=vSlip;
        }

        /*
        void explicitMethod_DGType(int edge, int edgeGrp)
        {
            //Dieu kien bien Maxwell-Smoluchowski
            int element(0), tempE(0), localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
            std::tie(element,tempE)=auxUlti::getMasterServantOfEdge(edge);

            //Tinh cac gia tri can thiet-----------------------------------------------------
            //Lay cac gia tri da luu
            double dux(0.0), duy(0.0), dvx(0.0), dvy(0.0), lambda(0.0),
                    TVal(SurfaceBCFields::TBc[localEdgeId]),
                    muVal(0.0), nx, ny, dTx, dTy, dTn;

            double TWall(bcValues::TBCFixed[edgeGrp - 1]), TJump(0.0), TC,
                    rhoC(rho[element][0]),
                    rhouC(rhou[element][0]),
                    rhovC(rhov[element][0]),
                    rhoEC(rhoE[element][0]),
                    uWall(bcValues::uBCFixed[edgeGrp - 1]), vWall(bcValues::vBCFixed[edgeGrp - 1]), a, b,
                    delta(meshVar::distanceFromCentroidToBCEdge[localEdgeId]);
            nx = auxUlti::getNormVectorComp(element, edge, 1);
            ny = auxUlti::getNormVectorComp(element, edge, 2);

            //Tinh gia tri tai normal projection cua center xuong BC
            a=meshVar::normProjectionOfCenterToBCEdge_standardSysCoor[localEdgeId][0];
            b=meshVar::normProjectionOfCenterToBCEdge_standardSysCoor[localEdgeId][1];

            double rhoBC(math::pointValue(element,a,b,1,2)),
                    rhouBC(math::pointValue(element,a,b,2,2)),
                    rhovBC(math::pointValue(element,a,b,3,2)),
                    dRhoBCx(math::pointAuxValue(element,a,b,1,1)),
                    dRhoBCy(math::pointAuxValue(element,a,b,1,2)),
                    dRhouBCx(math::pointAuxValue(element,a,b,2,1)),
                    dRhouBCy(math::pointAuxValue(element,a,b,2,2)),
                    dRhovBCx(math::pointAuxValue(element,a,b,3,1)),
                    dRhovBCy(math::pointAuxValue(element,a,b,3,2));

            //Tinh cac thanh phan dao ham theo 2 phuong
            //dTx=math::calcDivTFromAuxVariable(dRhoEBCx,dRhoBCx,rhoEBC,rhoBC);
            //dTy=math::calcDivTFromAuxVariable(dRhoEBCy,dRhoBCy,rhoEBC,rhoBC);

            dux=math::calcRhouvEDeriv(dRhouBCx,dRhoBCx,rhouBC,rhoBC);
            duy=math::calcRhouvEDeriv(dRhouBCy,dRhoBCy,rhouBC,rhoBC);

            dvx=math::calcRhouvEDeriv(dRhovBCx,dRhoBCx,rhovBC,rhoBC);
            dvy=math::calcRhouvEDeriv(dRhovBCy,dRhoBCy,rhovBC,rhoBC);

            //1. Smoluchowski temperature jump-----------------------------------------------------------------------------
            //Gia tri trung binh (tai centroid)
            TC=math::CalcTFromConsvVar(rhoC,rhouC,rhovC,rhoEC);

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
            SurfaceBCFields::TBc[localEdgeId]=TJump;

            //Tinh lai gia tri muVal tu TJump moi cap nhat
            muVal=math::CalcVisCoef(TJump);

            //Gia tri normal gradient
            dTn=(TJump-TC)*muVal/delta;

            //Tinh cac thanh phan dao ham theo 2 phuong
            dTx=dTn*nx;
            dTy=dTn*ny;

            //2. Maxwell slip velocity--------------------------------------------------------
            //Tinh lai gia tri lambda tu TJump moi cap nhat
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
                    uSlip=uWall-c*lambda*SdotndotPImc[0]-(3/4)*SdotDivT[0]/(rhoBC*TJump)-c*lambda*div_n_Sdotu[0],
                    vSlip=vWall-c*lambda*SdotndotPImc[1]-(3/4)*SdotDivT[1]/(rhoBC*TJump)-c*lambda*div_n_Sdotu[1];

            //Limit uSlip, vSlip
            if (fabs(uSlip)>fabs(iniValues::uIni))
            {
                uSlip=0.0;
            }

            if (fabs(vSlip)>fabs(iniValues::uIni))
            {
                vSlip=0.0;
            }

            //Update
            SurfaceBCFields::uBc[localEdgeId]=uSlip;
            SurfaceBCFields::vBc[localEdgeId]=vSlip;
        }*/

        void implicit2ndOrderMethod(int edge, int edgeGrp)
        {
            //Dieu kien bien Maxwell-Smoluchowski

            int element(0), tempE(0), localEdgeId(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
            std::tie(element,tempE)=auxUlti::getMasterServantOfEdge(edge);
            //Tinh cac gia tri can thiet
            double lambda(0.0),
                    TC, uC, vC,
                    rhoC(rho[element][0]),
                    rhouC(rhou[element][0]),
                    rhovC(rhov[element][0]),
                    rhoEC(rhoE[element][0]),
                    //TVal(SurfaceBCFields::TBc[localEdgeId]),
                    delta(meshVar::distanceFromCentroidToBCEdge[localEdgeId]),
                    muVal(0.0), nx, ny, dTx, dTy,dTn;

            double TWall(bcValues::TBCFixed[edgeGrp - 1]), TJump(0.0),
            uWall(bcValues::uBCFixed[edgeGrp - 1]), vWall(bcValues::vBCFixed[edgeGrp - 1]), a, b;

            //Tinh gia tri tai normal projection cua center xuong BC
            a=meshVar::normProjectionOfCenterToBCEdge_standardSysCoor[localEdgeId][0];
            b=meshVar::normProjectionOfCenterToBCEdge_standardSysCoor[localEdgeId][1];

            double rhoBC(math::pointValue(element,a,b,1,2));
            nx = auxUlti::getNormVectorComp(element, edge, 1);
            ny = auxUlti::getNormVectorComp(element, edge, 2);

            //Gia tri trung binh (tai centroid)
            TC=math::CalcTFromConsvVar(rhoC,rhouC,rhovC,rhoEC);
            uC=rhouC/rhoC;
            vC=rhovC/rhoC;

            //1. Smoluchowski temperature jump-----------------------------------------------------------------------------
            double AT(material::viscosityCoeff::Sutherland::As*pow(3.1416/(2*material::R),0.5)/rhoBC),
                    BT(2*material::gamma*(2-bcValues::sigmaT)/(bcValues::sigmaT*material::Pr*(material::gamma+1)));
            double ATT(delta+AT*BT),
                    BTT(material::viscosityCoeff::Sutherland::Ts*delta-TWall*delta-AT*BT*TC),
                    CTT(-TWall*material::viscosityCoeff::Sutherland::Ts*delta), T1, T2;

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
                //TJump=TVal;
            }
            //Update to surfaceBCfields
            SurfaceBCFields::TBc[localEdgeId]=TJump;

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
            double D1(muVal*(nx*nx/(3*delta)+ny*ny/delta)),
                    D2(muVal*(-2*nx*ny/(3*delta))),
                    D3(muVal*(-uC*ny*ny-uC*nx*nx/3+2*nx*ny*vC/3)),
                    D4(muVal*(-2*nx*ny/(3*delta))),
                    D5(muVal*(ny*ny/(3*delta)+nx*nx/delta)),
                    D6(muVal*(-vC*nx*nx-vC*ny*ny/3+2*nx*ny*uC/3));

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
            SurfaceBCFields::uBc[localEdgeId]=uSlip;
            SurfaceBCFields::vBc[localEdgeId]=vSlip;
        }
    }
}
