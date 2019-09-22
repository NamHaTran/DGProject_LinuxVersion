#include "DGProcLib.h"
#include "DGMath.h"
#include "VarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "dynamicVarDeclaration.h"
#include <tuple>
#include "DGBCsLib.h"
#include <algorithm>
#include "DGIOLib.h"
#include <iostream>
#include "DGLimiterLib.h"
#include <math.h>
#include <mpi.h>

namespace meshParam
{
	void GaussParam()
	{
		math::Gauss(mathVar::nGauss);
		math::GaussLobatto(mathVar::nGauss);  //run GaussLobatto for applying limiter
        for (int na = 0; na <= mathVar::nGauss; na++)  //nGauss is started from 0
		{
            for (int nb = 0; nb <= mathVar::nGauss; nb++)
			{
				mathVar::GaussPts[na][nb][0] = mathVar::xGauss[na];
				mathVar::GaussPts[na][nb][1] = mathVar::xGauss[nb];
				mathVar::GaussLobattoPts[na][nb][0] = mathVar::xGaussLobatto[na];
				mathVar::GaussLobattoPts[na][nb][1] = mathVar::xGaussLobatto[nb];

				mathVar::wGaussPts[na][nb][0] = mathVar::wGauss[na];
				mathVar::wGaussPts[na][nb][1] = mathVar::wGauss[nb];
				mathVar::wGaussLobattoPts[na][nb][0] = mathVar::wGaussLobatto[na];
				mathVar::wGaussLobattoPts[na][nb][1] = mathVar::wGaussLobatto[nb];
			}
		}
	}

	void JacobianParam()
	{
		/*2D Jacobi*/
		double a(0.0), b(0.0);
        for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
            for (int na = 0; na <= mathVar::nGauss; na++)  //nGauss is started from 0
			{
                for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					a = mathVar::GaussPts[na][nb][0];
					b = mathVar::GaussPts[na][nb][1];
					meshVar::J2D[ielem][na][nb] = math::J2DCal(ielem, a, b);
				}
			}
		}
	}

	void basisFcParam()
	{
		double a(0.0), b(0.0);
        for (int na = 0; na <= mathVar::nGauss; na++)
		{
            for (int nb = 0; nb <= mathVar::nGauss; nb++)
			{
				//Triangle
				a = mathVar::GaussPts[na][nb][0];
				b = mathVar::GaussPts[na][nb][1];
				math::basisFc(a, b, 3);
				math::dBasisFc(a, b, 3);
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					mathVar::BPts_Tri[order][na][nb] = mathVar::B[order];

					mathVar::dBaPts_Tri[order][na][nb] = mathVar::dBa[order];
					mathVar::dBbPts_Tri[order][na][nb] = mathVar::dBb[order];
				}

				//Quadrilateral
				math::basisFc(a, b, 4);
				math::dBasisFc(a, b, 4);
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					mathVar::BPts_Quad[order][na][nb] = mathVar::B[order];

					mathVar::dBaPts_Quad[order][na][nb] = mathVar::dBa[order];
					mathVar::dBbPts_Quad[order][na][nb] = mathVar::dBb[order];
				}
			}
		}
	}

	void derivCoordinates()
	{
		double dxa(0.0), dxb(0.0), dya(0.0), dyb(0.0);
		double a(0.0), b(0.0);
        for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
            for (int na = 0; na <= mathVar::nGauss; na++)
			{
                for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					a = mathVar::GaussPts[na][nb][0];
					b = mathVar::GaussPts[na][nb][1];
					std::tie(dxa, dxb, dya, dyb) = math::Calc_dxydab(ielem, a, b);
					meshVar::dxa[ielem][na][nb] = dxa;
					meshVar::dxb[ielem][na][nb] = dxb;
					meshVar::dya[ielem][na][nb] = dya;
					meshVar::dyb[ielem][na][nb] = dyb;
				}
			}
		}
	}

	void calcCellMetrics()
	{
		int elemType(0);
		double xCG(0.0), yCG(0.0), cellArea(0.0);
        for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
		{
			elemType = auxUlti::checkType(nelem);
			std::vector<double> xCoor(elemType, 0.0),
				yCoor(elemType, 0.0),
				xSubTriCoor(elemType, 0.0),
				ySubTriCoor(elemType, 0.0);
            for (int i = 0; i < elemType; i++)
			{
				std::tie(xCoor[i], yCoor[i]) = auxUlti::getElemCornerCoord(nelem, i);
			}
            //1. Calculate cell size
            if (elemType==3)
            {
                std::tie(meshVar::cellSize[nelem], cellArea) = math::geometricOp::calROfInscribedCircleOfTriElement(xCoor, yCoor);
            }
            else if (elemType==4)
            {
                std::tie(meshVar::cellSize[nelem], cellArea) = math::geometricOp::calSizeOfQuadElement(xCoor, yCoor);
            }

			//2. Compute geometric center of cell
			std::tie(xCG, yCG) = math::geometricOp::calcGeoCenter(xCoor, yCoor, elemType);

			//3. Compute centroid of cell
			if (elemType == 3)
			{
				meshVar::geoCenter[nelem][0] = xCG;
				meshVar::geoCenter[nelem][1] = yCG;
			}
			else  //quad element
			{
				//3. Compute geometric center of sub-triangles which creates polygon
				std::tie(meshVar::geoCenter[nelem][0], meshVar::geoCenter[nelem][1]) = math::geometricOp::calcQuadCentroid(nelem, xCG, yCG, cellArea);
			}

			//4. Calculate local cell size
			meshVar::localCellSize[nelem] = math::geometricOp::calLocalCellSize(nelem, cellArea);
		}
	}

    void calcEdgeLength()
    {
        //Jacobi 1D = edgeLength/2
        int pt1(-1), pt2(-1);
        for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
        {
            pt1=meshVar::inpoed[iedge][0];
            pt2=meshVar::inpoed[iedge][1];
            meshVar::J1D[iedge]=0.5*math::geometricOp::calDistBetween2Points(meshVar::Points[pt1][0], meshVar::Points[pt1][1], meshVar::Points[pt2][0], meshVar::Points[pt2][1]);
        }
    }

	void calcStiffMatrixCoeffs()
	{
        std::vector<std::vector<double>>StiffMatrix(mathVar::orderElem + 1,std::vector<double>(mathVar::orderElem + 1, 0.0));
        for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
		{
            for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
			{
				stiffMatrixCoeffs[nelem][iorder] = process::calculateStiffMatrixElement(nelem, iorder, iorder);
                //for (int iorder2 = 0; iorder2 <= mathVar::orderElem; ++iorder2) {
                    //StiffMatrix[iorder][iorder2] = process::calculateStiffMatrixElement(nelem, iorder, iorder2);
                //}
			}
		}
	}
}

namespace process
{
	void setIniValues()
	{
		iniValues::rhoIni = iniValues::pIni / (material::R*iniValues::TIni);
		iniValues::eIni = material::Cv*iniValues::TIni; //+0.5*(pow(iniValues::uIni, 2) + pow(iniValues::vIni, 2) + pow(iniValues::wIni, 2));
		iniValues::muIni = math::CalcVisCoef(iniValues::TIni);
		std::vector<double> iniRho(mathVar::orderElem + 1),
			iniRhou(mathVar::orderElem + 1),
			iniRhov(mathVar::orderElem + 1),
			iniRhoE(mathVar::orderElem + 1);

        for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
		{
			iniRho = process::calcIniValues(iniValues::rhoIni, nelement);
			iniRhou = process::calcIniValues(iniValues::rhoIni*iniValues::uIni, nelement);
			iniRhov = process::calcIniValues(iniValues::rhoIni*iniValues::vIni, nelement);
			iniRhoE = process::calcIniValues(iniValues::rhoIni*(iniValues::eIni + 0.5*(pow(iniValues::uIni, 2) + pow(iniValues::vIni, 2))), nelement);
            for (int i = 0; i <= mathVar::orderElem; i++)
			{
				rho[nelement][i] = iniRho[i];
				rhou[nelement][i] = iniRhou[i];
				rhov[nelement][i] = iniRhov[i];
				rhoE[nelement][i] = iniRhoE[i];
			}
		}

		//Calculate limit of rhoE
		//limitVal::rhoEUp = material::Cv*limitVal::TUp*limitVal::rhoUp;
		//limitVal::rhoEDwn = material::Cv*limitVal::TDwn*limitVal::rhoDwn;
	}

	std::vector<double> calcIniValues(double iniVal, int element)
	{
		std::vector<double> RHS(mathVar::orderElem + 1, 0.0);
		std::vector<double> iniVector(mathVar::orderElem + 1, 0.0);

		RHS = process::calcIniValuesRHS(element, iniVal);
        for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
		{
			iniVector[iorder] = RHS[iorder] / stiffMatrixCoeffs[element][iorder];
		}
		return iniVector;
	}

	std::vector<double> calcIniValuesRHS(int element, double iniVal)
	{
		std::vector<double> Out(mathVar::orderElem + 1, 0.0);
		std::vector<std::vector<double>> matrix(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
		double a(0.0), b(0.0);
        for (int order = 0; order <= mathVar::orderElem; order++)
		{
            for (int na = 0; na <= mathVar::nGauss; na++)
			{
                for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
					math::basisFc(a, b, auxUlti::checkType(element));
					matrix[na][nb] = mathVar::B[order];
				}
			}
			Out[order] = math::volumeInte(matrix, element)*iniVal;
		}
		return Out;
	}

	void calcVolumeGaussValues()
	{
		double a(0.0), b(0.0);
        std::vector<double>UVol(4,0.0);
        for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
		{
            for (int na = 0; na <= mathVar::nGauss; na++)
			{
                for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
                    UVol=math::pointUVars(nelem,a,b);
                    volumeFields::rhoVolGauss[nelem][na][nb] = UVol[0];
                    volumeFields::rhouVolGauss[nelem][na][nb] = UVol[1];
                    volumeFields::rhovVolGauss[nelem][na][nb] = UVol[2];
                    volumeFields::rhoEVolGauss[nelem][na][nb] = UVol[3];
				}
			}
		}
	}

    void calcTGauss()
    {
        if (flowProperties::massDiffusion)
        {
            //Volume
            double a(0.0), b(0.0), rhoX(0.0), rhoY(0.0);
            if (systemVar::auxVariables==1)
            {
                for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
                {
                    for (int na = 0; na <= mathVar::nGauss; na++)
                    {
                        for (int nb = 0; nb <= mathVar::nGauss; nb++)
                        {
                            std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
                            rhoX=math::pointAuxValue(nelem,a,b,1,1);
                            rhoY=math::pointAuxValue(nelem,a,b,1,2);
                            volumeFields::T[nelem][na][nb] = math::CalcTFromConsvVar_massDiff(volumeFields::rhoVolGauss[nelem][na][nb],volumeFields::rhouVolGauss[nelem][na][nb],volumeFields::rhovVolGauss[nelem][na][nb],volumeFields::rhoEVolGauss[nelem][na][nb],rhoX,rhoY);
                        }
                    }
                }
            }
            else if (systemVar::auxVariables==2)
            {
                for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
                {
                    for (int na = 0; na <= mathVar::nGauss; na++)
                    {
                        for (int nb = 0; nb <= mathVar::nGauss; nb++)
                        {
                            std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
                            rhoX=math::BR2Fncs::pointAuxValue_vol(nelem,a,b,1,1);
                            rhoY=math::BR2Fncs::pointAuxValue_vol(nelem,a,b,1,2);
                            volumeFields::T[nelem][na][nb] = math::CalcTFromConsvVar_massDiff(volumeFields::rhoVolGauss[nelem][na][nb],volumeFields::rhouVolGauss[nelem][na][nb],volumeFields::rhovVolGauss[nelem][na][nb],volumeFields::rhoEVolGauss[nelem][na][nb],rhoX,rhoY);
                        }
                    }
                }
            }

            //Surface
            int masterCell(-1), slaveCell(-1), bcGrp(0);
            double rhoMaster(0.0), rhouMaster(0.0), rhovMaster(0.0), rhoEMaster(0.0),
                    rhoSlave(0.0), rhouSlave(0.0), rhovSlave(0.0), rhoESlave(0.0),
                    rhoXMaster(0.0), rhoYMaster(0.0), rhoXSlave(0.0), rhoYSlave(0.0);
            for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
            {
                bcGrp = auxUlti::getBCType(iedge);
                if (bcGrp == 0)
                {
                    std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        std::tie(rhoMaster,rhoSlave)=auxUlti::getUAtInterfaces(iedge,masterCell,nG,1);
                        std::tie(rhouMaster,rhouSlave)=auxUlti::getUAtInterfaces(iedge,masterCell,nG,2);
                        std::tie(rhovMaster,rhovSlave)=auxUlti::getUAtInterfaces(iedge,masterCell,nG,3);
                        std::tie(rhoEMaster,rhoESlave)=auxUlti::getUAtInterfaces(iedge,masterCell,nG,4);
                        std::tie(rhoXMaster,rhoXSlave)=math::internalSurfaceDerivativeValue(iedge,masterCell,nG,1,1);
                        std::tie(rhoYMaster,rhoYSlave)=math::internalSurfaceDerivativeValue(iedge,masterCell,nG,1,2);
                        surfaceFields::T[iedge][nG]=math::CalcTFromConsvVar_massDiff(rhoMaster,rhouMaster,rhovMaster,rhoEMaster,rhoXMaster,rhoYMaster);
                        surfaceFields::T[iedge][nG + mathVar::nGauss + 1]=math::CalcTFromConsvVar_massDiff(rhoSlave,rhouSlave,rhovSlave,rhoESlave,rhoXSlave,rhoYSlave);
                    }
                }
            }
        }
        else {
            //Volume
            double a(0.0), b(0.0);
            for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
            {
                for (int na = 0; na <= mathVar::nGauss; na++)
                {
                    for (int nb = 0; nb <= mathVar::nGauss; nb++)
                    {
                        std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
                        volumeFields::T[nelem][na][nb] = math::pointValue(nelem, a, b, 6, 1);
                    }
                }
            }

            //Surface
            int masterCell(-1), slaveCell(-1), bcGrp(0);
            double rhoMaster(0.0), rhouMaster(0.0), rhovMaster(0.0), rhoEMaster(0.0),
                    rhoSlave(0.0), rhouSlave(0.0), rhovSlave(0.0), rhoESlave(0.0);
            for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
            {
                bcGrp = auxUlti::getBCType(iedge);
                if (bcGrp == 0)
                {
                    std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
                    {
                        std::tie(rhoMaster,rhoSlave)=auxUlti::getUAtInterfaces(iedge,masterCell,nG,1);
                        std::tie(rhouMaster,rhouSlave)=auxUlti::getUAtInterfaces(iedge,masterCell,nG,2);
                        std::tie(rhovMaster,rhovSlave)=auxUlti::getUAtInterfaces(iedge,masterCell,nG,3);
                        std::tie(rhoEMaster,rhoESlave)=auxUlti::getUAtInterfaces(iedge,masterCell,nG,4);
                        surfaceFields::T[iedge][nG]=math::CalcTFromConsvVar(rhoMaster,rhouMaster,rhovMaster,rhoEMaster);
                        surfaceFields::T[iedge][nG + mathVar::nGauss + 1]=math::CalcTFromConsvVar(rhoSlave,rhouSlave,rhovSlave,rhoESlave);
                    }
                }
            }
        }
    }

    void calcValuesAtInterface()
    {
        int masterCell(-1), slaveCell(-1), bcGrp(0);
        double tempMaster(0.0), tempSlave(0.0);
        for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
        {
            bcGrp = auxUlti::getBCType(iedge);
            if (bcGrp == 0)
            {
                std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    //rho
                    std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 1, 2);
                    surfaceFields::rho[iedge][nG] = tempMaster;
                    surfaceFields::rho[iedge][nG + mathVar::nGauss + 1] = tempSlave;

                    //rhou
                    std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 2, 2);
                    surfaceFields::rhou[iedge][nG] = tempMaster;
                    surfaceFields::rhou[iedge][nG + mathVar::nGauss + 1] = tempSlave;

                    //rhov
                    std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 3, 2);
                    surfaceFields::rhov[iedge][nG] = tempMaster;
                    surfaceFields::rhov[iedge][nG + mathVar::nGauss + 1] = tempSlave;

                    //rhoE
                    std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 4, 2);
                    surfaceFields::rhoE[iedge][nG] = tempMaster;
                    surfaceFields::rhoE[iedge][nG + mathVar::nGauss + 1] = tempSlave;
                }
            }
        }
    }

	namespace auxEq
    {
    void solveAuxEquation()
    {
        if (systemVar::auxVariables==1)
        {
            //std::vector<std::vector<double>> StiffMatrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));
            std::vector<double> rhoRHSTermOxDir(mathVar::orderElem + 1, 0.0),
                rhoRHSTermOyDir(mathVar::orderElem + 1, 0.0),
                rhouRHSTermOxDir(mathVar::orderElem + 1, 0.0),
                rhouRHSTermOyDir(mathVar::orderElem + 1, 0.0),
                rhovRHSTermOxDir(mathVar::orderElem + 1, 0.0),
                rhovRHSTermOyDir(mathVar::orderElem + 1, 0.0),
                rhoERHSTermOxDir(mathVar::orderElem + 1, 0.0),
                rhoERHSTermOyDir(mathVar::orderElem + 1, 0.0);

            std::vector<double> rhoxVector(mathVar::orderElem + 1, 0.0),
                rhoyVector(mathVar::orderElem + 1, 0.0),
                rhouxVector(mathVar::orderElem + 1, 0.0),
                rhouyVector(mathVar::orderElem + 1, 0.0),
                rhovxVector(mathVar::orderElem + 1, 0.0),
                rhovyVector(mathVar::orderElem + 1, 0.0),
                rhoExVector(mathVar::orderElem + 1, 0.0),
                rhoEyVector(mathVar::orderElem + 1, 0.0);

            for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
            {
                //1) Calculate Right hand side terms
                process::auxEq::BR1::CalcRHSTerm(nelement, rhoRHSTermOxDir, rhoRHSTermOyDir, rhouRHSTermOxDir, rhouRHSTermOyDir, rhovRHSTermOxDir, rhovRHSTermOyDir, rhoERHSTermOxDir, rhoERHSTermOyDir);

                //2) Solve for auxilary variables
                if (!flowProperties::massDiffusion)
                {
                    for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                    {
                        //Ox direction
                        BR1Vars::rhoX[nelement][iorder] = rhoRHSTermOxDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
                        //Oy direction
                        BR1Vars::rhoY[nelement][iorder] = rhoRHSTermOyDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
                    }
                }
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    //Ox direction
                    BR1Vars::rhouX[nelement][iorder] = rhouRHSTermOxDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
                    BR1Vars::rhovX[nelement][iorder] = rhovRHSTermOxDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
                    BR1Vars::rhoEX[nelement][iorder] = rhoERHSTermOxDir[iorder] / stiffMatrixCoeffs[nelement][iorder];

                    //Oy direction
                    BR1Vars::rhouY[nelement][iorder] = rhouRHSTermOyDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
                    BR1Vars::rhovY[nelement][iorder] = rhovRHSTermOyDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
                    BR1Vars::rhoEY[nelement][iorder] = rhoERHSTermOyDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
                }
            }
        }
        else if (systemVar::auxVariables==2)
        {
            for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
            {
                //1) Calculate volume divU
                process::auxEq::BR2::calcVolDivU(nelement);

                //2) Calculate suface divU
                process::auxEq::BR2::calcSurDivU(nelement);
            }
        }

    }
    namespace BR2 {
    void calcVolDivU(int element)
    {
        std::vector<std::vector<double>> gaussMX(mathVar::nGauss+1,std::vector<double>(mathVar::nGauss+1,0.0)),
                gaussMY(mathVar::nGauss+1,std::vector<double>(mathVar::nGauss+1,0.0));
        std::vector<double>rhoC(mathVar::orderElem+1,0.0),
                rhouC(mathVar::orderElem+1,0.0),
                rhovC(mathVar::orderElem+1,0.0),
                rhoEC(mathVar::orderElem+1,0.0);
        int elemType(auxUlti::checkType(element));
        double Bn(0.0),dBmX(0.0),dBmY(0.0),
                GammaRhoX(0.0), GammaRhoY(0.0),
                GammaRhouX(0.0), GammaRhouY(0.0),
                GammaRhovX(0.0), GammaRhovY(0.0),
                GammaRhoEX(0.0), GammaRhoEY(0.0),
                GammaX(0.0), GammaY(0.0), theta1(1.0), theta2(1.0);
        if (!flowProperties::massDiffusion)
        {
            rhoC=auxUlti::getElementConserValuesOfOrder(element,1);
        }
        rhouC=auxUlti::getElementConserValuesOfOrder(element,2);
        rhovC=auxUlti::getElementConserValuesOfOrder(element,3);
        rhoEC=auxUlti::getElementConserValuesOfOrder(element,4);
        for (int n=0;n<=mathVar::orderElem;n++) {
            GammaRhoX=0; GammaRhoY=0;
            GammaRhouX=0; GammaRhouY=0;
            GammaRhovX=0; GammaRhovY=0;
            GammaRhoEX=0; GammaRhoEY=0;

            for (int m=0;m<=mathVar::orderElem;m++) {
                if (m==0)
                {
                    theta1=1.0;
                    theta2=1.0;
                }
                else {
                    theta1=theta1Arr[element];
                    theta2=theta2Arr[element];
                }

                if (elemType==4)
                {
                    for (int na=0;na<=mathVar::nGauss;na++) {
                        for (int nb=0;nb<=mathVar::nGauss;nb++) {
                            Bn=mathVar::BPts_Quad[n][na][nb];
                            dBmX=math::Calc_dBxdBy(element, m, na, nb, 1);
                            dBmY=math::Calc_dBxdBy(element, m, na, nb, 2);
                            gaussMX[na][nb]=Bn*dBmX;
                            gaussMY[na][nb]=Bn*dBmY;
                        }
                    }
                }
                else if (elemType==3)
                {
                    for (int na=0;na<=mathVar::nGauss;na++) {
                        for (int nb=0;nb<=mathVar::nGauss;nb++) {
                            Bn=mathVar::BPts_Tri[n][na][nb];
                            dBmX=math::Calc_dBxdBy(element, m, na, nb, 1);
                            dBmY=math::Calc_dBxdBy(element, m, na, nb, 2);
                            gaussMX[na][nb]=Bn*dBmX;
                            gaussMY[na][nb]=Bn*dBmY;
                        }
                    }
                }

                GammaX=math::volumeInte(gaussMX,element);
                GammaY=math::volumeInte(gaussMY,element);

                if (!flowProperties::massDiffusion)
                {
                    GammaRhoX+=GammaX*rhoC[m]*theta1*theta2;
                    GammaRhoY+=GammaY*rhoC[m]*theta1*theta2;
                }

                GammaRhouX+=GammaX*rhouC[m]*theta2;
                GammaRhouY+=GammaY*rhouC[m]*theta2;

                GammaRhovX+=GammaX*rhovC[m]*theta2;
                GammaRhovY+=GammaY*rhovC[m]*theta2;

                GammaRhoEX+=GammaX*rhoEC[m]*theta2;
                GammaRhoEY+=GammaY*rhoEC[m]*theta2;
            }
            double stiffCoeff(stiffMatrixCoeffs[element][n]);
            if (!flowProperties::massDiffusion)
            {
                BR2Vars::rhoXVol[element][n]=GammaRhoX/stiffCoeff;
                BR2Vars::rhoYVol[element][n]=GammaRhoY/stiffCoeff;
            }

            BR2Vars::rhouXVol[element][n]=GammaRhouX/stiffCoeff;
            BR2Vars::rhouYVol[element][n]=GammaRhouY/stiffCoeff;

            BR2Vars::rhovXVol[element][n]=GammaRhovX/stiffCoeff;
            BR2Vars::rhovYVol[element][n]=GammaRhovY/stiffCoeff;

            BR2Vars::rhoEXVol[element][n]=GammaRhoEX/stiffCoeff;
            BR2Vars::rhoEYVol[element][n]=GammaRhoEY/stiffCoeff;
        }
    }

    void calcSurDivU(int element)
    {
        int elemType(auxUlti::checkType(element)), edgeId(0), BCType(-1);
        double a(0.0), b(0.0),BVar(0.0), nx(0.0), ny(0.0), UP(0.0), UM(0.0), BRConst(4.0),
                dRhoXVar(0.0), dRhoYVar(0.0),
                dRhouXVar(0.0), dRhouYVar(0.0),
                dRhovXVar(0.0), dRhovYVar(0.0),
                dRhoEXVar(0.0), dRhoEYVar(0.0);
        bool isMaster;
        std::vector<std::vector<double>> gaussVector(4, std::vector<double>(2, 0.0));
        std::vector<double> vectorRhoX(mathVar::nGauss+1,0.0),
                vectorRhoY(mathVar::nGauss+1,0.0),
                vectorRhouX(mathVar::nGauss+1,0.0),
                vectorRhouY(mathVar::nGauss+1,0.0),
                vectorRhovX(mathVar::nGauss+1,0.0),
                vectorRhovY(mathVar::nGauss+1,0.0),
                vectorRhoEX(mathVar::nGauss+1,0.0),
                vectorRhoEY(mathVar::nGauss+1,0.0),

                vectordRhoXVol_org(mathVar::orderElem+1,0.0),
                vectordRhoYVol_org(mathVar::orderElem+1,0.0),
                vectordRhouXVol_org(mathVar::orderElem+1,0.0),
                vectordRhouYVol_org(mathVar::orderElem+1,0.0),
                vectordRhovXVol_org(mathVar::orderElem+1,0.0),
                vectordRhovYVol_org(mathVar::orderElem+1,0.0),
                vectordRhoEXVol_org(mathVar::orderElem+1,0.0),
                vectordRhoEYVol_org(mathVar::orderElem+1,0.0);

        if (!flowProperties::massDiffusion)
        {
            for (int iorder=0;iorder<=mathVar::orderElem;iorder++) {
                vectordRhoXVol_org[iorder]=BR2Vars::rhoXVol[element][iorder];
                vectordRhoYVol_org[iorder]=BR2Vars::rhoYVol[element][iorder];
            }
        }
        for (int iorder=0;iorder<=mathVar::orderElem;iorder++) {
            vectordRhouXVol_org[iorder]=BR2Vars::rhouXVol[element][iorder];
            vectordRhouYVol_org[iorder]=BR2Vars::rhouYVol[element][iorder];
            vectordRhovXVol_org[iorder]=BR2Vars::rhovXVol[element][iorder];
            vectordRhovYVol_org[iorder]=BR2Vars::rhovYVol[element][iorder];
            vectordRhoEXVol_org[iorder]=BR2Vars::rhoEXVol[element][iorder];
            vectordRhoEYVol_org[iorder]=BR2Vars::rhoEYVol[element][iorder];
        }

        for (int iedge=0;iedge<elemType;iedge++) {
            edgeId=meshVar::inedel[element][iedge];
            BCType=auxUlti::getBCType(edgeId);
            isMaster=auxUlti::checkMaster(element,edgeId);
            nx=auxUlti::getNormVectorComp(element,edgeId,1);
            ny=auxUlti::getNormVectorComp(element,edgeId,2);

            if (!flowProperties::massDiffusion)
            {
                if (BCType!=0)
                {
                    for (int nG=0;nG<=mathVar::nGauss;nG++) {
                        std::tie(a,b)=auxUlti::getGaussSurfCoor(edgeId,element,nG);
                        gaussVector = auxEqBCsImplement(element, edgeId, nG);

                        //Rho
                        vectorRhoX[nG]=0.5*(gaussVector[0][1]-gaussVector[0][0])*nx;
                        vectorRhoY[nG]=0.5*(gaussVector[0][1]-gaussVector[0][0])*ny;
                    }
                }
                else {
                    for (int nG=0;nG<=mathVar::nGauss;nG++) {
                        std::tie(a,b)=auxUlti::getGaussSurfCoor(edgeId,element,nG);

                        //Rho
                        std::tie(UP,UM)=process::getInternalValuesFromCalculatedArrays(edgeId,element,nG,1);
                        vectorRhoX[nG]=0.5*(UM-UP)*nx;
                        vectorRhoY[nG]=0.5*(UM-UP)*ny;
                    }
                }
            }

            if (BCType!=0)
            {
                for (int nG=0;nG<=mathVar::nGauss;nG++) {
                    std::tie(a,b)=auxUlti::getGaussSurfCoor(edgeId,element,nG);
                    gaussVector = auxEqBCsImplement(element, edgeId, nG);

                    //Rhou
                    vectorRhouX[nG]=0.5*(gaussVector[1][1]-gaussVector[1][0])*nx;
                    vectorRhouY[nG]=0.5*(gaussVector[1][1]-gaussVector[1][0])*ny;

                    //Rhov
                    vectorRhovX[nG]=0.5*(gaussVector[2][1]-gaussVector[2][0])*nx;
                    vectorRhovY[nG]=0.5*(gaussVector[2][1]-gaussVector[2][0])*ny;

                    //RhoE
                    vectorRhoEX[nG]=0.5*(gaussVector[3][1]-gaussVector[3][0])*nx;
                    vectorRhoEY[nG]=0.5*(gaussVector[3][1]-gaussVector[3][0])*ny;
                }
            }
            else {
                for (int nG=0;nG<=mathVar::nGauss;nG++) {
                    std::tie(a,b)=auxUlti::getGaussSurfCoor(edgeId,element,nG);

                    //Rhou
                    std::tie(UP,UM)=process::getInternalValuesFromCalculatedArrays(edgeId,element,nG,2);
                    vectorRhouX[nG]=0.5*(UM-UP)*nx;
                    vectorRhouY[nG]=0.5*(UM-UP)*ny;

                    //Rhov
                    std::tie(UP,UM)=process::getInternalValuesFromCalculatedArrays(edgeId,element,nG,3);
                    vectorRhovX[nG]=0.5*(UM-UP)*nx;
                    vectorRhovY[nG]=0.5*(UM-UP)*ny;

                    //RhoE
                    std::tie(UP,UM)=process::getInternalValuesFromCalculatedArrays(edgeId,element,nG,4);
                    vectorRhoEX[nG]=0.5*(UM-UP)*nx;
                    vectorRhoEY[nG]=0.5*(UM-UP)*ny;
                }
            }

            for (int iorder=0;iorder<=mathVar::orderElem;iorder++) {
                double stiffCoeff(stiffMatrixCoeffs[element][iorder]);
                for (int nG=0;nG<=mathVar::nGauss;nG++) {
                    std::tie(a,b)=auxUlti::getGaussSurfCoor(edgeId,element,nG);
                    math::basisFc(a,b,elemType);
                    BVar=mathVar::B[iorder];
                    //Rho
                    vectorRhoX[nG]*=BVar;
                    vectorRhoY[nG]*=BVar;

                    //Rhou
                    vectorRhouX[nG]*=BVar;
                    vectorRhouY[nG]*=BVar;

                    //Rhov
                    vectorRhovX[nG]*=BVar;
                    vectorRhovY[nG]*=BVar;

                    //RhoE
                    vectorRhoEX[nG]*=BVar;
                    vectorRhoEY[nG]*=BVar;
                }

                if (!flowProperties::massDiffusion)
                {
                    dRhoXVar=math::surfaceInte(vectorRhoX,edgeId)/stiffCoeff;
                    dRhoYVar=math::surfaceInte(vectorRhoY,edgeId)/stiffCoeff;
                }

                dRhouXVar=math::surfaceInte(vectorRhouX,edgeId)/stiffCoeff;
                dRhouYVar=math::surfaceInte(vectorRhouY,edgeId)/stiffCoeff;

                dRhovXVar=math::surfaceInte(vectorRhovX,edgeId)/stiffCoeff;
                dRhovYVar=math::surfaceInte(vectorRhovY,edgeId)/stiffCoeff;

                dRhoEXVar=math::surfaceInte(vectorRhoEX,edgeId)/stiffCoeff;
                dRhoEYVar=math::surfaceInte(vectorRhoEY,edgeId)/stiffCoeff;

                if (isMaster)
                {
                    if (!flowProperties::massDiffusion)
                    {
                        BR2Vars::rhoXSurMaster[edgeId][iorder]=BRConst*dRhoXVar+vectordRhoXVol_org[iorder];
                        BR2Vars::rhoYSurMaster[edgeId][iorder]=BRConst*dRhoYVar+vectordRhoYVol_org[iorder];
                    }

                    BR2Vars::rhouXSurMaster[edgeId][iorder]=BRConst*dRhouXVar+vectordRhouXVol_org[iorder];
                    BR2Vars::rhouYSurMaster[edgeId][iorder]=BRConst*dRhouYVar+vectordRhouYVol_org[iorder];

                    BR2Vars::rhovXSurMaster[edgeId][iorder]=BRConst*dRhovXVar+vectordRhovXVol_org[iorder];
                    BR2Vars::rhovYSurMaster[edgeId][iorder]=BRConst*dRhovYVar+vectordRhovYVol_org[iorder];

                    BR2Vars::rhoEXSurMaster[edgeId][iorder]=BRConst*dRhoEXVar+vectordRhoEXVol_org[iorder];
                    BR2Vars::rhoEYSurMaster[edgeId][iorder]=BRConst*dRhoEYVar+vectordRhoEYVol_org[iorder];
                }
                else {
                    if (!flowProperties::massDiffusion)
                    {
                        BR2Vars::rhoXSurSlave[edgeId][iorder]=BRConst*dRhoXVar+vectordRhoXVol_org[iorder];
                        BR2Vars::rhoYSurSlave[edgeId][iorder]=BRConst*dRhoYVar+vectordRhoYVol_org[iorder];
                    }

                    BR2Vars::rhouXSurSlave[edgeId][iorder]=BRConst*dRhouXVar+vectordRhouXVol_org[iorder];
                    BR2Vars::rhouYSurSlave[edgeId][iorder]=BRConst*dRhouYVar+vectordRhouYVol_org[iorder];

                    BR2Vars::rhovXSurSlave[edgeId][iorder]=BRConst*dRhovXVar+vectordRhovXVol_org[iorder];
                    BR2Vars::rhovYSurSlave[edgeId][iorder]=BRConst*dRhovYVar+vectordRhovYVol_org[iorder];

                    BR2Vars::rhoEXSurSlave[edgeId][iorder]=BRConst*dRhoEXVar+vectordRhoEXVol_org[iorder];
                    BR2Vars::rhoEYSurSlave[edgeId][iorder]=BRConst*dRhoEYVar+vectordRhoEYVol_org[iorder];
                }

                if (!flowProperties::massDiffusion){
                    BR2Vars::rhoXVol[element][iorder]+=dRhoXVar;
                    BR2Vars::rhoYVol[element][iorder]+=dRhoYVar;
                }

                BR2Vars::rhouXVol[element][iorder]+=dRhouXVar;
                BR2Vars::rhouYVol[element][iorder]+=dRhouYVar;

                BR2Vars::rhovXVol[element][iorder]+=dRhovXVar;
                BR2Vars::rhovYVol[element][iorder]+=dRhovYVar;

                BR2Vars::rhoEXVol[element][iorder]+=dRhoEXVar;
                BR2Vars::rhoEYVol[element][iorder]+=dRhoEYVar;
            }
        }
    }


    }
    namespace BR1 {
    void calcValuesAtInterface()
    {
        //mu included
        int masterCell(-1), slaveCell(-1), bcGrp(0);
        double muMaster(0.0), muSlave(0.0), tempMaster(0.0), tempSlave(0.0),TMaster(0.0), TSlave(0.0);
        for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
        {
            bcGrp = auxUlti::getBCType(iedge);
            if (bcGrp == 0)
            {
                std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    std::tie(TMaster, TSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 6, 1);
                    muMaster=math::CalcVisCoef(TMaster);
                    muSlave=math::CalcVisCoef(TSlave);
                    surfaceFields::T[iedge][nG]=TMaster;
                    surfaceFields::T[iedge][nG+ mathVar::nGauss + 1]=TSlave;

                    //rho*mu
                    std::tie(tempMaster, tempSlave) = auxUlti::getUAtInterfaces(iedge, masterCell, nG, 1);
                    surfaceFields::aux_rho[iedge][nG] = tempMaster * muMaster;
                    surfaceFields::aux_rho[iedge][nG + mathVar::nGauss + 1] = tempSlave * muSlave;

                    //rhou*mu
                    std::tie(tempMaster, tempSlave) = auxUlti::getUAtInterfaces(iedge, masterCell, nG, 2);
                    surfaceFields::aux_rhou[iedge][nG] = tempMaster * muMaster;
                    surfaceFields::aux_rhou[iedge][nG + mathVar::nGauss + 1] = tempSlave * muSlave;

                    //rhov*mu
                    std::tie(tempMaster, tempSlave) = auxUlti::getUAtInterfaces(iedge, masterCell, nG, 3);
                    surfaceFields::aux_rhov[iedge][nG] = tempMaster * muMaster;
                    surfaceFields::aux_rhov[iedge][nG + mathVar::nGauss + 1] = tempSlave * muSlave;

                    //rhoE*mu
                    std::tie(tempMaster, tempSlave) = auxUlti::getUAtInterfaces(iedge, masterCell, nG, 4);
                    surfaceFields::aux_rhoE[iedge][nG] = tempMaster * muMaster;
                    surfaceFields::aux_rhoE[iedge][nG + mathVar::nGauss + 1] = tempSlave * muSlave;
                }
            }
        }
    }

    std::tuple<double, double> getInternalValuesFromCalculatedArrays(int edge, int element, int nG, int valType)
    {
        //valType = 1(rho*mu), 2(rhou*mu), 3(rhov*mu), 4(rhoE*mu)
        bool isMaster(auxUlti::checkMaster(element, edge));
        int locationPlus(-1), locationMinus(-1);
        double valPlus(0.0), valMinus(0.0);
        if (isMaster)
        {
            locationPlus = nG;
            locationMinus = nG + mathVar::nGauss + 1;
        }
        else
        {
            locationPlus = nG + mathVar::nGauss + 1;
            locationMinus = nG;
        }

        switch (valType)
        {
        case 1:
        {
            valPlus = surfaceFields::aux_rho[edge][locationPlus];
            valMinus = surfaceFields::aux_rho[edge][locationMinus];
        }
        break;
        case 2:
        {
            valPlus = surfaceFields::aux_rhou[edge][locationPlus];
            valMinus = surfaceFields::aux_rhou[edge][locationMinus];
        }
        break;
        case 3:
        {
            valPlus = surfaceFields::aux_rhov[edge][locationPlus];
            valMinus = surfaceFields::aux_rhov[edge][locationMinus];
        }
        break;
        case 4:
        {
            valPlus = surfaceFields::aux_rhoE[edge][locationPlus];
            valMinus = surfaceFields::aux_rhoE[edge][locationMinus];
        }
        break;
        default:
            break;
        }
        return std::make_tuple(valPlus, valMinus);
    }

    /*Function calculates right hand side terms of all conservative variables at all order in all directions*/
    void CalcRHSTerm(int element, std::vector<double> &rhoRHSOx, std::vector<double> &rhoRHSOy, std::vector<double> &rhouRHSOx, std::vector<double> &rhouRHSOy, std::vector<double> &rhovRHSOx, std::vector<double> &rhovRHSOy, std::vector<double> &rhoERHSOx, std::vector<double> &rhoERHSOy)
    {
        std::vector<double>
            //vectors for volume integrals
            rhoVolIntOx(mathVar::orderElem + 1, 0.0),
            rhoVolIntOy(mathVar::orderElem + 1, 0.0),
            rhouVolIntOx(mathVar::orderElem + 1, 0.0),
            rhouVolIntOy(mathVar::orderElem + 1, 0.0),
            rhovVolIntOx(mathVar::orderElem + 1, 0.0),
            rhovVolIntOy(mathVar::orderElem + 1, 0.0),
            rhoEVolIntOx(mathVar::orderElem + 1, 0.0),
            rhoEVolIntOy(mathVar::orderElem + 1, 0.0),

            //vectors for surface integrals
            rhoSurfIntOx(mathVar::orderElem + 1, 0.0),
            rhoSurfIntOy(mathVar::orderElem + 1, 0.0),
            rhouSurfIntOx(mathVar::orderElem + 1, 0.0),
            rhouSurfIntOy(mathVar::orderElem + 1, 0.0),
            rhovSurfIntOx(mathVar::orderElem + 1, 0.0),
            rhovSurfIntOy(mathVar::orderElem + 1, 0.0),
            rhoESurfIntOx(mathVar::orderElem + 1, 0.0),
            rhoESurfIntOy(mathVar::orderElem + 1, 0.0);

        /*1. Calculate volume integral term*/
        process::auxEq::BR1::calcVolumeIntegralTerms(element, rhoVolIntOx, rhouVolIntOx, rhovVolIntOx, rhoEVolIntOx, rhoVolIntOy, rhouVolIntOy, rhovVolIntOy, rhoEVolIntOy);

        /*2. Calculate surface integral term*/
        process::auxEq::BR1::calcSurfaceIntegralTerms(element, rhoSurfIntOx, rhouSurfIntOx, rhovSurfIntOx, rhoESurfIntOx, rhoSurfIntOy, rhouSurfIntOy, rhovSurfIntOy, rhoESurfIntOy);

        for (int order = 0; order <= mathVar::orderElem; order++)
        {
            rhoRHSOx[order] = -rhoVolIntOx[order] + rhoSurfIntOx[order];
            rhoRHSOy[order] = -rhoVolIntOy[order] + rhoSurfIntOy[order];

            rhouRHSOx[order] = -rhouVolIntOx[order] + rhouSurfIntOx[order];
            rhouRHSOy[order] = -rhouVolIntOy[order] + rhouSurfIntOy[order];

            rhovRHSOx[order] = -rhovVolIntOx[order] + rhovSurfIntOx[order];
            rhovRHSOy[order] = -rhovVolIntOy[order] + rhovSurfIntOy[order];

            rhoERHSOx[order] = -rhoEVolIntOx[order] + rhoESurfIntOx[order];
            rhoERHSOy[order] = -rhoEVolIntOy[order] + rhoESurfIntOy[order];
        }
    }

    void calcVolumeIntegralTerms(int element, std::vector<double> &rhoVolIntX, std::vector<double> &rhouVolIntX, std::vector<double> &rhovVolIntX, std::vector<double> &rhoEVolIntX, std::vector<double> &rhoVolIntY, std::vector<double> &rhouVolIntY, std::vector<double> &rhovVolIntY, std::vector<double> &rhoEVolIntY)
    {
        std::vector<std::vector<double>> rhoGsVol(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
            rhouGsVol(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
            rhovGsVol(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
            rhoEGsVol(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));

        //Calculates Gauss matrix
        if (!flowProperties::massDiffusion)
        {
            //rho -------------------------------------------------------------------------------------------
            rhoGsVol = process::auxEq::BR1::getGaussMatrixOfConserVar(element, 1);
            for (int order = 1; order <= mathVar::orderElem; order++)
            {
                rhoVolIntX[order] = process::volumeInte(element, rhoGsVol, order, 1);
                rhoVolIntY[order] = process::volumeInte(element, rhoGsVol, order, 2);
            }
            rhoVolIntX[0] = 0;
            rhoVolIntY[0] = 0;
        }
        //rhou ------------------------------------------------------------------------------------------
        rhouGsVol = process::auxEq::BR1::getGaussMatrixOfConserVar(element, 2);
        //rhov ------------------------------------------------------------------------------------------
        rhovGsVol = process::auxEq::BR1::getGaussMatrixOfConserVar(element, 3);
        //rhoE ------------------------------------------------------------------------------------------
        rhoEGsVol = process::auxEq::BR1::getGaussMatrixOfConserVar(element, 4);

        for (int order = 1; order <= mathVar::orderElem; order++)
        {
            rhouVolIntX[order] = process::volumeInte(element, rhouGsVol, order, 1);
            rhovVolIntX[order] = process::volumeInte(element, rhovGsVol, order, 1);
            rhoEVolIntX[order] = process::volumeInte(element, rhoEGsVol, order, 1);

            rhouVolIntY[order] = process::volumeInte(element, rhouGsVol, order, 2);
            rhovVolIntY[order] = process::volumeInte(element, rhovGsVol, order, 2);
            rhoEVolIntY[order] = process::volumeInte(element, rhoEGsVol, order, 2);
        }
        rhouVolIntX[0] = 0;
        rhovVolIntX[0] = 0;
        rhoEVolIntX[0] = 0;

        rhouVolIntY[0] = 0;
        rhovVolIntY[0] = 0;
        rhoEVolIntY[0] = 0;
    }

    void calcSurfaceIntegralTerms(int element, std::vector<double> &rhoSurfIntX, std::vector<double> &rhouSurfIntX, std::vector<double> &rhovSurfIntX, std::vector<double> &rhoESurfIntX, std::vector<double> &rhoSurfIntY, std::vector<double> &rhouSurfIntY, std::vector<double> &rhovSurfIntY, std::vector<double> &rhoESurfIntY)
    {
        /*User's guide:
        Input array rhoRHSTerm, rhouRHSTerm, rhovRHSTerm, rhoERHSTerm have following form:
        - number of row: orderElem + 1*/
        int elemType(auxUlti::checkType(element)), edgeName(0);
        std::vector<std::vector<double>> rhoFluxX(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
            rhouFluxX(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
            rhovFluxX(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
            rhoEFluxX(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
            rhoFluxY(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
            rhouFluxY(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
            rhovFluxY(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
            rhoEFluxY(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0));

        std::vector<double> rhoFluxXTemp(mathVar::nGauss + 1, 0.0),
            rhouFluxXTemp(mathVar::nGauss + 1, 0.0),
            rhovFluxXTemp(mathVar::nGauss + 1, 0.0),
            rhoEFluxXTemp(mathVar::nGauss + 1, 0.0),
            rhoFluxYTemp(mathVar::nGauss + 1, 0.0),
            rhouFluxYTemp(mathVar::nGauss + 1, 0.0),
            rhovFluxYTemp(mathVar::nGauss + 1, 0.0),
            rhoEFluxYTemp(mathVar::nGauss + 1, 0.0);

        /*1. Calculate flux of conservative variables at all Gauss points on all faces of element*/
        process::auxEq::BR1::getGaussVectorOfConserVar(element, rhoFluxX, rhouFluxX, rhovFluxX, rhoEFluxX, rhoFluxY, rhouFluxY, rhovFluxY, rhoEFluxY);

        /*2. Calculates surface integrals of all conservative variables at all order*/
        for (int nface = 0; nface < elemType; nface++)
        {
            edgeName = meshVar::inedel[element][nface];
            for (int nG = 0; nG <= mathVar::nGauss; nG++)
            {
                rhoFluxXTemp[nG] = rhoFluxX[nG][nface];
                rhoFluxYTemp[nG] = rhoFluxY[nG][nface];

                rhouFluxXTemp[nG] = rhouFluxX[nG][nface];
                rhovFluxXTemp[nG] = rhovFluxX[nG][nface];
                rhoEFluxXTemp[nG] = rhoEFluxX[nG][nface];

                rhouFluxYTemp[nG] = rhouFluxY[nG][nface];
                rhovFluxYTemp[nG] = rhovFluxY[nG][nface];
                rhoEFluxYTemp[nG] = rhoEFluxY[nG][nface];
            }

            for (int order = 0; order <= mathVar::orderElem; order++)
            {
                rhoSurfIntX[order] += process::surfaceInte(element, edgeName, rhoFluxXTemp, order);
                rhoSurfIntY[order] += process::surfaceInte(element, edgeName, rhoFluxYTemp, order);

                rhouSurfIntX[order] += process::surfaceInte(element, edgeName, rhouFluxXTemp, order);
                rhovSurfIntX[order] += process::surfaceInte(element, edgeName, rhovFluxXTemp, order);
                rhoESurfIntX[order] += process::surfaceInte(element, edgeName, rhoEFluxXTemp, order);

                rhouSurfIntY[order] += process::surfaceInte(element, edgeName, rhouFluxYTemp, order);
                rhovSurfIntY[order] += process::surfaceInte(element, edgeName, rhovFluxYTemp, order);
                rhoESurfIntY[order] += process::surfaceInte(element, edgeName, rhoEFluxYTemp, order);
            }
        }
    }

    std::vector<std::vector<double>> getGaussMatrixOfConserVar(int element, int valType)
    {
        //rhoE used for solving auxilary equation is advective rhoE!!
        std::vector<std::vector<double>> GaussMatrix(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
        for (int na = 0; na <= mathVar::nGauss; na++)
        {
            for (int nb = 0; nb <= mathVar::nGauss; nb++)
            {
                switch (valType)
                {
                case 1:
                {
                    GaussMatrix[na][nb] = volumeFields::rhoVolGauss[element][na][nb];
                }
                break;
                case 2:
                {
                    GaussMatrix[na][nb] = volumeFields::rhouVolGauss[element][na][nb];
                }
                break;
                case 3:
                {
                    GaussMatrix[na][nb] = volumeFields::rhovVolGauss[element][na][nb];
                }
                break;
                case 4:
                {
                    if (!flowProperties::massDiffusion)
                    {
                        GaussMatrix[na][nb] = volumeFields::rhoEVolGauss[element][na][nb];
                    }
                    else {
                        GaussMatrix[na][nb]=volumeFields::rhoVolGauss[element][na][nb]*material::Cv*volumeFields::T[element][na][nb]+0.5*(pow(volumeFields::rhouVolGauss[element][na][nb],2)+pow(volumeFields::rhovVolGauss[element][na][nb],2))/volumeFields::rhoVolGauss[element][na][nb];
                    }
                }
                break;
                default:
                    break;
                }
            }
        }
        return GaussMatrix;
    }

    std::vector<std::vector<double>> getVectorOfConserVarFluxesAtInternal(int edge, int element, int nG, double nx, double ny)
    {
        std::vector<std::vector<double>> gaussVector(4, std::vector<double>(2, 0.0)); //columns 0, 1 are Ox, Oy values
        double rhoM(0.0), rhoP(0.0), rhouM(0.0), rhouP(0.0), rhovM(0.0), rhovP(0.0), TP(0.0), TM(0.0), rhoEM(0.0), rhoEP(0.0);
        std::tie(rhoP, rhoM) = auxUlti::getUAtInterfaces(edge, element, nG, 1);
        std::tie(rhouP, rhouM) = auxUlti::getUAtInterfaces(edge, element, nG, 2);
        std::tie(rhovP, rhovM) = auxUlti::getUAtInterfaces(edge, element, nG, 3);
        if (flowProperties::massDiffusion)
        {
            std::tie(TP, TM)=auxUlti::getTAtInterfaces(edge,element,nG);
            rhoEP=rhoP*material::Cv*TP+0.5*(rhouP*rhouP+rhovP*rhovP)/rhoP;
            rhoEM=rhoM*material::Cv*TM+0.5*(rhouM*rhouM+rhovM*rhovM)/rhoM;
        }
        else {
            std::tie(rhoEP, rhoEM) = auxUlti::getUAtInterfaces(edge, element, nG, 4);
            gaussVector[0][0] = math::numericalFluxes::auxFlux(rhoM, rhoP, nx);
            gaussVector[0][1] = math::numericalFluxes::auxFlux(rhoM, rhoP, ny);
        }

        gaussVector[1][0] = math::numericalFluxes::auxFlux(rhouM, rhouP, nx);
        gaussVector[1][1] = math::numericalFluxes::auxFlux(rhouM, rhouP, ny);
        gaussVector[2][0] = math::numericalFluxes::auxFlux(rhovM, rhovP, nx);
        gaussVector[2][1] = math::numericalFluxes::auxFlux(rhovM, rhovP, ny);
        gaussVector[3][0] = math::numericalFluxes::auxFlux(rhoEM, rhoEP, nx);
        gaussVector[3][1] = math::numericalFluxes::auxFlux(rhoEM, rhoEP, ny);


        return gaussVector;
    }

    void getGaussVectorOfConserVar(int element, std::vector<std::vector<double>> &rhoFluxX, std::vector<std::vector<double>> &rhouFluxX, std::vector<std::vector<double>> &rhovFluxX, std::vector<std::vector<double>> &rhoEFluxX, std::vector<std::vector<double>> &rhoFluxY, std::vector<std::vector<double>> &rhouFluxY, std::vector<std::vector<double>> &rhovFluxY, std::vector<std::vector<double>> &rhoEFluxY)
    {
        /*User's guide:
        Input array rhoFlux, rhouFlux, rhovFlux, rhoEFlux have following form:
        - number of row: nGauss + 1*/

        int elemType(auxUlti::checkType(element)), edgeName(0);
        int faceBcType(0);
        std::vector<std::vector<double>> gaussVector(4, std::vector<double>(2, 0.0));

        for (int nface = 0; nface < elemType; nface++)
        {
            edgeName = meshVar::inedel[element][nface];
            faceBcType = auxUlti::getBCType(edgeName);
            double nx(auxUlti::getNormVectorComp(element, edgeName, 1)), ny(auxUlti::getNormVectorComp(element, edgeName, 2));
            if (faceBcType == 0)  //internal edge
            {
                for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
                {
                    gaussVector = process::auxEq::BR1::getVectorOfConserVarFluxesAtInternal(edgeName, element, nGauss, nx, ny);
                    rhoFluxX[nGauss][nface] = gaussVector[0][0];
                    rhouFluxX[nGauss][nface] = gaussVector[1][0];
                    rhovFluxX[nGauss][nface] = gaussVector[2][0];
                    rhoEFluxX[nGauss][nface] = gaussVector[3][0];

                    rhoFluxY[nGauss][nface] = gaussVector[0][1];
                    rhouFluxY[nGauss][nface] = gaussVector[1][1];
                    rhovFluxY[nGauss][nface] = gaussVector[2][1];
                    rhoEFluxY[nGauss][nface] = gaussVector[3][1];
                }
            }
            else  //boundary edge
            {
                for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
                {
                    gaussVector = auxEqBCsImplement(element, edgeName, nGauss);
                    rhoFluxX[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[0][0],gaussVector[0][1],nx);
                    rhouFluxX[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[1][0],gaussVector[1][1],nx);
                    rhovFluxX[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[2][0],gaussVector[2][1],nx);
                    rhoEFluxX[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[3][0],gaussVector[3][1],nx);

                    rhoFluxY[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[0][0],gaussVector[0][1],ny);
                    rhouFluxY[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[1][0],gaussVector[1][1],ny);
                    rhovFluxY[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[2][0],gaussVector[2][1],ny);
                    rhoEFluxY[nGauss][nface] = math::numericalFluxes::auxFlux(gaussVector[3][0],gaussVector[3][1],ny);
                }
            }
        }
    }
    }

    namespace massDiffusion
    {
        namespace BR1 {
        void solveDivRho()
        {
            std::vector<double> rhoRHSTermOxDir(mathVar::orderElem + 1, 0.0),
                rhoRHSTermOyDir(mathVar::orderElem + 1, 0.0);

            for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
            {
                //2) Calculate Right hand side terms
                process::auxEq::massDiffusion::BR1::CalcRHSTerm(nelement, rhoRHSTermOxDir, rhoRHSTermOyDir);

                //3) Solve for div(rho)
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    //Ox direction
                    BR1Vars::rhoX[nelement][iorder] = rhoRHSTermOxDir[iorder] / stiffMatrixCoeffs[nelement][iorder];

                    //Oy direction
                    BR1Vars::rhoY[nelement][iorder] = rhoRHSTermOyDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
                }
            }
        }

        void CalcRHSTerm(int element, std::vector<double> &rhoRHSOx, std::vector<double> &rhoRHSOy)
        {
            std::vector<double>
                //vectors for volume integrals
                rhoVolIntOx(mathVar::orderElem + 1, 0.0),
                rhoVolIntOy(mathVar::orderElem + 1, 0.0),

                //vectors for surface integrals
                rhoSurfIntOx(mathVar::orderElem + 1, 0.0),
                rhoSurfIntOy(mathVar::orderElem + 1, 0.0);

            /*1. Calculate volume integral term*/
            process::auxEq::massDiffusion::BR1::calcVolumeIntegralTerms(element, rhoVolIntOx, rhoVolIntOy);

            /*2. Calculate surface integral term*/
            process::auxEq::massDiffusion::BR1::calcSurfaceIntegralTerms(element, rhoSurfIntOx, rhoSurfIntOy);

            for (int order = 0; order <= mathVar::orderElem; order++)
            {
                rhoRHSOx[order] = -rhoVolIntOx[order] + rhoSurfIntOx[order];
                rhoRHSOy[order] = -rhoVolIntOy[order] + rhoSurfIntOy[order];
            }
        }

        void calcVolumeIntegralTerms(int element, std::vector<double> &rhoVolIntX, std::vector<double> &rhoVolIntY)
        {
            std::vector<std::vector<double>> rhoGsVol(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));

            //Calculates Gauss matrix
            //rho -------------------------------------------------------------------------------------------
            for (int na = 0; na <= mathVar::nGauss; na++)
            {
                for (int nb = 0; nb <= mathVar::nGauss; nb++)
                {
                    rhoGsVol[na][nb]=volumeFields::rhoVolGauss[element][na][nb];
                }
            }

            for (int order = 1; order <= mathVar::orderElem; order++)
            {
                rhoVolIntX[order] = process::volumeInte(element, rhoGsVol, order, 1);
                rhoVolIntY[order] = process::volumeInte(element, rhoGsVol, order, 2);
            }
            rhoVolIntX[0] = 0;
            rhoVolIntY[0] = 0;
        }

        void calcSurfaceIntegralTerms(int element, std::vector<double> &rhoSurfIntX, std::vector<double> &rhoSurfIntY)
        {
            int elemType(auxUlti::checkType(element)), edgeName(0);
            std::vector<std::vector<double>> rhoFluxX(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
                rhoFluxY(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0));

            std::vector<double> rhoFluxXTemp(mathVar::nGauss + 1, 0.0),
                rhoFluxYTemp(mathVar::nGauss + 1, 0.0);

            /*1. Calculate flux of rho at all Gauss points on all faces of element*/
            process::auxEq::massDiffusion::BR1::getGaussVectorOfRho(element, rhoFluxX, rhoFluxY);

            /*2. Calculates surface integrals of rho at all order*/
            for (int nface = 0; nface < elemType; nface++)
            {
                edgeName = meshVar::inedel[element][nface];
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    rhoFluxXTemp[nG] = rhoFluxX[nG][nface];
                    rhoFluxYTemp[nG] = rhoFluxY[nG][nface];
                }

                for (int order = 0; order <= mathVar::orderElem; order++)
                {
                    rhoSurfIntX[order] += process::surfaceInte(element, edgeName, rhoFluxXTemp, order);
                    rhoSurfIntY[order] += process::surfaceInte(element, edgeName, rhoFluxYTemp, order);
                }
            }
        }

        void getGaussVectorOfRho(int element, std::vector<std::vector<double>> &rhoFluxX, std::vector<std::vector<double>> &rhoFluxY)
        {
            int elemType(auxUlti::checkType(element)), edgeName(0);
            int faceBcType(0);
            std::vector<std::vector<double>> gaussVector(4, std::vector<double>(2, 0.0));

            for (int nface = 0; nface < elemType; nface++)
            {
                edgeName = meshVar::inedel[element][nface];
                faceBcType = auxUlti::getBCType(edgeName);
                double nx(auxUlti::getNormVectorComp(element, edgeName, 1)), ny(auxUlti::getNormVectorComp(element, edgeName, 2)), rhoP(0.0), rhoM(0.0);
                if (faceBcType == 0)  //internal edge
                {
                    for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
                    {
                        std::tie(rhoP, rhoM) = auxUlti::getUAtInterfaces(edgeName, element, nGauss, 1);
                        rhoFluxX[nGauss][nface] = math::numericalFluxes::auxFlux(rhoM, rhoP, nx);
                        rhoFluxY[nGauss][nface] = math::numericalFluxes::auxFlux(rhoM, rhoP, ny);
                    }
                }
                else  //boundary edge
                {
                    for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
                    {
                        std::tie(rhoP, rhoM)=rhoBCsImplement(element, edgeName, nGauss);
                        rhoFluxX[nGauss][nface]=math::numericalFluxes::auxFlux(rhoP, rhoM, nx);
                        rhoFluxY[nGauss][nface]=math::numericalFluxes::auxFlux(rhoP, rhoM, ny);
                    }
                }
            }
        }
        }

        namespace BR2 {
        void solveDivRho()
        {
            for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
            {
                //1) Calculate volume divU
                process::auxEq::massDiffusion::BR2::calcVolDivRho(nelement);

                //2) Calculate suface divU
                process::auxEq::massDiffusion::BR2::calcSurDivRho(nelement);
            }
        }

        void calcVolDivRho(int element)
        {
            std::vector<std::vector<double>> gaussMX(mathVar::nGauss+1,std::vector<double>(mathVar::nGauss+1,0.0)),
                    gaussMY(mathVar::nGauss+1,std::vector<double>(mathVar::nGauss+1,0.0));
            std::vector<double>rhoC(mathVar::orderElem+1,0.0);
            int elemType(auxUlti::checkType(element));
            double Bn(0.0),dBmX(0.0),dBmY(0.0),
                    GammaRhoX(0.0), GammaRhoY(0.0),
                    GammaX(0.0), GammaY(0.0), theta1(1.0), theta2(1.0);
            rhoC=auxUlti::getElementConserValuesOfOrder(element,1);
            for (int n=0;n<=mathVar::orderElem;n++) {
                GammaRhoX=0; GammaRhoY=0;

                for (int m=0;m<=mathVar::orderElem;m++) {
                    if (m==0)
                    {
                        theta1=1.0;
                        theta2=1.0;
                    }
                    else {
                        theta1=theta1Arr[element];
                        theta2=theta2Arr[element];
                    }

                    if (elemType==4)
                    {
                        for (int na=0;na<=mathVar::nGauss;na++) {
                            for (int nb=0;nb<=mathVar::nGauss;nb++) {
                                Bn=mathVar::BPts_Quad[n][na][nb];
                                dBmX=math::Calc_dBxdBy(element, m, na, nb, 1);
                                dBmY=math::Calc_dBxdBy(element, m, na, nb, 2);
                                gaussMX[na][nb]=Bn*dBmX;
                                gaussMY[na][nb]=Bn*dBmY;
                            }
                        }
                    }
                    else if (elemType==3)
                    {
                        for (int na=0;na<=mathVar::nGauss;na++) {
                            for (int nb=0;nb<=mathVar::nGauss;nb++) {
                                Bn=mathVar::BPts_Tri[n][na][nb];
                                dBmX=math::Calc_dBxdBy(element, m, na, nb, 1);
                                dBmY=math::Calc_dBxdBy(element, m, na, nb, 2);
                                gaussMX[na][nb]=Bn*dBmX;
                                gaussMY[na][nb]=Bn*dBmY;
                            }
                        }
                    }

                    GammaX=math::volumeInte(gaussMX,element);
                    GammaY=math::volumeInte(gaussMY,element);

                    GammaRhoX+=GammaX*rhoC[m]*theta1*theta2;
                    GammaRhoY+=GammaY*rhoC[m]*theta1*theta2;
                }
                double stiffCoeff(stiffMatrixCoeffs[element][n]);
                BR2Vars::rhoXVol[element][n]=GammaRhoX/stiffCoeff;
                BR2Vars::rhoYVol[element][n]=GammaRhoY/stiffCoeff;
            }
        }

        void calcSurDivRho(int element)
        {
            int elemType(auxUlti::checkType(element)), edgeId(0), BCType(-1);
            double a(0.0), b(0.0),BVar(0.0), nx(0.0), ny(0.0), UP(0.0), UM(0.0), BRConst(4.0),
                    dRhoXVar(0.0), dRhoYVar(0.0), rhoP(0.0), rhoM(0.0);
            bool isMaster;
            std::vector<double> vectorRhoX(mathVar::nGauss+1,0.0),
                    vectorRhoY(mathVar::nGauss+1,0.0),

                    vectordRhoXVol_org(mathVar::orderElem+1,0.0),
                    vectordRhoYVol_org(mathVar::orderElem+1,0.0);

            for (int iorder=0;iorder<=mathVar::orderElem;iorder++) {
                vectordRhoXVol_org[iorder]=BR2Vars::rhoXVol[element][iorder];
                vectordRhoYVol_org[iorder]=BR2Vars::rhoYVol[element][iorder];
            }

            for (int iedge=0;iedge<elemType;iedge++) {
                edgeId=meshVar::inedel[element][iedge];
                BCType=auxUlti::getBCType(edgeId);
                isMaster=auxUlti::checkMaster(element,edgeId);
                nx=auxUlti::getNormVectorComp(element,edgeId,1);
                ny=auxUlti::getNormVectorComp(element,edgeId,2);
                if (BCType!=0)
                {
                    for (int nG=0;nG<=mathVar::nGauss;nG++) {
                        std::tie(a,b)=auxUlti::getGaussSurfCoor(edgeId,element,nG);
                        std::tie(rhoP,rhoM)=rhoBCsImplement(element, edgeId, nG);

                        //Rho
                        vectorRhoX[nG]=0.5*(rhoM-rhoP)*nx;
                        vectorRhoY[nG]=0.5*(rhoM-rhoP)*ny;
                    }
                }
                else {
                    for (int nG=0;nG<=mathVar::nGauss;nG++) {
                        std::tie(a,b)=auxUlti::getGaussSurfCoor(edgeId,element,nG);

                        //Rho
                        std::tie(UP,UM)=process::getInternalValuesFromCalculatedArrays(edgeId,element,nG,1);
                        vectorRhoX[nG]=0.5*(UM-UP)*nx;
                        vectorRhoY[nG]=0.5*(UM-UP)*ny;
                    }
                }

                for (int iorder=0;iorder<=mathVar::orderElem;iorder++) {
                    double stiffCoeff(stiffMatrixCoeffs[element][iorder]);
                    for (int nG=0;nG<=mathVar::nGauss;nG++) {
                        std::tie(a,b)=auxUlti::getGaussSurfCoor(edgeId,element,nG);
                        math::basisFc(a,b,elemType);
                        BVar=mathVar::B[iorder];
                        //Rho
                        vectorRhoX[nG]*=BVar;
                        vectorRhoY[nG]*=BVar;
                    }
                    dRhoXVar=math::surfaceInte(vectorRhoX,edgeId)/stiffCoeff;
                    dRhoYVar=math::surfaceInte(vectorRhoY,edgeId)/stiffCoeff;

                    if (isMaster)
                    {
                        BR2Vars::rhoXSurMaster[edgeId][iorder]=BRConst*dRhoXVar+vectordRhoXVol_org[iorder];
                        BR2Vars::rhoYSurMaster[edgeId][iorder]=BRConst*dRhoYVar+vectordRhoYVol_org[iorder];
                    }
                    else {
                        BR2Vars::rhoXSurSlave[edgeId][iorder]=BRConst*dRhoXVar+vectordRhoXVol_org[iorder];
                        BR2Vars::rhoYSurSlave[edgeId][iorder]=BRConst*dRhoYVar+vectordRhoYVol_org[iorder];
                    }

                    BR2Vars::rhoXVol[element][iorder]+=dRhoXVar;
                    BR2Vars::rhoYVol[element][iorder]+=dRhoYVar;
                }
            }
        }
        }
    }

	}//end namespace auxEq

	namespace NSFEq
	{
        void calcFinvFvisAtInterface()
		{
			//mu included
			int masterCell(-1), slaveCell(-1), bcGrp(0);
			double uMaster(0.0), vMaster(0.0), totalEMaster(0.0), TMaster(0.0), pMaster(0.0),
				uSlave(0.0), vSlave(0.0), totalESlave(0.0), TSlave(0.0), pSlave(0.0), eMaster(0.0), eSlave(0.0),
                uMagM(0.0), uMagP(0.0), aM(0.0), aP(0.0),
                    umMaster(0.0), umSlave(0.0), vmMaster(0.0), vmSlave(0.0), muMaster(0.0), muSlave(0.0);
			std::vector<double> UMaster(4, 0.0), dUXMaster(4, 0.0), dUYMaster(4, 0.0),
				USlave(4, 0.0), dUXSlave(4, 0.0), dUYSlave(4, 0.0),
                CArray(mathVar::nGauss + 1, 0.0), vectorn(2, 0.0); // BetaArray(mathVar::nGauss + 1, 0.0),
			/*StressHeat matrix has form:
			[tauXx		tauXy		Qx]
			[tauYx		tauYy		Qy]
			*/
			std::vector<std::vector<double>> StressHeatP(2, std::vector<double>(3, 0.0));
			std::vector<std::vector<double>> StressHeatM(2, std::vector<double>(3, 0.0));
            for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
			{
				bcGrp = auxUlti::getBCType(iedge);
				if (bcGrp == 0)
				{
					std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
					vectorn[0] = auxUlti::getNormVectorComp(masterCell, iedge, 1);
					vectorn[1] = auxUlti::getNormVectorComp(masterCell, iedge, 2);
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
					{
                        std::tie(TMaster,TSlave)=auxUlti::getTAtInterfaces(iedge,masterCell,nG);
                        muMaster=math::CalcVisCoef(TMaster);
                        muSlave=math::CalcVisCoef(TSlave);
                        for (int i = 0; i < 4; i++)
						{
                            //std::tie(UMaster[i], USlave[i]) = math::internalSurfaceValue(iedge, masterCell, nG, i + 1, 2);
                            std::tie(UMaster[i], USlave[i]) = auxUlti::getUAtInterfaces(iedge, masterCell, nG, i + 1);
                            //d(U)/dx, d(U)/dy
                            std::tie(dUXMaster[i], dUXSlave[i]) = math::internalSurfaceDerivativeValue(iedge, masterCell, nG, i+1, 1);
                            std::tie(dUYMaster[i], dUYSlave[i]) = math::internalSurfaceDerivativeValue(iedge, masterCell, nG, i+1, 2);
                            dUXMaster[i]*=muMaster;
                            dUYMaster[i]*=muMaster;
                            dUXSlave[i]*=muSlave;
                            dUYSlave[i]*=muSlave;
						}

						uMaster = (UMaster[1] / UMaster[0]);
						uSlave = (USlave[1] / USlave[0]);

						vMaster = (UMaster[2] / UMaster[0]);
						vSlave = (USlave[2] / USlave[0]);

						totalEMaster = (UMaster[3] / UMaster[0]);
						totalESlave = (USlave[3] / USlave[0]);

                        /*calculate P*/
						pMaster = math::CalcP(TMaster, UMaster[0]);
						pSlave = math::CalcP(TSlave, USlave[0]);
						eMaster = material::Cv*TMaster;
						eSlave = material::Cv*TSlave;

						/*INVISCID TERMS*/
                        /*- Calculate total velocity components u_m*/
                        if (flowProperties::massDiffusion)
                        {
                            umMaster=math::massDiffusionFncs::calcTotalVelocity(UMaster[0],uMaster,dUXMaster[0]);
                            vmMaster=math::massDiffusionFncs::calcTotalVelocity(UMaster[0],vMaster,dUYMaster[0]);
                            umSlave=math::massDiffusionFncs::calcTotalVelocity(USlave[0],uSlave,dUXSlave[0]);
                            vmSlave=math::massDiffusionFncs::calcTotalVelocity(USlave[0],vSlave,dUYSlave[0]);
                        }
                        else
                        {
                            umMaster=uMaster;
                            vmMaster=vMaster;
                            umSlave=uSlave;
                            vmSlave=vSlave;
                        }

						/*calculate velocity magnitude*/
                        uMagP = sqrt(pow(uMaster, 2) + pow(vMaster, 2));
                        uMagM = sqrt(pow(uSlave, 2) + pow(vSlave, 2));

						/*calculate speed of sound*/
						aP = math::CalcSpeedOfSound(TMaster);
						aM = math::CalcSpeedOfSound(TSlave);

						/*calculate constant for Lax-Friederich flux*/
						CArray[nG] = math::numericalFluxes::constantC(uMagP, uMagM, aP, aM); 

						/*calculate inviscid terms*/
                        std::tie(surfaceFields::invis_rhoX[iedge][nG], surfaceFields::invis_rhouX[iedge][nG], surfaceFields::invis_rhovX[iedge][nG], surfaceFields::invis_rhoEX[iedge][nG]) = math::inviscidTerms::calcInvisTermsFromPriVars(UMaster[0], uMaster, umMaster, vMaster, vmMaster, totalEMaster, pMaster, 1);
                        std::tie(surfaceFields::invis_rhoY[iedge][nG], surfaceFields::invis_rhouY[iedge][nG], surfaceFields::invis_rhovY[iedge][nG], surfaceFields::invis_rhoEY[iedge][nG]) = math::inviscidTerms::calcInvisTermsFromPriVars(UMaster[0], uMaster, umMaster, vMaster, vmMaster, totalEMaster, pMaster, 2);

                        std::tie(surfaceFields::invis_rhoX[iedge][nG + mathVar::nGauss + 1], surfaceFields::invis_rhouX[iedge][nG + mathVar::nGauss + 1], surfaceFields::invis_rhovX[iedge][nG + mathVar::nGauss + 1], surfaceFields::invis_rhoEX[iedge][nG + mathVar::nGauss + 1]) = math::inviscidTerms::calcInvisTermsFromPriVars(USlave[0], uSlave, umSlave, vSlave, vmSlave, totalESlave, pSlave, 1);
                        std::tie(surfaceFields::invis_rhoY[iedge][nG + mathVar::nGauss + 1], surfaceFields::invis_rhouY[iedge][nG + mathVar::nGauss + 1], surfaceFields::invis_rhovY[iedge][nG + mathVar::nGauss + 1], surfaceFields::invis_rhoEY[iedge][nG + mathVar::nGauss + 1]) = math::inviscidTerms::calcInvisTermsFromPriVars(USlave[0], uSlave, umSlave, vSlave, vmSlave, totalESlave, pSlave, 2);

						/*calculate viscous terms*/
                        StressHeatP = math::viscousTerms::calcStressTensorAndHeatFlux(UMaster, dUXMaster, dUYMaster, TMaster);
                        std::tie(surfaceFields::Vis_rhoX[iedge][nG], surfaceFields::Vis_rhouX[iedge][nG], surfaceFields::Vis_rhovX[iedge][nG], surfaceFields::Vis_rhoEX[iedge][nG]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uMaster, vMaster, dUXMaster[0]/UMaster[0], 1);
                        std::tie(surfaceFields::Vis_rhoY[iedge][nG], surfaceFields::Vis_rhouY[iedge][nG], surfaceFields::Vis_rhovY[iedge][nG], surfaceFields::Vis_rhoEY[iedge][nG]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uMaster, vMaster, dUYMaster[0]/UMaster[0], 2);

                        StressHeatM = math::viscousTerms::calcStressTensorAndHeatFlux(USlave, dUXSlave, dUYSlave, TSlave);
                        std::tie(surfaceFields::Vis_rhoX[iedge][nG + mathVar::nGauss + 1], surfaceFields::Vis_rhouX[iedge][nG + mathVar::nGauss + 1], surfaceFields::Vis_rhovX[iedge][nG + mathVar::nGauss + 1], surfaceFields::Vis_rhoEX[iedge][nG + mathVar::nGauss + 1]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uSlave, vSlave, dUXSlave[0]/USlave[0], 1);
                        std::tie(surfaceFields::Vis_rhoY[iedge][nG + mathVar::nGauss + 1], surfaceFields::Vis_rhouY[iedge][nG + mathVar::nGauss + 1], surfaceFields::Vis_rhovY[iedge][nG + mathVar::nGauss + 1], surfaceFields::Vis_rhoEY[iedge][nG + mathVar::nGauss + 1]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uSlave, vSlave, dUYSlave[0]/USlave[0], 2);
					
						/*calculate constant for diffusive flux*/
                        //BetaArray[nG] = math::numericalFluxes::constantBeta(uMagP, uMagM, UMaster[0], USlave[0], eMaster, eSlave, pMaster, pSlave, StressHeatP, StressHeatM, vectorn);
					}
					LxFConst[iedge] = *std::max_element(CArray.begin(), CArray.end());
                    //DiffusiveFluxConst[iedge] = *std::max_element(BetaArray.begin(), BetaArray.end());
				}
			}
		}

		void solveNSFEquation(int RKOrder)
		{
			std::vector<double> rhoError(meshVar::nelem2D, 1.0),
				rhouError(meshVar::nelem2D, 1.0),
				rhovError(meshVar::nelem2D, 1.0),
				rhoEError(meshVar::nelem2D, 1.0);

			//std::vector<std::vector<double>> StiffMatrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));
			std::vector<double>
				RHSTerm1(mathVar::orderElem + 1, 0.0),
				RHSTerm2(mathVar::orderElem + 1, 0.0),
				RHSTerm3(mathVar::orderElem + 1, 0.0),
				RHSTerm4(mathVar::orderElem + 1, 0.0),

				ddtRhoVector(mathVar::orderElem + 1, 0.0),
				ddtRhouVector(mathVar::orderElem + 1, 0.0),
				ddtRhovVector(mathVar::orderElem + 1, 0.0),
				ddtRhoEVector(mathVar::orderElem + 1, 0.0),

				rhoVectorN(mathVar::orderElem + 1, 0.0),
				rhouVectorN(mathVar::orderElem + 1, 0.0),
				rhovVectorN(mathVar::orderElem + 1, 0.0),
				rhoEVectorN(mathVar::orderElem + 1, 0.0),

				UnVector(mathVar::orderElem + 1, 0.0);

            for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
			{
				//2) Calculate Right hand side terms
				process::NSFEq::CalcRHSTerm(nelement, RHSTerm1, RHSTerm2, RHSTerm3, RHSTerm4);

                //3) Solve for time derivartives of conservative variables
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
				{
					ddtRhoVector[iorder] = RHSTerm1[iorder] / stiffMatrixCoeffs[nelement][iorder];
					ddtRhouVector[iorder] = RHSTerm2[iorder] / stiffMatrixCoeffs[nelement][iorder];
					ddtRhovVector[iorder] = RHSTerm3[iorder] / stiffMatrixCoeffs[nelement][iorder];
					ddtRhoEVector[iorder] = RHSTerm4[iorder] / stiffMatrixCoeffs[nelement][iorder];

					switch (RKOrder)
					{
					case 1:
					{
						rhoResArr[nelement][iorder] = (1.0 / 6.0)*ddtRhoVector[iorder];
						rhouResArr[nelement][iorder] = (1.0 / 6.0)*ddtRhouVector[iorder];
						rhovResArr[nelement][iorder] = (1.0 / 6.0)*ddtRhovVector[iorder];
						rhoEResArr[nelement][iorder] = (1.0 / 6.0)*ddtRhoEVector[iorder];
					}
					break;
					case 2:
					{
						rhoResArr[nelement][iorder] += (1.0 / 6.0)*ddtRhoVector[iorder];
						rhouResArr[nelement][iorder] += (1.0 / 6.0)*ddtRhouVector[iorder];
						rhovResArr[nelement][iorder] += (1.0 / 6.0)*ddtRhovVector[iorder];
						rhoEResArr[nelement][iorder] += (1.0 / 6.0)*ddtRhoEVector[iorder];
					}
					break;
					case 3:
					{
						rhoResArr[nelement][iorder] += (2.0 / 3.0)*ddtRhoVector[iorder];
						rhouResArr[nelement][iorder] += (2.0 / 3.0)*ddtRhouVector[iorder];
						rhovResArr[nelement][iorder] += (2.0 / 3.0)*ddtRhovVector[iorder];
						rhoEResArr[nelement][iorder] += (2.0 / 3.0)*ddtRhoEVector[iorder];
					}
					break;
					default:
						break;
					}
				}

				//4) Solve time marching
				//rho
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rho[nelement][order];
				}
				rhoVectorN = process::NSFEq::solveTimeMarching(nelement, ddtRhoVector, UnVector, RKOrder, 1);
				//rhou
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rhou[nelement][order];
				}
				rhouVectorN = process::NSFEq::solveTimeMarching(nelement, ddtRhouVector, UnVector, RKOrder, 2);
				//rhov
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rhov[nelement][order];
				}
				rhovVectorN = process::NSFEq::solveTimeMarching(nelement, ddtRhovVector, UnVector, RKOrder, 3);
				//rhoE
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rhoE[nelement][order];
				}
				rhoEVectorN = process::NSFEq::solveTimeMarching(nelement, ddtRhoEVector, UnVector, RKOrder, 4);

                //5) Save results to conservative variables arrays
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					rhoN[nelement][order] = rhoVectorN[order];
					rhouN[nelement][order] = rhouVectorN[order];
					rhovN[nelement][order] = rhovVectorN[order];
					rhoEN[nelement][order] = rhoEVectorN[order];
				}
			}
		}

		void updateVariables()
		{
            rho = rhoN;
            rhou = rhouN;
            rhov = rhovN;
            rhoE = rhoEN;
		}

		/*Function calculates right hand side terms of all conservative variables at ONLY one order*/
		void CalcRHSTerm(int element, std::vector<double> &term1RHS, std::vector<double> &term2RHS, std::vector<double> &term3RHS, std::vector<double> &term4RHS)
		{
			std::vector<double>
				VolIntTerm1(mathVar::orderElem + 1, 0.0),
				VolIntTerm2(mathVar::orderElem + 1, 0.0),
				VolIntTerm3(mathVar::orderElem + 1, 0.0),
				VolIntTerm4(mathVar::orderElem + 1, 0.0),

				SurfIntTerm1(mathVar::orderElem + 1, 0.0),
				SurfIntTerm2(mathVar::orderElem + 1, 0.0),
				SurfIntTerm3(mathVar::orderElem + 1, 0.0),
				SurfIntTerm4(mathVar::orderElem + 1, 0.0);

			/*Volume integral term===========================================================================*/
			process::NSFEq::calcVolumeIntegralTerms(element, VolIntTerm1, VolIntTerm2, VolIntTerm3, VolIntTerm4);
			
			/*Surface integral term===========================================================================*/
			process::NSFEq::calcSurfaceIntegralTerms(element, SurfIntTerm1, SurfIntTerm2, SurfIntTerm3, SurfIntTerm4);
			
            for (int order = 0; order <= mathVar::orderElem; order++)
			{
				term1RHS[order] = VolIntTerm1[order] - SurfIntTerm1[order];
				term2RHS[order] = VolIntTerm2[order] - SurfIntTerm2[order];
				term3RHS[order] = VolIntTerm3[order] - SurfIntTerm3[order];
				term4RHS[order] = VolIntTerm4[order] - SurfIntTerm4[order];
			}
		}

		/*Function calculates Inviscid terms at Gauss point (a, b)*/
		std::vector<std::vector<double>> calcGaussInviscidTerm(int element, int na, int nb)
		{
			/*InviscidTerm 2D array has 4 row 2 column:
			- column 1: Ox direction
			- column 2: Oy direction*/
			std::vector<std::vector<double>> InviscidTerm(4, std::vector<double>(2, 0.0));
			double
                rhoVal(volumeFields::rhoVolGauss[element][na][nb]),
                rhouVal(volumeFields::rhouVolGauss[element][na][nb]),
                rhovVal(volumeFields::rhovVolGauss[element][na][nb]),
                rhoEVal(volumeFields::rhoEVolGauss[element][na][nb]),
				uVal(0.0),
				vVal(0.0),
				pVal(0.0),
                umVal(0.0),
                vmVal(0.0),
                totalE(0.0),
                a(0.0),b(0.0),muVal(0.0);

            std::tie(a,b)=auxUlti::getGaussCoor(na,nb);
			uVal = rhouVal / rhoVal;
			vVal = rhovVal / rhoVal;
			totalE = rhoEVal / rhoVal;
            pVal = math::CalcP(volumeFields::T[element][na][nb], rhoVal);
            muVal=math::CalcVisCoef(volumeFields::T[element][na][nb]);

            if (flowProperties::massDiffusion)
            {
                if (systemVar::auxVariables==1)
                {
                    umVal=math::massDiffusionFncs::calcTotalVelocity(rhoVal,uVal,math::pointAuxValue(element,a,b,1,1)*muVal);
                    vmVal=math::massDiffusionFncs::calcTotalVelocity(rhoVal,vVal,math::pointAuxValue(element,a,b,1,2)*muVal);
                }
                else if (systemVar::auxVariables==2)
                {
                    umVal=math::massDiffusionFncs::calcTotalVelocity(rhoVal,uVal,math::BR2Fncs::pointAuxValue_vol(element,a,b,1,1)*muVal);
                    vmVal=math::massDiffusionFncs::calcTotalVelocity(rhoVal,vVal,math::BR2Fncs::pointAuxValue_vol(element,a,b,1,2)*muVal);
                }
            }
            else {
                umVal=uVal;
                vmVal=vVal;
            }

            /*1. Ox direction*/
            std::tie(InviscidTerm[0][0], InviscidTerm[1][0], InviscidTerm[2][0], InviscidTerm[3][0]) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoVal, uVal, umVal, vVal, vmVal, totalE, pVal, 1);

			/*2. Oy direction*/
            std::tie(InviscidTerm[0][1], InviscidTerm[1][1], InviscidTerm[2][1], InviscidTerm[3][1]) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoVal, uVal, umVal, vVal, vmVal, totalE, pVal, 2);

			return InviscidTerm;
		}

		/*Function calculates Viscous terms at Gauss point (a, b)*/
		std::vector<std::vector<double>> calcGaussViscousTerm(int element, int na, int nb)
		{
            double uVal(0.0), vVal(0.0), a(0.0), b(0.0), muVal(0.0);
			/*ViscousTerm 2D array has 4 row 2 column:
			- column 1: Ox direction
			- column 2: Oy direction*/
			std::vector<std::vector<double>> ViscousTerm(4, std::vector<double>(2, 0.0));
			std::vector<std::vector<double>> StressHeatFlux(2, std::vector<double>(3, 0.0));
			std::vector<double> vectorU(4, 0.0);
			std::vector<double> vectordUx(4, 0.0);
			std::vector<double> vectordUy(4, 0.0);
			std::tie(a, b) = auxUlti::getGaussCoor(na, nb);

			/*calculate conservative and derivative variables*/
            vectorU[0] = volumeFields::rhoVolGauss[element][na][nb];
            vectorU[1] = volumeFields::rhouVolGauss[element][na][nb];
            vectorU[2] = volumeFields::rhovVolGauss[element][na][nb];
            vectorU[3] = volumeFields::rhoEVolGauss[element][na][nb];

            uVal=vectorU[1]/vectorU[0];
            vVal=vectorU[2]/vectorU[0];
            muVal=math::CalcVisCoef(volumeFields::T[element][na][nb]);

            /*
            for (int i = 0; i < 4; i++)
			{
                vectordUx[i] = math::pointAuxValue(element, a, b, i + 1, 1)*muVal;
                vectordUy[i] = math::pointAuxValue(element, a, b, i + 1, 2)*muVal;
            }
            */
            vectordUx=math::pointSVars(0,element,a,b,1,1);
            vectordUy=math::pointSVars(0,element,a,b,2,1);
            for (int i = 0; i < 4; i++)
            {
                vectordUx[i] *= muVal;
                vectordUy[i] *= muVal;
            }

			/*calculate stresses and heat fluxes*/
            StressHeatFlux = math::viscousTerms::calcStressTensorAndHeatFlux(vectorU, vectordUx, vectordUy, volumeFields::T[element][na][nb]);
            std::tie(ViscousTerm[0][0], ViscousTerm[1][0], ViscousTerm[2][0], ViscousTerm[3][0]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatFlux, uVal, vVal, vectordUx[0]/vectorU[0], 1);
            std::tie(ViscousTerm[0][1], ViscousTerm[1][1], ViscousTerm[2][1], ViscousTerm[3][1]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatFlux, uVal, vVal, vectordUy[0]/vectorU[0], 2);
			return ViscousTerm;
		}

        std::tuple<double, double> getFinvFvisAtInterfaces(int edge, int element, int nG, int mod, int direction, int valType)
        {
            /*
            -mod: 1 for inviscid, 2 for viscous
            -direction: 1 for Ox, 2 for Oy
            */
            bool isMaster(auxUlti::checkMaster(element, edge));
            int locationPlus(-1), locationMinus(-1);
            double valPlus(0.0), valMinus(0.0);
            if (isMaster)
            {
                locationPlus = nG;
                locationMinus = nG + mathVar::nGauss + 1;
            }
            else
            {
                locationPlus = nG + mathVar::nGauss + 1;
                locationMinus = nG;
            }

            switch (mod)
            {
            case 1: //inviscid
            {
                switch (direction)
                {
                case 1:
                {
                    switch (valType)
                    {
                    case 1:
                    {
                        valPlus = surfaceFields::invis_rhoX[edge][locationPlus];
                        valMinus = surfaceFields::invis_rhoX[edge][locationMinus];
                    }
                    break;
                    case 2:
                    {
                        valPlus = surfaceFields::invis_rhouX[edge][locationPlus];
                        valMinus = surfaceFields::invis_rhouX[edge][locationMinus];
                    }
                    break;
                    case 3:
                    {
                        valPlus = surfaceFields::invis_rhovX[edge][locationPlus];
                        valMinus = surfaceFields::invis_rhovX[edge][locationMinus];
                    }
                    break;
                    case 4:
                    {
                        valPlus = surfaceFields::invis_rhoEX[edge][locationPlus];
                        valMinus = surfaceFields::invis_rhoEX[edge][locationMinus];
                    }
                    break;
                    default:
                        break;
                    }
                }
                break;
                case 2:
                {
                    switch (valType)
                    {
                    case 1:
                    {
                        valPlus = surfaceFields::invis_rhoY[edge][locationPlus];
                        valMinus = surfaceFields::invis_rhoY[edge][locationMinus];
                    }
                    break;
                    case 2:
                    {
                        valPlus = surfaceFields::invis_rhouY[edge][locationPlus];
                        valMinus = surfaceFields::invis_rhouY[edge][locationMinus];
                    }
                    break;
                    case 3:
                    {
                        valPlus = surfaceFields::invis_rhovY[edge][locationPlus];
                        valMinus = surfaceFields::invis_rhovY[edge][locationMinus];
                    }
                    break;
                    case 4:
                    {
                        valPlus = surfaceFields::invis_rhoEY[edge][locationPlus];
                        valMinus = surfaceFields::invis_rhoEY[edge][locationMinus];
                    }
                    break;
                    default:
                        break;
                    }
                }
                break;
                default:
                    break;
                }
            }
                break;
            case 2: //viscous
            {
                switch (direction)
                {
                case 1:
                {
                    switch (valType)
                    {
                    case 1:
                    {
                        valPlus = surfaceFields::Vis_rhoX[edge][locationPlus];
                        valMinus = surfaceFields::Vis_rhoX[edge][locationMinus];
                    }
                    break;
                    case 2:
                    {
                        valPlus = surfaceFields::Vis_rhouX[edge][locationPlus];
                        valMinus = surfaceFields::Vis_rhouX[edge][locationMinus];
                    }
                    break;
                    case 3:
                    {
                        valPlus = surfaceFields::Vis_rhovX[edge][locationPlus];
                        valMinus = surfaceFields::Vis_rhovX[edge][locationMinus];
                    }
                    break;
                    case 4:
                    {
                        valPlus = surfaceFields::Vis_rhoEX[edge][locationPlus];
                        valMinus = surfaceFields::Vis_rhoEX[edge][locationMinus];
                    }
                    break;
                    default:
                        break;
                    }
                }
                break;
                case 2:
                {
                    switch (valType)
                    {
                    case 1:
                    {
                        valPlus = surfaceFields::Vis_rhoY[edge][locationPlus];
                        valMinus = surfaceFields::Vis_rhoY[edge][locationMinus];
                    }
                    break;
                    case 2:
                    {
                        valPlus = surfaceFields::Vis_rhouY[edge][locationPlus];
                        valMinus = surfaceFields::Vis_rhouY[edge][locationMinus];
                    }
                    break;
                    case 3:
                    {
                        valPlus = surfaceFields::Vis_rhovY[edge][locationPlus];
                        valMinus = surfaceFields::Vis_rhovY[edge][locationMinus];
                    }
                    break;
                    case 4:
                    {
                        valPlus = surfaceFields::Vis_rhoEY[edge][locationPlus];
                        valMinus = surfaceFields::Vis_rhoEY[edge][locationMinus];
                    }
                    break;
                    default:
                        break;
                    }
                }
                break;
                default:
                    break;
                }
            }
                break;
            default:
                break;
            }
            return std::make_tuple(valPlus, valMinus);
        }

		void calcVolumeIntegralTerms(int element, std::vector<double> &VolIntTerm1, std::vector<double> &VolIntTerm2, std::vector<double> &VolIntTerm3, std::vector<double> &VolIntTerm4)
		{
			/*User's guide:
			All input array have form:
			- number of rows: orderElem*/

			std::vector<std::vector<double>> ViscousTerms(4, std::vector<double>(2, 0.0));
			std::vector<std::vector<double>> InviscidTerms(4, std::vector<double>(2, 0.0));
			//std::vector<double> VolInt(4, 0.0);

			std::vector<std::vector<double>>
                GsVolX1(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
                GsVolY1(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),

                GsVolX2(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
                GsVolY2(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),

                GsVolX3(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
                GsVolY3(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),

                GsVolX4(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
                GsVolY4(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));

            for (int na = 0; na <= mathVar::nGauss; na++)
			{
                for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					//std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
					/*A INVISCID TERMS*/
					InviscidTerms = NSFEq::calcGaussInviscidTerm(element, na, nb);
					/*A1. Inviscid term on Ox direction*/
                    GsVolX1[na][nb] = InviscidTerms[0][0];
                    GsVolX2[na][nb] = InviscidTerms[1][0];
                    GsVolX3[na][nb] = InviscidTerms[2][0];
                    GsVolX4[na][nb] = InviscidTerms[3][0];
					/*A2. Inviscid term on Oy direction*/
                    GsVolY1[na][nb] = InviscidTerms[0][1];
                    GsVolY2[na][nb] = InviscidTerms[1][1];
                    GsVolY3[na][nb] = InviscidTerms[2][1];
                    GsVolY4[na][nb] = InviscidTerms[3][1];

					/*B VISCOUS TERMS*/
					ViscousTerms = NSFEq::calcGaussViscousTerm(element, na, nb);
					/*B1. Viscous term on Ox direction*/
                    GsVolX1[na][nb] += ViscousTerms[0][0];
                    GsVolX2[na][nb] += ViscousTerms[1][0];
                    GsVolX3[na][nb] += ViscousTerms[2][0];
                    GsVolX4[na][nb] += ViscousTerms[3][0];
					/*B2. Viscous term on Oy direction*/
                    GsVolY1[na][nb] += ViscousTerms[0][1];
                    GsVolY2[na][nb] += ViscousTerms[1][1];
                    GsVolY3[na][nb] += ViscousTerms[2][1];
                    GsVolY4[na][nb] += ViscousTerms[3][1];
				}
			}

            for (int order = 1; order <= mathVar::orderElem; order++)
			{
				/*CALCULATE INTEGRALS*/
                VolIntTerm1[order] = process::volumeInte(element, GsVolX1, order, 1);
                VolIntTerm1[order] += process::volumeInte(element, GsVolY1, order, 2);

                VolIntTerm2[order] = process::volumeInte(element, GsVolX2, order, 1);
                VolIntTerm2[order] += process::volumeInte(element, GsVolY2, order, 2);

                VolIntTerm3[order] = process::volumeInte(element, GsVolX3, order, 1);
                VolIntTerm3[order] += process::volumeInte(element, GsVolY3, order, 2);

                VolIntTerm4[order] = process::volumeInte(element, GsVolX4, order, 1);
                VolIntTerm4[order] += process::volumeInte(element, GsVolY4, order, 2);
			}
            VolIntTerm1[0] = 0;
            VolIntTerm2[0] = 0;
            VolIntTerm3[0] = 0;
            VolIntTerm4[0] = 0;
		}

		void calcSurfaceIntegralTerms(int element, std::vector<double> &SurfIntTerm1, std::vector<double> &SurfIntTerm2, std::vector<double> &SurfIntTerm3, std::vector<double> &SurfIntTerm4)
		{
			/*User's guide:
			All input array have form:
			- number of rows: orderElem*/

			int elemType(auxUlti::checkType(element)), edgeName(0);
			int faceBcType(0);

            std::vector<double> Flux1Temp(mathVar::nGauss + 1, 0.0),
                Flux2Temp(mathVar::nGauss + 1, 0.0),
                Flux3Temp(mathVar::nGauss + 1, 0.0),
                Flux4Temp(mathVar::nGauss + 1, 0.0);

			std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));

			std::vector<std::vector<double>>
                Flux1(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
                Flux2(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
                Flux3(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
                Flux4(mathVar::nGauss + 1, std::vector<double>(4, 0.0));

			//std::vector<double> SurInt(4, 0.0);

            for (int nface = 0; nface < elemType; nface++)
			{
                edgeName = meshVar::inedel[element][nface];
				faceBcType = auxUlti::getBCType(edgeName);

				if (faceBcType == 0)  //internal edge
				{
                    for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
					{
						Fluxes = process::NSFEq::getGaussVectorOfConserVarFluxesAtInternal(edgeName, element, nGauss);
                        Flux1[nGauss][nface] = Fluxes[0][0] + Fluxes[0][1];
                        Flux2[nGauss][nface] = Fluxes[1][0] + Fluxes[1][1];
                        Flux3[nGauss][nface] = Fluxes[2][0] + Fluxes[2][1];
                        Flux4[nGauss][nface] = Fluxes[3][0] + Fluxes[3][1];
					}
				}
				else  //boundary edge
				{
                    for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
					{
						Fluxes = NSFEqBCsImplement(element, edgeName, nGauss);
                        Flux1[nGauss][nface] = Fluxes[0][0] + Fluxes[0][1];
                        Flux2[nGauss][nface] = Fluxes[1][0] + Fluxes[1][1];
                        Flux3[nGauss][nface] = Fluxes[2][0] + Fluxes[2][1];
                        Flux4[nGauss][nface] = Fluxes[3][0] + Fluxes[3][1];
					}
				}
			}

            for (int order = 0; order <= mathVar::orderElem; order++)
			{
                for (int nface = 0; nface < elemType; nface++)
				{
                    edgeName = meshVar::inedel[element][nface];  //A BIG BUG!!!!!!
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
					{
                        Flux1Temp[nG] = Flux1[nG][nface];
                        Flux2Temp[nG] = Flux2[nG][nface];
                        Flux3Temp[nG] = Flux3[nG][nface];
                        Flux4Temp[nG] = Flux4[nG][nface];
					}
                    SurfIntTerm1[order] += process::surfaceInte(element, edgeName, Flux1Temp, order);
                    SurfIntTerm2[order] += process::surfaceInte(element, edgeName, Flux2Temp, order);
                    SurfIntTerm3[order] += process::surfaceInte(element, edgeName, Flux3Temp, order);
                    SurfIntTerm4[order] += process::surfaceInte(element, edgeName, Flux4Temp, order);
				}
			}
			//return SurInt;
		}

		std::vector<std::vector<double>> getGaussVectorOfConserVarFluxesAtInternal(int edgeName, int element, int nGauss)
		{
			std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
			double nx(auxUlti::getNormVectorComp(element, edgeName, 1)), ny(auxUlti::getNormVectorComp(element, edgeName, 2));
			double
				termX1P(0.0), termX1M(0.0),  //(rho*u)					or 0
				termX2P(0.0), termX2M(0.0),  //(rho*u^2 + p)			or tauxx
				termX3P(0.0), termX3M(0.0),  //(rho*u*v)				or tauxy
				termX4P(0.0), termX4M(0.0),  //(rho*totalE + p)*u		or tauxx*u + tauxy*v + Qx

				termY1P(0.0), termY1M(0.0),  //(rho*v)					or 0
				termY2P(0.0), termY2M(0.0),  //(rho*u*v)				or tauxy
				termY3P(0.0), termY3M(0.0),  //(rho*v^2 + p)			or tauyy
				termY4P(0.0), termY4M(0.0);  //(rho*totalE + p)*v		or tauxy*u + tauyy*v + Qy
			double rhoPlus(0.0), rhouPlus(0.0), rhovPlus(0.0), rhoEPlus(0.0), rhoMinus(0.0), rhouMinus(0.0), rhovMinus(0.0), rhoEMinus(0.0);

			/*INVISCID TERM*/
			//Get value
            std::tie(termX1P, termX1M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 1, 1, 1);
            std::tie(termX2P, termX2M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 1, 1, 2);
            std::tie(termX3P, termX3M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 1, 1, 3);
            std::tie(termX4P, termX4M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 1, 1, 4);

            std::tie(termY1P, termY1M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 1, 2, 1);
            std::tie(termY2P, termY2M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 1, 2, 2);
            std::tie(termY3P, termY3M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 1, 2, 3);
            std::tie(termY4P, termY4M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 1, 2, 4);

			std::tie(rhoPlus, rhoMinus) = process::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1);
			std::tie(rhouPlus, rhouMinus) = process::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2);
			std::tie(rhovPlus, rhovMinus) = process::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 3);
			std::tie(rhoEPlus, rhoEMinus) = process::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 4);

			/*Calculate fluxes*/
			Fluxes[0][0] = math::numericalFluxes::advectiveFlux(termX1P, termX1M, rhoPlus, rhoMinus, LxFConst[edgeName], nx) + math::numericalFluxes::advectiveFlux(termY1P, termY1M, rhoPlus, rhoMinus, LxFConst[edgeName], ny);
			Fluxes[1][0] = math::numericalFluxes::advectiveFlux(termX2P, termX2M, rhouPlus, rhouMinus, LxFConst[edgeName], nx) + math::numericalFluxes::advectiveFlux(termY2P, termY2M, rhouPlus, rhouMinus, LxFConst[edgeName], ny);
			Fluxes[2][0] = math::numericalFluxes::advectiveFlux(termX3P, termX3M, rhovPlus, rhovMinus, LxFConst[edgeName], nx) + math::numericalFluxes::advectiveFlux(termY3P, termY3M, rhovPlus, rhovMinus, LxFConst[edgeName], ny);
			Fluxes[3][0] = math::numericalFluxes::advectiveFlux(termX4P, termX4M, rhoEPlus, rhoEMinus, LxFConst[edgeName], nx) + math::numericalFluxes::advectiveFlux(termY4P, termY4M, rhoEPlus, rhoEMinus, LxFConst[edgeName], ny);
			
			/*VISCOUS TERM*/
			//Get value
            std::tie(termX1P, termX1M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 2, 1, 1);
            std::tie(termX2P, termX2M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 2, 1, 2);
            std::tie(termX3P, termX3M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 2, 1, 3);
            std::tie(termX4P, termX4M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 2, 1, 4);

            std::tie(termY1P, termY1M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 2, 2, 1);
            std::tie(termY2P, termY2M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 2, 2, 2);
            std::tie(termY3P, termY3M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 2, 2, 3);
            std::tie(termY4P, termY4M) = process::NSFEq::getFinvFvisAtInterfaces(edgeName, element, nGauss, 2, 2, 4);

			/*Calculate fluxes*/
			Fluxes[0][1] = math::numericalFluxes::diffusiveFlux(termX1M, termX1P, rhoPlus, rhoMinus, DiffusiveFluxConst[edgeName], nx) + math::numericalFluxes::diffusiveFlux(termY1M, termY1P, rhoPlus, rhoMinus, DiffusiveFluxConst[edgeName], ny);
			Fluxes[1][1] = math::numericalFluxes::diffusiveFlux(termX2M, termX2P, rhouPlus, rhouMinus, DiffusiveFluxConst[edgeName], nx) + math::numericalFluxes::diffusiveFlux(termY2M, termY2P, rhouPlus, rhouMinus, DiffusiveFluxConst[edgeName], ny);
			Fluxes[2][1] = math::numericalFluxes::diffusiveFlux(termX3M, termX3P, rhovPlus, rhovMinus, DiffusiveFluxConst[edgeName], nx) + math::numericalFluxes::diffusiveFlux(termY3M, termY3P, rhovPlus, rhovMinus, DiffusiveFluxConst[edgeName], ny);
			Fluxes[3][1] = math::numericalFluxes::diffusiveFlux(termX4M, termX4P, rhoEPlus, rhoEMinus, DiffusiveFluxConst[edgeName], nx) + math::numericalFluxes::diffusiveFlux(termY4M, termY4P, rhoEPlus, rhoEMinus, DiffusiveFluxConst[edgeName], ny);

			return Fluxes;
		}

		std::vector<double> solveTimeMarching(int element, std::vector<double> &ddtArr, std::vector<double> &UnArr, int RKOrder, int varType)
		{
			std::vector<double> OutArr(mathVar::orderElem + 1, 0.0);
			
			switch (RKOrder)
			{
			case 1:
			{
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					OutArr[order] = dt * ddtArr[order] + UnArr[order];
				}
			}
			break;
			case 2:
			{
				switch (varType)
				{
				case 1://rho
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = 0.25*(dt * ddtArr[order] + UnArr[order]) + (3.0 / 4.0)*rho0[element][order];
					}
				}
				break;
				case 2://rhou
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = 0.25*(dt * ddtArr[order] + UnArr[order]) + (3.0 / 4.0)*rhou0[element][order];
					}
				}
				break;
				case 3://rhov
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = 0.25*(dt * ddtArr[order] + UnArr[order]) + (3.0 / 4.0)*rhov0[element][order];
					}
				}
				break;
				case 4://rhoE
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = 0.25*(dt * ddtArr[order] + UnArr[order]) + (3.0 / 4.0)*rhoE0[element][order];
					}
				}
				break;
				default:
					break;
				}
			}
			break;
			case 3:
			{
				switch (varType)
				{
				case 1://rho
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = (2.0 / 3.0)*(dt * ddtArr[order] + UnArr[order]) + (1.0 / 3.0)*rho0[element][order];
					}
				}
				break;
				case 2://rhou
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = (2.0 / 3.0)*(dt * ddtArr[order] + UnArr[order]) + (1.0 / 3.0)*rhou0[element][order];
					}
				}
				break;
				case 3://rhov
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = (2.0 / 3.0)*(dt * ddtArr[order] + UnArr[order]) + (1.0 / 3.0)*rhov0[element][order];
					}
				}
				break;
				case 4://rhoE
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						OutArr[order] = (2.0 / 3.0)*(dt * ddtArr[order] + UnArr[order]) + (1.0 / 3.0)*rhoE0[element][order];
					}
				}
				break;
				default:
					break;
				}
			}
			break;
			default:
				break;
			}
			
			return OutArr;
		}
	}

	namespace timeDiscretization
	{
		void calcGlobalTimeStep()
		{
			std::vector<double>timeStepArr(meshVar::nelem2D, 1.0);
            for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
			{
				timeStepArr[nelement] = process::timeDiscretization::localTimeStep(nelement);
				
			}
			runTime += dt;
			dt = *std::min_element(timeStepArr.begin(), timeStepArr.end());  //find min value of vector
		}

		double localTimeStep(int element)
		{
			double deltaT(0.0), uVal(0.0), vVal(0.0), velocity(0.0), TVal(0.0), aSound(0.0), muVal(0.0), size(meshVar::cellSize[element]);
			uVal = rhou[element][0] / rho[element][0];
			vVal = rhov[element][0] / rho[element][0];
			velocity = sqrt(pow(uVal, 2) + pow(vVal, 2));
			TVal = math::CalcTFromConsvVar(rho[element][0], rhou[element][0], rhov[element][0], rhoE[element][0]);
			aSound = math::CalcSpeedOfSound(TVal);
			muVal = math::CalcVisCoef(TVal);
			
			deltaT = fabs((1.0 / pow(mathVar::orderElem + 1, 2))*(size*systemVar::CFL) / (aSound + velocity + muVal/size));
			return deltaT;
		}

		double localErrorEstimate(int element, std::vector<double> &ddtArr)
		{
			double xC(-1.0 / 3.0), yC(1.0 / 3.0), errorVal(0.0);

			if (auxUlti::checkType(element) == 4)
			{
				xC = 0.0;
				yC = 0.0;
			}

			//Compute basis function
			math::basisFc(xC, yC, auxUlti::checkType(element));
            for (int order = 0; order <= mathVar::orderElem; order++)
			{
				errorVal += ddtArr[order] * mathVar::B[order];
			}

			return abs(errorVal);
		}

		void globalErrorEstimate()
		{
			double xC(0.0), yC(0.0), aC(0.0), bC(0.0);
			std::vector<double> rhoRes(meshVar::nelem2D, 0.0),
				rhouRes(meshVar::nelem2D, 0.0),
				rhovRes(meshVar::nelem2D, 0.0),
				rhoERes(meshVar::nelem2D, 0.0);
			double rhoResGlobal(1.0), rhouResGlobal(1.0), rhovResGlobal(1.0), rhoEResGlobal(1.0);

            for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
			{
				std::tie(xC, yC) = auxUlti::getCellCentroid(nelement);
				std::tie(aC, bC) = math::inverseMapping(nelement, xC, yC);
				rhoRes[nelement] = math::calcResidualFromResidualOfOrder(nelement, aC, bC, 1);
				rhouRes[nelement] = math::calcResidualFromResidualOfOrder(nelement, aC, bC, 2);
				rhovRes[nelement] = math::calcResidualFromResidualOfOrder(nelement, aC, bC, 3);
				rhoERes[nelement] = math::calcResidualFromResidualOfOrder(nelement, aC, bC, 4);
			}
			rhoResGlobal = *std::max_element(rhoRes.begin(), rhoRes.end());
			rhouResGlobal = *std::max_element(rhouRes.begin(), rhouRes.end());
			rhovResGlobal = *std::max_element(rhovRes.begin(), rhovRes.end());
			rhoEResGlobal = *std::max_element(rhoERes.begin(), rhoERes.end());
			
			IO::residualOutput(rhoResGlobal, rhouResGlobal, rhovResGlobal, rhoEResGlobal);
		}

		void TVDRK_1step(int RKOrder)
		{
            //COMPUTE GAUSS VALUES AND INTERFACES VALUES
			process::calcVolumeGaussValues();
            process::calcValuesAtInterface();

            if (flowProperties::massDiffusion)
            {
                if (systemVar::auxVariables==1)
                {
                    process::auxEq::massDiffusion::BR1::solveDivRho();
                    auxUlti::functionsOfParallelComputing::sendReceivedRho();
                }
                else if (systemVar::auxVariables==2)
                {
                    process::auxEq::massDiffusion::BR2::solveDivRho();
                    //Chua sua cho BR2
                }
            }

            //CALCULATE T
            process::calcTGauss();

            //SOLVE AUXILARY EQUATION
			process::auxEq::solveAuxEquation();
            auxUlti::functionsOfParallelComputing::sendReceivedU();

			//SOLVE NSF EQUATION
            process::NSFEq::calcFinvFvisAtInterface();
			process::NSFEq::solveNSFEquation(RKOrder);
		}

		void TVDRK3()
		{
            rho0 = rho;
            rhou0 = rhou;
            rhov0 = rhov;
            rhoE0 = rhoE;

            for (int iRKOrder = 1; iRKOrder <= 3; iRKOrder++)
			{
				process::timeDiscretization::TVDRK_1step(iRKOrder);
                rho = rhoN;
                rhou = rhouN;
                rhov = rhovN;
                rhoE = rhoEN;
				limiter::limiter();
                auxUlti::functionsOfParallelComputing::sendReceiveU();
			}
		}

        namespace parallel {
        void minTimeStep()
        {
            MPI_Barrier(MPI_COMM_WORLD);
            if (systemVar::currentProc==0)
            {
                std::vector<double>vectorTimeStep(systemVar::totalProc,1.0);
                vectorTimeStep[0]=dt;
                for (int inode=1;inode<systemVar::totalProc;inode++)
                {
                    MPI_Recv(&vectorTimeStep[inode], 1, MPI_DOUBLE, inode, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                dt=*std::min_element(vectorTimeStep.begin(), vectorTimeStep.end());
            }
            else
            {
                MPI_Send(&dt, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }

            //Sau khi tim duoc dt min, send gia tri dt lai ve cac process
            if (systemVar::currentProc==0)
            {
                for (int inode=1;inode<systemVar::totalProc;inode++)
                {
                    MPI_Send(&dt, 1, MPI_DOUBLE, inode, 0, MPI_COMM_WORLD);
                }
            }
            else
            {
                MPI_Recv(&dt, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        }
	}

	std::tuple<double, double> getInternalValuesFromCalculatedArrays(int edge, int element, int nG, int valType)
	{
		//valType = 1(rho*mu), 2(rhou*mu), 3(rhov*mu), 4(rhoE*mu)
		bool isMaster(auxUlti::checkMaster(element, edge));
		int locationPlus(-1), locationMinus(-1);
		double valPlus(0.0), valMinus(0.0);
		if (isMaster)
		{
			locationPlus = nG;
			locationMinus = nG + mathVar::nGauss + 1;
		}
		else
		{
			locationPlus = nG + mathVar::nGauss + 1;
			locationMinus = nG;
		}

		switch (valType)
		{
		case 1:
		{
            valPlus = surfaceFields::rho[edge][locationPlus];
            valMinus = surfaceFields::rho[edge][locationMinus];
		}
		break;
		case 2:
		{
            valPlus = surfaceFields::rhou[edge][locationPlus];
            valMinus = surfaceFields::rhou[edge][locationMinus];
		}
		break;
		case 3:
		{
            valPlus = surfaceFields::rhov[edge][locationPlus];
            valMinus = surfaceFields::rhov[edge][locationMinus];
		}
		break;
		case 4:
		{
            valPlus = surfaceFields::rhoE[edge][locationPlus];
            valMinus = surfaceFields::rhoE[edge][locationMinus];
		}
		break;
		default:
			break;
		}
		return std::make_tuple(valPlus, valMinus);
	}

	double volumeInte(int elem, std::vector< std::vector<double> > &Ui, int order, int direction)
	{
		double Int(0.0);
		std::vector<std::vector<double>> A(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
		double dBi(0.0);

        for (int na = 0; na <= mathVar::nGauss; na++)
		{
            for (int nb = 0; nb <= mathVar::nGauss; nb++)
			{
				dBi = math::Calc_dBxdBy(elem, order, na, nb, direction);
                A[na][nb] = dBi * Ui[na][nb];
			}
		}
		Int = math::volumeInte(A, elem);

		return Int;
	}

	double surfaceInte(int elem, int edge, std::vector<double> &FluxVector, int order)
	{
		double Bi(0.0), a(0.0), b(0.0), inte(0.0);
		std::vector<double> F(mathVar::nGauss + 1, 0.0);

		//inte = 0.0;  //value of integral is reset for each order
        for (int nG = 0; nG <= mathVar::nGauss; nG++)
		{
			std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, elem, nG);
			math::basisFc(a, b, auxUlti::checkType(elem));  //mathVar::B is changed after this command line is excuted
			Bi = mathVar::B[order];
			F[nG] = Bi * FluxVector[nG];
		}
        inte = math::surfaceInte(F, edge);

		return inte;
	}

	std::vector<std::vector<double>> calculateStiffMatrix(int element)
	{
		std::vector<std::vector<double>> StiffMatrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));

        for (int order1 = 0; order1 <= mathVar::orderElem; order1++)
		{
            for (int order2 = 0; order2 <= mathVar::orderElem; order2++)
			{
				StiffMatrix[order1][order2] = process::calculateStiffMatrixElement(element, order1, order2);
			}
		}
		return StiffMatrix;
	}

	double calculateStiffMatrixElement(int element, int order1, int order2)
	{
		double B1(0.0), B2(0.0), Inte(0.0);
		std::vector<std::vector<double>> FMatrix(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
		int elemType(auxUlti::checkType(element));

		switch (elemType)
		{
		case 3:
		{
            for (int na = 0; na <= mathVar::nGauss; na++)
			{
                for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					B1 = mathVar::BPts_Tri[order1][na][nb];
					B2 = mathVar::BPts_Tri[order2][na][nb];
					FMatrix[na][nb] = B1 * B2;
				}
			}
		}
		break;
		case 4:
		{
            for (int na = 0; na <= mathVar::nGauss; na++)
			{
                for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					B1 = mathVar::BPts_Quad[order1][na][nb];
					B2 = mathVar::BPts_Quad[order2][na][nb];
					FMatrix[na][nb] = B1 * B2;
				}
			}
		}
		break;
		default:
			break;
		}
		Inte = math::volumeInte(FMatrix, element);
		return Inte;
	}

	bool checkRunningCond()
	{
		bool run(true);
		if (runTime >= systemVar::Ttime)
		{
			run = false;  //stop running if runtime bigger than Ttime (total time)
		}
		return run;
	}
}
