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
            pt1=meshVar::inpoed[0][iedge];
            pt2=meshVar::inpoed[1][iedge];
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
        for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
		{
            for (int na = 0; na <= mathVar::nGauss; na++)
			{
                for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
					rhoVolGauss[nelem][na][nb] = math::pointValue(nelem, a, b, 1, 2);
					rhouVolGauss[nelem][na][nb] = math::pointValue(nelem, a, b, 2, 2);
					rhovVolGauss[nelem][na][nb] = math::pointValue(nelem, a, b, 3, 2);
					rhoEVolGauss[nelem][na][nb] = math::pointValue(nelem, a, b, 4, 2);
				}
			}
		}
	}

	namespace auxEq
	{
		void calcValuesAtInterface()
		{
			//mu included
			int masterCell(-1), slaveCell(-1), bcGrp(0);
			double muMaster(0.0), muSlave(0.0), tempMaster(0.0), tempSlave(0.0);
            for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
			{
				bcGrp = auxUlti::getBCType(iedge);
				if (bcGrp == 0)
				{
					std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
                    for (int nG = 0; nG <= mathVar::nGauss; nG++)
					{
						std::tie(muMaster, muSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 7, 1);

						//rho*mu
						std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 1, 2);
						aux_interface_rho[iedge][nG] = tempMaster * muMaster;
						aux_interface_rho[iedge][nG + mathVar::nGauss + 1] = tempSlave * muSlave;
						interface_rho[iedge][nG] = tempMaster;
						interface_rho[iedge][nG + mathVar::nGauss + 1] = tempSlave;

						//rhou*mu
						std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 2, 2);
						aux_interface_rhou[iedge][nG] = tempMaster * muMaster;
						aux_interface_rhou[iedge][nG + mathVar::nGauss + 1] = tempSlave * muSlave;
						interface_rhou[iedge][nG] = tempMaster;
						interface_rhou[iedge][nG + mathVar::nGauss + 1] = tempSlave;

						//rhov*mu
						std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 3, 2);
						aux_interface_rhov[iedge][nG] = tempMaster * muMaster;
						aux_interface_rhov[iedge][nG + mathVar::nGauss + 1] = tempSlave * muSlave;
						interface_rhov[iedge][nG] = tempMaster;
						interface_rhov[iedge][nG + mathVar::nGauss + 1] = tempSlave;

						//rhoE*mu
						std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 4, 2);
						aux_interface_rhoE[iedge][nG] = tempMaster * muMaster;
						aux_interface_rhoE[iedge][nG + mathVar::nGauss + 1] = tempSlave * muSlave;
						interface_rhoE[iedge][nG] = tempMaster;
						interface_rhoE[iedge][nG + mathVar::nGauss + 1] = tempSlave;
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
				valPlus = aux_interface_rho[edge][locationPlus];
				valMinus = aux_interface_rho[edge][locationMinus];
			}
			break;
			case 2:
			{
				valPlus = aux_interface_rhou[edge][locationPlus];
				valMinus = aux_interface_rhou[edge][locationMinus];
			}
			break;
			case 3:
			{
				valPlus = aux_interface_rhov[edge][locationPlus];
				valMinus = aux_interface_rhov[edge][locationMinus];
			}
			break;
			case 4:
			{
				valPlus = aux_interface_rhoE[edge][locationPlus];
				valMinus = aux_interface_rhoE[edge][locationMinus];
			}
			break;
			default:
				break;
			}
			return std::make_tuple(valPlus, valMinus);
		}

		void solveAuxEquation()
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
				//1) Calculate Stiff matrix
				//StiffMatrix = process::calculateStiffMatrix(nelement);
				
				//2) Calculate Right hand side terms
				process::auxEq::CalcRHSTerm(nelement, rhoRHSTermOxDir, rhoRHSTermOyDir, rhouRHSTermOxDir, rhouRHSTermOyDir, rhovRHSTermOxDir, rhovRHSTermOyDir, rhoERHSTermOxDir, rhoERHSTermOyDir);

				//3) Solve for auxilary variables
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
				{
					//Ox direction
					rhoX[nelement][iorder] = rhoRHSTermOxDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
					rhouX[nelement][iorder] = rhouRHSTermOxDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
					rhovX[nelement][iorder] = rhovRHSTermOxDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
					rhoEX[nelement][iorder] = rhoERHSTermOxDir[iorder] / stiffMatrixCoeffs[nelement][iorder];

					//Oy direction
					rhoY[nelement][iorder] = rhoRHSTermOyDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
					rhouY[nelement][iorder] = rhouRHSTermOyDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
					rhovY[nelement][iorder] = rhovRHSTermOyDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
					rhoEY[nelement][iorder] = rhoERHSTermOyDir[iorder] / stiffMatrixCoeffs[nelement][iorder];
				}
			}
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
			process::auxEq::calcVolumeIntegralTerms(element, rhoVolIntOx, rhouVolIntOx, rhovVolIntOx, rhoEVolIntOx, rhoVolIntOy, rhouVolIntOy, rhovVolIntOy, rhoEVolIntOy);

			/*2. Calculate surface integral term*/
			process::auxEq::calcSurfaceIntegralTerms(element, rhoSurfIntOx, rhouSurfIntOx, rhovSurfIntOx, rhoESurfIntOx, rhoSurfIntOy, rhouSurfIntOy, rhovSurfIntOy, rhoESurfIntOy);

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
			//rho -------------------------------------------------------------------------------------------
			rhoGsVol = process::auxEq::getGaussMatrixOfConserVar(element, 1);
			//rhou ------------------------------------------------------------------------------------------
			rhouGsVol = process::auxEq::getGaussMatrixOfConserVar(element, 2);
			//rhov ------------------------------------------------------------------------------------------
			rhovGsVol = process::auxEq::getGaussMatrixOfConserVar(element, 3);
			//rhou ------------------------------------------------------------------------------------------
			rhoEGsVol = process::auxEq::getGaussMatrixOfConserVar(element, 4);

            for (int order = 1; order <= mathVar::orderElem; order++)
			{
				rhoVolIntX[order] = process::volumeInte(element, rhoGsVol, order, 1);
				rhouVolIntX[order] = process::volumeInte(element, rhouGsVol, order, 1);
				rhovVolIntX[order] = process::volumeInte(element, rhovGsVol, order, 1);
				rhoEVolIntX[order] = process::volumeInte(element, rhoEGsVol, order, 1);

				rhoVolIntY[order] = process::volumeInte(element, rhoGsVol, order, 2);
				rhouVolIntY[order] = process::volumeInte(element, rhouGsVol, order, 2);
				rhovVolIntY[order] = process::volumeInte(element, rhovGsVol, order, 2);
				rhoEVolIntY[order] = process::volumeInte(element, rhoEGsVol, order, 2);
			}
            rhoVolIntX[0] = 0;
            rhouVolIntX[0] = 0;
            rhovVolIntX[0] = 0;
            rhoEVolIntX[0] = 0;

            rhoVolIntY[0] = 0;
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
			process::auxEq::getGaussVectorOfConserVar(element, rhoFluxX, rhouFluxX, rhovFluxX, rhoEFluxX, rhoFluxY, rhouFluxY, rhovFluxY, rhoEFluxY);

			/*2. Calculates surface integrals of all conservative variables at all order*/
			
            for (int nface = 0; nface < elemType; nface++)
			{
				edgeName = meshVar::inedel[nface][element];
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
				{
					rhoFluxXTemp[nG] = rhoFluxX[nG][nface];
					rhouFluxXTemp[nG] = rhouFluxX[nG][nface];
					rhovFluxXTemp[nG] = rhovFluxX[nG][nface];
					rhoEFluxXTemp[nG] = rhoEFluxX[nG][nface];

					rhoFluxYTemp[nG] = rhoFluxY[nG][nface];
					rhouFluxYTemp[nG] = rhouFluxY[nG][nface];
					rhovFluxYTemp[nG] = rhovFluxY[nG][nface];
					rhoEFluxYTemp[nG] = rhoEFluxY[nG][nface];
				}

                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					rhoSurfIntX[order] += process::surfaceInte(element, edgeName, rhoFluxXTemp, order);
					rhouSurfIntX[order] += process::surfaceInte(element, edgeName, rhouFluxXTemp, order);
					rhovSurfIntX[order] += process::surfaceInte(element, edgeName, rhovFluxXTemp, order);
					rhoESurfIntX[order] += process::surfaceInte(element, edgeName, rhoEFluxXTemp, order);

					rhoSurfIntY[order] += process::surfaceInte(element, edgeName, rhoFluxYTemp, order);
					rhouSurfIntY[order] += process::surfaceInte(element, edgeName, rhouFluxYTemp, order);
					rhovSurfIntY[order] += process::surfaceInte(element, edgeName, rhovFluxYTemp, order);
					rhoESurfIntY[order] += process::surfaceInte(element, edgeName, rhoEFluxYTemp, order);
				}
			}
		}

		std::vector<std::vector<double>> getGaussMatrixOfConserVar(int element, int valType)
		{
			std::vector<std::vector<double>> GaussMatrix(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
			double muGs(0.0), a(0.0), b(0.0), val(0.0);
            for (int na = 0; na <= mathVar::nGauss; na++)
			{
                for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
					muGs = math::pointValue(element, a, b, 7, 1);
					switch (valType)
					{
					case 1:
					{
						val = rhoVolGauss[element][na][nb];
					}
					break;
					case 2:
					{
						val = rhouVolGauss[element][na][nb];
					}
					break;
					case 3:
					{
						val = rhovVolGauss[element][na][nb];
					}
					break;
					case 4:
					{
						val = rhoEVolGauss[element][na][nb];
					}
					break;
					default:
						break;
					}
					GaussMatrix[na][nb] = val*muGs;
				}
			}
			return GaussMatrix;
		}

		std::vector<std::vector<double>> getVectorOfConserVarFluxesAtInternal(int edge, int element, int nG, double nx, double ny)
		{
			std::vector<std::vector<double>> gaussVector(4, std::vector<double>(2, 0.0)); //columns 0, 1 are Ox, Oy values
			double tempP(0.0), tempM(0.0);
            for (int i = 0; i < 4; i++)
			{
				std::tie(tempP, tempM) = process::auxEq::getInternalValuesFromCalculatedArrays(edge, element, nG, i+1);
				gaussVector[i][0] = math::numericalFluxes::auxFlux(tempM, tempP, nx);
				gaussVector[i][1] = math::numericalFluxes::auxFlux(tempM, tempP, ny);
			}

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
				edgeName = meshVar::inedel[nface][element];
				faceBcType = auxUlti::getBCType(edgeName);
				double nx(auxUlti::getNormVectorComp(element, edgeName, 1)), ny(auxUlti::getNormVectorComp(element, edgeName, 2));
				if (faceBcType == 0)  //internal edge
				{
                    for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
					{
						gaussVector = process::auxEq::getVectorOfConserVarFluxesAtInternal(edgeName, element, nGauss, nx, ny);
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
			}
		}
	}//end namespace auxEq

	namespace NSFEq
	{
		void calcValuesAtInterface()
		{
			//mu included
			int masterCell(-1), slaveCell(-1), bcGrp(0);
			double uMaster(0.0), vMaster(0.0), totalEMaster(0.0), TMaster(0.0), pMaster(0.0),
				uSlave(0.0), vSlave(0.0), totalESlave(0.0), TSlave(0.0), pSlave(0.0), eMaster(0.0), eSlave(0.0),
                uMagM(0.0), uMagP(0.0), aM(0.0), aP(0.0);
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
                        for (int i = 0; i < 4; i++)
						{
							/*INVISCID TERMS*/
							std::tie(UMaster[i], USlave[i]) = math::internalSurfaceValue(iedge, masterCell, nG, i + 1, 2);
							/*VISCOUS TERMS*/
							std::tie(dUXMaster[i], dUXSlave[i]) = math::internalSurfaceDerivativeValue(iedge, masterCell, nG, i + 1, 1);  //dUx
							std::tie(dUYMaster[i], dUYSlave[i]) = math::internalSurfaceDerivativeValue(iedge, masterCell, nG, i + 1, 2);  //dUy
						}

						uMaster = (UMaster[1] / UMaster[0]);
						uSlave = (USlave[1] / USlave[0]);

						vMaster = (UMaster[2] / UMaster[0]);
						vSlave = (USlave[2] / USlave[0]);

						totalEMaster = (UMaster[3] / UMaster[0]);
						totalESlave = (USlave[3] / USlave[0]);

						/*calculate T and P*/
						TMaster = math::CalcTFromConsvVar(UMaster[0], UMaster[1], UMaster[2], UMaster[3]);
						TSlave = math::CalcTFromConsvVar(USlave[0], USlave[1], USlave[2], USlave[3]);
						pMaster = math::CalcP(TMaster, UMaster[0]);
						pSlave = math::CalcP(TSlave, USlave[0]);
						eMaster = material::Cv*TMaster;
						eSlave = material::Cv*TSlave;

						/*INVISCID TERMS*/
						/*calculate velocity magnitude*/
						uMagP = sqrt(pow(uMaster, 2) + pow(vMaster, 2));
						uMagM = sqrt(pow(uSlave, 2) + pow(vSlave, 2));

						/*calculate speed of sound*/
						aP = math::CalcSpeedOfSound(TMaster);
						aM = math::CalcSpeedOfSound(TSlave);

						/*calculate constant for Lax-Friederich flux*/
						CArray[nG] = math::numericalFluxes::constantC(uMagP, uMagM, aP, aM); 

						/*calculate inviscid terms*/
						
						std::tie(invis_interface_rhoX[iedge][nG], invis_interface_rhouX[iedge][nG], invis_interface_rhovX[iedge][nG], invis_interface_rhoEX[iedge][nG]) = math::inviscidTerms::calcInvisTermsFromPriVars(UMaster[0], uMaster, vMaster, totalEMaster, pMaster, 1);
						std::tie(invis_interface_rhoY[iedge][nG], invis_interface_rhouY[iedge][nG], invis_interface_rhovY[iedge][nG], invis_interface_rhoEY[iedge][nG]) = math::inviscidTerms::calcInvisTermsFromPriVars(UMaster[0], uMaster, vMaster, totalEMaster, pMaster, 2);

						std::tie(invis_interface_rhoX[iedge][nG + mathVar::nGauss + 1], invis_interface_rhouX[iedge][nG + mathVar::nGauss + 1], invis_interface_rhovX[iedge][nG + mathVar::nGauss + 1], invis_interface_rhoEX[iedge][nG + mathVar::nGauss + 1]) = math::inviscidTerms::calcInvisTermsFromPriVars(USlave[0], uSlave, vSlave, totalESlave, pSlave, 1);
						std::tie(invis_interface_rhoY[iedge][nG + mathVar::nGauss + 1], invis_interface_rhouY[iedge][nG + mathVar::nGauss + 1], invis_interface_rhovY[iedge][nG + mathVar::nGauss + 1], invis_interface_rhoEY[iedge][nG + mathVar::nGauss + 1]) = math::inviscidTerms::calcInvisTermsFromPriVars(USlave[0], uSlave, vSlave, totalESlave, pSlave, 2);

						/*calculate viscous terms*/
						StressHeatP = math::viscousTerms::calcStressTensorAndHeatFlux(UMaster, dUXMaster, dUYMaster);
						std::tie(Vis_interface_rhoX[iedge][nG], Vis_interface_rhouX[iedge][nG], Vis_interface_rhovX[iedge][nG], Vis_interface_rhoEX[iedge][nG]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uMaster, vMaster, 1);
						std::tie(Vis_interface_rhoY[iedge][nG], Vis_interface_rhouY[iedge][nG], Vis_interface_rhovY[iedge][nG], Vis_interface_rhoEY[iedge][nG]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uMaster, vMaster, 2);

						StressHeatM = math::viscousTerms::calcStressTensorAndHeatFlux(USlave, dUXSlave, dUYSlave);
						std::tie(Vis_interface_rhoX[iedge][nG + mathVar::nGauss + 1], Vis_interface_rhouX[iedge][nG + mathVar::nGauss + 1], Vis_interface_rhovX[iedge][nG + mathVar::nGauss + 1], Vis_interface_rhoEX[iedge][nG + mathVar::nGauss + 1]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uSlave, vSlave, 1);
						std::tie(Vis_interface_rhoY[iedge][nG + mathVar::nGauss + 1], Vis_interface_rhouY[iedge][nG + mathVar::nGauss + 1], Vis_interface_rhovY[iedge][nG + mathVar::nGauss + 1], Vis_interface_rhoEY[iedge][nG + mathVar::nGauss + 1]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uSlave, vSlave, 2);
					
						/*calculate constant for diffusive flux*/
						//BetaArray[nG] = math::numericalFluxes::constantBeta(uMagP, uMagM, UMaster[0], USlave[0], eMaster, eSlave, pMaster, pSlave, StressHeatP, StressHeatM, vectorn);
					}
					LxFConst[iedge] = *std::max_element(CArray.begin(), CArray.end());
					//DiffusiveFluxConst[iedge] = *std::max_element(BetaArray.begin(), BetaArray.end());
				}
			}
		}

		std::tuple<double, double> getInternalValuesFromCalculatedArrays(int edge, int element, int nG, int mod, int direction, int valType)
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
						valPlus = invis_interface_rhoX[edge][locationPlus];
						valMinus = invis_interface_rhoX[edge][locationMinus];
					}
					break;
					case 2:
					{
						valPlus = invis_interface_rhouX[edge][locationPlus];
						valMinus = invis_interface_rhouX[edge][locationMinus];
					}
					break;
					case 3:
					{
						valPlus = invis_interface_rhovX[edge][locationPlus];
						valMinus = invis_interface_rhovX[edge][locationMinus];
					}
					break;
					case 4:
					{
						valPlus = invis_interface_rhoEX[edge][locationPlus];
						valMinus = invis_interface_rhoEX[edge][locationMinus];
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
						valPlus = invis_interface_rhoY[edge][locationPlus];
						valMinus = invis_interface_rhoY[edge][locationMinus];
					}
					break;
					case 2:
					{
						valPlus = invis_interface_rhouY[edge][locationPlus];
						valMinus = invis_interface_rhouY[edge][locationMinus];
					}
					break;
					case 3:
					{
						valPlus = invis_interface_rhovY[edge][locationPlus];
						valMinus = invis_interface_rhovY[edge][locationMinus];
					}
					break;
					case 4:
					{
						valPlus = invis_interface_rhoEY[edge][locationPlus];
						valMinus = invis_interface_rhoEY[edge][locationMinus];
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
						valPlus = Vis_interface_rhoX[edge][locationPlus];
						valMinus = Vis_interface_rhoX[edge][locationMinus];
					}
					break;
					case 2:
					{
						valPlus = Vis_interface_rhouX[edge][locationPlus];
						valMinus = Vis_interface_rhouX[edge][locationMinus];
					}
					break;
					case 3:
					{
						valPlus = Vis_interface_rhovX[edge][locationPlus];
						valMinus = Vis_interface_rhovX[edge][locationMinus];
					}
					break;
					case 4:
					{
						valPlus = Vis_interface_rhoEX[edge][locationPlus];
						valMinus = Vis_interface_rhoEX[edge][locationMinus];
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
						valPlus = Vis_interface_rhoY[edge][locationPlus];
						valMinus = Vis_interface_rhoY[edge][locationMinus];
					}
					break;
					case 2:
					{
						valPlus = Vis_interface_rhouY[edge][locationPlus];
						valMinus = Vis_interface_rhouY[edge][locationMinus];
					}
					break;
					case 3:
					{
						valPlus = Vis_interface_rhovY[edge][locationPlus];
						valMinus = Vis_interface_rhovY[edge][locationMinus];
					}
					break;
					case 4:
					{
						valPlus = Vis_interface_rhoEY[edge][locationPlus];
						valMinus = Vis_interface_rhoEY[edge][locationMinus];
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

				//5) Save results to conservative variables array
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
            for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
			{
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					rho[nelement][order] = rhoN[nelement][order];
					rhou[nelement][order] = rhouN[nelement][order];
					rhov[nelement][order] = rhovN[nelement][order];
					rhoE[nelement][order] = rhoEN[nelement][order];
				}
			}
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
				rhoVal(rhoVolGauss[element][na][nb]),
				rhouVal(rhouVolGauss[element][na][nb]),
				rhovVal(rhovVolGauss[element][na][nb]),
				rhoEVal(rhoEVolGauss[element][na][nb]),
				uVal(0.0),
				vVal(0.0),
				pVal(0.0),
				totalE(0.0);

			uVal = rhouVal / rhoVal;
			vVal = rhovVal / rhoVal;
			totalE = rhoEVal / rhoVal;
			pVal = math::CalcP(math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal), rhoVal);

			/*1. Ox direction*/
			std::tie(InviscidTerm[0][0], InviscidTerm[1][0], InviscidTerm[2][0], InviscidTerm[3][0]) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoVal, uVal, vVal, totalE, pVal, 1);

			/*2. Oy direction*/
			std::tie(InviscidTerm[0][1], InviscidTerm[1][1], InviscidTerm[2][1], InviscidTerm[3][1]) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoVal, uVal, vVal, totalE, pVal, 2);

			return InviscidTerm;
		}

		/*Function calculates Viscous terms at Gauss point (a, b)*/
		std::vector<std::vector<double>> calcGaussViscousTerm(int element, int na, int nb)
		{
			double uVal(rhouVolGauss[element][na][nb] / rhoVolGauss[element][na][nb]), vVal(rhovVolGauss[element][na][nb] / rhoVolGauss[element][na][nb]), a(0.0), b(0.0);
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
			vectorU[0] = rhoVolGauss[element][na][nb];
			vectorU[1] = rhouVolGauss[element][na][nb];
			vectorU[2] = rhovVolGauss[element][na][nb];
			vectorU[3] = rhoEVolGauss[element][na][nb];
            for (int i = 0; i < 4; i++)
			{
				vectordUx[i] = math::pointAuxValue(element, a, b, i + 1, 1);
				vectordUy[i] = math::pointAuxValue(element, a, b, i + 1, 2);
			}

			/*calculate stresses and heat fluxes*/
			StressHeatFlux = math::viscousTerms::calcStressTensorAndHeatFlux(vectorU, vectordUx, vectordUy);
			std::tie(ViscousTerm[0][0], ViscousTerm[1][0], ViscousTerm[2][0], ViscousTerm[3][0]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatFlux, uVal, vVal, 1);
			std::tie(ViscousTerm[0][1], ViscousTerm[1][1], ViscousTerm[2][1], ViscousTerm[3][1]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatFlux, uVal, vVal, 2);
			return ViscousTerm;
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
					//ViscGsVolX1[na][nb] = ViscousTerms[0][0];
                    GsVolX2[na][nb] += ViscousTerms[1][0];
                    GsVolX3[na][nb] += ViscousTerms[2][0];
                    GsVolX4[na][nb] += ViscousTerms[3][0];
					/*B2. Viscous term on Oy direction*/
					//ViscGsVolY1[na][nb] = ViscousTerms[0][1];
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
				edgeName = meshVar::inedel[nface][element];
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
                    edgeName = meshVar::inedel[nface][element];  //A BIG BUG!!!!!!
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
			std::tie(termX1P, termX1M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 1, 1);
			std::tie(termX2P, termX2M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 1, 2);
			std::tie(termX3P, termX3M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 1, 3);
			std::tie(termX4P, termX4M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 1, 4);

			std::tie(termY1P, termY1M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 2, 1);
			std::tie(termY2P, termY2M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 2, 2);
			std::tie(termY3P, termY3M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 2, 3);
			std::tie(termY4P, termY4M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 2, 4);

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
			std::tie(termX1P, termX1M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 1, 1);
			std::tie(termX2P, termX2M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 1, 2);
			std::tie(termX3P, termX3M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 1, 3);
			std::tie(termX4P, termX4M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 1, 4);

			std::tie(termY1P, termY1M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 2, 1);
			std::tie(termY2P, termY2M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 2, 2);
			std::tie(termY3P, termY3M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 2, 3);
			std::tie(termY4P, termY4M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 2, 4);

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
			//COMPUTE GAUSS VALUES
			process::calcVolumeGaussValues();

			//SOLVE AUXILARY EQUATION
			process::auxEq::calcValuesAtInterface();
			process::auxEq::solveAuxEquation();

			//SOLVE NSF EQUATION
			process::NSFEq::calcValuesAtInterface();
			process::NSFEq::solveNSFEquation(RKOrder);
		}

		void TVDRK3()
		{
            for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
			{
                for (int order = 0; order <= mathVar::orderElem; order++)
				{
					rho0[nelement][order] = rho[nelement][order];
					rhou0[nelement][order] = rhou[nelement][order];
					rhov0[nelement][order] = rhov[nelement][order];
					rhoE0[nelement][order] = rhoE[nelement][order];
				}
			}

            for (int iRKOrder = 1; iRKOrder <= 3; iRKOrder++)
			{
				process::timeDiscretization::TVDRK_1step(iRKOrder);
                for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
				{
                    for (int order = 0; order <= mathVar::orderElem; order++)
					{
						rho[nelement][order] = rhoN[nelement][order];
						rhou[nelement][order] = rhouN[nelement][order];
						rhov[nelement][order] = rhovN[nelement][order];
						rhoE[nelement][order] = rhoEN[nelement][order];
					}
				}
				limiter::limiter();
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
			valPlus = interface_rho[edge][locationPlus];
			valMinus = interface_rho[edge][locationMinus];
		}
		break;
		case 2:
		{
			valPlus = interface_rhou[edge][locationPlus];
			valMinus = interface_rhou[edge][locationMinus];
		}
		break;
		case 3:
		{
			valPlus = interface_rhov[edge][locationPlus];
			valMinus = interface_rhov[edge][locationMinus];
		}
		break;
		case 4:
		{
			valPlus = interface_rhoE[edge][locationPlus];
			valMinus = interface_rhoE[edge][locationMinus];
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
		inte = math::surfaceInte(F, edge, elem);

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
