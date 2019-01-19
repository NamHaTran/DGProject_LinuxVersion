#include "DGProcLib.h"
#include "DGMath.h"
#include "varDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "dynamicVarDeclaration.h"
#include <tuple>
#include "DGBCsLib.h"
#include <algorithm>
#include "DGIOLib.h"
#include <iostream>

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

		/*1D Jacobi*/
		for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
		{
			std::tie(meshVar::J1D[iedge][0], meshVar::J1D[iedge][1]) = math::J1DCal(iedge);
		}
	}

	void basisFcParam()
	{
		double a(0.0), b(0.0);
		for (int na = 0; na <= mathVar::nGauss; na++)
		{
			for (int nb = 0; nb <= mathVar::nGauss; nb++)
			{
				a = mathVar::GaussPts[na][nb][0];
				b = mathVar::GaussPts[na][nb][1];
				math::basisFc(a, b);
				math::dBasisFc(a, b);
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					mathVar::BPts[order][na][nb] = mathVar::B[order];
					mathVar::dBaPts[order][na][nb] = mathVar::dBa[order];
					mathVar::dBbPts[order][na][nb] = mathVar::dBb[order];
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
		double xCG(0.0), yCG(0.0);
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
			//1. Calculate cell area (cell size)
			meshVar::cellSize[nelem] = math::geometricOp::calcPolygonArea(xCoor, yCoor, elemType);

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
				std::tie(meshVar::geoCenter[nelem][0], meshVar::geoCenter[nelem][1]) = math::geometricOp::calcQuadCentroid(nelem, xCG, yCG, meshVar::cellSize[nelem]);
			}
		}
	}
}

namespace process
{
	void setIniValues()
	{
		iniValues::rhoIni = iniValues::pIni / (material::R*iniValues::TIni);
		material::Cp = material::R*material::gamma / (material::gamma - 1);
		material::Cv = material::Cp - material::R;
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
	}

	std::vector<double> calcIniValues(double iniVal, int element)
	{
		std::vector<std::vector<double>> matrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));
		std::vector<double> RHS(mathVar::orderElem + 1, 0.0);
		std::vector<double> iniVector(mathVar::orderElem + 1, 0.0);

		matrix = process::calculateStiffMatrix(element);
		RHS = process::calcIniValuesRHS(element, iniVal);
		iniVector = math::SolveSysEqs(matrix, RHS);
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
					math::basisFc(a, b);
					matrix[na][nb] = mathVar::B[order];
				}
			}
			Out[order] = math::volumeInte(matrix, element)*iniVal;
		}
		return Out;
	}

	namespace auxEq
	{
		void solveAuxEquation()
		{
			std::vector<std::vector<double>> StiffMatrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));
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
				StiffMatrix = process::calculateStiffMatrix(nelement);
				
				//2) Calculate Right hand side terms
				process::auxEq::CalcRHSTerm(nelement, rhoRHSTermOxDir, rhoRHSTermOyDir, rhouRHSTermOxDir, rhouRHSTermOyDir, rhovRHSTermOxDir, rhovRHSTermOyDir, rhoERHSTermOxDir, rhoERHSTermOyDir);

				//3) Solve for auxilary variables
				//Ox direction
				rhoxVector = math::SolveSysEqs(StiffMatrix, rhoRHSTermOxDir);
				rhouxVector = math::SolveSysEqs(StiffMatrix, rhouRHSTermOxDir);
				rhovxVector = math::SolveSysEqs(StiffMatrix, rhovRHSTermOxDir);
				rhoExVector = math::SolveSysEqs(StiffMatrix, rhoERHSTermOxDir);

				//Oy direction
				rhoyVector = math::SolveSysEqs(StiffMatrix, rhoRHSTermOyDir);
				rhouyVector = math::SolveSysEqs(StiffMatrix, rhouRHSTermOyDir);
				rhovyVector = math::SolveSysEqs(StiffMatrix, rhovRHSTermOyDir);
				rhoEyVector = math::SolveSysEqs(StiffMatrix, rhoERHSTermOyDir);

				//4) Save results to auxilary variables array

				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					rhoX[nelement][order] = rhoxVector[order];
					rhoY[nelement][order] = rhoyVector[order];
					rhouX[nelement][order] = rhouxVector[order];
					rhouY[nelement][order] = rhouyVector[order];
					rhovX[nelement][order] = rhovxVector[order];
					rhovY[nelement][order] = rhovyVector[order];
					rhoEX[nelement][order] = rhoExVector[order];
					rhoEY[nelement][order] = rhoEyVector[order];
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
			//Ox direction
			process::auxEq::calcVolumeIntegralTerms(element, rhoVolIntOx, rhouVolIntOx, rhovVolIntOx, rhoEVolIntOx, 1);
			//Oy direction
			process::auxEq::calcVolumeIntegralTerms(element, rhoVolIntOy, rhouVolIntOy, rhovVolIntOy, rhoEVolIntOy, 2);

			/*2. Calculate surface integral term*/
			//Ox direction
			process::auxEq::calcSurfaceIntegralTerms(element, rhoSurfIntOx, rhouSurfIntOx, rhovSurfIntOx, rhoESurfIntOx, 1);
			//Oy direction
			process::auxEq::calcSurfaceIntegralTerms(element, rhoSurfIntOy, rhouSurfIntOy, rhovSurfIntOy, rhoESurfIntOy, 2);

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

		void calcVolumeIntegralTerms(int element, std::vector<double> &rhoVolInt, std::vector<double> &rhouVolInt, std::vector<double> &rhovVolInt, std::vector<double> &rhoEVolInt, int dir)
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

			for (int order = 0; order <= mathVar::orderElem; order++)
			{
				rhoVolInt[order] = process::volumeInte(element, rhoGsVol, order, dir);
				rhouVolInt[order] = process::volumeInte(element, rhouGsVol, order, dir);
				rhovVolInt[order] = process::volumeInte(element, rhovGsVol, order, dir);
				rhoEVolInt[order] = process::volumeInte(element, rhoEGsVol, order, dir);
			}
		}

		void calcSurfaceIntegralTerms(int element, std::vector<double> &rhoSurfInt, std::vector<double> &rhouSurfInt, std::vector<double> &rhovSurfInt, std::vector<double> &rhoESurfInt, int dir)
		{
			/*User's guide:
			Input array rhoRHSTerm, rhouRHSTerm, rhovRHSTerm, rhoERHSTerm have following form:
			- number of row: orderElem + 1*/

			std::vector<std::vector<double>> rhoFlux(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
				rhouFlux(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
				rhovFlux(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
				rhoEFlux(mathVar::nGauss + 1, std::vector<double>(4, 0.0));

			int elemType(auxUlti::checkType(element)), edgeName(0);
			std::vector<double> rhoFluxTemp(mathVar::nGauss + 1, 0.0),
				rhouFluxTemp(mathVar::nGauss + 1, 0.0),
				rhovFluxTemp(mathVar::nGauss + 1, 0.0),
				rhoEFluxTemp(mathVar::nGauss + 1, 0.0);

			/*1. Calculate flux of conservative variables at all Gauss points on all faces of element*/
			process::auxEq::getGaussVectorOfConserVar(element, rhoFlux, rhouFlux, rhovFlux, rhoEFlux, dir);

			/*2. Calculates surface integrals of all conservative variables at all order*/
			for (int order = 0; order <= mathVar::orderElem; order++)
			{
				for (int nface = 0; nface < elemType; nface++)
				{
					for (int nG = 0; nG <= mathVar::nGauss; nG++)
					{
						rhoFluxTemp[nG] = rhoFlux[nG][nface];
						rhouFluxTemp[nG] = rhouFlux[nG][nface];
						rhovFluxTemp[nG] = rhovFlux[nG][nface];
						rhoEFluxTemp[nG] = rhoEFlux[nG][nface];	
					}
					rhoSurfInt[order] += process::surfaceInte(element, edgeName, rhoFluxTemp, order);
					rhouSurfInt[order] += process::surfaceInte(element, edgeName, rhouFluxTemp, order);
					rhovSurfInt[order] += process::surfaceInte(element, edgeName, rhovFluxTemp, order);
					rhoESurfInt[order] += process::surfaceInte(element, edgeName, rhoEFluxTemp, order);
				}
			}
		}

		std::vector<std::vector<double>> getGaussMatrixOfConserVar(int element, int valType)
		{
			std::vector<std::vector<double>> GaussMatrix(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
			double muGs(0.0), a(0.0), b(0.0);
			for (int na = 0; na <= mathVar::nGauss; na++)
			{
				for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
					muGs = math::pointValue(element, a, b, 7, 1);
					GaussMatrix[na][nb] = math::pointValue(element, a, b, valType, 2)*muGs;
				}
			}
			return GaussMatrix;
		}

		std::vector<double> getVectorOfConserVarFluxesAtInternal(int edge, int element, int nG, double nVectorComp)
		{
			std::vector<double> Flux(4, 0.0);
			double muGsP(0.0), muGsM(0.0), valP(0.0), valM(0.0);
			std::tie(muGsP, muGsM) = math::internalSurfaceValue(edge, element, nG, 7, 1);

			//rho flux
			std::tie(valP, valM) = math::internalSurfaceValue(edge, element, nG, 1, 2);
			Flux[0] = math::numericalFluxes::auxFlux(muGsM*valM, muGsP*valP, nVectorComp);

			//rhou flux
			std::tie(valP, valM) = math::internalSurfaceValue(edge, element, nG, 2, 2);
			Flux[1] = math::numericalFluxes::auxFlux(muGsM*valM, muGsP*valP, nVectorComp);

			//rhov flux
			std::tie(valP, valM) = math::internalSurfaceValue(edge, element, nG, 3, 2);
			Flux[2] = math::numericalFluxes::auxFlux(muGsM*valM, muGsP*valP, nVectorComp);

			//rhoE flux
			std::tie(valP, valM) = math::internalSurfaceValue(edge, element, nG, 4, 2);
			Flux[3] = math::numericalFluxes::auxFlux(muGsM*valM, muGsP*valP, nVectorComp);

			return Flux;
		}

		void getGaussVectorOfConserVar(int element, std::vector<std::vector<double>> &rhoFlux, std::vector<std::vector<double>> &rhouFlux, std::vector<std::vector<double>> &rhovFlux, std::vector<std::vector<double>> &rhoEFlux, int dir)
		{
			/*User's guide:
			Input array rhoFlux, rhouFlux, rhovFlux, rhoEFlux have following form:
			- number of row: nGauss + 1
			- number of column: default is 4 (4 faces per element)*/

			int elemType(auxUlti::checkType(element)), edgeName(0);
			int faceBcType(0);
			double nVectorComp(0.0);
			std::vector<double> Flux(4, 0.0);

			for (int nface = 0; nface < elemType; nface++)
			{
				edgeName = meshVar::inedel[nface][element];

				faceBcType = auxUlti::getBCType(edgeName);
				nVectorComp = (auxUlti::getNormVectorComp(element, edgeName, dir));

				if (faceBcType == 0)  //internal edge
				{
					for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
					{
						Flux = process::auxEq::getVectorOfConserVarFluxesAtInternal(edgeName, element, nGauss, nVectorComp);
						rhoFlux[nGauss][nface] = Flux[0];
						rhouFlux[nGauss][nface] = Flux[1];
						rhovFlux[nGauss][nface] = Flux[2];
						rhoEFlux[nGauss][nface] = Flux[3];
					}
				}
				else  //boundary edge
				{
					for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
					{
						Flux = auxEqBCsImplement(element, edgeName, nGauss, nVectorComp);
						rhoFlux[nGauss][nface] = Flux[0];
						rhouFlux[nGauss][nface] = Flux[1];
						rhovFlux[nGauss][nface] = Flux[2];
						rhoEFlux[nGauss][nface] = Flux[3];
					}
				}
			}
		}
	}//end namespace auxEq

	namespace NSFEq
	{
		void solveNSFEquation()
		{
			std::vector<double> rhoError(meshVar::nelem2D, 1.0),
				rhouError(meshVar::nelem2D, 1.0),
				rhovError(meshVar::nelem2D, 1.0),
				rhoEError(meshVar::nelem2D, 1.0);

			std::vector<std::vector<double>> StiffMatrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));
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
				
				UnVector(mathVar::orderElem + 1, 0.0),
				
				timeStepArr(meshVar::nelem2D,1.0);

			double rhoRes(1.0), rhouRes(1.0), rhovRes(1.0), rhoERes(1.0);

			for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
			{
				//1) Calculate Stiff matrix
				StiffMatrix = process::calculateStiffMatrix(nelement);

				//2) Calculate Right hand side terms
				process::NSFEq::CalcRHSTerm(nelement, RHSTerm1, RHSTerm2, RHSTerm3, RHSTerm4);

				//3) Solve for derivartives of conservative variables
				ddtRhoVector = math::SolveSysEqs(StiffMatrix, RHSTerm1);
				ddtRhouVector = math::SolveSysEqs(StiffMatrix, RHSTerm2);
				ddtRhovVector = math::SolveSysEqs(StiffMatrix, RHSTerm3);
				ddtRhoEVector = math::SolveSysEqs(StiffMatrix, RHSTerm4);

				//4) Solve time marching
				//rho
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rho[nelement][order];
				}
				rhoVectorN = process::NSFEq::solveTimeMarching(ddtRhoVector, UnVector);
				//rhou
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rhou[nelement][order];
				}
				rhouVectorN = process::NSFEq::solveTimeMarching(ddtRhouVector, UnVector);
				//rhov
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rhov[nelement][order];
				}
				rhovVectorN = process::NSFEq::solveTimeMarching(ddtRhovVector, UnVector);
				//rhoE
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rhoE[nelement][order];
				}
				rhoEVectorN = process::NSFEq::solveTimeMarching(ddtRhoEVector, UnVector);

				//5) Save results to conservative variables array
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					rhoN[nelement][order] = rhoVectorN[order];
					rhouN[nelement][order] = rhouVectorN[order];
					rhovN[nelement][order] = rhovVectorN[order];
					rhoEN[nelement][order] = rhoEVectorN[order];
				}

				//6) Estimate Residuals
				//
				rhoError[nelement] = (process::Euler::localErrorEstimate(nelement, ddtRhoVector));
				rhouError[nelement] = (process::Euler::localErrorEstimate(nelement, ddtRhouVector));
				rhovError[nelement] = (process::Euler::localErrorEstimate(nelement, ddtRhovVector));
				rhoEError[nelement] = (process::Euler::localErrorEstimate(nelement, ddtRhoEVector));

				//7) Compute local time step
				timeStepArr[nelement] = process::Euler::localTimeStep(nelement);
			}
			runTime += dt;
			dt = *std::min_element(timeStepArr.begin(), timeStepArr.end());  //find min value of vector

			std::tie(rhoRes, rhouRes, rhovRes, rhoERes) = process::Euler::globalErrorEstimate(rhoError, rhouError, rhovError, rhoEError);
			IO::residualOutput(rhoRes, rhouRes, rhovRes, rhoERes);
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
		std::vector<std::vector<double>> calcGaussInviscidTerm(int element, double a, double b)
		{
			/*InviscidTerm 2D array has 4 row 2 column:
			- column 1: Ox direction
			- column 2: Oy direction*/
			std::vector<std::vector<double>> InviscidTerm(4, std::vector<double>(2, 0.0));
			double
				rhoVal(math::pointValue(element, a, b, 1, 1)),
				uVal(math::pointValue(element, a, b, 2, 1)),
				vVal(math::pointValue(element, a, b, 3, 1)),
				eVal(math::pointValue(element, a, b, 4, 1)),
				pVal(math::pointValue(element, a, b, 5, 1));

			double totalE(eVal + 0.5*(uVal*uVal + vVal * vVal));

			/*1. Ox direction*/
			std::tie(InviscidTerm[0][0], InviscidTerm[1][0], InviscidTerm[2][0], InviscidTerm[3][0]) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoVal, uVal, vVal, totalE, pVal, 1);

			/*2. Oy direction*/
			std::tie(InviscidTerm[0][1], InviscidTerm[1][1], InviscidTerm[2][1], InviscidTerm[3][1]) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoVal, uVal, vVal, totalE, pVal, 2);

			return InviscidTerm;
		}

		/*Function calculates Viscous terms at Gauss point (a, b)*/
		std::vector<std::vector<double>> calcGaussViscousTerm(int element, double a, double b)
		{
			double uVal(math::pointValue(element, a, b, 2, 1)), vVal(math::pointValue(element, a, b, 3, 1)), muVal(math::pointValue(element, a, b, 7, 1));
			/*ViscousTerm 2D array has 4 row 2 column:
			- column 1: Ox direction
			- column 2: Oy direction*/
			std::vector<std::vector<double>> ViscousTerm(4, std::vector<double>(2, 0.0));
			std::vector<std::vector<double>> StressHeatFlux(2, std::vector<double>(3, 0.0));
			std::vector<double> vectorU(4, 0.0);
			std::vector<double> vectordUx(4, 0.0);
			std::vector<double> vectordUy(4, 0.0);

			/*calculate conservative and derivative variables*/

			vectorU[0] = math::pointValue(element, a, b, 1, 2);
			vectorU[1] = math::pointValue(element, a, b, 2, 2);
			vectorU[2] = math::pointValue(element, a, b, 3, 2);
			vectorU[3] = math::pointValue(element, a, b, 4, 2);

			vectordUx[0] = math::pointAuxValue(element, a, b, 1, 1);
			vectordUy[0] = math::pointAuxValue(element, a, b, 1, 2);

			vectordUx[1] = math::pointAuxValue(element, a, b, 2, 1);
			vectordUy[1] = math::pointAuxValue(element, a, b, 2, 2);

			vectordUx[2] = math::pointAuxValue(element, a, b, 3, 1);
			vectordUy[2] = math::pointAuxValue(element, a, b, 3, 2);

			vectordUx[3] = math::pointAuxValue(element, a, b, 4, 1);
			vectordUy[3] = math::pointAuxValue(element, a, b, 4, 2);

			/*calculate stresses and heat fluxes*/
			StressHeatFlux = math::viscousTerms::calcStressTensorAndHeatFlux(muVal, vectorU, vectordUx, vectordUy);
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
				InvisGsVolX1(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				InvisGsVolY1(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),

				InvisGsVolX2(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				InvisGsVolY2(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),

				InvisGsVolX3(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				InvisGsVolY3(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),

				InvisGsVolX4(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				InvisGsVolY4(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));

			std::vector<std::vector<double>>
				ViscGsVolX2(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				ViscGsVolY2(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),

				ViscGsVolX3(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				ViscGsVolY3(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),

				ViscGsVolX4(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				ViscGsVolY4(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));

			double a(0.0), b(0.0);

			for (int na = 0; na <= mathVar::nGauss; na++)
			{
				for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
					/*A INVISCID TERMS*/
					InviscidTerms = NSFEq::calcGaussInviscidTerm(element, a, b);
					/*A1. Inviscid term on Ox direction*/
					InvisGsVolX1[na][nb] = InviscidTerms[0][0];
					InvisGsVolX2[na][nb] = InviscidTerms[1][0];
					InvisGsVolX3[na][nb] = InviscidTerms[2][0];
					InvisGsVolX4[na][nb] = InviscidTerms[3][0];
					/*A2. Inviscid term on Oy direction*/
					InvisGsVolY1[na][nb] = InviscidTerms[0][1];
					InvisGsVolY2[na][nb] = InviscidTerms[1][1];
					InvisGsVolY3[na][nb] = InviscidTerms[2][1];
					InvisGsVolY4[na][nb] = InviscidTerms[3][1];

					/*B VISCOUS TERMS*/
					ViscousTerms = NSFEq::calcGaussViscousTerm(element, a, b);
					/*B1. Viscous term on Ox direction*/
					//ViscGsVolX1[na][nb] = ViscousTerms[0][0];
					ViscGsVolX2[na][nb] = ViscousTerms[1][0];
					ViscGsVolX3[na][nb] = ViscousTerms[2][0];
					ViscGsVolX4[na][nb] = ViscousTerms[3][0];
					/*B2. Viscous term on Oy direction*/
					//ViscGsVolY1[na][nb] = ViscousTerms[0][1];
					ViscGsVolY2[na][nb] = ViscousTerms[1][1];
					ViscGsVolY3[na][nb] = ViscousTerms[2][1];
					ViscGsVolY4[na][nb] = ViscousTerms[3][1];
				}
			}

			for (int order = 0; order <= mathVar::orderElem; order++)
			{
				/*CALCULATE INTEGRALS*/
				/*A INVISCID TERMS*/
				/*- integral term of (dB/dx)*(rho*u)*/
				VolIntTerm1[order] += process::volumeInte(element, InvisGsVolX1, order, 1);
				/*- integral term of (dB/dy)*(rho*v)*/
				VolIntTerm1[order] += process::volumeInte(element, InvisGsVolY1, order, 2);

				/*- integral term of (dB/dx)*(rho*u^2 + p)*/
				VolIntTerm2[order] += process::volumeInte(element, InvisGsVolX2, order, 1);
				/*- integral term of (dB/dy)*(rho*u*v)*/
				VolIntTerm2[order] += process::volumeInte(element, InvisGsVolY2, order, 2);

				/*- integral term of (dB/dx)*(rho*u*v)*/
				VolIntTerm3[order] += process::volumeInte(element, InvisGsVolX3, order, 1);
				/*- integral term of (dB/dy)*(rho*v^2 + p)*/
				VolIntTerm3[order] += process::volumeInte(element, InvisGsVolY3, order, 2);

				/*- integral term of (dB/dx)*(rho*totalE + p)*u*/
				VolIntTerm4[order] += process::volumeInte(element, InvisGsVolX4, order, 1);
				/*- integral term of (dB/dy)*(rho*totalE + p)*v*/
				VolIntTerm4[order] += process::volumeInte(element, InvisGsVolY4, order, 2);

				/*B VISCOUS TERMS*/
				/*- integral term of (dB/dx)*(0.0)*/
				VolIntTerm1[order] += 0.0;
				/*- integral term of (dB/dy)*(0.0)*/
				VolIntTerm1[order] += 0.0;

				/*- integral term of (dB/dx)*(tau_xx)*/
				VolIntTerm2[order] += process::volumeInte(element, ViscGsVolX2, order, 1);
				/*- integral term of (dB/dy)*(tau_xy)*/
				VolIntTerm2[order] += process::volumeInte(element, ViscGsVolY2, order, 2);

				/*- integral term of (dB/dx)*(tau_xy)*/
				VolIntTerm3[order] += process::volumeInte(element, ViscGsVolX3, order, 1);
				/*- integral term of (dB/dy)*(tau_yy)*/
				VolIntTerm3[order] += process::volumeInte(element, ViscGsVolY3, order, 2);

				/*- integral term of (dB/dx)*(u*tau_xx + v*tau_xy + Qx)*/
				VolIntTerm4[order] += process::volumeInte(element, ViscGsVolX4, order, 1);
				/*- integral term of (dB/dy)*(u*tau_xy + v*tau_yy + Qy)*/
				VolIntTerm4[order] += process::volumeInte(element, ViscGsVolY4, order, 2);
				/*End volume terms===========================================================================*/
			}
			//return VolInt;
		}

		void calcSurfaceIntegralTerms(int element, std::vector<double> &SurfIntTerm1, std::vector<double> &SurfIntTerm2, std::vector<double> &SurfIntTerm3, std::vector<double> &SurfIntTerm4)
		{
			/*User's guide:
			All input array have form:
			- number of rows: orderElem*/

			int elemType(auxUlti::checkType(element)), edgeName(0);
			int faceBcType(0);

			std::vector<double> inviscFlux1Temp(mathVar::nGauss + 1, 0.0),
				inviscFlux2Temp(mathVar::nGauss + 1, 0.0),
				inviscFlux3Temp(mathVar::nGauss + 1, 0.0),
				inviscFlux4Temp(mathVar::nGauss + 1, 0.0),
				
				ViscFlux1Temp(mathVar::nGauss + 1, 0.0),
				ViscFlux2Temp(mathVar::nGauss + 1, 0.0),
				ViscFlux3Temp(mathVar::nGauss + 1, 0.0),
				ViscFlux4Temp(mathVar::nGauss + 1, 0.0);

			std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));

			std::vector<std::vector<double>>
				inviscFlux1(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
				inviscFlux2(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
				inviscFlux3(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
				inviscFlux4(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),

				ViscFlux1(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
				ViscFlux2(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
				ViscFlux3(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
				ViscFlux4(mathVar::nGauss + 1, std::vector<double>(4, 0.0));

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
						inviscFlux1[nGauss][nface] = Fluxes[0][0];
						inviscFlux2[nGauss][nface] = Fluxes[1][0];
						inviscFlux3[nGauss][nface] = Fluxes[2][0];
						inviscFlux4[nGauss][nface] = Fluxes[3][0];

						ViscFlux1[nGauss][nface] = Fluxes[0][1];
						ViscFlux2[nGauss][nface] = Fluxes[1][1];
						ViscFlux3[nGauss][nface] = Fluxes[2][1];
						ViscFlux4[nGauss][nface] = Fluxes[3][1];
					}
				}
				else  //boundary edge
				{
					for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
					{
						Fluxes = NSFEqBCsImplement(element, edgeName, nGauss);
						inviscFlux1[nGauss][nface] = Fluxes[0][0];
						inviscFlux2[nGauss][nface] = Fluxes[1][0];
						inviscFlux3[nGauss][nface] = Fluxes[2][0];
						inviscFlux4[nGauss][nface] = Fluxes[3][0];

						ViscFlux1[nGauss][nface] = Fluxes[0][1];
						ViscFlux2[nGauss][nface] = Fluxes[1][1];
						ViscFlux3[nGauss][nface] = Fluxes[2][1];
						ViscFlux4[nGauss][nface] = Fluxes[3][1];
					}
				}
			}

			for (int order = 0; order <= mathVar::orderElem; order++)
			{
				for (int nface = 0; nface < elemType; nface++)
				{
					for (int nG = 0; nG <= mathVar::nGauss; nG++)
					{
						inviscFlux1Temp[nG] = inviscFlux1[nG][nface];
						inviscFlux2Temp[nG] = inviscFlux2[nG][nface];
						inviscFlux3Temp[nG] = inviscFlux3[nG][nface];
						inviscFlux4Temp[nG] = inviscFlux4[nG][nface];

						ViscFlux1Temp[nG] = ViscFlux1[nG][nface];
						ViscFlux2Temp[nG] = ViscFlux1[nG][nface];
						ViscFlux3Temp[nG] = ViscFlux1[nG][nface];
						ViscFlux4Temp[nG] = ViscFlux1[nG][nface];
					}
					SurfIntTerm1[order] += process::surfaceInte(element, edgeName, inviscFlux1Temp, order) + process::surfaceInte(element, edgeName, ViscFlux1Temp, order);
					SurfIntTerm2[order] += process::surfaceInte(element, edgeName, inviscFlux2Temp, order) + process::surfaceInte(element, edgeName, ViscFlux2Temp, order);
					SurfIntTerm3[order] += process::surfaceInte(element, edgeName, inviscFlux3Temp, order) + process::surfaceInte(element, edgeName, ViscFlux3Temp, order);
					SurfIntTerm4[order] += process::surfaceInte(element, edgeName, inviscFlux4Temp, order) + process::surfaceInte(element, edgeName, ViscFlux4Temp, order);
				}
			}
			//return SurInt;
		}

		std::vector<std::vector<double>> getGaussVectorOfConserVarFluxesAtInternal(int edgeName, int element, int nGauss)
		{
			std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
			std::vector<double> UMinus(4, 0.0),
				UPlus(4, 0.0),
				dUXMinus(4, 0.0),
				dUXPlus(4, 0.0),
				dUYMinus(4, 0.0),
				dUYPlus(4, 0.0),
				normalVector(2, 0.0);

			/*Normal vector*/
			normalVector[0] = auxUlti::getNormVectorComp(element, edgeName, 1);
			normalVector[1] = auxUlti::getNormVectorComp(element, edgeName, 2);

			/*INVISCID TERMS*/
			std::tie(UPlus[0], UMinus[0]) = math::internalSurfaceValue(edgeName, element, nGauss, 1, 2);  //rho
			std::tie(UPlus[1], UMinus[1]) = math::internalSurfaceValue(edgeName, element, nGauss, 2, 2);  //rhou
			std::tie(UPlus[2], UMinus[2]) = math::internalSurfaceValue(edgeName, element, nGauss, 3, 2);  //rhov
			std::tie(UPlus[3], UMinus[3]) = math::internalSurfaceValue(edgeName, element, nGauss, 4, 2);  //rhoE
			
			/*VISCOUS TERMS*/
			std::tie(dUXPlus[0], dUXMinus[0]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 1, 1);  //drhox
			std::tie(dUYPlus[0], dUYMinus[0]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 1, 2);  //drhoy

			std::tie(dUXPlus[1], dUXMinus[1]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 2, 1);  //drhoux
			std::tie(dUYPlus[1], dUYMinus[1]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 2, 2);  //drhoux

			std::tie(dUXPlus[2], dUXMinus[2]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 3, 1);  //drhovx
			std::tie(dUYPlus[2], dUYMinus[2]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 3, 2);  //drhovx

			std::tie(dUXPlus[3], dUXMinus[3]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 4, 1);  //drhoEx
			std::tie(dUYPlus[3], dUYMinus[3]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 4, 2);  //drhoEx
			
			/*Calculate fluxes*/
			Fluxes = math::numericalFluxes::NSFEqFluxFromConserVars(UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, normalVector);
			return Fluxes;
		}

		std::vector<double> solveTimeMarching(std::vector<double> &ddtArr, std::vector<double> &UnArr)
		{
			std::vector<double> OutArr(mathVar::orderElem + 1, 0.0);

			if (systemVar::ddtScheme==1) //Euler scheme
			{
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					OutArr[order] = dt * ddtArr[order] + UnArr[order];
				}
			}
			return OutArr;
		}
	}

	namespace Euler
	{
		double localTimeStep(int element)
		{
			//Here, xC and yC are coordinates of center of standard elements
			double deltaT(0.0), xC(-1.0 / 3.0), yC(1.0 / 3.0), size(0.0);

			if (auxUlti::checkType(element) == 4)
			{
				xC = 0.0;
				yC = 0.0;
			}

			//std::tie(xC, yC, size) = auxUlti::getCellMetrics(element);
			double uVal(math::pointValue(element, xC, yC, 2, 1)),
				vVal(math::pointValue(element, xC, yC, 3, 1)), velocity(0.0),
				TVal(math::pointValue(element, xC, yC, 6, 1)), aSound(0.0), LocalMach(0.0);
			if (TVal<=0 || TVal != TVal)
			{
				std::cout << "Negative T is detected at element " << element + meshVar::nelem1D + 1 << std::endl;
				TVal = iniValues::TIni;
				system("pause");
				double rhoVal(math::pointValue(element, xC, yC, 1, 2)),
					rhouVal(math::pointValue(element, xC, yC, 2, 2)),
					rhovVal(math::pointValue(element, xC, yC, 3, 2)),
					rhoEVal(math::pointValue(element, xC, yC, 4, 2));
				std::cout << "value rho, rhou, rhov, rhoe = " << rhoVal << ", " << rhouVal << ", " << rhovVal << ", " << rhoEVal << std::endl;
			}

			velocity = sqrt(pow(uVal, 2) + pow(vVal, 2));
			aSound = math::CalcSpeedOfSound(TVal);
			LocalMach = velocity / aSound;
			size = meshVar::cellSize[element];
			double muVal(math::CalcVisCoef(TVal));

			deltaT = (1.0 / pow((mathVar::orderElem + 1 + 1), 2))*(size*systemVar::CFL) / (fabs(velocity) + (aSound / LocalMach) + (muVal / size));

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
			math::basisFc(xC, yC);
			for (int order = 0; order <= mathVar::orderElem; order++)
			{
				errorVal += ddtArr[order] * mathVar::B[order];
			}

			return errorVal;
		}

		std::tuple <double, double, double, double> globalErrorEstimate(std::vector<double> &RhoError, std::vector<double> &RhouError, std::vector<double> &RhovError, std::vector<double> &RhoEError)
		{
			double rhoRes(*std::max_element(RhoError.begin(), RhoError.end())),
				rhouRes(*std::max_element(RhouError.begin(), RhouError.end())),
				rhovRes(*std::max_element(RhovError.begin(), RhovError.end())),
				rhoERes(*std::max_element(RhoEError.begin(), RhoEError.end()));
			return std::make_tuple(rhoRes, rhouRes, rhovRes, rhoERes);
		}
	}

	namespace limiter
	{
		void limiter()
		{
			double theta1(0.0), theta2(0.0);
			if (systemVar::limiter==1)  //positivity preserving
			{
				for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
				{
					std::tie(theta1, theta2) = limiter::Pp::calcPpLimiterCoef(nelem);
					theta1Arr[nelem] = theta1;
					theta2Arr[nelem] = theta2;
				}
			}
			else if (systemVar::limiter == 0)  //No limiter
			{
				for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
				{
					theta1Arr[nelem] = 1.0;
					theta2Arr[nelem] = 1.0;
				}
			}
		}

		namespace Pp
		{
			//Function calculates coefficients of positivity preserving limiter
			std::tuple<double, double> calcPpLimiterCoef(int element)
			{
				double meanRho(0.0), minRho(0.0), theta1(0.0), theta2(0.0), omega(0.0);
				int elemType(auxUlti::checkType(element));

				double meanRhou(0.0), meanRhov(0.0), meanRhoE(0.0);

				/*Note: according to Kontzialis et al, positivity preserving limiter for quadrilateral element, which is presented on Zhang's paper,
				shown a very good effect on results. Because of that, Zhang's limiter is used in this code for both triangular and quadrilateral elements*/

				//Find theta1
				minRho = math::limiter::calcMinRhoQuad(element);
				meanRho = math::limiter::calcMeanConsvVarQuad(element, 1);

				meanRhou = math::limiter::calcMeanConsvVarQuad(element, 2);
				meanRhov = math::limiter::calcMeanConsvVarQuad(element, 3);
				meanRhoE = math::limiter::calcMeanConsvVarQuad(element, 4);
				double meanT(math::CalcTFromConsvVar(meanRho, meanRhou, meanRhov, meanRhoE));
				double meanP(math::CalcP(meanT, meanRho));
				
				//Compute theta1
				std::tie(theta1, omega) = math::limiter::calcTheta1Coeff(meanRho, minRho, meanP);

				//Find theta2
				std::vector<double> vectort(2 * (mathVar::nGauss + 1) * (mathVar::nGauss + 1), 0.0);
				int index(0);

				//meanRhou = math::limiter::calcMeanConsvVarQuad(element, 2);
				//meanRhov = math::limiter::calcMeanConsvVarQuad(element, 3);
				//meanRhoE = math::limiter::calcMeanConsvVarQuad(element, 4);

				//Save mean values to arrays
				meanVals[element][0] = meanRho;
				meanVals[element][1] = meanRhou;
				meanVals[element][2] = meanRhov;
				meanVals[element][3] = meanRhoE;

				for (int na = 0; na <= mathVar::nGauss; na++)
				{
					for (int nb = 0; nb <= mathVar::nGauss; nb++)
					{
						vectort[index] = math::limiter::calcTheta2Coeff(element, na, nb, theta1, omega, meanRho, meanRhou, meanRhov, meanRhoE, 1);
						index++;

						vectort[index] = math::limiter::calcTheta2Coeff(element, na, nb, theta1, omega, meanRho, meanRhou, meanRhov, meanRhoE, 2);
						index++;
					}
				}
				theta2 = *std::min_element(vectort.begin(), vectort.end());  //find min value of vector

				return std::make_tuple(theta1, theta2);
			}
		}

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
			math::basisFc(a, b);  //mathVar::B is changed after this command line is excuted
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
		for (int na = 0; na <= mathVar::nGauss; na++)
		{
			for (int nb = 0; nb <= mathVar::nGauss; nb++)
			{
				B1 = mathVar::BPts[order1][na][nb];
				B2 = mathVar::BPts[order2][na][nb];
				FMatrix[na][nb] = B1 * B2;
			}
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