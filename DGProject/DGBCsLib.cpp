#include "DGBCsLib.h"
#include <vector>
#include "DGMath.h"
#include "VarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <tuple>
#include "DGMessagesLib.h"

/*NOTES:
- All boundary functions must return numerical fluxes (rho, rhou, rhov, rhoE fluxes) at surface Gauss points!
- Method index: 1: weak weak Riemann, 2: weak Prescribed

Boundary conditions compatibility
|U					|T					|p					||advectionBC		|diffusionBC		|auxilaryBC			|
+-------------------+-------------------+-------------------++------------------+-------------------+-------------------+
|1. inOutFlow		|1. inOutFlow		|1. inOutFlow		||inOutFlow			|interiorExtrapolate|inOutFlow			|
|	value u v w		|	value T			|	value p			||					|					|					|
+-------------------+-------------------+-------------------++------------------+-------------------+-------------------+
|2. noSlip			|2. WallIsothermal	|2. zeroGradient	||noSlipIsoThermal	|interiorExtrapolate|noSlipIsoThermal	|
|					|	value T			|					||					|					|					|
+-------------------+-------------------+-------------------++------------------+-------------------+-------------------+
|2. noSlip			|3. WallAdiabatic	|2. zeroGradient	||noSlipAdiabatic	|interiorExtrapolate|noSlipAdiabatic	|
+-------------------+-------------------+-------------------++------------------+-------------------+-------------------+
|3.	fixedValue		|4. fixedValue		|3. fixedValue		||fixedValues		|interiorExtrapolate|fixedValues		|
|	value u v w		|	value T			|	value p			||					|					|					|
+-------------------+-------------------+-------------------++------------------+-------------------+-------------------+
|4.	zeroGradient	|5. zeroGradient	|2. zeroGradient	||zeroGradient		|interiorExtrapolate|zeroGradient		|
|	value u v w		|	value T			|	value p			||					|					|					|
+-------------------+-------------------+-------------------++------------------+-------------------+-------------------+
*/

std::vector <double> auxEqBCsImplement(int element, int edge, int nG, double n)
{
	std::vector<double> Fluxes(4, 0.0);
	std::vector<double> MinusVal(4, 0.0);
	double muPlus(0.0), muMinus(0.0),
		rhoPlus(0.0),
		rhoMinus(0.0),
		rhouPlus(0.0),
		rhouMinus(0.0),
		rhovPlus(0.0),
		rhovMinus(0.0),
		rhoEPlus(0.0),
		rhoEMinus(0.0);
	double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));

	std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
	rhoPlus = math::pointValue(element, a, b, 1, 2);
	rhouPlus = math::pointValue(element, a, b, 2, 2);
	rhovPlus = math::pointValue(element, a, b, 3, 2);
	rhoEPlus = math::pointValue(element, a, b, 4, 2);
	muPlus = math::pointValue(element, a, b, 7, 1);
	/*NOTE: can use math::SurfaceValueFromMaster to calculate plus values, but it will decrease perfomance*/

	/*Apply boundary condition*/
	MinusVal = BCSupportFncs::boundaryMinusValsCalculator(3, element, edge, rhoPlus, rhouPlus, rhovPlus, rhoEPlus, nx, ny);
	rhoMinus = MinusVal[0];
	rhouMinus = MinusVal[1];
	rhovMinus = MinusVal[2];
	rhoEMinus = MinusVal[3];
	double TMinus(math::CalcTFromConsvVar(rhoMinus, rhouMinus, rhovMinus, rhoEMinus));
	muMinus = math::CalcVisCoef(TMinus);

	/*Calculate fluxes*/
	Fluxes[0] = math::numericalFluxes::auxFlux(rhoMinus*muMinus, rhoPlus*muPlus, n);
	Fluxes[1] = math::numericalFluxes::auxFlux(rhouMinus*muMinus, rhouPlus*muPlus, n);
	Fluxes[2] = math::numericalFluxes::auxFlux(rhovMinus*muMinus, rhovPlus*muPlus, n);
	Fluxes[3] = math::numericalFluxes::auxFlux(rhoEMinus*muMinus, rhoEPlus*muPlus, n);
	return Fluxes;
}

std::vector<std::vector<double>> NSFEqBCsImplement(int element, int edge, int nG)
{
	/*Fluxes array has the following form:
	- column 0: advective fluxes
	- column 1: diffusive fluxes*/
	std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
	std::vector<double> UMinus(4, 0.0),
		UPlus(4, 0.0),
		dUXMinus(4, 0.0),
		dUXPlus(4, 0.0),
		dUYMinus(4, 0.0),
		dUYPlus(4, 0.0),
		normalVector(2, 0.0);

	double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));

	normalVector[0] = nx;
	normalVector[1] = ny;

	std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
	UPlus[0] = math::pointValue(element, a, b, 1, 2);
	UPlus[1] = math::pointValue(element, a, b, 2, 2);
	UPlus[2] = math::pointValue(element, a, b, 3, 2);
	UPlus[3] = math::pointValue(element, a, b, 4, 2);

	dUXPlus[0] = math::pointAuxValue(element, a, b, 1, 1);
	dUXPlus[1] = math::pointAuxValue(element, a, b, 2, 1);
	dUXPlus[2] = math::pointAuxValue(element, a, b, 3, 1);
	dUXPlus[3] = math::pointAuxValue(element, a, b, 4, 1);

	dUYPlus[0] = math::pointAuxValue(element, a, b, 1, 2);
	dUYPlus[1] = math::pointAuxValue(element, a, b, 2, 2);
	dUYPlus[2] = math::pointAuxValue(element, a, b, 3, 2);
	dUYPlus[3] = math::pointAuxValue(element, a, b, 4, 2);
	/*NOTE: can use math::SurfaceValueFromMaster to calculate plus values, but it will decrease perfomance*/

	/*Apply boundary condition*/
	UMinus = BCSupportFncs::boundaryMinusValsCalculator(1, element, edge, UPlus[0], UPlus[1], UPlus[2], UPlus[3], nx, ny);

	dUXMinus = BCSupportFncs::boundaryMinusValsCalculator(2, element, edge, dUXPlus[0], dUXPlus[1], dUXPlus[2], dUXPlus[3], 0.0, 0.0);
	dUYMinus = BCSupportFncs::boundaryMinusValsCalculator(2, element, edge, dUYPlus[0], dUYPlus[1], dUYPlus[2], dUYPlus[3], 0.0, 0.0);

	/*Calculate fluxes*/
	Fluxes = math::numericalFluxes::NSFEqFluxFromConserVars(UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, normalVector);
	return Fluxes;
}

namespace BCSupportFncs
{
	std::vector <double> boundaryMinusValsCalculator(int option, int element, int edge, double rhoP, double rhouP, double rhovP, double rhoEP, double nx, double ny)
	{
		/*If option = 2 (diffusion boundary condition), input values (rhoP, rhouP, rhovP, rhoEP) are derivative type.*/
		int edgeGrp(auxUlti::getGrpOfEdge(edge));
		std::vector<double> MinusVal(4, 0.0);
		int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]), pType(bcValues::pBcType[edgeGrp - 1]), method(meshVar::BoundaryType[edgeGrp - 1][2]);

		if (option==1 || option == 3)  //advection boundary condition
		{
			if (UType == 1 && TType == 1 && pType == 1)
			{
				MinusVal = advectionBCs::patch::inOutFlow(method, edgeGrp, rhoP, rhouP, rhovP, rhoEP, nx, ny);  //working
			}
			else if (UType == 2 && TType == 2 && pType == 2)
			{
				MinusVal = advectionBCs::wall::noSlipIsoThermal(method, edgeGrp, rhoP, rhouP, rhovP, rhoEP);
			}
			else if (UType == 2 && TType == 3 && pType == 2)
			{
				MinusVal = advectionBCs::wall::noSlipAdiabatic(method, element, edgeGrp, rhoP, rhouP, rhovP, rhoEP);  //working
			}
			else if (UType == 3 && TType == 4 && pType == 3)
			{
				MinusVal = strongBCs::fixedValues(rhoP, rhouP, rhovP, rhoEP, edgeGrp);  //working
			}
			else if (UType == 4 && TType == 5 && pType == 2)
			{
				MinusVal = advectionBCs::patch::zeroGradient(element, rhoP, rhouP, rhovP, rhoEP);  //working
			}
			else
			{
				std::string errorStr = message::BcCompatibleError(edgeGrp);
				message::writeLog(systemVar::pwd, systemVar::caseName, errorStr);
			}
		}
		else if (option==2)  //diffusion boundary condition
		{
			MinusVal = diffusionBCs::interiorExtrapolate(rhoP, rhouP, rhovP, rhoEP);
		}
		/*else if (option==3)  //auxilary boundary condition
		{
			if (UType==1 && TType==1 && pType==1)
			{
				//MinusVal = auxilaryBCs::patch::inOutFlow(rhoP, rhouP, rhovP, rhoEP);
				MinusVal = advectionBCs::patch::inOutFlow(method, edgeGrp, rhoP, rhouP, rhovP, rhoEP, nx, ny);
			}
			else if (UType == 2 && TType == 2 && pType == 2)
			{
				MinusVal = advectionBCs::wall::noSlipIsoThermal(method, edgeGrp, rhoP, rhouP, rhovP, rhoEP);
			}
			else if (UType == 2 && TType == 3 && pType == 2)
			{
				MinusVal = advectionBCs::wall::noSlipAdiabatic(method, edgeGrp, rhoP, rhouP, rhovP, rhoEP);
			}
			else if (UType == 3 && TType == 4 && pType == 3)
			{
				MinusVal = strongBCs::fixedValues(rhoP, rhouP, rhovP, rhoEP, edgeGrp);
			}
			else if (UType == 4 && TType == 5 && pType == 2)
			{
				MinusVal = advectionBCs::patch::zeroGradient(rhoP, rhouP, rhovP, rhoEP);
			}
			else
			{
				std::string errorStr = message::BcCompatibleError(edgeGrp);
				message::writeLog(systemVar::pwd, systemVar::caseName, errorStr);
			}
		}
		*/
		return MinusVal;
	}

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

		normUMag = math::vectorDotProduct(U,normVector);

		if (normUMag>0)
		{
			inflow = false;
		}
		return inflow;
	}

	std::tuple<double, double> calcInOutFlowVelocity(double vExternx, double vExterny, double nx, double ny, double vBCMag)
	{
		double vBCx(0.0), vBCy(0.0);
		std::vector<double> vExternVector(2, 0.0), normVector(2, 0.0);
		normVector[0] = nx;
		normVector[1] = ny;
		vExternVector[0] = vExternx;
		vExternVector[1] = vExterny;
		double vExternMag(math::vectorDotProduct(vExternVector, normVector));
		vBCx = vExternx + (vBCMag - vExternMag)*nx;
		vBCy = vExterny + (vBCMag - vExternMag)*ny;
		return std::make_tuple(vBCx, vBCy);
	}
}

namespace advectionBCs
{
	namespace wall
	{
		std::vector <double> noSlipIsoThermal(int method, int edgeGrp, double rhoP, double rhouP, double rhovP, double rhoEP)
		{
			std::vector<double> MinusVal(4, 0.0);

			if (method == 1)  //weak Riemann
			{
				double rhoEBC(rhoP*(material::Cv*bcValues::TBC[edgeGrp - 1]));
				MinusVal[0] = rhoP;
				MinusVal[1] = -rhouP;
				MinusVal[2] = -rhovP;
				//MinusVal[3] = rhoEP;
				MinusVal[2] = rhoEBC;
			}
			else if (method == 2)  //weak Prescribed
			{
				double rhoBC(rhoP), rhouBC(0.0), rhovBC(0.0), rhoEBC(rhoP*(material::Cv*bcValues::TBC[edgeGrp - 1]));
				MinusVal[0] = 2 * rhoBC - rhoP;
				MinusVal[1] = 2 * rhouBC - rhouP;
				MinusVal[2] = 2 * rhovBC - rhovP;
				MinusVal[3] = 2 * rhoEBC - rhoEP;
			}
			return MinusVal;
		}

		std::vector <double> noSlipAdiabatic(int method, int element, int edgeGrp, double rhoP, double rhouP, double rhovP, double rhoEP)
		{
			std::vector<double> MinusVal(4, 0.0);
			double rhoBC(rhoP),
				rhouBC(0.0),
				rhovBC(0.0),
				rhoEBC(rhoEP);

			if (method == 1)  //weak Riemann
			{
				MinusVal[0] = rhoP;
				MinusVal[1] = -rhouP;
				MinusVal[2] = -rhovP;
				MinusVal[3] = rhoEP;
			}
			else if (method == 2)  //weak Prescribed
			{
				//double rhoBC(rhoP), rhouBC(0.0), rhovBC(0.0), rhoEBC(rhoEP);  //value of rhoEBC = rhoP*(material::Cv*bcValues::TBC[edgeGrp - 1]), it's usage is not verified!!
				MinusVal[0] = 2 * rhoBC - rhoP;
				MinusVal[1] = 2 * rhouBC - rhouP;
				MinusVal[2] = 2 * rhovBC - rhovP;
				MinusVal[3] = 2 * rhoEBC - rhoEP;
			}
			return MinusVal;
		}
	}

	namespace patch
	{
		std::vector <double> inOutFlow(int method, int edgeGrp, double rhoP, double rhouP, double rhovP, double rhoEP, double nx, double ny)
		{
			std::vector<double> MinusVal(4, 0.0);
			double uPlus(rhouP / rhoP), vPlus(rhovP / rhoP);
			bool inflow(BCSupportFncs::checkInflow(uPlus, vPlus, nx, ny));
			double TInternal(math::CalcTFromConsvVar(rhoP, rhouP, rhovP, rhoEP));
			//bool subsonic(auxUlti::checkSubSonicLocally(TInternal, uPlus, vPlus));// this line checks subsonic locally, not generally

			if (method==1)  //weak Riemann
			{
				double pInternal(0);
				pInternal = rhoP * material::R*TInternal;
				if (inflow==true)
				{
					//Apply weak Riemann infinite value
					MinusVal[0] = bcValues::pBC[edgeGrp - 1] / (material::R*bcValues::TBC[edgeGrp - 1]);
					MinusVal[1] = MinusVal[0] * bcValues::uBC[edgeGrp - 1];
					MinusVal[2] = MinusVal[0] * bcValues::vBC[edgeGrp - 1];
					MinusVal[3] = MinusVal[0] * (bcValues::TBC[edgeGrp - 1] * material::Cv + 0.5*(pow(bcValues::uBC[edgeGrp - 1],2) + pow(bcValues::vBC[edgeGrp - 1], 2)));
				}
				else if (inflow == false)
				{
					//Apply Partially non-reflective pressure outflow
					if (refValues::subsonic == true)
					{
						MinusVal[0] = rhoP;
						MinusVal[1] = rhouP;
						MinusVal[2] = rhovP;
						MinusVal[3] = (2 * bcValues::pBC[edgeGrp - 1] - pInternal) / (material::gamma - 1) + 0.5*rhoP*(pow(uPlus, 2) + pow(vPlus, 2));
					}
					else if (refValues::subsonic == false)
					{
						MinusVal[0] = rhoP;
						MinusVal[1] = rhouP;
						MinusVal[2] = rhovP;
						MinusVal[3] = rhoEP;
					}
				}
			}
			else if (method==2)  //weal Prescribed
			{
				double RPlus(0.0), RMinus(0.0), rhoBc(0.0), rhouBc(0.0), rhovBc(0.0), rhoEBc(0.0),
					pBc(0.0), TBc(0.0), uBc(0.0), vBc(0.0);
				double vInternalNormMag(0.0), vExternalNormMag(0.0), cInternal(0.0), SBc(0.0);
				std::vector<double> vInternalVector(2, 0.0), normVector(2, 0.0);
				double uExternal(bcValues::uBC[edgeGrp - 1]), vExternal(bcValues::vBC[edgeGrp - 1]), cExternal(0.0);
				std::vector<double> vExternalVector(2, 0.0);
				double VBc(0.0), cBc(0.0), rhoExternal(0.0);
				
				normVector[0] = nx;
				normVector[1] = ny;

				rhoExternal = bcValues::pBC[edgeGrp - 1] / (material::R*bcValues::TBC[edgeGrp - 1]);
				vExternalVector[0] = uExternal;
				vExternalVector[1] = vExternal;
				vExternalNormMag = math::vectorDotProduct(vExternalVector, normVector);
				cExternal = math::CalcSpeedOfSound(bcValues::TBC[edgeGrp - 1]);

				vInternalVector[0] = uPlus;
				vInternalVector[1] = vPlus;
				vInternalNormMag = math::vectorDotProduct(vInternalVector, normVector);
				cInternal = math::CalcSpeedOfSound(TInternal);

				if (inflow == true)  //inflow boundary condtion
				{
					if (refValues::subsonic == true)
					{
						RPlus = vInternalNormMag + (2 * cInternal) / (material::gamma - 1);
					}
					else if (refValues::subsonic == false)
					{
						RPlus = vExternalNormMag + (2 * cExternal) / (material::gamma - 1);
					}
					RMinus = vExternalNormMag - (2 * cExternal) / (material::gamma - 1);

					VBc = 0.5*(RPlus + RMinus);
					cBc = 0.25*(material::gamma - 1)*(RPlus - RMinus);
					SBc = pow(cExternal, 2) / (material::gamma*pow(rhoExternal, material::gamma - 1));
					rhoBc = (pow(cBc, 2)) / (material::gamma*SBc);
					pBc = rhoBc * pow(cBc, 2) / material::gamma;
					TBc = pBc / (rhoBc*material::R);
					std::tie(uBc, vBc) = BCSupportFncs::calcInOutFlowVelocity(uExternal, vExternal, nx, ny, VBc);
				}
				else if (inflow == false)  //outflow boundary condition
				{
					RPlus = vInternalNormMag + (2 * cInternal) / (material::gamma - 1);
					if (refValues::subsonic == true)
					{
						RMinus = vExternalNormMag - (2 * cExternal) / (material::gamma - 1);
					}
					else if (refValues::subsonic == false)
					{
						RMinus = vInternalNormMag - (2 * cInternal) / (material::gamma - 1);
					}
					VBc = 0.5*(RPlus + RMinus);
					cBc = 0.25*(material::gamma - 1)*(RPlus - RMinus);
					SBc = pow(cInternal, 2) / (material::gamma*pow(rhoP, material::gamma - 1));
					rhoBc = (pow(cBc, 2)) / (material::gamma*SBc);
					pBc = rhoBc * pow(cBc, 2) / material::gamma;
					TBc = pBc / (rhoBc*material::R);
					std::tie(uBc, vBc) = BCSupportFncs::calcInOutFlowVelocity(uPlus, vPlus, nx, ny, VBc);
				}
				rhouBc = rhoBc * uBc;
				rhovBc = rhoBc * vBc;
				rhoEBc = rhoBc * (material::Cv*TBc + 0.5*(pow(uBc,2) + pow(vBc,2)));

				MinusVal[0] = 2 * rhoBc - rhoP;
				MinusVal[1] = 2 * rhouBc - rhouP;
				MinusVal[2] = 2 * rhovBc - rhovP;
				MinusVal[3] = 2 * rhoEBc - rhoEP;
			}
			else if (method == 3)  //strong BCS
			{
				if (inflow == true)
				{
					/*
					MinusVal[0] = bcValues::pBC[edgeGrp - 1] / (material::R*bcValues::TBC[edgeGrp - 1]);
					MinusVal[1] = MinusVal[0] * bcValues::uBC[edgeGrp - 1];
					MinusVal[2] = MinusVal[0] * bcValues::vBC[edgeGrp - 1];
					MinusVal[3] = MinusVal[0] * (bcValues::TBC[edgeGrp - 1] * material::Cv + 0.5*(pow(bcValues::uBC[edgeGrp - 1], 2) + pow(bcValues::vBC[edgeGrp - 1], 2)));*/

					double rhoBC(bcValues::pBC[edgeGrp - 1] / (material::R*bcValues::TBC[edgeGrp - 1]));
					double rhouBC(rhoBC * bcValues::uBC[edgeGrp - 1]), rhovBC(rhoBC * bcValues::vBC[edgeGrp - 1]),
						rhoEBC(rhoBC * (bcValues::TBC[edgeGrp - 1] * material::Cv + 0.5*(pow(bcValues::uBC[edgeGrp - 1], 2) + pow(bcValues::vBC[edgeGrp - 1], 2))));  //value of rhoEBC = rhoP*(material::Cv*bcValues::TBC[edgeGrp - 1]), it's usage is not verified!!
					MinusVal[0] = 2 * rhoBC - rhoP;
					MinusVal[1] = 2 * rhouBC - rhouP;
					MinusVal[2] = 2 * rhovBC - rhovP;
					MinusVal[3] = 2 * rhoEBC - rhoEP;
				}
				else if (inflow == false)
				{
					//Apply Partially non-reflective pressure outflow
					MinusVal[0] = rhoP;
					MinusVal[1] = rhouP;
					MinusVal[2] = rhovP;
					MinusVal[3] = rhoEP;
				}
			}
			return MinusVal;
		}

		
		std::vector <double> zeroGradient(int element, double rhoP, double rhouP, double rhovP, double rhoEP)
		{
			std::vector<double> MinusVal(4, 0.0);
			MinusVal[0] = rhoP;
			MinusVal[1] = rhouP;
			MinusVal[2] = rhovP;
			MinusVal[3] = rhoEP;
			return MinusVal;
		}
		

		/*
		std::vector <double> zeroGradient(int element, double rhoP, double rhouP, double rhovP, double rhoEP)
		{
			std::vector<double> MinusVal(4, 0.0);
			double rhoBC(math::centerValue(element, 1, 2)),
				rhouBC(math::centerValue(element, 2, 2)),
				rhovBC(math::centerValue(element, 3, 2)),
				rhoEBC(math::centerValue(element, 4, 2));

			MinusVal[0] = 2 * rhoBC - rhoP;
			MinusVal[1] = 2 * rhouBC - rhouP;
			MinusVal[2] = 2 * rhovBC - rhovP;
			MinusVal[3] = 2 * rhoEBC - rhoEP;
			return MinusVal;
		}
		*/
	}
}

namespace diffusionBCs
{
	std::vector <double> interiorExtrapolate(double drhoP, double drhouP, double drhovP, double drhoEP)
	{
		std::vector<double> MinusVal(4, 0.0);
		MinusVal[0] = drhoP;
		MinusVal[1] = drhouP;
		MinusVal[2] = drhovP;
		MinusVal[3] = drhoEP;
		return MinusVal;
	}
}

namespace auxilaryBCs
{
	namespace wall
	{
		//for auxilary equation, use the same approach as advective term
		std::vector <double> noSlipAdiabatic(double rhoP, double rhouP, double rhovP, double rhoEP)
		{
			std::vector<double> MinusVal(4, 0.0);
			double rhoBC(rhoP), rhouBC(0.0), rhovBC(0.0), rhoEBC(rhoEP);
			MinusVal[0] = 2 * rhoBC - rhoP;
			MinusVal[1] = 2 * rhouBC - rhouP;
			MinusVal[2] = 2 * rhovBC - rhovP;
			MinusVal[3] = 2 * rhoEBC - rhoEP;
			return MinusVal;
		}

		std::vector <double> noSlipIsoThermal(double rhoP, double rhouP, double rhovP, double rhoEP, int edgeGrp)
		{
			std::vector<double> MinusVal(4, 0.0);
			double rhoBC(rhoP), rhouBC(0.0), rhovBC(0.0), rhoEBC(rhoP*(material::Cv*bcValues::TBC[edgeGrp - 1]));
			MinusVal[0] = 2 * rhoBC - rhoP;
			MinusVal[1] = 2 * rhouBC - rhouP;
			MinusVal[2] = 2 * rhovBC - rhovP;
			MinusVal[3] = 2 * rhoEBC - rhoEP;
			return MinusVal;
		}
	}

	namespace patch
	{
		std::vector <double> inOutFlow(double rhoP, double rhouP, double rhovP, double rhoEP)
		{
			std::vector<double> MinusVal(4, 0.0);
			double rhoBC(rhoP), rhouBC(rhouP), rhovBC(rhovP), rhoEBC(rhoEP);
			MinusVal[0] = 2 * rhoBC - rhoP;
			MinusVal[1] = 2 * rhouBC - rhouP;
			MinusVal[2] = 2 * rhovBC - rhovP;
			MinusVal[3] = 2 * rhoEBC - rhoEP;
			return MinusVal;
		}
	}
}

namespace strongBCs
{
	std::vector <double> fixedValues(double rhoP, double rhouP, double rhovP, double rhoEP, int edgeGrp)
	{
		std::vector<double> MinusVal(4, 0.0);
		double rhoBC(bcValues::pBC[edgeGrp - 1] / (material::R*bcValues::TBC[edgeGrp - 1]));
		double rhouBC(rhoBC*bcValues::uBC[edgeGrp - 1]),
			rhovBC(rhoBC*bcValues::vBC[edgeGrp - 1]),
			rhoEBC(rhoBC*(bcValues::TBC[edgeGrp - 1] * material::Cv + 0.5*(pow(bcValues::uBC[edgeGrp - 1], 2) + pow(bcValues::vBC[edgeGrp - 1], 2))));
		MinusVal[0] = 2 * rhoBC - rhoP;
		MinusVal[1] = 2 * rhouBC - rhouP;
		MinusVal[2] = 2 * rhovBC - rhovP;
		MinusVal[3] = 2 * rhoEBC - rhoEP;

		return MinusVal;
	}
}