#include "DGBCsLib.h"
#include <vector>
#include "DGMath.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <tuple>
#include "DGMessagesLib.h"
#include <math.h>

/*NOTES:
- All boundary functions must return numerical fluxes (rho, rhou, rhov, rhoE fluxes) at surface Gauss points!
- Method index: 1: weak weak Riemann, 2: weak Prescribed

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
        U:
        + 3:
        movingWall
        velocity        u v w
*/

std::vector<std::vector<double>> NSFEqBCsImplement(int element, int edge, int nG)
{
	/*Fluxes array has the following form:
	- column 0: advective fluxes
	- column 1: diffusive fluxes*/
	std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
	int edgeGrp(auxUlti::getGrpOfEdge(edge));
	int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]), pType(bcValues::pBcType[edgeGrp - 1]), method(meshVar::BoundaryType[edgeGrp - 1][2]);
	if (UType == 1 && TType == 1 && pType == 1)
	{
		switch (method)
		{
		case 1:  //weak Riemann method
			Fluxes = NSFEqBCs::weakRiemann::patch::inFlow(element, edge, edgeGrp, nG);
			break;
		case 2:  //weak Prescribed
			Fluxes = NSFEqBCs::weakPrescribed::patch::inOutFlow(element, edge, nG);
			break;
		default:
			break;
		}
	}
	else if (UType == 4 && TType == 4 && pType == 4)
	{
		switch (method)
		{
		case 1:  //weak Riemann method
			Fluxes = NSFEqBCs::weakRiemann::patch::outFlow(element, edge, edgeGrp, nG);
			break;
		case 2:  //weak Prescribed
			Fluxes = NSFEqBCs::weakPrescribed::patch::inOutFlow(element, edge, nG);
			break;
		default:
			break;
		}
	}
    else if ((UType == 2 || UType == 3) && TType == 2 && pType == 2)
	{
		switch (method)
		{
		case 1:  //weak Riemann method
            Fluxes = NSFEqBCs::weakRiemann::wall::wallIsoThermal(element, edge, edgeGrp, nG);
			break;
		case 2:  //weak Prescribed
            Fluxes = NSFEqBCs::weakPrescribed::wall::wallIsoThermal(element, edge, nG);
			break;
		default:
			break;
		}
	}
    else if ((UType == 2 || UType == 3) && TType == 3 && pType == 2)
	{
		switch (method)
		{
		case 1:  //weak Riemann method
            Fluxes = NSFEqBCs::weakRiemann::wall::wallAdiabatic(element, edge, edgeGrp, nG);
			break;
		case 2:  //weak Prescribed
            Fluxes = NSFEqBCs::weakPrescribed::wall::wallAdiabatic(element, edge, nG);
			break;
		default:
			break;
		}
	}
	else if (UType == 7 && TType == 7 && pType == 7)
	{
		Fluxes = NSFEqBCs::weakRiemann::Symmetry(element, edge, nG);
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
	std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
	int edgeGrp(auxUlti::getGrpOfEdge(edge));
	int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]), pType(bcValues::pBcType[edgeGrp - 1]), method(meshVar::BoundaryType[edgeGrp - 1][2]);

	if (UType == 1 && TType == 1 && pType == 1)
	{
		switch (method)
		{
		case 1:  //weak Riemann method
			Fluxes = auxilaryBCs::weakRiemann::patch::inFlow(element, edge, edgeGrp, nG);
			break;
		case 2:  //weak Prescribed
			BCSupportFncs::weakPrescribedFluxes::calcInFlowBCVals(element, edge, edgeGrp, nG);
            Fluxes = auxilaryBCs::weakPrescribed::auxFluxesAtBC(element, edge, edgeGrp, nG);
			break;
		default:
			break;
		}
	}
	else if (UType == 4 && TType == 4 && pType == 4)
	{
		switch (method)
		{
		case 1:  //weak Riemann method
			Fluxes = auxilaryBCs::weakRiemann::patch::outFlow(element, edge, edgeGrp, nG);
			break;
		case 2:  //weak Prescribed
			BCSupportFncs::weakPrescribedFluxes::calcOutFlowBCVals(element, edge, edgeGrp, nG);
            Fluxes = auxilaryBCs::weakPrescribed::auxFluxesAtBC(element, edge, edgeGrp, nG);
			break;
		default:
			break;
		}
	}
    else if ((UType == 2 || UType == 3) && TType == 2 && pType == 2)
	{
		switch (method)
		{
		case 1:  //weak Riemann method
            Fluxes = auxilaryBCs::weakRiemann::wall::wallIsoThermal(element, edge, edgeGrp, nG);
			break;
		case 2:  //weak Prescribed
			BCSupportFncs::weakPrescribedFluxes::calcWallIsothermalBCVals(element, edge, edgeGrp, nG);
            Fluxes = auxilaryBCs::weakPrescribed::auxFluxesAtBC(element, edge, edgeGrp, nG);
			break;
		default:
			break;
		}
	}
    else if ((UType == 2 || UType == 3) && TType == 3 && pType == 2)
	{
		switch (method)
		{
		case 1:  //weak Riemann method
            Fluxes = auxilaryBCs::weakRiemann::wall::wallAdiabatic(element, edge, edgeGrp, nG);
			break;
		case 2:  //weak Prescribed
            BCSupportFncs::weakPrescribedFluxes::calcWallAdiabaticBCVals(element, edge, edgeGrp, nG);
            Fluxes = auxilaryBCs::weakPrescribed::auxFluxesAtBC(element, edge, edgeGrp, nG);
			break;
		default:
			break;
		}
	}
	else if (UType == 7 && TType == 7 && pType == 7)
	{
		Fluxes = auxilaryBCs::Symmetry(element, edge, nG);
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
    double rhoP(0.0), rhoM(0.0), rhoFluxX(0.0), rhoFluxY(0.0), a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
    std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
    rhoP = math::pointValue(element, a, b, 1, 2);

    if (UType == 1 && TType == 1 && pType == 1)
    {
        rhoM = (bcValues::pBC[edgeGrp - 1] / (material::R*bcValues::TBC[edgeGrp - 1]));
    }
    else if ((UType == 4 && TType == 4 && pType == 4) || ((UType == 2 || UType == 3) && (TType == 2 || TType == 3) && pType == 2) || (UType == 7 && TType == 7 && pType == 7))
    {
        rhoM = rhoP;
    }
    else
    {
        std::string errorStr = message::BcCompatibleError(edgeGrp);
        message::writeLog(systemVar::pwd, systemVar::caseName, errorStr);
    }
    rhoFluxX = math::numericalFluxes::auxFlux(rhoP, rhoM, nx);
    rhoFluxY = math::numericalFluxes::auxFlux(rhoP, rhoM, ny);
    return std::make_tuple(rhoFluxX,rhoFluxY);
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

	
	namespace weakPrescribedFluxes
	{
		std::vector<std::vector<double>> NSFEqFluxes_Wall(std::vector<double> &UBc, std::vector<double> &dUXBc, std::vector<double> &dUYBc, std::vector<double> &norm)
		{
			std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
			double rhoBc(UBc[0]), rhouBc(UBc[1]), rhovBc(UBc[2]), rhoEBc(UBc[3]);

			double
				uBc(rhouBc / rhoBc),
				vBc(rhovBc / rhoBc),
				totalEBc(rhoEBc / rhoBc),
				TBc(math::CalcTFromConsvVar(rhoBc, rhouBc, rhovBc, rhoEBc)),
				pBc(math::CalcP(TBc, rhoBc));
				//muBc(math::CalcVisCoef(TBc));

			double
				termX1Bc(0.0),  //(rho*u)					or 0
				termX2Bc(0.0),  //(rho*u^2 + p)				or tauxx
				termX3Bc(0.0),  //(rho*u*v)					or tauxy
				termX4Bc(0.0),  //(rho*totalE + p)*u		or tauxx*u + tauxy*v + Qx

				termY1Bc(0.0),  //(rho*v)					or 0
				termY2Bc(0.0),  //(rho*u*v)					or tauxy
				termY3Bc(0.0),  //(rho*v^2 + p)				or tauyy
				termY4Bc(0.0);  //(rho*totalE + p)*v		or tauxy*u + tauyy*v + Qy

			//calculate advective terms
			std::tie(termX1Bc, termX2Bc, termX3Bc, termX4Bc) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoBc, uBc, vBc, totalEBc, pBc, 1);
			std::tie(termY1Bc, termY2Bc, termY3Bc, termY4Bc) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoBc, uBc, vBc, totalEBc, pBc, 2);
			Fluxes[0][0] = termX1Bc * norm[0] + termY1Bc * norm[1];
			Fluxes[1][0] = termX2Bc * norm[0] + termY2Bc * norm[1];
			Fluxes[2][0] = termX3Bc * norm[0] + termY3Bc * norm[1];
			Fluxes[3][0] = termX4Bc * norm[0] + termY4Bc * norm[1];

			//calculate diffusive terms
			std::vector<std::vector<double>> StressHeatBc(2, std::vector<double>(3, 0.0));

			StressHeatBc = math::viscousTerms::calcStressTensorAndHeatFlux(UBc, dUXBc, dUYBc);
			std::tie(termX1Bc, termX2Bc, termX3Bc, termX4Bc) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatBc, uBc, vBc, 1);
			std::tie(termY1Bc, termY2Bc, termY3Bc, termY4Bc) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatBc, uBc, vBc, 2);
			Fluxes[0][1] = termX1Bc * norm[0] + termY1Bc * norm[1];
			Fluxes[1][1] = termX2Bc * norm[0] + termY2Bc * norm[1];
			Fluxes[2][1] = termX3Bc * norm[0] + termY3Bc * norm[1];
			Fluxes[3][1] = termX4Bc * norm[0] + termY4Bc * norm[1];
			return Fluxes;
		}

		void calcInFlowBCVals(int element, int edge, int edgeGrp, int nG)
		{
			std::vector<double>	UPlus(4, 0.0),
				norm(2, 0.0),
				vInternalVector(2, 0.0);
            double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2)),
                RPlus(0.0), RMinus(0.0), rhoBc(0.0), pBc(0.0), TBc(0.0), uBc(0.0), vBc(0.0);;
			std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
			norm[0] = nx;
			norm[1] = ny;

			double vInternalNormMag(0.0), vExternalNormMag(0.0), cInternal(0.0), SBc(0.0);
			double uExternal(bcValues::uBC[edgeGrp - 1]), vExternal(bcValues::vBC[edgeGrp - 1]), cExternal(0.0);
			std::vector<double> vExternalVector(2, 0.0);
			double VBc(0.0), cBc(0.0), rhoExternal(0.0);

			rhoExternal = bcValues::pBC[edgeGrp - 1] / (material::R*bcValues::TBC[edgeGrp - 1]);
			vExternalVector[0] = uExternal;
			vExternalVector[1] = vExternal;
			vExternalNormMag = math::vectorDotProduct(vExternalVector, norm);
			cExternal = math::CalcSpeedOfSound(bcValues::TBC[edgeGrp - 1]);

			for (int i = 0; i < 4; i++)
			{
				UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
			}

			double uPlus(UPlus[1] / UPlus[0]), vPlus(UPlus[2] / UPlus[0]);
			double TInternal(math::CalcTFromConsvVar(UPlus[0], UPlus[1], UPlus[2], UPlus[3]));

			vInternalVector[0] = uPlus;
			vInternalVector[1] = vPlus;
			vInternalNormMag = math::vectorDotProduct(vInternalVector, norm);
			cInternal = math::CalcSpeedOfSound(TInternal);
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
			rhoBc = pow((pow(cBc, 2)) / (material::gamma*SBc), 1 / (material::gamma - 1));
			pBc = rhoBc * pow(cBc, 2) / material::gamma;
			TBc = pBc / (rhoBc*material::R);
			std::tie(uBc, vBc) = BCSupportFncs::calcInOutFlowVelocity(uExternal, vExternal, nx, ny, VBc);

			int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
			SurfaceBCFields::rhoBc[nG][loc] = rhoBc;
			SurfaceBCFields::rhouBc[nG][loc] = rhoBc * uBc;
			SurfaceBCFields::rhovBc[nG][loc] = rhoBc * vBc;
			SurfaceBCFields::rhoEBc[nG][loc] = rhoBc * (material::Cv*TBc + 0.5*(pow(uBc, 2) + pow(vBc, 2)));
		}

		void calcOutFlowBCVals(int element, int edge, int edgeGrp, int nG)
		{
			std::vector<double>	UPlus(4, 0.0),
				norm(2, 0.0),
				vInternalVector(2, 0.0);
            double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2)),
                RPlus(0.0), RMinus(0.0), rhoBc(0.0), pBc(0.0), TBc(0.0), uBc(0.0), vBc(0.0);;
			std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
			norm[0] = nx;
			norm[1] = ny;

			double vInternalNormMag(0.0), vExternalNormMag(0.0), cInternal(0.0), SBc(0.0);
			double uExternal(bcValues::uBC[edgeGrp - 1]), vExternal(bcValues::vBC[edgeGrp - 1]), cExternal(0.0);
			std::vector<double> vExternalVector(2, 0.0);
			double VBc(0.0), cBc(0.0), rhoExternal(0.0);

			rhoExternal = bcValues::pBC[edgeGrp - 1] / (material::R*bcValues::TBC[edgeGrp - 1]);
			vExternalVector[0] = uExternal;
			vExternalVector[1] = vExternal;
			vExternalNormMag = math::vectorDotProduct(vExternalVector, norm);
			cExternal = math::CalcSpeedOfSound(bcValues::TBC[edgeGrp - 1]);

			for (int i = 0; i < 4; i++)
			{
				UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
			}

			double uPlus(UPlus[1] / UPlus[0]), vPlus(UPlus[2] / UPlus[0]);
			double TInternal(math::CalcTFromConsvVar(UPlus[0], UPlus[1], UPlus[2], UPlus[3]));

			vInternalVector[0] = uPlus;
			vInternalVector[1] = vPlus;
			vInternalNormMag = math::vectorDotProduct(vInternalVector, norm);
			cInternal = math::CalcSpeedOfSound(TInternal);

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
			SBc = pow(cInternal, 2) / (material::gamma*pow(UPlus[0], material::gamma - 1));
			rhoBc = pow((pow(cBc, 2)) / (material::gamma*SBc), 1 / (material::gamma - 1));
			pBc = rhoBc * pow(cBc, 2) / material::gamma;
			TBc = pBc / (rhoBc*material::R);
			std::tie(uBc, vBc) = BCSupportFncs::calcInOutFlowVelocity(uPlus, vPlus, nx, ny, VBc);

			int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
			SurfaceBCFields::rhoBc[nG][loc] = rhoBc;
			SurfaceBCFields::rhouBc[nG][loc] = rhoBc * uBc;
			SurfaceBCFields::rhovBc[nG][loc] = rhoBc * vBc;
			SurfaceBCFields::rhoEBc[nG][loc] = rhoBc * (material::Cv*TBc + 0.5*(pow(uBc, 2) + pow(vBc, 2)));
		}

		void calcWallIsothermalBCVals(int element, int edge, int edgeGrp, int nG)
		{
			double a(0.0), b(0.0);
			std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
			int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
			SurfaceBCFields::rhoBc[nG][loc] = math::pointValue(element, a, b, 1, 2);
            if (bcValues::UBcType[edgeGrp - 1] == 2)
            {
                SurfaceBCFields::rhouBc[nG][loc] = 0;
                SurfaceBCFields::rhovBc[nG][loc] = 0;
            }
            else if (bcValues::UBcType[edgeGrp - 1] == 3)
            {
                SurfaceBCFields::rhouBc[nG][loc] = bcValues::uBC[edgeGrp - 1]*SurfaceBCFields::rhoBc[nG][loc];
                SurfaceBCFields::rhovBc[nG][loc] = bcValues::vBC[edgeGrp - 1]*SurfaceBCFields::rhoBc[nG][loc];
            }
			SurfaceBCFields::rhoEBc[nG][loc] = SurfaceBCFields::rhoBc[nG][loc] * (material::R*bcValues::TBC[edgeGrp - 1]) / (material::gamma - 1);
		}

        void calcWallAdiabaticBCVals(int element, int edge, int edgeGrp, int nG)
		{
			std::vector<double> UBc(4, 0.0),
				UPlus(4, 0.0);

			double a(0.0), b(0.0);
			std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

			for (int i = 0; i < 4; i++)
			{
				UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
			}
			double TP(math::CalcTFromConsvVar(UPlus[0], UPlus[1], UPlus[2], UPlus[3]));
			int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
			SurfaceBCFields::rhoBc[nG][loc] = math::pointValue(element, a, b, 1, 2);
            if (bcValues::UBcType[edgeGrp - 1] == 2)
            {
                SurfaceBCFields::rhouBc[nG][loc] = 0;
                SurfaceBCFields::rhovBc[nG][loc] = 0;
            }
            else if (bcValues::UBcType[edgeGrp - 1] == 3)
            {
                SurfaceBCFields::rhouBc[nG][loc] = bcValues::uBC[edgeGrp - 1]*SurfaceBCFields::rhoBc[nG][loc];
                SurfaceBCFields::rhovBc[nG][loc] = bcValues::vBC[edgeGrp - 1]*SurfaceBCFields::rhoBc[nG][loc];
            }
			SurfaceBCFields::rhoEBc[nG][loc] = SurfaceBCFields::rhoBc[nG][loc] * material::R*TP / (material::gamma - 1);
		}

		std::vector<double> distributeBCValsToArray(int nG, int edge)
		{
			std::vector<double> Arr(4, 0.0);
			int loc(auxUlti::getAdressOfBCEdgesOnBCValsArray(edge));
			Arr[0] = SurfaceBCFields::rhoBc[nG][loc];
			Arr[1] = SurfaceBCFields::rhouBc[nG][loc];
			Arr[2] = SurfaceBCFields::rhovBc[nG][loc];
			Arr[3] = SurfaceBCFields::rhoEBc[nG][loc];
			return Arr;
		}
	}
}

namespace NSFEqBCs
{
	namespace weakRiemann
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
					norm(2, 0.0);
                double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
				std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
				norm[0] = nx;
				norm[1] = ny;

				//A1 approach
				for (int i = 0; i < 4; i++)
				{
					UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
					dUXPlus[i] = math::pointAuxValue(element, a, b, i + 1, 1);
					dUYPlus[i] = math::pointAuxValue(element, a, b, i + 1, 2);
				}

				UMinus[0] = UPlus[0];
                if (bcValues::UBcType[edgeGrp - 1] == 2)
                {
                    UMinus[1] = 0.0;
                    UMinus[2] = 0.0;
                }
                else if (bcValues::UBcType[edgeGrp - 1] == 3) {
                    UMinus[1] = UMinus[0]*bcValues::uBC[edgeGrp - 1]*UMinus[0];
                    UMinus[2] = UMinus[0]*bcValues::vBC[edgeGrp - 1]*UMinus[0];
                }
				UMinus[3] = UPlus[0]*material::Cv*bcValues::TBC[edgeGrp - 1];
				//with isothermal BC, dUXMinus = dUXPlus, dUYMinus = dUYPlus
				Fluxes = math::numericalFluxes::NSFEqAdvDiffFluxFromConserVars(edge, UPlus, UMinus, dUXPlus, dUXPlus, dUYPlus, dUYPlus, norm);
				return Fluxes;
			}

            std::vector <std::vector<double>> wallAdiabatic(int element, int edge, int edgeGrp, int nG)
			{
				std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
				std::vector<double> UMinus(4, 0.0),
					UPlus(4, 0.0),
					dUXPlus(4, 0.0), dUXMinus(4, 0.0),
					dUYPlus(4, 0.0), dUYMinus(4, 0.0),
					norm(2, 0.0);
                double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
				std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
				norm[0] = nx;
				norm[1] = ny;

				//A1 approach
				for (int i = 0; i < 4; i++)
				{
					UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
					dUXPlus[i] = math::pointAuxValue(element, a, b, i + 1, 1);
					dUYPlus[i] = math::pointAuxValue(element, a, b, i + 1, 2);
					dUXMinus[i] = dUXPlus[i];
					dUYMinus[i] = dUYPlus[i];
				}
				//zero normal temperature gradient (sua lai theo cach tong quat)
				dUXMinus[3] = 0;
				dUYMinus[3] = 0;

				UMinus[0] = UPlus[0];
                if (bcValues::UBcType[edgeGrp - 1] == 2)
                {
                    UMinus[1] = 0.0;
                    UMinus[2] = 0.0;
                }
                else if (bcValues::UBcType[edgeGrp - 1] == 3) {
                    UMinus[1] = UMinus[0]*bcValues::uBC[edgeGrp - 1]*UMinus[0];
                    UMinus[2] = UMinus[0]*bcValues::vBC[edgeGrp - 1]*UMinus[0];
                }
				UMinus[3] = UPlus[0]*material::Cv*math::CalcTFromConsvVar(UPlus[0], UPlus[1], UPlus[2], UPlus[3]);

				Fluxes = math::numericalFluxes::NSFEqAdvDiffFluxFromConserVars(edge, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
				return Fluxes;
			}
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
					norm(2, 0.0);
                double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
				std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
				norm[0] = nx;
				norm[1] = ny;

				for (int i = 0; i < 4; i++)
				{
					UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
					dUXPlus[i] = math::pointAuxValue(element, a, b, i + 1, 1);
					dUYPlus[i] = math::pointAuxValue(element, a, b, i + 1, 2);
				}

				//Apply weak Riemann infinite value
				UMinus[0] = (bcValues::pBC[edgeGrp - 1] / (material::R*bcValues::TBC[edgeGrp - 1]));
				UMinus[1] = UMinus[0] * bcValues::uBC[edgeGrp - 1];
				UMinus[2] = UMinus[0] * bcValues::vBC[edgeGrp - 1];
				UMinus[3] = UMinus[0] * (bcValues::TBC[edgeGrp - 1] * material::Cv + 0.5*(pow(bcValues::uBC[edgeGrp - 1], 2) + pow(bcValues::vBC[edgeGrp - 1], 2)));

				Fluxes = math::numericalFluxes::NSFEqAdvDiffFluxFromConserVars(edge, UPlus, UMinus, dUXPlus, dUXPlus, dUYPlus, dUYPlus, norm);
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
					norm(2, 0.0);
                double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
				std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
				norm[0] = nx;
				norm[1] = ny;

				for (int i = 0; i < 4; i++)
				{
					UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
					dUXPlus[i] = math::pointAuxValue(element, a, b, i + 1, 1);
					dUYPlus[i] = math::pointAuxValue(element, a, b, i + 1, 2);
				}

				double uPlus(UPlus[1] / UPlus[0]), vPlus(UPlus[2] / UPlus[0]);

				//Apply PNR (2), R (1)
				int implementation(1);
				switch (implementation)
				{
				case 1: //R
				{
					if (refValues::subsonic)
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
						UMinus[3] = bcValues::pBC[edgeGrp - 1] / (material::gamma - 1) + 0.5*UPlus[0] * (pow(uPlus, 2) + pow(vPlus, 2));
					}
					else
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
						UMinus[3] = UPlus[3];
					}
				}
				break;
				case 2: //PNR
				{
					double TInternal(math::CalcTFromConsvVar(UPlus[0], UPlus[1], UPlus[2], UPlus[3]));
					double pInternal(0);
					pInternal = UPlus[0] * material::R*TInternal;
					if (refValues::subsonic)
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
						UMinus[3] = (2 * bcValues::pBC[edgeGrp - 1] - pInternal) / (material::gamma - 1) + 0.5*UPlus[0] * (pow(uPlus, 2) + pow(vPlus, 2));
					}
					else
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
						UMinus[3] = UPlus[3];
					}
				}
				break;
				default:
					break;
				}

				Fluxes = math::numericalFluxes::NSFEqAdvDiffFluxFromConserVars(edge, UPlus, UMinus, dUXPlus, dUXPlus, dUYPlus, dUYPlus, norm);
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
			double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
			std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
			norm[0] = nx;
			norm[1] = ny;

			for (int i = 0; i < 4; i++)
			{
				UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
				//UMinus[i] = UPlus[i];
				dUXPlus[i] = math::pointAuxValue(element, a, b, i + 1, 1);
				dUYPlus[i] = math::pointAuxValue(element, a, b, i + 1, 2);
				dUXMinus[i] = dUXPlus[i] - 2 * (dUXPlus[i] * nx + dUYPlus[i] * ny)*nx;
				dUYMinus[i] = dUYPlus[i] - 2 * (dUXPlus[i] * nx + dUYPlus[i] * ny)*ny;
				//dUXMinus[i] = dUXPlus[i];
				//dUYMinus[i] = dUYPlus[i];
			}
			UMinus[0] = UPlus[0];
			UMinus[1] = UPlus[1] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*nx;
			UMinus[2] = UPlus[2] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*ny;
			UMinus[3] = UPlus[3];
			
			Fluxes = math::numericalFluxes::NSFEqAdvDiffFluxFromConserVars(edge, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
			return Fluxes;
		}
	}

	namespace weakPrescribed
	{
		namespace wall
		{
            std::vector <std::vector<double>> wallIsoThermal(int element, int edge, int nG)
			{
				std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
				std::vector<double> UBc(4, 0.0),
					dUXBc(4, 0.0),
					dUYBc(4, 0.0),
					norm(2, 0.0);
                double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
				std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
				norm[0] = nx;
				norm[1] = ny;
				UBc = BCSupportFncs::weakPrescribedFluxes::distributeBCValsToArray(nG, edge);

				//with isothermal BC, dUXMinus = dUXPlus, dUYMinus = dUYPlus
				for (int i = 0; i < 4; i++)
				{
					dUXBc[i] = math::pointAuxValue(element, a, b, i + 1, 1);
					dUYBc[i] = math::pointAuxValue(element, a, b, i + 1, 2);
				}

				Fluxes = BCSupportFncs::weakPrescribedFluxes::NSFEqFluxes_Wall(UBc, dUXBc, dUYBc, norm);
				return Fluxes;
			}

            std::vector <std::vector<double>> wallAdiabatic(int element, int edge, int nG)
			{
				std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
				std::vector<double> UBc(4, 0.0),
					UPlus(4, 0.0),
					dUXBc(4, 0.0),
					dUYBc(4, 0.0),
					norm(2, 0.0);
                double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
				std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
				norm[0] = nx;
				norm[1] = ny;
				for (int i = 0; i < 4; i++)
				{
					dUXBc[i] = math::pointAuxValue(element, a, b, i + 1, 1);
					dUYBc[i] = math::pointAuxValue(element, a, b, i + 1, 2);
				}
				UBc = BCSupportFncs::weakPrescribedFluxes::distributeBCValsToArray(nG, edge);

				//zero normal temperature gradient
				dUXBc[3] = 0;
				dUYBc[3] = 0;

				Fluxes = BCSupportFncs::weakPrescribedFluxes::NSFEqFluxes_Wall(UBc, dUXBc, dUYBc, norm);
				return Fluxes;
			}
		}

		namespace patch
		{
			std::vector <std::vector<double>> inOutFlow(int element, int edge, int nG)
			{
				std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
				std::vector<double>
					UPlus(4, 0.0),
					UBc(4, 0.0),
					dUXPlus(4, 0.0),
					dUYPlus(4, 0.0),
					norm(2, 0.0);
				double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
				std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
				norm[0] = nx;
				norm[1] = ny;

				for (int i = 0; i < 4; i++)
				{
					UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
					dUXPlus[i] = math::pointAuxValue(element, a, b, i + 1, 1);
					dUYPlus[i] = math::pointAuxValue(element, a, b, i + 1, 2);
				}
				UBc = BCSupportFncs::weakPrescribedFluxes::distributeBCValsToArray(nG, edge);

				Fluxes = math::numericalFluxes::NSFEqAdvDiffFluxFromConserVars(edge, UPlus, UBc, dUXPlus, dUXPlus, dUYPlus, dUYPlus, norm);
				return Fluxes;
			}
		}
	}
}

namespace auxilaryBCs
{
	namespace weakPrescribed
	{
        std::vector <std::vector<double>> auxFluxesAtBC(int element, int edge,  int edgeGrp, int nG)
		{
			//General formulae: hS_BC = UBc*n
			std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
			std::vector<double> UBc(4, 0.0);
			double TBc(0.0), muBc(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            UBc = BCSupportFncs::weakPrescribedFluxes::distributeBCValsToArray(nG, edge);
            if (bcValues::TBcType[edgeGrp - 1] == 2)
            {
                TBc = bcValues::TBC[edgeGrp - 1];
            }
            else {
                TBc = math::CalcTFromConsvVar(UBc[0], UBc[1], UBc[2], UBc[3]);
            }
			muBc = math::CalcVisCoef(TBc);
			for (int i = 0; i < 4; i++)
			{
				Fluxes[i][0] = UBc[i] * nx*muBc;
				Fluxes[i][1] = UBc[i] * ny*muBc;
			}
			return Fluxes;
		}
	}

	namespace weakRiemann
	{
		namespace wall
		{
            std::vector <std::vector<double>> wallIsoThermal(int element, int edge, int edgeGrp, int nG)
			{
				//columns 0, 1 are plus, minus values
				std::vector<std::vector<double>> gaussVector(4, std::vector<double>(2, 0.0)), Fluxes(4, std::vector<double>(2, 0.0));
				double a(0.0), b(0.0);
				std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
				double muM(math::CalcVisCoef(bcValues::TBC[edgeGrp - 1])), muP(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
				for (int i = 0; i < 4; i++)
				{
					gaussVector[i][0] = math::pointValue(element, a, b, i + 1, 2);
				}
				muP = math::CalcVisCoef(math::CalcTFromConsvVar(gaussVector[0][0], gaussVector[1][0], gaussVector[2][0], gaussVector[3][0]));

				gaussVector[0][1] = gaussVector[0][0];
                if (bcValues::UBcType[edgeGrp - 1] == 2)
                {
                    gaussVector[1][1] = 0;
                    gaussVector[2][1] = 0;
                }
                else if (bcValues::UBcType[edgeGrp - 1] == 3)
                {
                    gaussVector[1][1] = bcValues::uBC[edgeGrp - 1]*gaussVector[0][1];
                    gaussVector[2][1] = bcValues::vBC[edgeGrp - 1]*gaussVector[0][1];
                }
				gaussVector[3][1] = gaussVector[0][1]*material::Cv*bcValues::TBC[edgeGrp - 1];

				for (int i = 0; i < 4; i++)
				{
					gaussVector[i][0] = gaussVector[i][0] * muP;
					gaussVector[i][1] = gaussVector[i][1] * muM;
					Fluxes[i][0] = math::numericalFluxes::auxFlux(gaussVector[i][1], gaussVector[i][0], nx);
					Fluxes[i][1] = math::numericalFluxes::auxFlux(gaussVector[i][1], gaussVector[i][0], ny);
				}
				return Fluxes;
			}

            std::vector <std::vector<double>> wallAdiabatic(int element, int edge, int edgeGrp, int nG)
			{
				//columns 0, 1 are plus, minus values
				std::vector<std::vector<double>> gaussVector(4, std::vector<double>(2, 0.0)), Fluxes(4, std::vector<double>(2, 0.0));
				double a(0.0), b(0.0), muP(0.0), muM(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
				std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

				for (int i = 0; i < 4; i++)
				{
					gaussVector[i][0] = math::pointValue(element, a, b, i + 1, 2);
				}
				muP = math::CalcVisCoef(math::CalcTFromConsvVar(gaussVector[0][0], gaussVector[1][0], gaussVector[2][0], gaussVector[3][0]));
				muM = muP;
				gaussVector[0][1] = gaussVector[0][0];
                if (bcValues::UBcType[edgeGrp - 1] == 2)
                {
                    gaussVector[1][1] = 0;
                    gaussVector[2][1] = 0;
                }
                else if (bcValues::UBcType[edgeGrp - 1] == 3)
                {
                    gaussVector[1][1] = bcValues::uBC[edgeGrp - 1]*gaussVector[0][1];
                    gaussVector[2][1] = bcValues::vBC[edgeGrp - 1]*gaussVector[0][1];
                }
				gaussVector[3][1] = gaussVector[0][1] * material::Cv*math::CalcTFromConsvVar(gaussVector[0][0], gaussVector[1][0], gaussVector[2][0], gaussVector[3][0]);

				for (int i = 0; i < 4; i++)
				{
					gaussVector[i][0] = gaussVector[i][0] * muP;
					gaussVector[i][1] = gaussVector[i][1] * muM;
					Fluxes[i][0] = math::numericalFluxes::auxFlux(gaussVector[i][1], gaussVector[i][0], nx);
					Fluxes[i][1] = math::numericalFluxes::auxFlux(gaussVector[i][1], gaussVector[i][0], ny);
				}
				return Fluxes;
			}
		}

		namespace patch
		{
			std::vector <std::vector<double>> inFlow(int element, int edge, int edgeGrp, int nG)
			{
				std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
				std::vector<double>
					UPlus(4, 0.0),
                    UMinus(4, 0.0);
				double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2)), rhoEBC(0), muP(0.0), muM(0.0);
				std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

				for (int i = 0; i < 4; i++)
				{
					UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
				}

				//Apply weak Riemann infinite value
				UMinus[0] = (bcValues::pBC[edgeGrp - 1] / (material::R*bcValues::TBC[edgeGrp - 1]));
				UMinus[1] = UMinus[0] * bcValues::uBC[edgeGrp - 1];
				UMinus[2] = UMinus[0] * bcValues::vBC[edgeGrp - 1];
				UMinus[3] = UMinus[0] * (bcValues::TBC[edgeGrp - 1] * material::Cv + 0.5*(pow(bcValues::uBC[edgeGrp - 1], 2) + pow(bcValues::vBC[edgeGrp - 1], 2)));

				muP = (math::CalcVisCoef(math::CalcTFromConsvVar(UPlus[0], UPlus[1], UPlus[2], UPlus[3])));
				muM = (math::CalcVisCoef(math::CalcTFromConsvVar(UMinus[0], UMinus[1], UMinus[2], UMinus[3])));

				for (int i = 0; i < 4; i++)
				{
					Fluxes[i][0] = math::numericalFluxes::auxFlux(UMinus[i] * muM, UPlus[i] * muP, nx);
					Fluxes[i][1] = math::numericalFluxes::auxFlux(UMinus[i] * muM, UPlus[i] * muP, ny);
				}
				return Fluxes;
			}

			std::vector <std::vector<double>> outFlow(int element, int edge, int edgeGrp, int nG)
			{
				std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
				std::vector<double>
					UPlus(4, 0.0),
					UMinus(4, 0.0),
					norm(2, 0.0);
				double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2)), rhoEBC(0), muP(0.0), muM(0.0);
				std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
				norm[0] = nx;
				norm[1] = ny;

				for (int i = 0; i < 4; i++)
				{
					UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
				}

				double uPlus(UPlus[1] / UPlus[0]), vPlus(UPlus[2] / UPlus[0]);

				//Apply PNR (2), R (1)
				int implementation(1);
				switch (implementation)
				{
				case 1: //R
				{
					if (refValues::subsonic)
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
						UMinus[3] = bcValues::pBC[edgeGrp - 1] / (material::gamma - 1) + 0.5*UPlus[0] * (pow(uPlus, 2) + pow(vPlus, 2));
					}
					else
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
						UMinus[3] = UPlus[3];

					}
				}
				break;
				case 2: //PNR
				{
					double TInternal(math::CalcTFromConsvVar(UPlus[0], UPlus[1], UPlus[2], UPlus[3]));
					double pInternal(0);
					pInternal = UPlus[0] * material::R*TInternal;
					if (refValues::subsonic)
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
						UMinus[3] = (2 * bcValues::pBC[edgeGrp - 1] - pInternal) / (material::gamma - 1) + 0.5*UPlus[0] * (pow(uPlus, 2) + pow(vPlus, 2));
					}
					else
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
						UMinus[3] = UPlus[3];
					}
				}
				break;
				default:
					break;
				}
				muP = (math::CalcVisCoef(math::CalcTFromConsvVar(UPlus[0], UPlus[1], UPlus[2], UPlus[3])));
				muM = (math::CalcVisCoef(math::CalcTFromConsvVar(UMinus[0], UMinus[1], UMinus[2], UMinus[3])));

				for (int i = 0; i < 4; i++)
				{
					Fluxes[i][0] = math::numericalFluxes::auxFlux(UMinus[i] * muM, UPlus[i] * muP, nx);
					Fluxes[i][1] = math::numericalFluxes::auxFlux(UMinus[i] * muM, UPlus[i] * muP, ny);
				}
				return Fluxes;
			}
		}
	}

	std::vector <std::vector<double>> Symmetry(int element, int edge, int nG)
	{
		std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
		std::vector<double>
			UPlus(4, 0.0),
			UMinus(4, 0.0),
			norm(2, 0.0);
		double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
		std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
		norm[0] = nx;
		norm[1] = ny;

		for (int i = 0; i < 4; i++)
		{
			UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
			UMinus[i] = UPlus[i];
		}
		double muP(math::CalcVisCoef(math::CalcTFromConsvVar(UPlus[0], UPlus[1], UPlus[2], UPlus[3])));

		UMinus[0] = UPlus[0];
		UMinus[1] = UPlus[1] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*nx;
		UMinus[2] = UPlus[2] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*ny;
		UMinus[3] = UPlus[3];

		UMinus[0] = UPlus[0];
		for (int i = 0; i < 4; i++)
		{
			Fluxes[i][0] = math::numericalFluxes::auxFlux(UMinus[i] * muP, UPlus[i] * muP, nx);
			Fluxes[i][1] = math::numericalFluxes::auxFlux(UMinus[i] * muP, UPlus[i] * muP, ny);
		}
		return Fluxes;
	}
}
