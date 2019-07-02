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

    namespace auxilaryBCs {
        void calcUPlus(int element, int edge, int nG, std::vector<double> &UPlus)
        {
            double a(0.0), b(0.0), TP(0.0);
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            UPlus=math::pointUVars(element, a, b);

            if (flowProperties::massDiffusion)
            {
                double rhoXP(math::pointAuxValue(element,a,b,1,1)), rhoYP(math::pointAuxValue(element,a,b,1,1));
                TP=math::CalcTFromConsvVar_massDiff(UPlus[0],UPlus[1],UPlus[2],UPlus[3],rhoXP,rhoYP);
                UPlus[3]=UPlus[0]*material::Cv*TP+0.5*(UPlus[1]*UPlus[1]+UPlus[2]*UPlus[2])/UPlus[0];
            }
            else {
                TP=math::CalcTFromConsvVar(UPlus[0],UPlus[1],UPlus[2],UPlus[3]);
            }
            surfaceFields::T[edge][nG]=TP;
        }
    }

    namespace NSFEqBCs {
    void calcdUPlus(int element, double a, double b, std::vector<double> &dUXPlus, std::vector<double> &dUYPlus)
    {
        //Compute dU+
        /*
        for (int i = 0; i < 4; i++)
        {
            dUXPlus[i] = math::pointAuxValue(element, a, b, i + 1, 1);
            dUYPlus[i] = math::pointAuxValue(element, a, b, i + 1, 2);
        }
        */
        dUXPlus=math::pointSVars(element,a,b,1);
        dUYPlus=math::pointSVars(element,a,b,2);
    }

    std::tuple<double, double, double, double> calcTotalVelocity(int BCType, double rhoP, double rhoM, double uP, double uM, double vP, double vM, double mudRhoXP, double mudRhoXM, double mudRhoYP, double mudRhoYM)
    {
        double umP(0.0), vmP(0.0), umM(0.0), vmM(0.0);
        umP=math::massDiffusionFncs::calcTotalVelocity(rhoP,uP,mudRhoXP);
        vmP=math::massDiffusionFncs::calcTotalVelocity(rhoP,vP,mudRhoYP);
        if (BCType==1) //type wall
        {
            umM=0;
            vmM=0;
        }
        else {
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
        uMagP = sqrt(pow(uPlus, 2) + pow(vPlus, 2));
        uMagM = sqrt(pow(uMinus, 2) + pow(vMinus, 2));

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
        Fluxes[0][0] = math::numericalFluxes::advectiveFlux(termX1P, termX1M, rhoPlus, rhoMinus, C, nx) + math::numericalFluxes::advectiveFlux(termY1P, termY1M, rhoPlus, rhoMinus, C, ny);
        Fluxes[1][0] = math::numericalFluxes::advectiveFlux(termX2P, termX2M, rhouPlus, rhouMinus, C, nx) + math::numericalFluxes::advectiveFlux(termY2P, termY2M, rhouPlus, rhouMinus, C, ny);
        Fluxes[2][0] = math::numericalFluxes::advectiveFlux(termX3P, termX3M, rhovPlus, rhovMinus, C, nx) + math::numericalFluxes::advectiveFlux(termY3P, termY3M, rhovPlus, rhovMinus, C, ny);
        Fluxes[3][0] = math::numericalFluxes::advectiveFlux(termX4P, termX4M, rhoEPlus, rhoEMinus, C, nx) + math::numericalFluxes::advectiveFlux(termY4P, termY4M, rhoEPlus, rhoEMinus, C, ny);

        /*VISCOUS TERMS*/
        /*calculate viscous terms*/
        StressHeatP = math::viscousTerms::calcStressTensorAndHeatFlux(UPlus, dUXPlus, dUYPlus, TPlus);
        StressHeatM = math::viscousTerms::calcStressTensorAndHeatFlux(UMinus, dUXMinus, dUYMinus, TMinus);
        /*
        if (auxUlti::getBCType(edge) == 0)
        {
            Beta=DiffusiveFluxConst[edge];
        }
        else
        {
            double eP(material::Cv*TPlus), eM(material::Cv*TMinus);
            std::vector<double> nP(2, 0.0);
            nP[0] = nx;
            nP[0] = ny;
            Beta = math::numericalFluxes::constantBeta(uMagP, uMagM, UPlus[0], UMinus[0], eP, eM, pPlus, pMinus, StressHeatP, StressHeatM, nP);
        }
        */
        std::tie(termX1P, termX2P, termX3P, termX4P) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uPlus, vPlus, dUXPlus[0]/rhoPlus, 1);
        std::tie(termY1P, termY2P, termY3P, termY4P) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uPlus, vPlus, dUYPlus[0]/rhoPlus, 2);
        std::tie(termX1M, termX2M, termX3M, termX4M) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uMinus, vMinus, dUXMinus[0]/rhoMinus, 1);
        std::tie(termY1M, termY2M, termY3M, termY4M) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uMinus, vMinus, dUYMinus[0]/rhoMinus, 2);

        /*Calculate fluxes*/
        Fluxes[0][1] = math::numericalFluxes::diffusiveFlux(termX1M, termX1P, rhoPlus, rhoMinus, Beta, nx) + math::numericalFluxes::diffusiveFlux(termY1M, termY1P, rhoPlus, rhoMinus, Beta, ny);
        Fluxes[1][1] = math::numericalFluxes::diffusiveFlux(termX2M, termX2P, rhouPlus, rhouMinus, Beta, nx) + math::numericalFluxes::diffusiveFlux(termY2M, termY2P, rhouPlus, rhouMinus, Beta, ny);
        Fluxes[2][1] = math::numericalFluxes::diffusiveFlux(termX3M, termX3P, rhovPlus, rhovMinus, Beta, nx) + math::numericalFluxes::diffusiveFlux(termY3M, termY3P, rhovPlus, rhovMinus, Beta, ny);
        Fluxes[3][1] = math::numericalFluxes::diffusiveFlux(termX4M, termX4P, rhoEPlus, rhoEMinus, Beta, nx) + math::numericalFluxes::diffusiveFlux(termY4M, termY4P, rhoEPlus, rhoEMinus, Beta, ny);

        return Fluxes;
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
            double a(0.0), b(0.0), TPlus(0.0), muPlus(0.0),TMinus(bcValues::TBC[edgeGrp - 1]), muMinus(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            norm[0] = nx;
            norm[1] = ny;

            //Compute U+
            for (int i = 0; i < 4; i++)
            {
                std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
            }
            TPlus=surfaceFields::T[edge][nG];
            muPlus=math::CalcVisCoef(TPlus);
            muMinus=math::CalcVisCoef(TMinus);

            //Compute dU+
            BCSupportFncs::NSFEqBCs::calcdUPlus(element,a,b,dUXPlus,dUYPlus);
            dUXMinus=dUXPlus;
            dUYMinus=dUYPlus;
            for (int i = 0; i < 4; i++)
            {
                dUXPlus[i]*=muPlus;
                dUYPlus[i]*=muPlus;
                dUXMinus[i]*=muMinus;
                dUYMinus[i]*=muMinus;
            }

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
            double a(0.0), b(0.0), TPlus(0.0), muPlus(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            norm[0] = nx;
            norm[1] = ny;

            //Compute U+
            for (int i = 0; i < 4; i++)
            {
                std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
            }
            TPlus=surfaceFields::T[edge][nG];
            muPlus=math::CalcVisCoef(TPlus);

            //Compute dU+
            BCSupportFncs::NSFEqBCs::calcdUPlus(element,a,b,dUXPlus,dUYPlus);
            for (int i = 0; i < 4; i++)
            {
                dUXPlus[i]*=muPlus;
                dUYPlus[i]*=muPlus;
            }

            //zero normal temperature gradient (sua lai theo cach tong quat)
            dUXMinus = dUXPlus;
            dUYMinus = dUYPlus;
            dUXMinus[3] = 0;
            dUYMinus[3] = 0;

            Fluxes = BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge,1,TPlus, TPlus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
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
                dUXMinus(4, 0.0),
                dUYMinus(4, 0.0),
                norm(2, 0.0);
            double a(0.0), b(0.0), TPlus(0.0), muPlus(0.0), TMinus(bcValues::TBC[edgeGrp - 1]), muMinus(math::CalcVisCoef(bcValues::TBC[edgeGrp - 1])), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            norm[0] = nx;
            norm[1] = ny;

            //Compute U+
            for (int i = 0; i < 4; i++)
            {
                std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
            }
            TPlus=surfaceFields::T[edge][nG];
            muPlus=math::CalcVisCoef(TPlus);

            //Compute dU+
            BCSupportFncs::NSFEqBCs::calcdUPlus(element,a,b,dUXPlus,dUYPlus);
            dUXMinus=dUXPlus;
            dUYMinus=dUYPlus;
            for (int i = 0; i < 4; i++)
            {
                dUXPlus[i]*=muPlus;
                dUYPlus[i]*=muPlus;
                dUXMinus[i]*=muMinus;
                dUYMinus[i]*=muMinus;
            }

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
            double a(0.0), b(0.0), TPlus(0.0), muPlus(0.0), TMinus(0.0), muMinus(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            norm[0] = nx;
            norm[1] = ny;

            //Compute U+
            for (int i = 0; i < 4; i++)
            {
                std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
            }
            TPlus=surfaceFields::T[edge][nG];
            muPlus=math::CalcVisCoef(TPlus);
            if(refValues::subsonic)
            {
                TMinus=bcValues::TBC[edgeGrp-1];
            }
            else {
                TMinus=TPlus;
            }
            muMinus=math::CalcVisCoef(TMinus);

            //Compute dU+
            BCSupportFncs::NSFEqBCs::calcdUPlus(element,a,b,dUXPlus,dUYPlus);
            dUXMinus=dUXPlus;
            dUYMinus=dUYPlus;
            for (int i = 0; i < 4; i++)
            {
                dUXPlus[i]*=muPlus;
                dUYPlus[i]*=muPlus;
                dUXMinus[i]*=muMinus;
                dUYMinus[i]*=muMinus;
            }

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
        double a(0.0), b(0.0), TPlus(0.0), muPlus(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
        std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
        norm[0] = nx;
        norm[1] = ny;

        //Compute U+
        for (int i = 0; i < 4; i++)
        {
            std::tie(UPlus[i],UMinus[i])=auxUlti::getUAtInterfaces(edge,element,nG,i+1);
        }
        TPlus=surfaceFields::T[edge][nG];
        muPlus=math::CalcVisCoef(TPlus);

        //Compute dU+/-
        BCSupportFncs::NSFEqBCs::calcdUPlus(element,a,b,dUXPlus,dUYPlus);
        for (int i = 0; i < 4; i++)
        {
            dUXPlus[i]*=muPlus;
            dUYPlus[i]*=muPlus;
            dUXMinus[i] = dUXPlus[i] - 2 * (dUXPlus[i] * nx + dUYPlus[i] * ny)*nx;
            dUYMinus[i] = dUYPlus[i] - 2 * (dUXPlus[i] * nx + dUYPlus[i] * ny)*ny;
        }

        Fluxes = BCSupportFncs::NSFEqBCs::NSFEqFluxes(edge,3, TPlus, TPlus, UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, norm);
        return Fluxes;
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
            double a(0.0), b(0.0);
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            double nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));

            //Compute plus values--------------------------------------------------------
            BCSupportFncs::auxilaryBCs::calcUPlus(element,edge,nG,UPlus);

            //Compute minus values--------------------------------------------------------
            UMinus[0] = UPlus[0];
            //UMinus[0]=UPlus[0]*surfaceFields::T[edge][nG]/bcValues::TBC[edgeGrp - 1];
            if (bcValues::UBcType[edgeGrp - 1] == 2)
            {
                UMinus[1] = 0;
                UMinus[2] = 0;
                UMinus[3] = UMinus[0]*material::Cv*bcValues::TBC[edgeGrp - 1];
            }
            else if (bcValues::UBcType[edgeGrp - 1] == 3)
            {
                UMinus[1] = bcValues::uBC[edgeGrp - 1]*UMinus[0];
                UMinus[2] = bcValues::vBC[edgeGrp - 1]*UMinus[0];
                UMinus[3] = UMinus[0]*(material::Cv*bcValues::TBC[edgeGrp - 1] + 0.5*(pow(bcValues::uBC[edgeGrp - 1],2)+pow(bcValues::vBC[edgeGrp - 1],2)));
            }
            //-----------------------------------------------------------------------------

            for (int i = 0; i < 4; i++)
            {
                Fluxes[i][0] = math::numericalFluxes::auxFlux(UMinus[i], UPlus[i], nx);
                Fluxes[i][1] = math::numericalFluxes::auxFlux(UMinus[i], UPlus[i], ny);
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

        std::vector <std::vector<double>> wallAdiabatic(int element, int edge, int edgeGrp, int nG)
        {
            //columns 0, 1 are plus, minus values
            std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
            std::vector<double>
                UPlus(4, 0.0),
                UMinus(4, 0.0);
            double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

            //Compute plus values--------------------------------------------------------
            BCSupportFncs::auxilaryBCs::calcUPlus(element,edge,nG,UPlus);

            //Compute minus values--------------------------------------------------------
            UMinus[0] = UPlus[0];
            if (bcValues::UBcType[edgeGrp - 1] == 2)
            {
                UMinus[1] = 0;
                UMinus[2] = 0;
                UMinus[3] = UPlus[3]-0.5*(pow(UPlus[1],2)+pow(UPlus[2],2))/UPlus[0];
            }
            else if (bcValues::UBcType[edgeGrp - 1] == 3)
            {
                UMinus[1] = bcValues::uBC[edgeGrp - 1]*UMinus[0];
                UMinus[2] = bcValues::vBC[edgeGrp - 1]*UMinus[0];
                UMinus[3] = UMinus[0]*((UPlus[3]-0.5*(pow(UPlus[1],2)+pow(UPlus[2],2))/UPlus[0]) + 0.5*(pow(bcValues::uBC[edgeGrp - 1],2)+pow(bcValues::vBC[edgeGrp - 1],2)));
            }
            //-----------------------------------------------------------------------------

            for (int i = 0; i < 4; i++)
            {
                Fluxes[i][0] = math::numericalFluxes::auxFlux(UMinus[i], UPlus[i], nx);
                Fluxes[i][1] = math::numericalFluxes::auxFlux(UMinus[i], UPlus[i], ny);
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
    }

    namespace patch
    {
        std::vector <std::vector<double>> inFlow(int element, int edge, int edgeGrp, int nG)
        {
            std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
            std::vector<double>
                UPlus(4, 0.0),
                UMinus(4, 0.0);
            double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);

            //Compute plus values--------------------------------------------------------
            BCSupportFncs::auxilaryBCs::calcUPlus(element,edge,nG,UPlus);

            //Compute minus values--------------------------------------------------------
            //Apply weak Riemann infinite value
            UMinus[0] = (bcValues::pBC[edgeGrp - 1] / (material::R*bcValues::TBC[edgeGrp - 1]));
            UMinus[1] = UMinus[0] * bcValues::uBC[edgeGrp - 1];
            UMinus[2] = UMinus[0] * bcValues::vBC[edgeGrp - 1];
            UMinus[3] = UMinus[0] * (bcValues::TBC[edgeGrp - 1] * material::Cv + 0.5*(pow(bcValues::uBC[edgeGrp - 1], 2) + pow(bcValues::vBC[edgeGrp - 1], 2)));
            //-----------------------------------------------------------------------------

            for (int i = 0; i < 4; i++)
            {
                Fluxes[i][0] = math::numericalFluxes::auxFlux(UMinus[i], UPlus[i], nx);
                Fluxes[i][1] = math::numericalFluxes::auxFlux(UMinus[i], UPlus[i], ny);
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
            std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
            std::vector<double>
                UPlus(4, 0.0),
                UMinus(4, 0.0),
                norm(2, 0.0);
            double a(0.0), b(0.0), nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
            std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, element, nG);
            norm[0] = nx;
            norm[1] = ny;

            //Compute plus values--------------------------------------------------------
            BCSupportFncs::auxilaryBCs::calcUPlus(element,edge,nG,UPlus);

            //Compute minus values--------------------------------------------------------
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
                double pInternal(0);
                pInternal = UPlus[0] * material::R*math::CalcTFromConsvVar(UPlus[0],UPlus[1],UPlus[2],UPlus[3]);
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
            //---------------------------------------------------------------------------

            for (int i = 0; i < 4; i++)
            {
                Fluxes[i][0] = math::numericalFluxes::auxFlux(UMinus[i], UPlus[i], nx);
                Fluxes[i][1] = math::numericalFluxes::auxFlux(UMinus[i], UPlus[i], ny);
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

        //Compute plus values--------------------------------------------------------
        BCSupportFncs::auxilaryBCs::calcUPlus(element,edge,nG,UPlus);

        //Compute minus values--------------------------------------------------------
		UMinus[0] = UPlus[0];
		UMinus[1] = UPlus[1] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*nx;
		UMinus[2] = UPlus[2] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*ny;
		UMinus[3] = UPlus[3];
        //----------------------------------------------------------------------------

		for (int i = 0; i < 4; i++)
		{
            Fluxes[i][0] = math::numericalFluxes::auxFlux(UMinus[i], UPlus[i], nx);
            Fluxes[i][1] = math::numericalFluxes::auxFlux(UMinus[i], UPlus[i], ny);
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
}
