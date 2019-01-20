#ifndef DGBCSLIB_H_INCLUDED
#define DGBCSLIB_H_INCLUDED
#include <vector>
/*NOTES:
- All boundary functions must return numerical fluxes (rho, rhou, rhov, rhoE fluxes) at surface Gauss points!
- Method index: 1: weak weak Riemann, 2: weak Prescribed

Boundary conditions compatibility
		|U					|T					|p					|
		+-------------------+-------------------+-------------------+
		|1. inOutFlow		|1. inOutFlow		|1. inOutFlow		|
		|	Value u v w		|	Value T			|	Value p			|
		+-------------------+-------------------+-------------------+
		|2. noSlip			|2. WallIsothermal	|2. zeroGradient	|
		|					|	Value T			|					|
		+-------------------+-------------------+-------------------+
		|2. noSlip			|3. WallAdiabatic	|2. zeroGradient	|
		+-------------------+-------------------+-------------------+
		|7.	symmetry		|7. symmetry		|7. symmetry		|
		+-------------------+-------------------+-------------------+
*/

/*Function applies advective and diffusive boundary conditions to edges in NSF equation, it returns values of fluxes at Gauss point*/
std::vector<std::vector<double>> NSFEqBCsImplement(int element, int edge, int nG);

/*Function applies auxilary boundary conditions to edges in auxilary equation, it returns values of auxiraly flux at Gauss point*/
std::vector<std::vector<double>> auxEqBCsImplement(int element, int edge, int nG);

namespace NSFEqBCs
{
	namespace weakRiemann
	{
		namespace wall
		{
			/*Function computes numerical flux at isothermal wall by using weakRiemann approach*/
			std::vector <std::vector<double>> noSlipIsoThermal(int element, int edge, int edgeGrp, int nG);

			/*Function computes numerical flux at adiabatic wall by using weakRiemann approach*/
			std::vector <std::vector<double>> noSlipAdiabatic(int element, int edge, int nG);
		}

		namespace patch
		{
			/*Function computes numerical flux at inflow/outflow by using weakRiemann approach*/
			std::vector <std::vector<double>> inOutFlow(int element, int edge, int edgeGrp, int nG);
		}

		/*Function computes numerical flux at symmetry BC by using weakRiemann approach*/
		std::vector <std::vector<double>> Symmetry(int element, int edge, int nG);
	}
	
	namespace weakPrescribed
	{
		namespace wall
		{
			//Function computes numerical flux at isothermal wall by using weakPrescribed approach
			std::vector <std::vector<double>> noSlipIsoThermal(int element, int edge, int nG);

			//Function computes numerical flux at adiabatic wall by using weakPrescribed approach
			std::vector <std::vector<double>> noSlipAdiabatic(int element, int edge, int nG);
		}

		namespace patch
		{
			//Function computes numerical flux at inflow/outflow by using weakPrescribed approach
			std::vector <std::vector<double>> inOutFlow(int element, int edge, int nG);
		}
	}
}

namespace auxilaryBCs
{
	namespace weakPrescribed
	{
		/*Function computes numerical flux of auxilary variables at in/out flow using weakPrescribed approach*/
		std::vector <std::vector<double>> auxFluxesAtBC(int element, int edge, int nG);
	}
	
	namespace weakRiemann
	{
		namespace wall
		{
			std::vector <std::vector<double>> noslipIsoThermal(int element, int edge, int edgeGrp, int nG);

			std::vector <std::vector<double>> noslipAdiabatic(int element, int edge, int nG);
		}

		namespace patch {
			std::vector <std::vector<double>> inOutFlow(int element, int edge, int edgeGrp, int nG);
		}
	}

	/*Function computes numerical flux of auxilary variables at symmetry BC by using weakRiemann approach*/
	std::vector <std::vector<double>> Symmetry(int element, int edge, int nG);
}

namespace BCSupportFncs
{
	/*Function returns true if at considering face (edge), flow is going into computational domain, returns false if flow is going out of computational domain.
	- nx, ny are surface unit normal vector components.*/
	bool checkInflow(double u, double v, double nx, double ny);

	/*Function supports for inflow/outflow boundary condition, it calculates inflow/outflow velocity components at surface*/
	std::tuple<double, double> calcInOutFlowVelocity(double vExternx, double vExterny, double nx, double ny, double vBCMag);

	namespace weakPrescribedFluxes
	{
		/*Function computes numerical fluxes at wall by using weakPrescribed approach*/
		std::vector<std::vector<double>> NSFEqFluxes_Wall(std::vector<double> &UBc, std::vector<double> &dUXBc, std::vector<double> &dUYBc, std::vector<double> &norm);

		/*Function computes BC values at inflow/outflow*/
		void calcInOutFlowBCVals(int element, int edge, int edgeGrp, int nG);

		/*Function computes BC values at isothermal wall*/
		void calcWallIsothermalBCVals(int element, int edge, int edgeGrp, int nG);

		/*Function computes BC values at adiabatic wall*/
		void calcWallAdiabaticBCVals(int element, int edge, int nG);

		/*Function distributes BC values to output Array*/
		std::vector<double> distributeBCValsToArray(int nG, int edge);
	}
}

#endif // DGBCSLIB_H_INCLUDED
