#ifndef DGBCSLIB_H_INCLUDED
#define DGBCSLIB_H_INCLUDED
#include <vector>
/*NOTES:
- All boundary functions must return numerical fluxes (rho, rhou, rhov, rhoE fluxes) at surface Gauss points!
- Method index: 1: weak weak Riemann, 2: weak Prescribed

Boundary conditions compatibility
|U					|T					|p					||advectionBC		|diffusionBC		|auxilaryBC			|
|-------------------+-------------------+-------------------||------------------+-------------------+-------------------+
|1. inOutFlow		|1. inOutFlow		|1. inOutFlow		||inOutFlow			|interiorExtrapolate|inOutFlow			|
|	value u v w		|	value T			|	value p			||					|					|					|
+-------------------+-------------------+-------------------++------------------+-------------------+-------------------+
|2. noSlip			|2. WallIsothermal	|2. zeroGradient	||noSlipIsoThermal	|interiorExtrapolate|noSlipIsoThermal	|
					|	value T			|					||					|					|					|
+-------------------+-------------------+-------------------++------------------+-------------------+-------------------+
|2. noSlip			|3. WallAdiabatic	|2. zeroGradient	||noSlipAdiabatic	|interiorExtrapolate|noSlipAdiabatic	|
+-------------------+-------------------+-------------------++------------------+-------------------+-------------------+
|3.	fixedValue		|4. fixedValue		|3. fixedValue		||fixedValues		|interiorExtrapolate|fixedValues		|
|	value u v w		|	value T			|	value p			||					|					|					|
+-------------------+-------------------+-------------------++------------------+-------------------+-------------------+
*/

/*Function applies auxilary boundary conditions to edges in auxilary equation, it returns values of auxiraly flux at Gauss point*/
std::vector <double> auxEqBCsImplement(int element, int edge, int nG, double n);

/*Function applies advective and diffusive boundary conditions to edges in NSF equation, it returns values of fluxes at Gauss point*/
std::vector<std::vector<double>> NSFEqBCsImplement(int element, int edge, int nG);

namespace advectionBCs
{
	namespace wall
	{
		/*no slip, isothermal wall boundary condition:
		For velocity filed: this BC evaluates the interior velocity field at the boundary V_interior and then calculates external velocity field such that
			V_minus = -V_interior = -V_plus
		For density and internal energy, these values are extrapolated from the interior
			rho_minus = rho_plus
			e_minus = e_plus*/
		std::vector <double> noSlipIsoThermal(int method, int edgeGrp, double rhoP, double rhouP, double rhovP, double rhoEP);

		/*no slip, adiabatic wall boundary condition:
		The advection term in this case is exactly equivalent to the isothermal no-slip wall BCs*/
		std::vector <double> noSlipAdiabatic(int method, int element, int edgeGrp, double rhoP, double rhouP, double rhovP, double rhoEP);
	}

	namespace patch
	{
		/*inflow-outflow boundary condition. This boundary condition uses Riemman invariants to calculates boundary values.
		It also checks flow direction automatically, a surface with negative surface normal velocity is inlet, a surface with positive surface normal velocity is outlet*/
		std::vector <double> inOutFlow(int method, int edgeGrp, double rhoP, double rhouP, double rhovP, double rhoEP, double nx, double ny);

		/*zeroGradient boundary condition*/
		//std::vector <double> zeroGradient(double rhoP, double rhouP, double rhovP, double rhoEP);
		std::vector <double> zeroGradient(int element, double rhoP, double rhouP, double rhovP, double rhoEP);
	}

	//symmetry
}

namespace diffusionBCs
{
	std::vector <double> interiorExtrapolate(double drhoP, double drhouP, double drhovP, double drhoEP);
}

namespace auxilaryBCs
{
	namespace wall
	{
		/*no slip, adiabatic wall boundary condition*/
		std::vector <double> noSlipAdiabatic(double rhoP, double rhouP, double rhovP, double rhoEP);

		/*no slip, isothermal wall boundary condition*/
		std::vector <double> noSlipIsoThermal(double rhoP, double rhouP, double rhovP, double rhoEP, int edgeGrp);
	}

	namespace patch
	{
		/*inflow/outflow boundary condition*/
		std::vector <double> inOutFlow(double rhoP, double rhouP, double rhovP, double rhoEP);
	}
	//symmetry
	//std::vector <double> symmetry(double rhoP, double rhouP, double rhovP, double rhoEP);
}

namespace strongBCs
{
	/*Fixed values boundary condition: All the variables AT THE BOUNDARY are specified from pre-deffined conditions (strong boundary condition)*/
	std::vector <double> fixedValues(double rhoP, double rhouP, double rhovP, double rhoEP, int edgeGrp);
}

namespace BCSupportFncs
{
	/*Function returns minus values at Gauss point, available options:
	- option == 1: advection boundary condition
	- option == 2: diffusion boundary condition
	- option == 3: auxilary boundary condition*/
	std::vector <double> boundaryMinusValsCalculator(int option, int element,int edge, double rhoP, double rhouP, double rhovP, double rhoEP, double nx, double ny);

	/*Function returns true if at considering face (edge), flow is going into computational domain, returns false if flow is going out of computational domain.
	- nx, ny are surface unit normal vector components.*/
	bool checkInflow(double u, double v, double nx, double ny);

	/*Function supports for inflow/outflow boundary condition, it calculates inflow/outflow velocity components at surface*/
	std::tuple<double, double> calcInOutFlowVelocity(double vExternx, double vExterny, double nx, double ny, double vBCMag);
}

#endif // DGBCSLIB_H_INCLUDED
