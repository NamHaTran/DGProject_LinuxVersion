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

//Implement bondary condition of Rho (use when massDiffusion is on)
//Method weakRiemann is used
std::tuple<double, double> rhoBCsImplement(int element, int edge, int nG);

namespace NSFEqBCs
{
    namespace wall
    {
        /*Function computes numerical flux at isothermal wall by using weakRiemann approach*/
        std::vector <std::vector<double>> wallIsoThermal(int element, int edge, int edgeGrp, int nG);

        /*Function computes numerical flux at adiabatic wall by using weakRiemann approach*/
        std::vector <std::vector<double>> wallAdiabatic(int element, int edge, int edgeGrp, int nG);
    }

    namespace patch
    {
        /*Function computes numerical flux at inflow/outflow by using weakRiemann approach*/
        std::vector <std::vector<double>> inFlow(int element, int edge, int edgeGrp, int nG);

        /*Function computes numerical flux at inflow/outflow by using weakRiemann approach*/
        std::vector <std::vector<double>> outFlow(int element, int edge, int edgeGrp, int nG);
    }

    /*Function computes numerical flux at symmetry BC by using weakRiemann approach*/
    std::vector <std::vector<double>> Symmetry(int element, int edge, int nG);
}

namespace auxilaryBCs
{	
    namespace wall
    {
        std::vector <std::vector<double>> wallIsoThermal(int element, int edge, int edgeGrp, int nG);

        std::vector <std::vector<double>> wallAdiabatic(int element, int edge, int edgeGrp, int nG);
    }

    namespace patch {
        std::vector <std::vector<double>> inFlow(int element, int edge, int edgeGrp, int nG);

        std::vector <std::vector<double>> outFlow(int element, int edge, int edgeGrp, int nG);
    }

	/*Function computes numerical flux of auxilary variables at symmetry BC by using weakRiemann approach*/
	std::vector <std::vector<double>> Symmetry(int element, int edge, int nG);
}

namespace BCSupportFncs
{
	/*Function returns true if at considering face (edge), flow is going into computational domain, returns false if flow is going out of computational domain.
	- nx, ny are surface unit normal vector components.*/
	bool checkInflow(double u, double v, double nx, double ny);

    namespace auxilaryBCs {
        void calcUPlus(int element, int edge, int nG, std::vector<double> &UPlus);
    }

    namespace NSFEqBCs {
        void calcdUPlus(int edge, int element, double a, double b, std::vector<double> &dUXPlus, std::vector<double> &dUYPlus);

        /*Function calculates fluxes of all advective and diffusive terms of NSF equation from conservative variables*/
        std::vector<std::vector<double>> NSFEqFluxes(int edge, int BCType, double TPlus, double TMinus, std::vector<double> &UPlus, std::vector<double> &UMinus, std::vector<double> &dUXPlus, std::vector<double> &dUXMinus, std::vector<double> &dUYPlus, std::vector<double> &dUYMinus, std::vector<double> &normVector);

        std::tuple<double, double, double, double> calcTotalVelocity(int BCType, double rhoP, double rhoM, double uP, double uM, double vP, double vM, double mudRhoXP, double mudRhoXM, double mudRhoYP, double mudRhoYM);
    }
}

#endif // DGBCSLIB_H_INCLUDED
