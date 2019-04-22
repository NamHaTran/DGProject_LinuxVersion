#ifndef DYNAMICVARDECLARATION_H_INCLUDED
#define DYNAMICVARDECLARATION_H_INCLUDED
#include "VarDeclaration.h"
#include <vector>

namespace meshVar
{
	extern std::vector<std::vector<double>> Points, normalVector;
	extern std::vector<std::vector<int>> Elements1D,
		Elements2D,
		BoundaryType; //BoundaryType: column 0 is boundary group (from 0), column 1 is boundary type (1, 2, 3, 4), column 3 is boundary method
	extern std::vector<int>  markPointsAtBC, MasterElemOfEdge; //NOTE: boundary point Id (BCPtsId) saved in this vector is (real BCPtsId + 1)
	
	/*Gauss points on edges*/
	extern std::vector<std::vector<double>> edgeGaussPoints_a, edgeGaussPoints_b;
	
	/*Vector contents BC edges name and location of them on BC values arrays*/
	extern std::vector<int>adressOfBCVals;
	extern std::vector<std::vector<int>>neighboringElements;

    extern std::vector<std::vector<double>> geoCenter;
    extern std::vector<double> cellSize, localCellSize, cellArea;

    /*derivatives dx/da, dx/db, dy/da, dy/db*/
    extern std::vector<std::vector<std::vector<double>>>dxa, dxb, dya, dyb;

    /*Jacobian*/
    extern std::vector<std::vector<std::vector<double>>> J2D;
    extern std::vector<double> J1D;
}

namespace mathVar {
    extern std::vector<double> wGauss, xGauss, wGaussLobatto, xGaussLobatto, B, dBa, dBb;
    extern std::vector<std::vector<std::vector<double>>> GaussPts, wGaussPts, GaussLobattoPts, wGaussLobattoPts,
        BPts_Quad, dBaPts_Quad, dBbPts_Quad, BPts_Tri, dBaPts_Tri, dBbPts_Tri;
}

/*Conservative variables declaration*/
extern std::vector<std::vector<double>> rho, rhou, rhov, rhoE, rhoN, rhouN, rhovN, rhoEN, rho0, rhou0, rhov0, rhoE0, rhoResArr, rhouResArr, rhovResArr, rhoEResArr;

/*Primary variables declaration
extern std::vector<std::vector<double>>u, v, e, p, T, mu;*/

/*Auxilary variables
//X direction*/
extern std::vector<std::vector<double>> rhoX, rhouX, rhovX, rhoEX;

/*Y direction*/
extern std::vector<std::vector<double>> rhoY, rhouY, rhovY, rhoEY;

/*Interface values*/
//conservative variable
extern std::vector<std::vector<double>> interface_rho, interface_rhou, interface_rhov, interface_rhoE;

//auxilary equaiton
extern std::vector<std::vector<double>> aux_interface_rho, aux_interface_rhou, aux_interface_rhov, aux_interface_rhoE;

//X direction*/
extern std::vector<std::vector<double>> invis_interface_rhoX, invis_interface_rhouX, invis_interface_rhovX, invis_interface_rhoEX,
Vis_interface_rhoX, Vis_interface_rhouX, Vis_interface_rhovX, Vis_interface_rhoEX;

/*Y direction*/
extern std::vector<std::vector<double>> invis_interface_rhoY, invis_interface_rhouY, invis_interface_rhovY, invis_interface_rhoEY,
Vis_interface_rhoY, Vis_interface_rhouY, Vis_interface_rhovY, Vis_interface_rhoEY;

//Lax-Friedrich constant
extern std::vector<double> LxFConst;

//Diffusion flux constant
extern std::vector<double> DiffusiveFluxConst;

//time step
extern double dt, runTime;

//Limiting coefficients
extern std::vector<double>
theta1Arr,
theta2Arr;

//StiffMatrixCoefficients
extern std::vector<std::vector<double>> stiffMatrixCoeffs;

//Volume values
extern std::vector<std::vector<std::vector<double>>> rhoVolGauss, rhouVolGauss, rhovVolGauss, rhoEVolGauss;

namespace SurfaceBCFields
{
	extern std::vector<std::vector<double>> rhoBc, rhouBc, rhovBc, rhoEBc;
	extern std::vector<std::vector<int>>BCPointsInfor;
}

//for debuging
namespace debug
{
	//extern std::vector<double> minRhoArr, minRhoeArr;
}
#endif // DYNAMICVARDECLARATION_H_INCLUDED
