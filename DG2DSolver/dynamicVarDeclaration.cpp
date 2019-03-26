﻿#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include <vector>

namespace meshVar
{
	std::vector<std::vector<double>> Points, normalVector;
	std::vector<std::vector<int>> Elements1D, Elements2D, BoundaryType; ////BoundaryType: column 0 is boundary group (from 0), column 1 is boundary type (1, 2, 3, 4), column 3 is boundary method
	std::vector<int>  markPointsAtBC, MasterElemOfEdge;

	/*Gauss points on edges*/
	std::vector<std::vector<double>>
		edgeGaussPoints_a(1, std::vector<double>(1, 0.0)),
		edgeGaussPoints_b(1, std::vector<double>(1, 0.0));

	/*Vector contents BC edges name and location of them on BC values arrays*/
	std::vector<int>adressOfBCVals(0, 0);
	std::vector<std::vector<int>>neighboringElements(1, std::vector<int>(1, 0));

    std::vector<std::vector<double>> geoCenter(0, std::vector<double>(2, 0.0));
    std::vector<double> cellSize(0, 0.0), cellArea(0, 0.0), localCellSize(0, 0.0);

    /*derivatives dx/da, dx/db, dy/da, dy/db*/
    std::vector<std::vector<std::vector<double>>>dxa, dxb, dya, dyb;

    /*Jacobian*/
    std::vector<std::vector<std::vector<double>>> J2D;
    std::vector<std::vector<double>> J1D;
}

/*Conservative variables declaration*/
std::vector<std::vector<double>>
rho(1, std::vector<double>(1, 0.0)),
rhou(1, std::vector<double>(1, 0.0)),
rhov(1, std::vector<double>(1, 0.0)),
rhoE(1, std::vector<double>(1, 0.0)),

rhoN(1, std::vector<double>(1, 0.0)),
rhouN(1, std::vector<double>(1, 0.0)),
rhovN(1, std::vector<double>(1, 0.0)),
rhoEN(1, std::vector<double>(1, 0.0)),

rho0(1, std::vector<double>(1, 0.0)),
rhou0(1, std::vector<double>(1, 0.0)),
rhov0(1, std::vector<double>(1, 0.0)),
rhoE0(1, std::vector<double>(1, 0.0)),

rhoResArr(1, std::vector<double>(1, 0.0)),
rhouResArr(1, std::vector<double>(1, 0.0)),
rhovResArr(1, std::vector<double>(1, 0.0)),
rhoEResArr(1, std::vector<double>(1, 0.0));

/*Primary variables declaration*/
/*
std::vector<std::vector<double>>
u(1, std::vector<double>(1, 0.0)),
v(1, std::vector<double>(1, 0.0)),
e(1, std::vector<double>(1, 0.0)),
p(1, std::vector<double>(1, 0.0)),
T(1, std::vector<double>(1, 0.0)),
mu(1, std::vector<double>(1, 0.0));
*/

/*Auxilary variables
//X direction*/
std::vector<std::vector<double>>
rhoX(1, std::vector<double>(1, 0.0)),
rhouX(1, std::vector<double>(1, 0.0)),
rhovX(1, std::vector<double>(1, 0.0)),
rhoEX(1, std::vector<double>(1, 0.0));

/*Y direction*/
std::vector<std::vector<double>>
rhoY(1, std::vector<double>(1, 0.0)),
rhouY(1, std::vector<double>(1, 0.0)),
rhovY(1, std::vector<double>(1, 0.0)),
rhoEY(1, std::vector<double>(1, 0.0));

/*Interface conservative variables*/
std::vector<std::vector<double>>
interface_rho(1, std::vector<double>(1, 0.0)),
interface_rhou(1, std::vector<double>(1, 0.0)),
interface_rhov(1, std::vector<double>(1, 0.0)),
interface_rhoE(1, std::vector<double>(1, 0.0));

/*Interface values*/
//Auxilary equation
std::vector<std::vector<double>>
aux_interface_rho(1, std::vector<double>(1, 0.0)),
aux_interface_rhou(1, std::vector<double>(1, 0.0)),
aux_interface_rhov(1, std::vector<double>(1, 0.0)),
aux_interface_rhoE(1, std::vector<double>(1, 0.0));

//NSF equation
//X direction
std::vector<std::vector<double>>
invis_interface_rhoX(1, std::vector<double>(1, 0.0)),
invis_interface_rhouX(1, std::vector<double>(1, 0.0)),
invis_interface_rhovX(1, std::vector<double>(1, 0.0)),
invis_interface_rhoEX(1, std::vector<double>(1, 0.0));

std::vector<std::vector<double>>
Vis_interface_rhoX(1, std::vector<double>(1, 0.0)),
Vis_interface_rhouX(1, std::vector<double>(1, 0.0)),
Vis_interface_rhovX(1, std::vector<double>(1, 0.0)),
Vis_interface_rhoEX(1, std::vector<double>(1, 0.0));

//Y direction
std::vector<std::vector<double>>
invis_interface_rhoY(1, std::vector<double>(1, 0.0)),
invis_interface_rhouY(1, std::vector<double>(1, 0.0)),
invis_interface_rhovY(1, std::vector<double>(1, 0.0)),
invis_interface_rhoEY(1, std::vector<double>(1, 0.0));

std::vector<std::vector<double>>
Vis_interface_rhoY(1, std::vector<double>(1, 0.0)),
Vis_interface_rhouY(1, std::vector<double>(1, 0.0)),
Vis_interface_rhovY(1, std::vector<double>(1, 0.0)),
Vis_interface_rhoEY(1, std::vector<double>(1, 0.0));

//Lax-Friedrich constant
std::vector<double> LxFConst(1, 0.0);
std::vector<double> DiffusiveFluxConst(1, 0.0);

//StiffMatrixCoefficients
std::vector<std::vector<double>>
stiffMatrixCoeffs(1, std::vector<double>(1, 0.0));

//Limiting coefficients
std::vector<double>
theta1Arr(1, 1.0),
theta2Arr(1, 1.0);

//Volume values
std::vector<std::vector<std::vector<double>>>
rhoVolGauss(1, std::vector<std::vector<double>>(1, std::vector<double>(1, 0.0))),
rhouVolGauss(1, std::vector<std::vector<double>>(1, std::vector<double>(1, 0.0))),
rhovVolGauss(1, std::vector<std::vector<double>>(1, std::vector<double>(1, 0.0))),
rhoEVolGauss(1, std::vector<std::vector<double>>(1, std::vector<double>(1, 0.0)));

namespace SurfaceBCFields
{
	std::vector<std::vector<double>> rhoBc(1, std::vector<double>(1, 0.0)),
		rhouBc(1, std::vector<double>(1, 0.0)),
		rhovBc(1, std::vector<double>(1, 0.0)),
		rhoEBc(1, std::vector<double>(1, 0.0));
	std::vector<std::vector<int>>BCPointsInfor;
}

//for debugging
namespace debug
{
	//std::vector<double>
		//minRhoArr(1, 1.0),
		//minRhoeArr(1, 1.0);
}
