#include "VarDeclaration.h"
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
    std::vector<double> J1D;

    /*Array needed for decomposing case*/
    std::vector<int> Elem2DlocalIdWithRank, rankOf2DElem;
    std::vector<std::vector<int>> PointslocalIdWithRank;
}

namespace mathVar {
    std::vector<double> wGauss, xGauss, wGaussLobatto, xGaussLobatto, B, dBa, dBb;
    std::vector<std::vector<std::vector<double>>> GaussPts, wGaussPts, GaussLobattoPts, wGaussLobattoPts,
        BPts_Quad, dBaPts_Quad, dBbPts_Quad, BPts_Tri, dBaPts_Tri, dBbPts_Tri;
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

/*Auxilary variables*/
namespace BR1Vars {
//X direction
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
}

namespace BR2Vars {
//X direction
std::vector<std::vector<double>>
rhoXVol(1, std::vector<double>(1, 0.0)),
rhouXVol(1, std::vector<double>(1, 0.0)),
rhovXVol(1, std::vector<double>(1, 0.0)),
rhoEXVol(1, std::vector<double>(1, 0.0));

/*Y direction*/
std::vector<std::vector<double>>
rhoYVol(1, std::vector<double>(1, 0.0)),
rhouYVol(1, std::vector<double>(1, 0.0)),
rhovYVol(1, std::vector<double>(1, 0.0)),
rhoEYVol(1, std::vector<double>(1, 0.0));

//X direction
std::vector<std::vector<double>>
rhoXSurMaster(1, std::vector<double>(1, 0.0)),
rhouXSurMaster(1, std::vector<double>(1, 0.0)),
rhovXSurMaster(1, std::vector<double>(1, 0.0)),
rhoEXSurMaster(1, std::vector<double>(1, 0.0)),
rhoXSurSlave(1, std::vector<double>(1, 0.0)),
rhouXSurSlave(1, std::vector<double>(1, 0.0)),
rhovXSurSlave(1, std::vector<double>(1, 0.0)),
rhoEXSurSlave(1, std::vector<double>(1, 0.0));

/*Y direction*/
std::vector<std::vector<double>>
rhoYSurMaster(1, std::vector<double>(1, 0.0)),
rhouYSurMaster(1, std::vector<double>(1, 0.0)),
rhovYSurMaster(1, std::vector<double>(1, 0.0)),
rhoEYSurMaster(1, std::vector<double>(1, 0.0)),
rhoYSurSlave(1, std::vector<double>(1, 0.0)),
rhouYSurSlave(1, std::vector<double>(1, 0.0)),
rhovYSurSlave(1, std::vector<double>(1, 0.0)),
rhoEYSurSlave(1, std::vector<double>(1, 0.0));
}

/*Interface conservative variables*/
namespace surfaceFields {
    std::vector<std::vector<double>>
    rho(1, std::vector<double>(1, 0.0)),
    rhou(1, std::vector<double>(1, 0.0)),
    rhov(1, std::vector<double>(1, 0.0)),
    rhoE(1, std::vector<double>(1, 0.0));

    /*Interface values*/
    //Auxilary equation
    std::vector<std::vector<double>>
    aux_rho(1, std::vector<double>(1, 0.0)),
    aux_rhou(1, std::vector<double>(1, 0.0)),
    aux_rhov(1, std::vector<double>(1, 0.0)),
    aux_rhoE(1, std::vector<double>(1, 0.0));

    //NSF equation
    //X direction
    std::vector<std::vector<double>>
    invis_rhoX(1, std::vector<double>(1, 0.0)),
    invis_rhouX(1, std::vector<double>(1, 0.0)),
    invis_rhovX(1, std::vector<double>(1, 0.0)),
    invis_rhoEX(1, std::vector<double>(1, 0.0));

    std::vector<std::vector<double>>
    Vis_rhoX(1, std::vector<double>(1, 0.0)),
    Vis_rhouX(1, std::vector<double>(1, 0.0)),
    Vis_rhovX(1, std::vector<double>(1, 0.0)),
    Vis_rhoEX(1, std::vector<double>(1, 0.0));

    //Y direction
    std::vector<std::vector<double>>
    invis_rhoY(1, std::vector<double>(1, 0.0)),
    invis_rhouY(1, std::vector<double>(1, 0.0)),
    invis_rhovY(1, std::vector<double>(1, 0.0)),
    invis_rhoEY(1, std::vector<double>(1, 0.0));

    std::vector<std::vector<double>>
    Vis_rhoY(1, std::vector<double>(1, 0.0)),
    Vis_rhouY(1, std::vector<double>(1, 0.0)),
    Vis_rhovY(1, std::vector<double>(1, 0.0)),
    Vis_rhoEY(1, std::vector<double>(1, 0.0));

    std::vector<std::vector<double>> T(1, std::vector<double>(1, 0.0));
}

namespace volumeFields {
    //Volume values
    std::vector<std::vector<std::vector<double>>>
    rhoVolGauss(1, std::vector<std::vector<double>>(1, std::vector<double>(1, 0.0))),
    rhouVolGauss(1, std::vector<std::vector<double>>(1, std::vector<double>(1, 0.0))),
    rhovVolGauss(1, std::vector<std::vector<double>>(1, std::vector<double>(1, 0.0))),
    rhoEVolGauss(1, std::vector<std::vector<double>>(1, std::vector<double>(1, 0.0))),
    drhoXVolGauss(1, std::vector<std::vector<double>>(1, std::vector<double>(1, 0.0))),
    drhoYVolGauss(1, std::vector<std::vector<double>>(1, std::vector<double>(1, 0.0))),
    T(1, std::vector<std::vector<double>>(1, std::vector<double>(1, 0.0)));
}

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
