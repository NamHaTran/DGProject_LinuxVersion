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
	
    /*Vector contents BC edges id and location of them on BC values arrays*/
    //extern std::vector<int>adressOfBCVals;
	extern std::vector<std::vector<int>>neighboringElements;

    extern std::vector<std::vector<double>> geoCenter;
    extern std::vector<double> cellSize, localCellSize, cellArea;

    /*derivatives dx/da, dx/db, dy/da, dy/db*/
    extern std::vector<std::vector<std::vector<double>>>dxa, dxb, dya, dyb;

    /*Jacobian*/
    extern std::vector<std::vector<std::vector<double>>> J2D;
    extern std::vector<double> J1D;

    /*Arrays for parallel computing*/
    extern std::vector<int> Elem2DlocalIdWithRank, rankOf2DElem;
    extern std::vector<std::vector<int>> PointslocalIdWithRank;

    /*Mesh connection*/
    extern std::vector<std::vector<int>> meshConnection;
}

namespace mathVar {
    extern std::vector<double> wGauss, xGauss, wGaussLobatto, xGaussLobatto, B, dBa, dBb;
    extern std::vector<std::vector<std::vector<double>>> GaussPts, wGaussPts, GaussLobattoPts, wGaussLobattoPts,
        BPts_Quad, dBaPts_Quad, dBbPts_Quad, BPts_Tri, dBaPts_Tri, dBbPts_Tri;
}

/*Conservative variables declaration
NOTE: in case of mass diffusion:
- rho, rhou, rhov are computed using advective velocity (u, v)
- rhoE is conputed using mass velocity which is defined as um=u+Jd/rho (Jd is mass diffusive flux Jd = -Dm*grad(rho))
At initial condition i assume that grad(rho) = 0 so (rhoE)initial is conputed using advective velocity.
*/
extern std::vector<std::vector<double>> rho, rhou, rhov, rhoE, rhoN, rhouN, rhovN, rhoEN, rho0, rhou0, rhov0, rhoE0, rhoResArr, rhouResArr, rhovResArr, rhoEResArr;

/*Primary variables declaration
extern std::vector<std::vector<double>>u, v, e, p, T, mu;*/

/*Auxilary variables
NOTE: in case of mass diffusion:
- S=mu*grad(U) with all components of vector U being computed by using advective velocity
(now i define rhoE_adv is total energy with advective kinetic energy)
*/
namespace BR1Vars {
//X direction
extern std::vector<std::vector<double>> rhoX, rhouX, rhovX, rhoEX;

/*Y direction*/
extern std::vector<std::vector<double>> rhoY, rhouY, rhovY, rhoEY;
}

namespace BR2Vars {
//X direction
extern std::vector<std::vector<double>>rhoXVol,rhouXVol,rhovXVol,rhoEXVol;

/*Y direction*/
extern std::vector<std::vector<double>>rhoYVol,rhouYVol,rhovYVol,rhoEYVol;

//X direction
extern std::vector<std::vector<double>>rhoXSurMaster,rhouXSurMaster,rhovXSurMaster,rhoEXSurMaster,rhoXSurSlave,rhouXSurSlave,rhovXSurSlave,rhoEXSurSlave;

/*Y direction*/
extern std::vector<std::vector<double>>rhoYSurMaster,rhouYSurMaster,rhovYSurMaster,rhoEYSurMaster,rhoYSurSlave,rhouYSurSlave,rhovYSurSlave,rhoEYSurSlave;
}

namespace surfaceFields {
    /*Interface values*/
    //conservative variable
    extern std::vector<std::vector<double>> rho, rhou, rhov, rhoE;

    //auxilary equaiton
    extern std::vector<std::vector<double>> aux_rho, aux_rhou, aux_rhov, aux_rhoE;

    //X direction*/
    extern std::vector<std::vector<double>> invis_rhoX, invis_rhouX, invis_rhovX, invis_rhoEX,
    Vis_rhoX, Vis_rhouX, Vis_rhovX, Vis_rhoEX;

    /*Y direction*/
    extern std::vector<std::vector<double>> invis_rhoY, invis_rhouY, invis_rhovY, invis_rhoEY,
    Vis_rhoY, Vis_rhouY, Vis_rhovY, Vis_rhoEY;

    extern std::vector<std::vector<double>> T;
}

namespace volumeFields {
    //Volume values
    extern std::vector<std::vector<std::vector<double>>> rhoVolGauss, rhouVolGauss, rhovVolGauss, rhoEVolGauss, drhoXVolGauss, drhoYVolGauss, T;
}

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

namespace SurfaceBCFields
{
	extern std::vector<std::vector<double>> rhoBc, rhouBc, rhovBc, rhoEBc;
	extern std::vector<std::vector<int>>BCPointsInfor;
}

namespace limitVal {
    extern std::vector<bool> troubleCellsMarker;
}

//for debuging
namespace debug
{
	//extern std::vector<double> minRhoArr, minRhoeArr;
}

namespace parallelBuffer {
    extern std::vector<std::vector<double>> rho,rhou,rhov,rhoE,
    drhoX,drhouX,drhovX,drhoEX,
    drhoY,drhouY,drhovY,drhoEY,
    xCoor, yCoor, aCoor, bCoor;

    extern std::vector<double> theta1, theta2;
    extern std::vector<int> elemType;
}
#endif // DYNAMICVARDECLARATION_H_INCLUDED
