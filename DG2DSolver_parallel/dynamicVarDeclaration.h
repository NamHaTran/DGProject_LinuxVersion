#ifndef DYNAMICVARDECLARATION_H_INCLUDED
#define DYNAMICVARDECLARATION_H_INCLUDED
#include "VarDeclaration.h"
#include <vector>

namespace meshVar
{
    /*Edges informations*/
    extern int **inpoed;
    /*Cot 5 chua edgeId cua edge tai bien ung voi Element1D array*/

    /*Elements surrounding element*/
    extern int **esuel;

    extern int **inpoel;

    /*Edges of element*/
    //number of row is 4 because of default quad element
    extern int **inedel;
    extern int **ineled;

    extern double **Points;
    extern double **normalVector;

    extern int **Elements1D;
    extern int **Elements2D;
    extern int **BoundaryType; //BoundaryType: column 0 is boundary group (from 0), column 1 is boundary type (1, 2, 3, 4), column 3 is boundary method
    extern int *markPointsAtBC;
    extern int *MasterElemOfEdge;//NOTE: boundary point Id (BCPtsId) saved in this vector is (real BCPtsId + 1)
	
	/*Gauss points on edges*/
    extern double **edgeGaussPoints_a;
    extern double **edgeGaussPoints_b;
	
    /*Vector contents BC edges id and location of them on BC values arrays*/
    //extern std::vector<int>adressOfBCVals;
    extern int **neighboringElements;

    extern double **geoCenter;
    extern double *cellSize;
    extern double *cellArea;
    extern double *localCellSize;

    /*derivatives dx/da, dx/db, dy/da, dy/db*/
    extern double **dxa;
    extern double **dxb;
    extern double **dya;
    extern double **dyb;
	
    /*Jacobian*/
    extern double **J2D;
    extern double *J1D;

    /*Arrays for parallel computing*/
    extern int *Elem2DlocalIdWithRank;
    extern int *rankOf2DElem;
    extern int **PointslocalIdWithRank;

    /*Mesh connection*/
    extern int **meshConnection;

    /*For Maxwell-Smoluchowsky BC*/
    extern double **normProjectionOfCenterToBCEdge_realSysCoor;
    extern double **normProjectionOfCenterToBCEdge_standardSysCoor;
    extern double *distanceFromCentroidToBCEdge;
}

namespace mathVar {
    extern double *wGauss;
    extern double *xGauss;
    extern double *wGaussLobatto;
    extern double *xGaussLobatto;
    extern double *B;
    extern double *dBa;
    extern double *dBb;
	
    extern double **GaussPts;
    extern double **wGaussPts;
    extern double **GaussLobattoPts;
    extern double **wGaussLobattoPts;
    extern double **BPts_Quad;
    extern double **dBaPts_Quad;
    extern double **dBbPts_Quad;
    extern double **BPts_Tri;
    extern double **dBaPts_Tri;
    extern double **dBbPts_Tri;
}

/*Conservative variables declaration
NOTE: in case of mass diffusion:
- rho, rhou, rhov are computed using advective velocity (u, v)
- rhoE is conputed using mass velocity which is defined as um=u+Jd/rho (Jd is mass diffusive flux Jd = -Dm*grad(rho))
At initial condition i assume that grad(rho) = 0 so (rhoE)initial is conputed using advective velocity.
*/
extern double **rho;
extern double **rhou;
extern double **rhov;
extern double **rhoE;

extern double **rhoN;
extern double **rhouN;
extern double **rhovN;
extern double **rhoEN;

extern double **rho0;
extern double **rhou0;
extern double **rhov0;
extern double **rhoE0;

extern double **rhoResArr;
extern double **rhouResArr;
extern double **rhovResArr;
extern double **rhoEResArr;

/*Primary variables declaration
extern std::vector<std::vector<double>>u, v, e, p, T, mu;*/

/*Auxilary variables
NOTE: in case of mass diffusion:
- S=mu*grad(U) with all components of vector U being computed by using advective velocity
(now i define rhoE_adv is total energy with advective kinetic energy)
*/
namespace BR1Vars {
//X direction
extern double **rhoX;
extern double **rhouX;
extern double **rhovX;
extern double **rhoEX;

/*Y direction*/
extern double **rhoY;
extern double **rhouY;
extern double **rhovY;
extern double **rhoEY;
}

namespace BR2Vars {
//X direction
extern double **rhoXVol;
extern double **rhouXVol;
extern double **rhovXVol;
extern double **rhoEXVol;

/*Y direction*/
extern double **rhoYVol;;
extern double **rhouYVol;
extern double **rhovYVol;
extern double **rhoEYVol;

//X direction
extern double **rhoXSurMaster;
extern double **rhouXSurMaster;
extern double **rhovXSurMaster;
extern double **rhoEXSurMaster;
extern double **rhoXSurSlave;
extern double **rhouXSurSlave;
extern double **rhovXSurSlave;
extern double **rhoEXSurSlave;

/*Y direction*/
extern double **rhoYSurMaster;
extern double **rhouYSurMaster;
extern double **rhovYSurMaster;
extern double **rhoEYSurMaster;
extern double **rhoYSurSlave;
extern double **rhouYSurSlave;
extern double **rhovYSurSlave;
extern double **rhoEYSurSlave;
}

namespace surfaceFields {
    /*Interface values*/
    extern double **rho;
    extern double **rhou;
    extern double **rhov;
    extern double **rhoE;

    /*Interface values*/
    //Auxilary equation
    extern double **aux_rho;
    extern double **aux_rhou;
    extern double **aux_rhov;
    extern double **aux_rhoE;

    //NSF equation
    //inviscid/viscous flux
    extern double **invis_rho;
    extern double **invis_rhou;
    extern double **invis_rhov;
    extern double **invis_rhoE;

    extern double **Vis_rho;
    extern double **Vis_rhou;
    extern double **Vis_rhov;
    extern double **Vis_rhoE;

    extern double **T;
}

namespace volumeFields {
    //Volume values
	//Cac array cua volumeFields la 3D array
    extern double **rhoVolGauss;
    extern double **rhouVolGauss;
    extern double **rhovVolGauss;
    extern double **rhoEVolGauss;
    extern double **drhoXVolGauss;
    extern double **drhoYVolGauss;
    extern double **T;
}

//Lax-Friedrich constant
extern double *LxFConst;
extern double *DiffusiveFluxConst;

//StiffMatrixCoefficients
extern double **stiffMatrixCoeffs;

//Limiting coefficients
extern double *theta1Arr;
extern double *theta2Arr;

namespace SurfaceBCFields
{
    //extern int **BCPointsInfor;
    extern std::vector<std::vector<int>> BCPointsInfor;
    /*
     * - 2 fields SurfaceBCFields::TBC, SurfaceBCFields::uBC va SurfaceBCFields::vBC dung de
     * luu gia tri TJump, uSlip va vSlip khi dung dieu kien bien temperatureJump va slip.
     * - gia tri tren cac field nay update theo thoi gian, khac gia tri o cac field bcValues::TBC, bcValues::uBC
     * va bcValues::vBC (la gia tri fixed doc tu cac file T, U ban dau).
    */
    extern double *TBc;
    extern double *uBc;
    extern double *vBc;
    extern int *localGlobalBCEdgesMatching;
}

namespace limitVal {
    extern bool *troubleCellsMarker;
}

//for debuging
namespace debug
{
	//extern std::vector<double> minRhoArr, minRhoeArr;
}

namespace parallelBuffer {
    extern double **rho;
    extern double **rhou;
    extern double **rhov;
    extern double **rhoE;

    //Neu div scheme la BR2, cac mang nay chua S_Surface
    extern double **drhoX;
    extern double **drhouX;
    extern double **drhovX;
    extern double **drhoEX;
    extern double **drhoY;
    extern double **drhouY;
    extern double **drhovY;
    extern double **drhoEY;

    extern double **xCoor;
    extern double **yCoor;
    extern double **aCoor;
    extern double **bCoor;

    extern double *theta1;
    extern double *theta2;
    extern int *elemType;
}
#endif // DYNAMICVARDECLARATION_H_INCLUDED
