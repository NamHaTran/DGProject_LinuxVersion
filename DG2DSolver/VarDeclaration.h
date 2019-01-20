#ifndef VARDECLARATION_H_INCLUDED
#define VARDECLARATION_H_INCLUDED
#include <string>
#include "ConstDeclaration.h"
#include <vector>

namespace systemVar
{
	extern std::string wD;
	extern std::string caseName;
	extern std::string pwd;
	extern std::string cmd;
	extern bool endKey;

	extern double CFL; //Courant number
	extern double Ttime; //Total time
	extern int wrtI; //write interval
	extern bool wrtLog, loadSavedCase;

	/*
	time discretization scheme	|keyWord	|index		|
	----------------------------|-----------|-----------|
	-Euler						|Euler		|1			|
	-Runge-Kutta 2 order		|RK2		|2			|
	-Runge-Kutta 3 order		|RK3		|3			|
	-Total Variation Diminishing|TVDRK2		|4			|
	Runge-Kutta 2 order			|			| 			|
	-Total Variation Diminishing|TVDRK2		|5			|
	Runge-Kutta 3 order			|			| 			|
	----------------------------|-----------|-----------|*/
	extern int ddtScheme;

	//constant for limiter
	extern double epsilon;

	extern int iterCount, savingCout;
	extern std::vector<double> rhoResNormVector,
		rhouResNormVector,
		rhovResNormVector,
		rhoEResNormVector;
	extern double rhoResNorm, rhouResNorm, rhovResNorm, rhoEResNorm;
}

namespace meshVar
{
	extern double Points[pointsArrSize][3];
	extern int Elements1D[pointsArrSize][3],
		Elements2D[elements2DArrSize][4],
		BoundaryType[bcSize][3];  //BoundaryType: column 0 is boundary group (from 0), column 1 is boundary type (1, 2, 3, 4), column 3 is boundary method
	extern int nelem1D, nelem2D, npoin, nBc;
	extern double normalVector[2][2 * elements2DArrSize];  //row 1 contents nx, row 2 contents ny, row 3 contents master element of edge
	extern int MasterElemOfEdge[2 * elements2DArrSize];  //array content master element of iedge, use it with normalVector to get information of normal vector of edge

	/*Default values*/
	//number of nodes per element
	extern const int nnode;
	//number of edges per element (default is 4)
	extern const int nedel;

	/*Elements surrounding points*/
	extern int esup1[4 * pointsArrSize], esup2[pointsArrSize + 1], inpoel[5][elements2DArrSize];

	/*Points surrounding points*/
	extern int psup1[5 * pointsArrSize], psup2[pointsArrSize + 1];

	/*Elements surrounding element*/
	extern int esuel[4][elements2DArrSize];

	/*Edges informations*/
	extern int inpoed[4][2 * elements2DArrSize];

	/*Edges of element*/
	extern int inedel[4][elements2DArrSize], //number of row is 4 because of default quad element
		ineled[3][2 * elements2DArrSize];

	/*Variables help to save mesh data*/
	extern int inpoedCount;  //can be used for normalVector, MasterElemOfEdge, ineled

	/*Jacobian*/
	extern double J2D[elements2DArrSize][maxGauss][maxGauss], J1D[pointsArrSize][2];

	/*derivatives dx/da, dx/db, dy/da, dy/db*/
	extern double dxa[elements2DArrSize][maxGauss][maxGauss],
		dxb[elements2DArrSize][maxGauss][maxGauss],
		dya[elements2DArrSize][maxGauss][maxGauss],
		dyb[elements2DArrSize][maxGauss][maxGauss];

	extern std::vector<std::vector<double>> geoCenter;
	extern std::vector<double> cellSize;
	extern std::vector<double> localCellSize;
	extern int numBCEdges;
}

namespace mathVar
{
	extern int nGauss, orderElem;
	extern double wGauss[maxGauss], xGauss[maxGauss], wGaussLobatto[maxGauss], xGaussLobatto[maxGauss];
	extern double B[maxOrder], dBa[maxOrder], dBb[maxOrder];
	extern double BPts_Tri[maxOrder][maxGauss][maxGauss], dBaPts_Tri[maxOrder][maxGauss][maxGauss], dBbPts_Tri[maxOrder][maxGauss][maxGauss],
		BPts_Quad[maxOrder][maxGauss][maxGauss], dBaPts_Quad[maxOrder][maxGauss][maxGauss], dBbPts_Quad[maxOrder][maxGauss][maxGauss];
	extern double GaussPts[maxGauss][maxGauss][2], //coordinate a is array (..,..,1), coordinate b is array (..,..,2)
		wGaussPts[maxGauss][maxGauss][2], //weights on a direction (w1) is array (..,..,1), weights on b direction (w2) is array (..,..,2)
		GaussLobattoPts[maxGauss][maxGauss][2],
		wGaussLobattoPts[maxGauss][maxGauss][2];
}

namespace material
{
	extern double gamma, R, Pr, As, Ts, Cp, Cv;
}

namespace iniValues
{
	extern double uIni, vIni, wIni, pIni, TIni, muIni, rhoIni, eIni;
}

namespace bcValues
{
	extern double uBC[bcSize],
		vBC[bcSize],
		wBC[bcSize],
		pBC[bcSize],
		TBC[bcSize];

	/*Values for Maxwell-Smoluchovsky condition*/
	extern double TWall[bcSize], uWall[bcSize], vWall[bcSize], wWall[bcSize];

	/*Variables help identify boundary condition at each group*/
	extern int 	UBcType[bcSize],
		pBcType[bcSize],
		TBcType[bcSize];
}

namespace refValues
{
	extern bool subsonic;
}

//time step
extern double dt;
extern double runTime;

namespace limitVal
{
	extern double TUp, TDwn;
	extern double rhoUp, rhoDwn;
	extern double rhoEUp, rhoEDwn;
	extern bool limitTOrNot, limitFlagLocal, limitFlagGlobal;
	extern int numOfLimitCell;

	namespace pAdaptive
	{
		extern bool limitFlagLocal, limitFlagGlobal;
		extern int numOfLimitCell;
	}

	extern std::vector<std::string> limiterName;
	extern bool PositivityPreserving, PAdaptive;
	namespace PositivityPreservingSettings
	{
		/*version:
		- 1: full
		- 2: simplified
		*/
		extern int version;
	}
}
#endif // VARDECLARATION_H_INCLUDED