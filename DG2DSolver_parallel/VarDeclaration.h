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
	extern bool wrtLog, loadSavedCase, initializedOrNot, runPreProcess;

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
	extern double rhoResNorm, rhouResNorm, rhovResNorm, rhoEResNorm;

    //1: BR1, 2: BR2
    extern int auxVariables;

    //For parallel computing
    extern int totalProc, currentProc;
    extern bool runDecomposeCaseFnc, parallelMode;
}

namespace meshVar
{
	extern int nelem1D, nelem2D, npoin, nBc;

	/*Default values*/
	//number of nodes per element
	extern const int nnode;
	//number of edges per element (default is 4)
	extern const int nedel;

	/*Elements surrounding points*/
    extern int esup1[4 * pointsArrSize], esup2[pointsArrSize + 1], inpoel[elements2DArrSize][5];

	/*Points surrounding points*/
	extern int psup1[5 * pointsArrSize], psup2[pointsArrSize + 1];

	/*Elements surrounding element*/
    extern int esuel[elements2DArrSize][4];

	/*Edges informations*/
    extern int inpoed[2 * elements2DArrSize][5];
    /*Cot 5 chua edgeId cua edge tai bien ung voi Element1D array*/

	/*Edges of element*/
    extern int inedel[elements2DArrSize][4], //number of row is 4 because of default quad element
        ineled[2 * elements2DArrSize][3];

	/*Variables help to save mesh data*/
	extern int inpoedCount;  //can be used for normalVector, MasterElemOfEdge, ineled

	extern int numBCEdges;

    //Number of BC groups
    extern int numBCGrp;
}

namespace mathVar
{
    extern int nGauss, orderElem, orderOfAccuracy;
}

namespace material
{
	extern double gamma, R, Pr, As, Ts, Cp, Cv;
    namespace massDiffusion {
    extern
    //coefficient of self-diffusion
    double DmCoeff;
    }
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

namespace flowProperties
{
    extern bool viscous, massDiffusion;
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

namespace controlFlag {
namespace sequence {
extern bool checkUnvReader,
checkBCsHelper,
checkUnvHelper,
reSubmit,
mappResults,
exportMeshToMetis,
testMeshPartitionResult,
debug_checkElement,
decomposeCase;
}

namespace parallel {
extern bool checkDGRun;
}
}
#endif // VARDECLARATION_H_INCLUDED
