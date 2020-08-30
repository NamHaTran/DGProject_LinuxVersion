#ifndef VARDECLARATION_H_INCLUDED
#define VARDECLARATION_H_INCLUDED
#include <string>
#include "ConstDeclaration.h"
#include <vector>
#include "DGMessagesLib.h"

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

    extern int iterCount, savingCout, loadConstCount;
	extern double rhoResNorm, rhouResNorm, rhovResNorm, rhoEResNorm;

    //1: BR1, 2: BR2
    extern int auxVariables;

    //For parallel computing
    extern int totalProc, currentProc;
    extern bool runDecomposeCaseFnc, parallelMode, runCheckParMesh;

    //Detect first iteration
    extern bool firstIter;
	
	extern std::string headerFile;

    extern int **sendRecvOrder;

    extern int sendRecvOrder_length;

    //Read/write mode
    extern std::string readWriteMode;
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
    extern std::vector<int> esup1, esup2;

	/*Points surrounding points*/
	extern std::vector<int> psup1, psup2;

	/*Variables help to save mesh data*/
	extern int inpoedCount;  //can be used for normalVector, MasterElemOfEdge, ineled

	extern int numBCEdges;

    //Number of BC groups
    extern int numBCGrp;
}

namespace mathVar
{
    extern int nGauss, orderElem, orderOfAccuracy;
    extern bool solveTFailed;
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
	extern double uBCFixed[bcSize],
		vBCFixed[bcSize],
		wBCFixed[bcSize],
		pBCFixed[bcSize],
		TBCFixed[bcSize];

	/*Variables help identify boundary condition at each group*/
	extern int 	UBcType[bcSize],
		pBcType[bcSize],
		TBcType[bcSize];

    /*Variables of slip & temperature jump boundary conditions*/
    extern double sigmaU, sigmaT;

    /*Flags of time varying BCs*/
    //U
    extern bool slipBCFlag;
    //T
    extern bool temperatureJump;
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
    extern bool PositivityPreserving, PAdaptive, massDiffusion, runningMassDiffLimiter//bien nay kiem tra xem mass diffusion limiter co dang chay hay khong
    ;
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
exportMeshToMetis,
testMeshPartitionResult,
debug_checkElement,
decomposeCase,
reconstructLatestTime,
checkPartitionedMesh;
}

namespace parallel {
extern bool checkDGRun,
mapResults;
}

/* Cac flag trong namespace nay check xem folder chua kq co file TSurface, uSurface ... hay khong
 * Neu khong co, cac gia tri cua cac field nay set ve gia tri mac dinh o dieu kien bien
*/
namespace fileAvailFlags {
extern bool fileTSurface,
    fileuSurface,
    filevSurface;
}
}

namespace warningFlag {
extern bool reversedFlowOccur;
}

namespace DGSchemes {
    //Dung cho truong hop mass diffusion, true neu giai T implicit, false neu giai T explicit
    extern bool solveTImplicitly;

    namespace fluxControl {
        //Tai version hien tai, flux cua bien phu va viscous flux la central flux, flux control chi apply cho
        //convective flux
        extern bool LxF, Roe, HLL, HLLC, central;
    }
}
#endif // VARDECLARATION_H_INCLUDED
