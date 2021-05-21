#ifndef VARDECLARATION_H_INCLUDED
#define VARDECLARATION_H_INCLUDED
#include <string>
#include "ConstDeclaration.h"
#include <vector>
#include "DGMessagesLib.h"

namespace systemVar
{
    //! Working directory
	extern std::string wD;
    //! Name of submitted case
	extern std::string caseName;
	extern std::string pwd;
    //! Inputed comamnd
	extern std::string cmd;
    //! Flag of 'end' command
	extern bool endKey;

    //! Courant number
    extern double CFL;
    //! Total time
    extern double Ttime;
    //! write interval
    extern int wrtI;
	extern bool wrtLog, loadSavedCase, initializedOrNot, runPreProcess;

    /*! Type of time discretization scheme
     * Time discretization scheme	|Key word	|index		|
        ----------------------------|-----------|-----------|
        Euler						|Euler		|1			|
        Runge-Kutta 2 order		|RK2		|2			|
        Runge-Kutta 3 order		|RK3		|3			|
        Total Variation Diminishing|TVDRK2		|4			|
        Runge-Kutta 2 order			|			| 			|
        Total Variation Diminishing|TVDRK2		|5			|
        Runge-Kutta 3 order			|			| 			|*/
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
    //! Number of 1D element
    extern int nelem1D,
    //! Number of 2D element
    nelem2D,
    //! Number of point
    npoin,
    //! Number of boundary (which can be applied boundary condition)
    nBc;

	/*Default values*/
    //! Number of nodes per element
    extern const int nnode;
    //! Number of edges per element (default is 4)
    extern const int nedel;

    //! Array supports for finding elements surrounding point
    extern std::vector<int> esup1,
    //! Array supports for finding elements surrounding point
    esup2;

    //! Array supports for finding points surrounding point
    extern std::vector<int> psup1,
    //! Array supports for finding points surrounding point
    psup2;

    //! Variables help to save mesh data, can be used for normalVector, MasterElemOfEdge, ineled
    extern int inpoedCount;

    //! Number of boundary edge
    extern int numBCEdges;

    //! Number of BC groups
    extern int numBCGrp;
}

namespace mathVar
{
    extern int
    //! Number of basis function mode. Smallest is 0
    orderElem,
    orderOfAccuracy;
    //! Flag to trigger T failed error
    extern bool solveTFailed;

    //! Number of volume Gauss point in 1 direction. So total number of Gauss points is nVolGauss*nVolGauss
    extern int nGauss;

    //! Number of Gauss points on 1 surface
    //nSurGauss;
}

namespace material
{
    //! Heat capacity ratio, is the ratio of the heat capacity at constant pressure (Cp) to heat capacity at constant volume (Cv)
    extern double gamma,
    //! Ideal gas constant
    R,
    //! Prandtl number, is the ratio of momentum diffusivity to thermal diffusivity
    Pr,
    //! Heat capacity at constant pressure
    Cp,
    //! Heat capacity at constant volume
    Cv;
    namespace massDiffusion {
        extern
        //! Coefficient of self-diffusion
        double DmCoeff;
    }

    namespace viscosityCoeff {
        namespace Sutherland {
            extern double
            //! Coefficient of Sutherland model
            As,
            //! Reference Temperature of Sutherland model
            Ts;
        }

        namespace powerLaw_VHS {
            extern double
            //! Mass of 1 mole of gas used in Power Law - VHS model (Argon 39.948g/mol)
            molMass,
            //! Boltzmann constant
            kBoltzmann,
            //! Temperature exponent used in Power Law - VHS model
            omega,
            //! Reference temperature used in Power Law - VHS model
            TRef,
            //! Molecule diameter used in Power Law - VHS model
            dRef;
        }

        namespace constant {
            //! Constant mu
            extern double mu;
        }
    }

    namespace viscousityModel {
        extern bool constant,
        power_VHS,
        sutherland;
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

namespace flowProperties
{
    extern bool viscous, massDiffusion;
    extern double Mach;
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

namespace debugVars {
extern int element;
}
#endif // VARDECLARATION_H_INCLUDED
