#include "VarDeclaration.h"
#include <string>
#include "ConstDeclaration.h"
#include <vector>

namespace systemVar
{
	std::string wD;
	std::string caseName;
	std::string pwd;
	std::string cmd;
	bool endKey(false);

	double CFL(0.5); //Courant number
	double Ttime(0.0); //Total time
	int wrtI(0); //write interval
    bool wrtLog(true), loadSavedCase(false);

    int ddtScheme(1);
    double epsilon(1e-13);

	int iterCount(0), savingCout(0);
	double rhoResNorm(1.0), rhouResNorm(1.0), rhovResNorm(1.0), rhoEResNorm(1.0);

	bool initializedOrNot(false), runPreProcess(false);

    //1: BR1, 2: BR2
    int auxVariables(1);

    //For parallel computing
    int totalProc(4), currentProc(0);
    bool runDecomposeCaseFnc(false), parallelMode(true);
}

namespace meshVar
{
    int nelem1D(0), nelem2D(0), npoin(0), nBc(0);

	/*Default value*/
	//number of nodes per element
	int const nnode(4);
	//number of edges per element (default is 4)
	int const nedel(4);  //change nedel value at file .h

	/*Elements surrounding point*/
    std::vector<int> esup1, esup2;
	std::vector<std::vector<int>> inpoel(elements2DArrSize,std::vector<int>(5));

	/*Points surrounding point*/
	std::vector<int> psup1(5 * pointsArrSize), psup2(pointsArrSize + 1);

	/*Elements surrounding element*/
    std::vector<std::vector<int>> esuel(elements2DArrSize,std::vector<int>(4)); //Default is 4 faces

	/*Edges informations*/
    std::vector<std::vector<int>> inpoed(2*elements2DArrSize,std::vector<int>(5));
	/*column 3 contents group which edge belongs to (group 0 is internal group),
	column 4 contents type of boundary (type 0 is internal edge)*/

	/*Edges of element*/
    //column index is element index, each row in column contents index of edge belong to element, number of row is 4 because of default quad element
    //column index is edge index, each row in column contents index of element which edge is belong to, row 3 contents pointer
    std::vector<std::vector<int>> inedel(elements2DArrSize,std::vector<int>(4)),
    	ineled(2*elements2DArrSize,std::vector<int>(3));

	/*Variables help to save mesh data*/
	int inpoedCount(0);  //can be used for normalVector, MasterElemOfEdge, ineled

	int numBCEdges(0);

    //Number of BC groups
    int numBCGrp(0);
}

namespace mathVar
{
    int nGauss(2), orderElem(0), orderOfAccuracy(0);
}

namespace material
{
	double gamma(1.4), R(287.0), Pr(0.72), As(0.001), Ts(110.4), Cp(0.0), Cv(0.0);
    namespace massDiffusion {
    double
    //coefficient of self-diffusion
    DmCoeff(0.0);
    }
}

namespace iniValues
{
	double uIni(0.0), vIni(0.0), wIni(0.0), pIni(0.0), TIni(0.0), muIni(0.0), rhoIni(0.0), eIni(0.0);
}

namespace bcValues
{
	double uBC[bcSize] = { 0.0 },
		vBC[bcSize] = { 0.0 },
		wBC[bcSize] = { 0.0 },
		pBC[bcSize] = { 0.0 },
		TBC[bcSize] = { 0.0 };

	/*Values for Maxwell-Smoluchovsky condition*/
	double TWall[bcSize] = {}, uWall[bcSize] = {}, vWall[bcSize] = {}, wWall[bcSize] = {};

	/*Variables help identify boundary condition at each group*/
	int UBcType[bcSize] = {};
	int pBcType[bcSize] = {};
	int TBcType[bcSize] = {};
}

namespace refValues
{
	bool subsonic(true);
}

namespace flowProperties
{
    bool viscous(true), massDiffusion(false);
}

//time step
double dt(1e-20);
double runTime(0.0);

namespace limitVal
{
	double TUp(5000), TDwn(10);
	double rhoUp(12.5), rhoDwn(0.0001);
	double rhoEUp(0.0), rhoEDwn(0.0);  //computed at initialValue function
	bool limitTOrNot(false), limitFlagLocal(false), limitFlagGlobal(false);
	int numOfLimitCell(0);

	namespace pAdaptive
	{
		bool limitFlagLocal(false), limitFlagGlobal(false);
		int numOfLimitCell(0);
	}

	std::vector<std::string> limiterName;

    //keys of limiter
    bool PositivityPreserving(false), PAdaptive(false), massDiffusion(false),
        runningMassDiffLimiter(false) //bien kiem tra xem mass diffusion limiter co dang chay hay khong
    ;
	namespace PositivityPreservingSettings
	{
		/*version:
		- 1: full
		- 2: simplified
		*/
		int version(2);
	}
}

namespace controlFlag {
namespace sequence {
bool checkUnvReader(false),
checkBCsHelper(false),
checkUnvHelper(false),
reSubmit(false),
exportMeshToMetis(false),
testMeshPartitionResult(false),
debug_checkElement(false),
decomposeCase(false);
}

namespace parallel {
bool checkDGRun(false),
mapResults(false);
}
}
