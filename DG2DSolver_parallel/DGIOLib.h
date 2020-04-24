#ifndef DGIOLIB_H_INCLUDED
#define DGIOLIB_H_INCLUDED
#include <string>
#include <vector>
namespace IO
{
	/*Function displays program's logo*/
	void dispLogo();

	/*Function gets working directory inputted from keyboard, case name, etc...*/
	void getCase();

	/*Function reads mesh data from files*/
    void loadMesh(std::string mode);

	/*Function saves mesh data to files*/
    void SaveMeshInfor(std::string mode);

	/*Function reads constants from folder constant*/
    void loadConstants(std::string mode);

	/*Function reads u, v, p, U values from folder 0*/
    void loadpTU(std::string mode);

	/*Function reads non scalar values from file*/
    void readNonScalar(std::string mode);

	/*Function reads scalar values from file*/
    void readScalar(std::string fileName, std::string mode);

    /*Function reads number of total processes in DGOPtions*/
    void readNumberOfCores();

	/*Function read informations of limiters*/
	void loadLimiterSettings();

	/*User's guide:
	This function returns datas of type double and type int read from files.
	Input arguments:
	- fileName: name of file you want to read (file name and extension)
	- direction: working diectory contents data file (with no "\\" characters at the end)
	- keyWordsDbl: array contents keyWords of double values listed in file
	- keyWordsInt: array contents keyWords of int values listed in file
	- outDbl: output array contents double values
	- outInt: output array contents int values*/
	void readDataFile(std::string fileName, std::string direction, std::string keyWordsDbl[], std::string keyWordsInt[], std::string keyWordsBool[], std::string keyWordsStr[], double *outDbl, int *outInt, bool *outBool, std::string *Str, int numParamDbl, int numParamInt, int numParamBool, int numParamStr);

    void readDecomposedMeshInfor();

    std::tuple<double**,int> read2DArray(int column, std::string location, std::string fileName);

    std::tuple<int**,int> read2DIntArray(int column, std::string location, std::string fileName);

    std::tuple<double*,int> read1DArray(std::string location, std::string fileName);

    std::tuple<int*,int> read1DIntArray(std::string location, std::string fileName);

	/*Function writes residuals on console*/
	void residualOutput(double rhoRes, double rhouRes, double rhovRes, double rhoERes);

    /*Function writes discreted fields to files*/
    void writeDiscretedFields(std::string Loc, std::string fileName, double **Var);

    /*Function reads discreted fields to files*/
    void readDiscretedFields(std::string Loc, std::string fileName, double **Var);

    void saveCase();

    void saveCase_reconstruct();

    void loadCase(std::string mode);

    void loadTime();

    void write2DIntArrayToFile(std::vector<std::vector<int>> &array, std::string loc, std::string name, int numRow, int numCol);

    void write2DDoubleArrayToFile(std::vector<std::vector<double>> &array, std::string loc, std::string name, int numRow, int numCol);

    void openFileToAppend(std::string Loc, std::string content);

    void write2DDoubleVectorToFile(std::string location, std::string fileName, std::vector<std::vector<double>> &vector);

    void write1DDoubleVectorToFile(std::string location, std::string fileName, double *vector, int length);

    void write1DIntVectorToFile(std::string location, std::string fileName, int *vector, int length);

    void writeResiduals(int iter, double rhoRes, double rhouRes, double rhovRes, double rhoERes);

	namespace importCase {

		void importResultsFromAnotherCase();

		void mappSourceToCurrent(std::string fileLoc, std::vector<std::vector<double>> &currentResult);
	}
}
#endif // DGIOLIB_H_INCLUDED
