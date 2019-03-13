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
	void loadMesh();

	/*Function saves mesh data to files*/
	void SaveMeshInfor();

	/*Function reads constants from folder constant*/
	void loadConstants();

	/*Function reads u, v, p, U values from folder 0*/
	void loadpTU();

	/*Function reads non scalar values from file*/
	void readNonScalar();

	/*Function reads scalar values from file*/
	void readScalar(std::string fileName);

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

	/*Function writes residuals on console*/
	void residualOutput(double rhoRes, double rhouRes, double rhovRes, double rhoERes);

	void saveCase();

	void loadCase();

	namespace importCase {

		void importResultsFromAnotherCase();

		void mappSourceToCurrent(std::string fileLoc, std::vector<std::vector<double>> &currentResult);
	}
}
#endif // DGIOLIB_H_INCLUDED
