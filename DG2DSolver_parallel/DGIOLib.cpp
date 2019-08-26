#include "DGIOLib.h"
#include "DGMath.h"
#include "DGMessagesLib.h"
#include "ConstDeclaration.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "DGProcLib.h"
#include <math.h>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <cstdlib> //include this library to create folder in linux
#include <mpi.h>

namespace IO
{
	void dispLogo()
	{
		std::string logoStr(" ");
		logoStr = message::headerFile();
		std::cout << logoStr << std::endl;
	}

	void getCase()
	{
        if (systemVar::currentProc==0)
        {
            std::cout << "***Getting case's information***\n";
        }
		systemVar::wD = auxUlti::workingdir();

		/*Get caseName from SubmitCase*/
		std::string submitLoc(" ");  //submitLoc contents location of submitingCase
        submitLoc = systemVar::wD + "/CASES/SubmitCase.txt";

		std::ifstream submitCaseFlux(submitLoc.c_str());
		if (submitCaseFlux)
		{
			std::string line, keyWord;
            int keyWFlag(0);
			while (std::getline(submitCaseFlux, line))  //Read submitCaseFlux line by line
			{
				std::istringstream line2str(line);
				std::vector<std::string> ptr;
				//Split <line2str> stringstream into array of words
				while ((line2str >> keyWord))
				{
					ptr.push_back(keyWord);
				}

                int numWrd = static_cast<int>(ptr.size());

				if (numWrd == 2)
				{
					std::string str1(ptr[0]), str2("caseName");
					if (str1.compare(str2) == 0)
					{
						systemVar::caseName = ptr[1];
						keyWFlag = 1;
                        if (systemVar::currentProc==0)
                        {
                            std::cout << "	Case " << systemVar::caseName << " has been submitted\n";
                        }
					}
				}
			}
			if (keyWFlag == 0)
			{
				std::cout << message::undfKeyW("caseName", submitLoc) << std::endl;
			}
		}
		else
		{
            message::writeLog((systemVar::wD + "/CASES/"), "", message::opFError("SubmitCase.txt", submitLoc));
		}
        systemVar::pwd = systemVar::wD + "/CASES/" + systemVar::caseName;
	}

    void loadMesh(std::string mode)
	{
        std::string  Elem1DLoc, ptLoc, Elem2DLoc, bcLoc;
		/*Declare loading locations*/
        if (mode.compare("p")!=0)
        {
            Elem1DLoc = systemVar::pwd + "/Constant/Mesh/Elements1D.txt";
            ptLoc = systemVar::pwd + "/Constant/Mesh/Points.txt";
            Elem2DLoc = systemVar::pwd + "/Constant/Mesh/Elements2D.txt";
            bcLoc = systemVar::pwd + "/Constant/boundaryPatch.txt";
        }
        else {
            Elem1DLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc) + "/Constant/Mesh/Elements1D.txt";
            ptLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc) + "/Constant/Mesh/Points.txt";
            Elem2DLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc) + "/Constant/Mesh/Elements2D.txt";
            bcLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc) + "/Constant/boundaryPatch.txt";
        }

		/*Load Points*/
		std::ifstream ptFlux(ptLoc.c_str());
		if (ptFlux)
		{
			std::string line(" ");
			double x, y, z;
			meshVar::npoin = 0;
			while (std::getline(ptFlux, line))
			{
				auxUlti::addRowTo2DDoubleArray(meshVar::Points, 3);
				std::istringstream ptData(line);
				ptData >> meshVar::npoin >> x >> y >> z;
				meshVar::Points[meshVar::npoin - 1][0] = x;
				meshVar::Points[meshVar::npoin - 1][1] = y;
				meshVar::Points[meshVar::npoin - 1][2] = z;
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("Points.txt", ptLoc));
		}

		/*Load 1D elements*/
		std::ifstream Elem1DFlux(Elem1DLoc.c_str());
		if (Elem1DFlux)
		{
			int node1(0), node2(0), bcGrp(0);  //nelem1D is number of 1D elements
			std::string line(" ");
			while (std::getline(Elem1DFlux, line))
			{
				auxUlti::addRowTo2DIntArray(meshVar::Elements1D, 3);
				std::istringstream elData(line);
				elData >> meshVar::nelem1D >> node1 >> node2 >> bcGrp;
				meshVar::Elements1D[meshVar::nelem1D - 1][0] = node1 - 1;
				meshVar::Elements1D[meshVar::nelem1D - 1][1] = node2 - 1;
				meshVar::Elements1D[meshVar::nelem1D - 1][2] = bcGrp;
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("Elements1D.txt", Elem1DLoc));
		}

		/*Load 2D elements*/
		std::ifstream Elem2DFlux(Elem2DLoc.c_str());
		if (Elem2DFlux)
		{
			int temp(0), node1(0), node2(0), node3(0), node4(0);  //nelem2D is number of 2D elements
			std::string line(" ");
			meshVar::nelem2D = 0;
			while (std::getline(Elem2DFlux, line))
			{
				auxUlti::addRowTo2DIntArray(meshVar::Elements2D, 4);
				meshVar::nelem2D++;
				std::istringstream elData(line);
				elData >> temp >> node1 >> node2 >> node3 >> node4;
				meshVar::Elements2D[meshVar::nelem2D - 1][0] = node1 - 1;
				meshVar::Elements2D[meshVar::nelem2D - 1][1] = node2 - 1;
				meshVar::Elements2D[meshVar::nelem2D - 1][2] = node3 - 1;
				meshVar::Elements2D[meshVar::nelem2D - 1][3] = node4 - 1;
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("Elements2D.txt", Elem2DLoc));
		}

		/*Load boundary conditions*/
		std::ifstream bcFlux(bcLoc.c_str());
		if (bcFlux)
		{
            int boundIndex(0);
			std::string line(" "), keyWord(" ");
			while (std::getline(bcFlux, line))
			{
				std::istringstream line2str(line);
				std::vector<std::string> ptr;
				//Split <line2str> stringstream into array of words
				while ((line2str >> keyWord))
				{
					if ((keyWord.compare(" ") != 0))
					{
						ptr.push_back(keyWord);
					}
				}

                int numWrd = static_cast<int>(ptr.size());
				if (numWrd != 0)
				{
					std::string str1(ptr[0]);
					if (str1.compare("{") == 0)  //0 means str1 is the same as str2
					{
						//start group
						boundIndex++;
					}

					if ((numWrd >= 2)&(boundIndex != 0))
					{
						if ((boundIndex <= bcSize))
						{
							auxUlti::addRowTo2DIntArray(meshVar::BoundaryType, 3);
							meshVar::BoundaryType[boundIndex - 1][0] = boundIndex;
							std::string str2 = ptr[1];
							if (str1.compare("Type") == 0)
							{
								if (str2.compare("wall") == 0)
								{
									meshVar::BoundaryType[boundIndex - 1][1] = 1;  //type 1 is wall
								}
								else if (str2.compare("patch") == 0)
								{
									meshVar::BoundaryType[boundIndex - 1][1] = 2;  //type 2 is patch
								}
								else if (str2.compare("symmetry") == 0)
								{
									meshVar::BoundaryType[boundIndex - 1][1] = 3;  //type 3 is symmetry
								}
                                else if (str2.compare("matched") == 0)
                                {
                                    meshVar::BoundaryType[boundIndex - 1][1] = 4;  //type 4 is matched (boundary condition type for parallel computing)
                                }
								else
								{
									std::string errorStr = "Boundary type <" + ptr[1] + R"(> is unknown, available boundary types in file boundaryPatch are:
	wall
	patch
    symmetry
    matched)";
									message::writeLog(systemVar::pwd, systemVar::caseName, errorStr);
								}
							}
						}
						else
						{
                            std::string errorStr = "Maximum number of boundaries DG2D solver supports is only " + std::to_string(bcSize) + ". Make sure the number of boundaries is less than " + std::to_string(bcSize);
							message::writeLog(systemVar::pwd, systemVar::caseName, errorStr);
						}
					}
				}
				meshVar::nBc = boundIndex;
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("boundaryPatch.txt", bcLoc));
		}

        if (mode.compare("p")==0)
        {
            std::string mshConnect = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc) + "/Constant/Mesh/meshConnection.txt";
            std::ifstream mshConnectFlux(mshConnect.c_str());
            if (mshConnectFlux)
            {
                std::string line("");
                int counter = 0;
                while (std::getline(mshConnectFlux, line))
                {
                    auxUlti::addRowTo2DIntArray(meshVar::meshConnection, 3);
                    std::istringstream Data(line);
                    Data >> meshVar::meshConnection[counter][0] >> meshVar::meshConnection[counter][1] >> meshVar::meshConnection[counter][2];
                    counter++;
                }
            }
            else
            {
                message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("meshConnection.txt", mshConnect));
            }
        }
	}

    void SaveMeshInfor(std::string mode)
	{
		/*List of arrays need to be saved:
		- inedel: array contents information of edges belong to elements, column index is element index, each row in column contents index of edge belong to element, number of row is 4 because of default quad element.
		- ineled: array contents information of elements sharing edges, column index is edge index, each row in column contents index of element which edge is belong to, row 3 contents pointer.
		- inpoed: array contents information of points belong to edges, column 3 contents group which edge belongs to (group 0 is internal group), column 4 contents type of boundary (type 0 is internal edge)
		- normalVector: array contents information of normal vector of edges
		- MasterElemOfEdge: array content master element of iedge, use it with normalVector to get information of normal vector of edge*/

        std::string headerFile(message::headerFile()), inedelLoc, ineledLoc, inpoedLoc, normVectorLoc, MasterElemOfEdgeLoc;

        std::string  Elem1DLoc, ptLoc, Elem2DLoc, bcLoc;
        /*Declare loading locations*/
        if (mode.compare("p")!=0)
        {
            inedelLoc = systemVar::pwd + "/Constant/Mesh/inedel.txt";
            ineledLoc = systemVar::pwd + "/Constant/Mesh/ineled.txt";
            inpoedLoc = systemVar::pwd + "/Constant/Mesh/inpoed.txt";
            normVectorLoc = systemVar::pwd + "/Constant/Mesh/normalVector.txt";
            MasterElemOfEdgeLoc = systemVar::pwd + "/Constant/Mesh/MasterElemOfEdge.txt";
        }
        else {
            inedelLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc)+ "/Constant/Mesh/inedel.txt";
            ineledLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc)+ "/Constant/Mesh/ineled.txt";
            inpoedLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc)+ "/Constant/Mesh/inpoed.txt";
            normVectorLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc)+ "/Constant/Mesh/normalVector.txt";
            MasterElemOfEdgeLoc = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc)+ "/Constant/Mesh/MasterElemOfEdge.txt";
        }

		/*inedel*/
		std::ofstream Fluxinedel(inedelLoc.c_str());
		if (Fluxinedel)
		{
			Fluxinedel << headerFile << std::endl << "	inedel array\n" << "\n";
            for (int i = 0; i < meshVar::nelem2D; i++)
			{
                for (int j = 0; j < 4; j++)
				{
					Fluxinedel << meshVar::inedel[i][j] << " ";
				}
				Fluxinedel << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("inedel.txt", inedelLoc));
		}

		/*ineled*/
		std::ofstream Fluxineled(ineledLoc.c_str());
		if (Fluxineled)
		{
			Fluxineled << headerFile << std::endl << "	ineled array\n" << "\n";
            for (int i = 0; i < meshVar::inpoedCount; i++)
			{
                for (int j = 0; j < 2; j++)
				{
					Fluxineled << meshVar::ineled[i][j] << " ";
				}
				Fluxineled << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("ineled.txt", ineledLoc));
		}

		/*inpoed*/
		std::ofstream Fluxinpoed(inpoedLoc.c_str());
		if (Fluxinpoed)
		{
			Fluxinpoed << headerFile << std::endl << "	inpoed array\n" << "\n";
            for (int i = 0; i < meshVar::inpoedCount; i++)
			{
                for (int j = 0; j < 4; j++)
				{
					Fluxinpoed << meshVar::inpoed[i][j] << " ";
				}
				Fluxinpoed << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("inpoed.txt", inpoedLoc));
		}

		/*normalVector*/
		std::ofstream FluxnormVector(normVectorLoc.c_str());
		if (FluxnormVector)
		{
			FluxnormVector << headerFile << std::endl << "	normalVector array\n" << "\n";
			for (int i = 0; i < meshVar::inpoedCount; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					FluxnormVector << meshVar::normalVector[i][j] << " ";
				}
				FluxnormVector << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("normalVector.txt", normVectorLoc));
		}

		/*inpoed*/
		std::ofstream FluxMasterElemOfEdge(MasterElemOfEdgeLoc.c_str());
		if (FluxMasterElemOfEdge)
		{
			FluxMasterElemOfEdge << headerFile << std::endl << "	MasterElemOfEdge array\n" << "\n";
			for (int i = 0; i < meshVar::inpoedCount; i++)
			{
				FluxMasterElemOfEdge << meshVar::MasterElemOfEdge[i] << " ";
				FluxMasterElemOfEdge << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("MasterElemOfEdge.txt", MasterElemOfEdgeLoc));
		}
	}

    void loadConstants(std::string mode)
	{
		/*Read DGOptions*/
		std::string DGOptfileName("DGOptions.txt");
        std::string DGOptLoc(systemVar::wD + "/CASES/" + systemVar::caseName + "/System");
        std::string DGOptkeyWordsDouble[2] = { "CourantNumber", "totalTime(s)" }, DGOptkeyWordsInt[4] = {"numberOfGaussPoints","orderOfAccuracy", "writeInterval","totalProcess"}, DGOptkeyWordsBool[2] = { "writeLog", "loadSavedCase"}, DGOptkeyWordsStr[2] = {"ddtScheme", "runningMode"};
        double DGOptoutDB[2] = {};
        int DGOptoutInt[4] = {};
		bool DGOptoutBool[2] = {};
		std::string DGOptoutStr[2] = {};

        readDataFile(DGOptfileName, DGOptLoc, DGOptkeyWordsDouble, DGOptkeyWordsInt, DGOptkeyWordsBool, DGOptkeyWordsStr, DGOptoutDB, DGOptoutInt, DGOptoutBool, DGOptoutStr, 2, 4, 2, 2);
		
		systemVar::CFL = DGOptoutDB[0];
		systemVar::Ttime = DGOptoutDB[1];
        mathVar::nGauss = DGOptoutInt[0];
        mathVar::orderElem = DGOptoutInt[1];
        systemVar::wrtI = DGOptoutInt[2];
        systemVar::totalProc = DGOptoutInt[3];
		systemVar::wrtLog = DGOptoutBool[0];
		systemVar::loadSavedCase = DGOptoutBool[1];

		if (DGOptoutStr[0].compare("Euler") == 0)
		{
			systemVar::ddtScheme = 1;
		}
		else if (DGOptoutStr[0].compare("TVDRK2") == 0)
		{
			systemVar::ddtScheme = 2;
		}
		else if (DGOptoutStr[0].compare("TVDRK3") == 0)
		{
			systemVar::ddtScheme = 3;
		}

        if (DGOptoutStr[1].compare("parallel") == 0)
        {
            systemVar::parallelMode = true;
        }
        else if (DGOptoutStr[1].compare("sequence") == 0)
        {
            systemVar::parallelMode = false;
            systemVar::totalProc=1;
        }
		
		/*Read Material*/
		std::string MatfileName("Material.txt");
        std::string MatLoc;
        if (mode.compare("p")==0)
        {
            MatLoc = systemVar::wD + "/CASES/" + systemVar::caseName + "/Processor"+std::to_string(systemVar::currentProc) + "/Constant";
        }
        else {
            MatLoc = systemVar::wD + "/CASES/" + systemVar::caseName + "/Constant";
        }
        std::string MatkeyWordsDouble[6] = { "gammaRatio", "gasConstant", "PrandtlNumber", "SutherlandAs", "SutherlandTs" , "DmCoef"}, MatkeyWordsInt[1] = {}, MatkeyWordsBool[1] = {}, MatkeyWordsStr[1] = {};
        double MatoutDB[6] = {};
		int MatoutInt[1] = {};
		bool MatoutBool[1] = {};
		std::string MatoutStr[1] = {};

        readDataFile(MatfileName, MatLoc, MatkeyWordsDouble, MatkeyWordsInt, MatkeyWordsBool, MatkeyWordsStr, MatoutDB, MatoutInt, MatoutBool, MatoutStr, 6, 0, 0, 0);
		material::gamma = MatoutDB[0];
		material::R = MatoutDB[1];
		material::Pr = MatoutDB[2];
		material::As = MatoutDB[3];
		material::Ts = MatoutDB[4];
        material::massDiffusion::DmCoeff = MatoutDB[5];

		material::Cp = material::R*material::gamma / (material::gamma - 1);
		material::Cv = material::Cp - material::R;

        /*Read FlowProperties*/
        std::string fileName("FlowProperties.txt");
        std::string Loc(systemVar::wD + "/CASES/" + systemVar::caseName + "/System");
        std::string keyWordsDouble[1] = {}, keyWordsInt[1] = {}, keyWordsBool[2] = {"Viscosity", "MassDiffusion"}, keyWordsStr[1] = {};
        double outDB[1] = {};
        int outInt[1] = {};
        bool outBool[2] = {};
        std::string outStr[1] = {};

        readDataFile(fileName, Loc, keyWordsDouble, keyWordsInt, keyWordsBool, keyWordsStr, outDB, outInt, outBool, outStr, 0, 0, 2, 0);
        flowProperties::viscous=outBool[0];
        flowProperties::massDiffusion=outBool[1];
        if (!flowProperties::massDiffusion)
        {
            material::massDiffusion::DmCoeff=0.0;
        }

        /*Read DGSchemes*/
        fileName=("DGSchemes.txt");
        Loc=(systemVar::wD + "/CASES/" + systemVar::caseName + "/System");
        std::string DGSchemeskeyWordsDouble[1] = {}, DGSchemeskeyWordsInt[1] = {}, DGSchemeskeyWordsBool[1] = {}, DGSchemeskeyWordsStr[1] = {"diffusionTermScheme"};
        double DGSchemesoutDB[1] = {};
        int DGSchemesoutInt[1] = {};
        bool DGSchemesoutBool[1] = {};
        std::string DGSchemesoutStr[1] = {};

        readDataFile(fileName, Loc, DGSchemeskeyWordsDouble, DGSchemeskeyWordsInt, DGSchemeskeyWordsBool, DGSchemeskeyWordsStr, DGSchemesoutDB, DGSchemesoutInt, DGSchemesoutBool, DGSchemesoutStr, 0, 0, 0, 1);
        if (DGSchemesoutStr[0].compare("BR1")==0)
        {
            systemVar::auxVariables=1;
        }
        else if (DGSchemesoutStr[0].compare("BR2")==0)
        {
            systemVar::auxVariables=2;
        }
        else {
            std::string str0("diffusionTermScheme '"+DGSchemesoutStr[0]+"' is not a diffusion scheme.");
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, str0);
        }
	}

	void loadLimiterSettings()
	{
        std::string FileDir(systemVar::wD + "/CASES/" + systemVar::caseName + "/System"), fileName("LimiterSettings.txt");
        std::string FileLoc(FileDir + "/" + fileName);
		std::ifstream FileFlux(FileLoc.c_str());
		if (FileFlux)
		{
			std::string line, Word;
			while (std::getline(FileFlux, line))
			{
			label:
				std::istringstream line2str(line);
				std::vector<std::string> ptr;
				//Split <line2str> stringstream into array of words
				while ((line2str >> Word))
				{
					ptr.push_back(Word);
				}

                int numWrd = static_cast<int>(ptr.size());
				if (numWrd >= 2)
				{
					std::istringstream strdata(ptr[1]);

					if (ptr[0].compare("limiter") == 0) //get selected limiter(s)
					{
						for (int i = 0; i < numWrd - 1; i++)
						{
							limitVal::limiterName.push_back(ptr[i + 1]);
						}
					}

					if (limitVal::limiterName.size() > 0)
					{
                        for (int ilimiter = 0; ilimiter < static_cast<int>(limitVal::limiterName.size()); ilimiter++)
						{
							if (limitVal::limiterName[ilimiter].compare("PositivityPreserving") == 0) //PositivityPreserving settings
							{
								limitVal::PositivityPreserving = true;
								std::getline(FileFlux, line); //jump
								while (std::getline(FileFlux, line))
								{
									std::istringstream line2str(line);
									while ((line2str >> Word))
									{
										if (Word.compare("version") == 0)
										{
											line2str >> Word;
											if (Word.compare("simplified") == 0)
											{
												limitVal::PositivityPreservingSettings::version = 2;
											}
											else if (Word.compare("full") == 0)
											{
												limitVal::PositivityPreservingSettings::version = 1;
											}
											else
											{
												std::cout << Word << " is not available version of Positivity Preserving limiter, version will be set to Simplified as a default\n";
												limitVal::PositivityPreservingSettings::version = 2;
											}
											line2str >> Word;
											goto label;
										}
									}
								}
							}
							if ((limitVal::limiterName[ilimiter].compare("PAdaptive") == 0) || (limitVal::limiterName[ilimiter].compare("pAdaptive") == 0)) //PositivityPreserving settings
							{
								limitVal::PAdaptive = true;
								std::getline(FileFlux, line); //jump
								goto label;
							}
						}
					}
				}
			}
		}
		else
		{
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, FileLoc));
		}
	}

    void loadpTU(std::string mode)
	{
		/*Read U*/  //U must be read first of all
        readNonScalar(mode);

		/*Read T*/
        readScalar("T", mode);

		/*Read p*/
        readScalar("p", mode);
	}

	void readDataFile(std::string fileName, std::string direction, std::string keyWordsDbl[], std::string keyWordsInt[], std::string keyWordsBool[], std::string keyWordsStr[], double *outDbl, int *outInt, bool *outBool, std::string *outStr, int numParamDbl, int numParamInt, int numParamBool, int numParamStr)  //Declaration of funciton which returns pointer of 1D array
	{
		/*User's guide:
		This function returns datas of type double and type int read from files.
		Input arguments:
		- fileName: name of file you want to read (file name and extension)
        - direction: working diectory contents data file (with no "/" characters at the end)
		- keyWordsDbl: array contents keyWords of double values listed in file
		- keyWordsInt: array contents keyWords of int values listed in file
		- keyWordsBool: array contents keyWords of bool values listed in file
		- keyWordsStr: array contents keyWords of string values listed in file
		- outDbl: output array contents double values
		- outInt: output array contents int values
		- numParamDbl, numParamInt, numParamBool, numParamStr: number of double, int, bool, string parameters*/

		double dataDbl(0.0);
        int dataInt(0);
		std::string dataStr("abc");
        if (systemVar::currentProc==0)
        {
            std::cout << "	Reading " << fileName <<"\n";
        }

        std::string FileLoc(direction + "/" + fileName);
		std::ifstream FileFlux(FileLoc.c_str());
		int indexDbl(0), indexInt(0), indexBool(0), indexStr(0);
		if (FileFlux)
		{
			std::string line, keyWord;
			while (std::getline(FileFlux, line))
			{
				std::istringstream line2str(line);
				std::vector<std::string> ptr;
				//Split <line2str> stringstream into array of words
				while ((line2str >> keyWord))
				{
					ptr.push_back(keyWord);
				}

                int numWrd = static_cast<int>(ptr.size());

				if (numWrd >= 2)
				{
					std::string str1(ptr[0]), str2(" "), str3(" "), str4(" "), str5(" ");
					std::istringstream strdata(ptr[1]);
					if (indexDbl<numParamDbl)
					{
						str2 = keyWordsDbl[indexDbl];
					}
					if (indexInt<numParamInt)
					{
						str3 = keyWordsInt[indexInt];
					}
					if (indexBool<numParamBool)
					{
						str4 = keyWordsBool[indexBool];
					}
					if (indexStr<numParamStr)
					{
						str5 = keyWordsStr[indexStr];
					}

					if (str1.compare(str2) == 0)  //double value
					{
						strdata >> dataDbl;
						outDbl[indexDbl]=dataDbl;
						indexDbl++;
					}
					else if ((str1.compare(str3) == 0))  //int value
					{
						strdata >> dataInt;
						outInt[indexInt] = dataInt;
						indexInt++;
					}
					else if ((str1.compare(str4) == 0))  //bool value
					{
						if ((ptr[1].compare("true") == 0) || (ptr[1].compare("yes") == 0))
						{
							outBool[indexBool] = true;
						}
						else if ((ptr[1].compare("false") == 0) || (ptr[1].compare("no") == 0))
						{
							outBool[indexBool] = false;
						}
						indexBool++;
					}
					else if ((str1.compare(str5) == 0))  //string value
					{
						outStr[indexStr] = ptr[1];
						indexStr++;
					}
				}
				if ((indexDbl>=numParamDbl) && (indexInt>=numParamInt) && (indexBool >= numParamBool) && (indexStr >= numParamStr))
				{
					break;
				}
			}
		}
		else
		{
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, FileLoc));
		}
	}

    void readNonScalar(std::string mode)
	{
		/*NOTES:
		Boundary conditions compatibility
		|U					|T					|p					|
		+-------------------+-------------------+-------------------+
		|1. inFlow			|1. inFlow			|1. inFlow			|
		|	Value u v w		|	Value T			|	Value p			|
		+-------------------+-------------------+-------------------+
		|2. noSlip			|2. WallIsothermal	|2. zeroGradient	|
		|					|	Value T			|					|
		+-------------------+-------------------+-------------------+
		|2. noSlip			|3. WallAdiabatic	|2. zeroGradient	|
		+-------------------+-------------------+-------------------+
		|7.	symmetry		|7. symmetry		|7. symmetry		|
		+-------------------+-------------------+-------------------+
		|4. outFlow			|4. outFlow			|4. outFlow			|
		|	Value u v w		|	Value T			|	Value p			|
		+-------------------+-------------------+-------------------+
        U:
        + 3:
        movingWall
        velocity        u v w
		*/

        std::string fileName("U.txt"), tempStr(""), Loc;
        if (mode.compare("p")==0)
        {
            Loc = (systemVar::wD + "/CASES/" + systemVar::caseName + "/Processor" + std::to_string(systemVar::currentProc) + "/0");
        }
        else {
            Loc = (systemVar::wD + "/CASES/" + systemVar::caseName + "/0");
        }

        if (systemVar::currentProc==0)
        {
            std::cout << "	Reading " << fileName << "\n";
        }

        std::string FileLoc(Loc + "/" + fileName);
		std::ifstream FileFlux(FileLoc.c_str());
        int bcGrp(0);

		if (FileFlux)
		{
			std::string line, keyWord;
			while (std::getline(FileFlux, line))
			{
				std::istringstream line2str(line);
				std::vector<std::string> ptr;
				//Split <line2str> stringstream into array of words
				while ((line2str >> keyWord))
				{
					ptr.push_back(keyWord);
				}

                int numWrd = static_cast<int>(ptr.size());
				if (numWrd != 0)
				{
					std::string str1(ptr[0]);
					if (str1.compare("initialValue") == 0)  //initial data
					{
						std::istringstream str_u(ptr[1]);
						std::istringstream str_v(ptr[2]);
						std::istringstream str_w(ptr[3]);
						str_u >> iniValues::uIni;
						str_v >> iniValues::vIni;
						str_w >> iniValues::wIni;
					}
					else if (str1.compare("Type") == 0)  //Bc Type
					{
						std::string str0(ptr[1]);

						if (meshVar::BoundaryType[bcGrp - 1][1] == 1)  //WALL
						{
							if ((str0.compare("noSlip") == 0))  //Type noSlip
							{
								bcValues::UBcType[bcGrp - 1] = 2;
								bcValues::uBC[bcGrp - 1] = 0.0;
								bcValues::vBC[bcGrp - 1] = 0.0;
								bcValues::wBC[bcGrp - 1] = 0.0;
							}
							else if ((str0.compare("slip") == 0))  //Type slip  BO SUNG SAU
							{
								//Use Maxwell-Smoluchovsky boundary condition
								bcValues::UBcType[bcGrp - 1] = 5;
								bcValues::uBC[bcGrp - 1] = 0.0;
								bcValues::vBC[bcGrp - 1] = 0.0;
								bcValues::vBC[bcGrp - 1] = 0.0;

								std::getline(FileFlux, line);
								std::istringstream fixedUStream(line);
								fixedUStream >> tempStr >> bcValues::uWall[bcGrp - 1] >> bcValues::vWall[bcGrp - 1] >> bcValues::wWall[bcGrp - 1];
							}
                            else if ((str0.compare("movingWall") == 0))  //Type movingWall
							{
								bcValues::UBcType[bcGrp - 1] = 3;
								std::getline(FileFlux, line);
								std::istringstream fixedUStream(line);
								fixedUStream >> tempStr >> bcValues::uBC[bcGrp - 1] >> bcValues::vBC[bcGrp - 1] >> bcValues::wBC[bcGrp - 1];
							}
							else
							{
                                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "U", "wall"));
							}
						}
						else if (meshVar::BoundaryType[bcGrp - 1][1] == 2)  //PATCH
						{
							if ((str0.compare("inFlow") == 0))  //Type inletOutlet
							{
								bcValues::UBcType[bcGrp - 1] = 1;
								std::getline(FileFlux, line);
								std::istringstream fixedUStream(line);
								fixedUStream >> tempStr >> bcValues::uBC[bcGrp - 1] >> bcValues::vBC[bcGrp - 1] >> bcValues::wBC[bcGrp - 1];
							}
							else if ((str0.compare("outFlow") == 0))  //Type inletOutlet
							{
								bcValues::UBcType[bcGrp - 1] = 4;
								std::getline(FileFlux, line);
								std::istringstream fixedUStream(line);
								fixedUStream >> tempStr >> bcValues::uBC[bcGrp - 1] >> bcValues::vBC[bcGrp - 1] >> bcValues::wBC[bcGrp - 1];
							}
							else
							{
                                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "U", "patch"));
							}
						}
						else if (meshVar::BoundaryType[bcGrp - 1][1] == 3)  //SYMMETRY
						{
							if ((str0.compare("symmetry") == 0))  //Type symmetry
							{
								bcValues::UBcType[bcGrp - 1] = 7;
								bcValues::uBC[bcGrp - 1] = 0.0;
								bcValues::vBC[bcGrp - 1] = 0.0;
								bcValues::wBC[bcGrp - 1] = 0.0;
							}
							else
							{
                                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "U", "symmetry"));
							}
						}
                        else if (meshVar::BoundaryType[bcGrp - 1][1] == 4)  //MATCHED
                        {
                            if ((str0.compare("matched") == 0))  //Type matched
                            {
                                bcValues::UBcType[bcGrp - 1] = 10;
                                bcValues::uBC[bcGrp - 1] = 0.0;
                                bcValues::vBC[bcGrp - 1] = 0.0;
                                bcValues::wBC[bcGrp - 1] = 0.0;
                            }
                            else
                            {
                                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "U", "matched"));
                            }
                        }
					}
					else if ((str1.compare("Group") == 0))  //Group
					{
						std::istringstream str_bcGrp(ptr[1]);
						str_bcGrp >> bcGrp;
					}
				}
			}
		}
		else
		{
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, FileLoc));
		}
	}

    void readScalar(std::string fileName, std::string mode)
	{
		/*NOTES:
		Boundary conditions compatibility
		|U					|T					|p					|
		+-------------------+-------------------+-------------------+
		|1. inFlow			|1. inFlow			|1. inFlow			|
		|	value u v w		|	value T			|	value p			|
		+-------------------+-------------------+-------------------+
		|2. noSlip			|2. WallIsothermal	|2. zeroGradient	|
		|					|	Value T			|					|
		+-------------------+-------------------+-------------------+
		|2. noSlip			|3. WallAdiabatic	|2. zeroGradient	|
		+-------------------+-------------------+-------------------+
		|4. outFlow			|4. outFlow			|4. outFlow			|
		|	value u v w		|	value T			|	value p			|
		+-------------------+-------------------+-------------------+
		|7.	symmetry		|7. symmetry		|7. symmetry		|
		+-------------------+-------------------+-------------------+
        |10.matched 		|10.matched 		|10.matched 		|
        +-------------------+-------------------+-------------------+
		*/

		fileName = fileName + ".txt";
		std::string tempStr("");
        std::string Loc;
        if (mode.compare("p")==0)
        {
            Loc=systemVar::wD + "/CASES/" + systemVar::caseName + "/Processor" + std::to_string(systemVar::currentProc) + "/0";
        }
        else {
            Loc=systemVar::wD + "/CASES/" + systemVar::caseName + "/0";
        }

        if (systemVar::currentProc==0)
        {
            std::cout << "	Reading " << fileName << "\n";
        }
        std::string FileLoc(Loc + "/" + fileName);
		std::ifstream FileFlux(FileLoc.c_str());
		
        int bcGrp(0);

		if (FileFlux)
		{
			std::string line, keyWord;
			if (fileName.compare("p.txt") == 0)
			{
				while (std::getline(FileFlux, line))
				{
					std::istringstream line2str(line);
					std::vector<std::string> ptr;
					//Split <line2str> stringstream into array of words
					while ((line2str >> keyWord))
					{
						ptr.push_back(keyWord);
					}

                    int numWrd = static_cast<int>(ptr.size());
					if (numWrd != 0)
					{
						std::string str1(ptr[0]);
						if (str1.compare("initialValue") == 0)  //initial data
						{
							std::istringstream str_val(ptr[1]);
							str_val >> iniValues::pIni;;
						}
						else if (str1.compare("Type") == 0)  //Bc Type
						{
							std::string str0(ptr[1]);

							if (meshVar::BoundaryType[bcGrp - 1][1] == 1)  //WALL
							{
								if ((str0.compare("zeroGradient") == 0))
								{
									bcValues::pBcType[bcGrp - 1] = 2;
									bcValues::pBC[bcGrp - 1] = iniValues::pIni;
								}
								else
								{
                                    message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "p", "wall"));
								}
							}
							else if (meshVar::BoundaryType[bcGrp - 1][1] == 2)  //PATCH
							{
								if ((str0.compare("inFlow") == 0))  //Type inletOutlet
								{
									bcValues::pBcType[bcGrp - 1] = 1;
									std::getline(FileFlux, line);
									std::istringstream Stream(line);
									Stream >> tempStr >> bcValues::pBC[bcGrp - 1];
								}
								else if ((str0.compare("zeroGradient") == 0))
								{
									bcValues::pBcType[bcGrp - 1] = 2;
									bcValues::pBC[bcGrp - 1] = iniValues::pIni;
								}
								else if ((str0.compare("outFlow") == 0))  //Type inletOutlet
								{
									bcValues::pBcType[bcGrp - 1] = 4;
									std::getline(FileFlux, line);
									std::istringstream Stream(line);
									Stream >> tempStr >> bcValues::pBC[bcGrp - 1];
								}
								else
								{
                                    message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "p", "patch"));
								}
							}
							else if (meshVar::BoundaryType[bcGrp - 1][1] == 3)
							{
								if ((str0.compare("symmetry") == 0))  //Type symmetry
								{
									bcValues::pBcType[bcGrp - 1] = 7;
									bcValues::pBC[bcGrp - 1] = iniValues::pIni;
								}
								else
								{
                                    message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "p", "symmetry"));
								}
							}
                            else if (meshVar::BoundaryType[bcGrp - 1][1] == 4)
                            {
                                if ((str0.compare("matched") == 0))  //Type matched
                                {
                                    bcValues::pBcType[bcGrp - 1] = 10;
                                    bcValues::pBC[bcGrp - 1] = iniValues::pIni;
                                }
                                else
                                {
                                    message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "p", "matched"));
                                }
                            }
						}
						else if ((str1.compare("Group") == 0))  //Group
						{
							std::istringstream str_bcGrp(ptr[1]);
							str_bcGrp >> bcGrp;
						}
					}
				}
			}
			else if (fileName.compare("T.txt") == 0)
			{
				while (std::getline(FileFlux, line))
				{
					std::istringstream line2str(line);
					std::vector<std::string> ptr;
					//Split <line2str> stringstream into array of words
					while ((line2str >> keyWord))
					{
						ptr.push_back(keyWord);
					}

                    int numWrd = static_cast<int>(ptr.size());
					if (numWrd != 0)
					{
						std::string str1(ptr[0]);
						if (str1.compare("initialValue") == 0)  //initial data
						{
							std::istringstream str_val(ptr[1]);
							str_val >> iniValues::TIni;;
						}
						else if (str1.compare("Type") == 0)  //Bc Type
						{
							std::string str0(ptr[1]);

							if (meshVar::BoundaryType[bcGrp - 1][1] == 1)  //WALL
							{
								if ((str0.compare("WallAdiabatic") == 0))
								{
									bcValues::TBcType[bcGrp - 1] = 3;
									bcValues::TBC[bcGrp - 1] = iniValues::TIni;
								}
								else if ((str0.compare("WallIsothermal") == 0))
								{
									bcValues::TBcType[bcGrp - 1] = 2;
									std::getline(FileFlux, line);
									std::istringstream Stream(line);
									Stream >> tempStr >> bcValues::TBC[bcGrp - 1];
								}
								else if ((str0.compare("temperatureJump") == 0))  //Type temperatureJump, BO SUNG SAU
								{
									bcValues::TBcType[bcGrp - 1] = 6;
									bcValues::TBC[bcGrp - 1] = iniValues::TIni;
									std::getline(FileFlux, line);
									std::istringstream Stream(line);
									Stream >> tempStr >> bcValues::TWall[bcGrp - 1];
								}
								else
								{
                                    message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "T", "wall"));
								}
							}
							else if (meshVar::BoundaryType[bcGrp - 1][1] == 2)  //PATCH
							{
								if ((str0.compare("inFlow") == 0))
								{
									bcValues::TBcType[bcGrp - 1] = 1;
									std::getline(FileFlux, line);
									std::istringstream Stream(line);
									Stream >> tempStr >> bcValues::TBC[bcGrp - 1];
								}
								else if ((str0.compare("outFlow") == 0))
								{
									bcValues::TBcType[bcGrp - 1] = 4;
									std::getline(FileFlux, line);
									std::istringstream Stream(line);
									Stream >> tempStr >> bcValues::TBC[bcGrp - 1];
								}
								else
								{
                                    message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "T", "patch"));
								}
							}
							else if (meshVar::BoundaryType[bcGrp - 1][1] == 3)
							{
								if ((str0.compare("symmetry") == 0))  //Type symmetry
								{
									bcValues::TBcType[bcGrp - 1] = 7;
									bcValues::TBC[bcGrp - 1] = iniValues::TIni;
								}
								else
								{
                                    message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "T", "symmetry"));
								}
							}
                            else if (meshVar::BoundaryType[bcGrp - 1][1] == 4)
                            {
                                if ((str0.compare("matched") == 0))  //Type matched
                                {
                                    bcValues::TBcType[bcGrp - 1] = 10;
                                    bcValues::TBC[bcGrp - 1] = iniValues::TIni;
                                }
                                else
                                {
                                    message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "T", "matched"));
                                }
                            }
						}
						else if ((str1.compare("Group") == 0))  //Group
						{
							std::istringstream str_bcGrp(ptr[1]);
							str_bcGrp >> bcGrp;
						}
					}
				}
			}

		}
		else
		{
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, FileLoc));
		}

		for (int i = 0; i < meshVar::nBc; i++)
		{
			if ((bcValues::UBcType[i]!= bcValues::TBcType[i])&&(bcValues::UBcType[i]==5))
			{
                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::SlipBcCompatibleError(i+1));
			}
		}
	}

	void residualOutput(double rhoResGlobal, double rhouResGlobal, double rhovResGlobal, double rhoEResGlobal)
	{
		if (systemVar::iterCount==1)
		{
			systemVar::rhoResNorm = rhoResGlobal;
			systemVar::rhouResNorm = rhouResGlobal;
			systemVar::rhovResNorm = rhovResGlobal;
			systemVar::rhoEResNorm = rhoEResGlobal;
			rhoResGlobal = 1.0;
			rhouResGlobal = 1.0;
			rhovResGlobal = 1.0;
			rhoEResGlobal = 1.0;
		}
		else
		{
			rhoResGlobal /= systemVar::rhoResNorm;
			rhouResGlobal /= systemVar::rhouResNorm;
			rhovResGlobal /= systemVar::rhovResNorm;
			rhoEResGlobal /= systemVar::rhoEResNorm;
		}
		std::cout << "Time step: " << dt << std::endl;
		std::cout << "Residuals: ddt(rho)=" << rhoResGlobal << ", ddt(rhou)=" << rhouResGlobal << ", ddt(rhov)=" << rhovResGlobal << ", ddt(rhoE)=" << rhoEResGlobal << std::endl;
		std::cout << std::endl;
	}

    void saveCase(std::string mode)
	{
		std::string iter_str = std::to_string(systemVar::iterCount);
        std::string fileName("rho.txt"), Loc;
        Loc = systemVar::wD + "/CASES/" + systemVar::caseName;
        if (mode.compare("p")==0)
        {
            Loc=Loc+"/Processor"+std::to_string(systemVar::currentProc) + "/" + iter_str;
        }
        else {
            Loc=Loc+"/" +iter_str;
        }
        auxUlti::createFolder(Loc);

		/*Conservative variables*/
        std::string fileLoc(Loc + "/" + fileName);
		std::ofstream fileFluxRho(fileLoc.c_str());
		for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
		{
			for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
			{
				fileFluxRho << rho[nelem][iorder] << " ";
			}
			fileFluxRho << std::endl;
		}

		fileName = "rhou.txt";
        fileLoc = (Loc + "/" + fileName);
		std::ofstream fileFluxRhou(fileLoc.c_str());
		for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
		{
			for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
			{
				fileFluxRhou << rhou[nelem][iorder] << " ";
			}
			fileFluxRhou << std::endl;
		}

		fileName = "rhov.txt";
        fileLoc = (Loc + "/" + fileName);
		std::ofstream fileFluxRhov(fileLoc.c_str());
		for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
		{
			for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
			{
				fileFluxRhov << rhov[nelem][iorder] << " ";
			}
			fileFluxRhov << std::endl;
		}

		fileName = "rhoE.txt";
        fileLoc = (Loc + "/" + fileName);
		std::ofstream fileFluxRhoE(fileLoc.c_str());
		for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
		{
			for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
			{
				fileFluxRhoE << rhoE[nelem][iorder] << " ";
			}
			fileFluxRhoE << std::endl;
		}
		/*end of saving conservative variables*/

		/*Residual normalized coeffs*/
		fileName = "ResidualNormCoeffs.txt";
        fileLoc = (Loc + "/" + fileName);
		std::ofstream fileFluxResNorm(fileLoc.c_str());
		fileFluxResNorm << systemVar::rhoResNorm << " " << systemVar::rhouResNorm << " " << systemVar::rhovResNorm << " " << systemVar::rhoEResNorm << std::endl;

		/*Time informations*/
        Loc = systemVar::wD + "/CASES/" + systemVar::caseName;
		fileName = "time.txt";
        fileLoc = (Loc + "/" + fileName);
		std::ofstream fileFluxTime(fileLoc.c_str());
		fileFluxTime << systemVar::iterCount << std::endl;
	}

    void loadCase(std::string mode)
	{
        //Read file time.txt
        std::string fileName("time.txt"), Loc(systemVar::wD + "/CASES/" + systemVar::caseName);
        std::string fileLoc(Loc + "/" + fileName);
		std::ifstream FileFluxTime(fileLoc.c_str());
		if (FileFluxTime)
		{
			std::string line;
			std::getline(FileFluxTime, line);
			std::istringstream lineflux(line);
			lineflux >> systemVar::iterCount;
		}
		else
		{
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, Loc));
		}

		std::string iter_str = std::to_string(systemVar::iterCount);
        Loc = systemVar::wD + "/CASES/" + systemVar::caseName;
        if (mode.compare("p")==0)
        {
            Loc=Loc+"/Processor"+std::to_string(systemVar::currentProc) + "/" + iter_str;
        }
        else {
            Loc=Loc+"/" +iter_str;
        }

		//Read rho
		fileName = "rho.txt";
        fileLoc = (Loc + "/" + fileName);
		std::ifstream FileFluxRho(fileLoc.c_str());
		if (FileFluxRho)
		{
			int nelement(0);
			std::string line;
			while (std::getline(FileFluxRho, line))
			{
				std::istringstream lineflux(line);
				for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
				{
					lineflux >> rho[nelement][iorder];
				}
				nelement++;
			}
		}
		else
		{
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, Loc));
		}

		//Read rhou
		fileName = "rhou.txt";
        fileLoc = (Loc + "/" + fileName);
		std::ifstream FileFluxRhou(fileLoc.c_str());
		if (FileFluxRhou)
		{
			int nelement(0);
			std::string line;
			while (std::getline(FileFluxRhou, line))
			{
				std::istringstream lineflux(line);
				for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
				{
					lineflux >> rhou[nelement][iorder];
				}
				nelement++;
			}
		}
		else
		{
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, Loc));
		}

		//Read rhov
		fileName = "rhov.txt";
        fileLoc = (Loc + "/" + fileName);
		std::ifstream FileFluxRhov(fileLoc.c_str());
		if (FileFluxRhov)
		{
			int nelement(0);
			std::string line;
			while (std::getline(FileFluxRhov, line))
			{
				std::istringstream lineflux(line);
				for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
				{
					lineflux >> rhov[nelement][iorder];
				}
				nelement++;
			}
		}
		else
		{
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, Loc));
		}

		//Read rhoE
		fileName = "rhoE.txt";
        fileLoc = (Loc + "/" + fileName);
		std::ifstream FileFluxRhoE(fileLoc.c_str());
		if (FileFluxRhoE)
		{
			int nelement(0);
			std::string line;
			while (std::getline(FileFluxRhoE, line))
			{
				std::istringstream lineflux(line);
				for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
				{
					lineflux >> rhoE[nelement][iorder];
				}
				nelement++;
			}
		}
		else
		{
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, Loc));
		}

		//Read residual norm coeffs
		fileName = "ResidualNormCoeffs.txt";
        fileLoc = (Loc + "/" + fileName);
		std::ifstream FileFluxResNorm(fileLoc.c_str());
		if (FileFluxResNorm)
		{
			std::string line;
			std::getline(FileFluxResNorm, line);
			std::istringstream lineflux(line);
			lineflux >> systemVar::rhoResNorm >> systemVar::rhouResNorm >> systemVar::rhovResNorm >> systemVar::rhoEResNorm;
		}
		else
		{
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, Loc));
		}
	}

	void write2DDoubleArrayToFile(std::vector<std::vector<double>> &array, std::string loc, int numRow, int numCol)
	{
		std::ofstream Flux(loc.c_str());
		if (Flux)
		{
            for (int i = 0; i < numRow; i++)
			{
                for (int j = 0; j < numCol; j++)
				{
					Flux << array[i][j] << " ";
				}
				Flux << "\n";
			}
		}
		else
		{
			std::cout<<"Can not open file "<<loc<<" to write.\n";
			std::cout << "DGSolver will exit after you hit return.\n";
			exit(EXIT_FAILURE);
		}
	}

	void write2DIntArrayToFile(std::vector<std::vector<int>> &array, std::string loc, int numRow, int numCol)
	{
		std::ofstream Flux(loc.c_str());
		if (Flux)
		{
            for (int i = 0; i < numRow; i++)
			{
                for (int j = 0; j < numCol; j++)
				{
					Flux << array[i][j] << " ";
				}
				Flux << "\n";
			}
		}
		else
		{
			std::cout<<"Can not open file "<<loc<<" to write.\n";
			std::cout << "DGSolver will exit after you hit return.\n";
			exit(EXIT_FAILURE);
		}
	}

    void openFileToAppend(std::string Loc, std::string content)
    {
        std::fstream fs;
          fs.open (Loc, std::fstream::app);

          fs << content;

          fs.close();
    }

	namespace importCase
	{
		void importResultsFromAnotherCase()
		{
            //SUA TAM THOI
            //import case parallel
            std::string sourceCaseName(" "), sourceLoc(" "), timeLoc(" "), fileName(" "), fileLoc(" "), processor("");

            if (systemVar::currentProc==0)
            {
                std::cout << "NOTE: source case must has the same mesh with current case.\n" << "Enter source name: ";
                std::cin>>sourceCaseName;
                //Send command to all processors
                for (int irank=1;irank<systemVar::totalProc;irank++) {
                    auxUlti::functionsOfParallelComputing::sendString(sourceCaseName,irank,2);
                }
            }
            else {
                sourceCaseName=auxUlti::functionsOfParallelComputing::receiveString(0,2);
            }

            sourceLoc = systemVar::wD + "/CASES/" + sourceCaseName;
            timeLoc = sourceLoc + "/time.txt";

            //read source DGOptions
            /*Read DGOptions*/
            std::string DGOptfileName("DGOptions.txt");
            std::string DGOptLoc(systemVar::wD + "/CASES/" + sourceCaseName + "/System");
            std::string DGOptkeyWordsDouble[1] = {}, DGOptkeyWordsInt[1] = {}, DGOptkeyWordsBool[1] = {}, DGOptkeyWordsStr[1] = {"runningMode"};
            double DGOptoutDB[1] = {};
            int DGOptoutInt[1] = {};
            bool DGOptoutBool[1] = {};
            std::string DGOptoutStr[1] = {};

            readDataFile(DGOptfileName, DGOptLoc, DGOptkeyWordsDouble, DGOptkeyWordsInt, DGOptkeyWordsBool, DGOptkeyWordsStr, DGOptoutDB, DGOptoutInt, DGOptoutBool, DGOptoutStr, 0, 0, 0, 1);
            if ((DGOptoutStr[0].compare("parallel") == 0)&& systemVar::parallelMode)
            {
                processor="Processor"+std::to_string(systemVar::currentProc);
            }
            else if (DGOptoutStr[0].compare("sequence") == 0 && !systemVar::parallelMode)
            {
                processor="";
            }
            else {
                std::cout << "runningMode option of source and destination cases are not the same.\n";
                std::cout << "DGSolver will exit after you hit return.\n";
                exit(EXIT_FAILURE);
            }

			//get time at source folder
			std::ifstream timeFlux(timeLoc.c_str());
			if (timeFlux)
			{
				int time(0);
				std::string line;
				std::getline(timeFlux, line);
				std::istringstream lineflux(line);
				lineflux >> time;

                if (systemVar::currentProc==0)
                {
                    std::cout << "Source time: " << std::to_string(time) <<std::endl;
                    std::cout << "Mapping fields:\n";
                }
				//read results from source folder
				fileName = "rho.txt";
                fileLoc = (sourceLoc + "/" + processor + "/" + std::to_string(time) + "/" + fileName);
				mappSourceToCurrent(fileLoc, rho);
                if (systemVar::currentProc==0)
                {
                    std::cout << "	rho\n";
                }

				fileName = "rhou.txt";
                fileLoc = (sourceLoc + "/" + processor + "/" + std::to_string(time) + "/" + fileName);
				mappSourceToCurrent(fileLoc, rhou);
                if (systemVar::currentProc==0)
                {
                    std::cout << "	rhou\n";
                }

				fileName = "rhov.txt";
                fileLoc = (sourceLoc + "/" + processor + "/" + std::to_string(time) + "/" + fileName);
				mappSourceToCurrent(fileLoc, rhov);
                if (systemVar::currentProc==0)
                {
                    std::cout << "	rhov\n";
                }

				fileName = "rhoE.txt";
                fileLoc = (sourceLoc + "/" + processor + "/" + std::to_string(time) + "/" + fileName);
				mappSourceToCurrent(fileLoc, rhoE);
                if (systemVar::currentProc==0)
                {
                    std::cout << "	rhoE\n";
                    std::cout << "DONE\n";
                }
			}
			else
			{
                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError("time.txt", sourceLoc));
			}
		}

		void mappSourceToCurrent(std::string fileLoc, std::vector<std::vector<double>> &currentResult)
		{
			std::ifstream FileFluxRho(fileLoc.c_str());
			std::vector<double> iniValue(mathVar::orderElem + 1, 0.0);
			if (FileFluxRho)
			{
				int nelement(0);
				double value(0.0);
				std::string line;
				while (std::getline(FileFluxRho, line))
				{
					std::istringstream lineflux(line);
					lineflux >> value;
					iniValue = process::calcIniValues(value, nelement);
					for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
					{
						currentResult[nelement][iorder] = iniValue[iorder];
					}
					nelement++;
				}
			}
			else
			{
                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileLoc, fileLoc));
			}
		}
	}
}
