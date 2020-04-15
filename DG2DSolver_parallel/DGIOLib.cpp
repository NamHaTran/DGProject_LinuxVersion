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
            std::string line(" "), checkStr;
			double x, y, z;
            int iPoint(0);
            bool startToRead(false);

            //Doc so luong points (line dau tien cua file)
            meshVar::npoin=auxUlti::lookForDataOfKeyword(ptLoc,"NumberOfEntities");

            //Resize mang meshVar::Points
            meshVar::Points = auxUlti::resize2DArray(meshVar::npoin,3,0.0);
			
            int temp;
			while (std::getline(ptFlux, line))
			{
                //auxUlti::addRowTo2DDoubleArray(meshVar::Points, 3);
                if (startToRead==false)
                {
                    std::istringstream ptData(line);
                    ptData >> checkStr;
                    if (checkStr.compare("{") == 0)
                    {
                        startToRead=true;
                    }
                }
                else
                {
					if (iPoint<meshVar::npoin)
					{
						std::istringstream ptData(line);
                        ptData >> temp >> x >> y >> z;
						meshVar::Points[iPoint][0] = x;
						meshVar::Points[iPoint][1] = y;
						meshVar::Points[iPoint][2] = z;
						iPoint++;
                    }
					else 
					{
						break;
					}
                }
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
            int node1(0), node2(0), bcGrp(0), i1D(0);  //nelem1D is number of 1D elements
            std::string line(" "), checkStr(" ");
            bool startToRead(false);

            //Doc so luong element 1D
            meshVar::nelem1D=auxUlti::lookForDataOfKeyword(Elem1DLoc,"NumberOfEntities");

            //Resize mang meshVar::Elements1D
            meshVar::Elements1D = auxUlti::resize2DIntArray(meshVar::nelem1D,3,-1);

            int temp;
			while (std::getline(Elem1DFlux, line))
			{
                //auxUlti::addRowTo2DIntArray(meshVar::Elements1D, 3);
                if (startToRead==false)
                {
                    std::istringstream elData(line);
                    elData >> checkStr;
                    if (checkStr.compare("{") == 0)
                    {
                        startToRead=true;
                    }
                }
                else
                {
					if (i1D<meshVar::nelem1D)
					{
						std::istringstream elData(line);
                        elData >> temp>> node1 >> node2 >> bcGrp;
						meshVar::Elements1D[i1D][0] = node1 - 1;
						meshVar::Elements1D[i1D][1] = node2 - 1;
						meshVar::Elements1D[i1D][2] = bcGrp;
						i1D++;
					}
					else
					{
						break;
					}
                }
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
            int node1(0), node2(0), node3(0), node4(0), i2D(0);  //nelem2D is number of 2D elements
            std::string line(" "), checkStr;
            bool startToRead(false);

            //Doc so luong element 2D
            meshVar::nelem2D=auxUlti::lookForDataOfKeyword(Elem2DLoc,"NumberOfEntities");

            //Resize mang meshVar::Elements2D
            meshVar::Elements2D = auxUlti::resize2DIntArray(meshVar::nelem2D,4,-1);

            int temp;
			while (std::getline(Elem2DFlux, line))
			{
                //auxUlti::addRowTo2DIntArray(meshVar::Elements2D, 4);
                if (startToRead==false)
                {
                    std::istringstream elData(line);
                    elData >> checkStr;
                    if (checkStr.compare("{") == 0)
                    {
                        startToRead=true;
                    }
                }
                else
                {
					if (i2D<meshVar::nelem2D)
					{
						std::istringstream elData(line);
                        elData >> temp>> node1 >> node2 >> node3 >> node4;
						meshVar::Elements2D[i2D][0] = node1 - 1;
						meshVar::Elements2D[i2D][1] = node2 - 1;
						meshVar::Elements2D[i2D][2] = node3 - 1;
						meshVar::Elements2D[i2D][3] = node4 - 1;
						i2D++;
					}
					else
					{
						break;
					}
                }
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("Elements2D.txt", Elem2DLoc));
		}

        /*Load boundary conditions from boundaryPatch file*/
		std::ifstream bcFlux(bcLoc.c_str());
		if (bcFlux)
		{
            int boundIndex(0);
			std::string line(" "), keyWord(" ");
            //Doc so luong boundary
            meshVar::nBc=auxUlti::lookForDataOfKeyword(bcLoc,"NumberOfEntities");
			
			//bien nBc se tang them 1 (vi them 1 bien matched) khi readingMode parallel
            if (mode.compare("p")==0)
            {
                meshVar::nBc++;
            }
			
            meshVar::BoundaryType = auxUlti::resize2DIntArray(meshVar::nBc,3,0);

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
                //meshVar::nBc = boundIndex;
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("boundaryPatch.txt", bcLoc));
		}

        if (mode.compare("p")==0)
        {
            int meshConnectionLength(0);
            bool startToRead(false);
            std::string mshConnect = systemVar::pwd + "/Processor" + std::to_string(systemVar::currentProc) + "/Constant/Mesh/meshConnection.txt",
                    checkStr;
            std::ifstream mshConnectFlux(mshConnect.c_str());
            if (mshConnectFlux)
            {
                std::string line("");
                int counter = 0;
                meshConnectionLength=auxUlti::lookForDataOfKeyword(mshConnect,"NumberOfEntities");
                meshVar::meshConnection = auxUlti::resize2DIntArray(meshConnectionLength,3,0);

                while (std::getline(mshConnectFlux, line))
                {
                    //auxUlti::addRowTo2DIntArray(meshVar::meshConnection, 3);
                    if (startToRead==false)
                    {
                        std::istringstream Data(line);
                        Data >> checkStr;
                        if (checkStr.compare("{") == 0)
                        {
                            startToRead=true;
                        }
                    }
                    else
                    {
						if (counter<meshConnectionLength)
						{
							std::istringstream Data(line);
							Data >> meshVar::meshConnection[counter][0] >> meshVar::meshConnection[counter][1] >> meshVar::meshConnection[counter][2];
							counter++;
						}
						else
						{
							break;
						}
                    }
                }
            }
            else
            {
                message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("meshConnection.txt", mshConnect));
            }
        }

        //Read sendRecvOrder
        if (mode.compare("p")==0)
        {
            bool startToRead(false);
            std::string sendRecvOrderLoc = systemVar::pwd + "/Constant/Mesh/sendRecvOrder.txt",
                    checkStr;
            std::ifstream sendRecvOrderFlux(sendRecvOrderLoc.c_str());
            if (sendRecvOrderFlux)
            {
                std::string line("");
                int counter = 0;
                systemVar::sendRecvOrder_length=auxUlti::lookForDataOfKeyword(sendRecvOrderLoc,"NumberOfEntities");
                systemVar::sendRecvOrder = auxUlti::resize2DIntArray(systemVar::sendRecvOrder_length,2,0);

                while (std::getline(sendRecvOrderFlux, line))
                {
                    //auxUlti::addRowTo2DIntArray(meshVar::meshConnection, 3);
                    if (startToRead==false)
                    {
                        std::istringstream Data(line);
                        Data >> checkStr;
                        if (checkStr.compare("{") == 0)
                        {
                            startToRead=true;
                        }
                    }
                    else
                    {
                        if (counter<systemVar::sendRecvOrder_length)
                        {
                            std::istringstream Data(line);
                            Data >> systemVar::sendRecvOrder[counter][0] >> systemVar::sendRecvOrder[counter][1];
                            counter++;
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }
            else
            {
                message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("sendRecvOrder.txt", sendRecvOrderLoc));
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

        std::string inedelLoc, ineledLoc, inpoedLoc, normVectorLoc, MasterElemOfEdgeLoc, Elem1DLoc, ptLoc, Elem2DLoc, bcLoc;
		
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
			Fluxinedel << systemVar::headerFile << std::endl << "	inedel array\n" << "\n";
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
			Fluxineled << systemVar::headerFile << std::endl << "	ineled array\n" << "\n";
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
			Fluxinpoed << systemVar::headerFile << std::endl << "	inpoed array\n" << "\n";
            for (int i = 0; i < meshVar::inpoedCount; i++)
			{
                for (int j = 0; j < 5; j++)
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
			FluxnormVector << systemVar::headerFile << std::endl << "	normalVector array\n" << "\n";
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
			FluxMasterElemOfEdge << systemVar::headerFile << std::endl << "	MasterElemOfEdge array\n" << "\n";
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

    void readNumberOfCores()
    {
        /*Read totalProcess in DGOptions*/
        std::string DGOptfileName("DGOptions.txt");
        std::string DGOptLoc(systemVar::wD + "/CASES/" + systemVar::caseName + "/System");
        std::string DGOptkeyWordsDouble[1] = {" "}, DGOptkeyWordsInt[1] = {"totalProcess"}, DGOptkeyWordsBool[1] = {" "}, DGOptkeyWordsStr[1] = {" "};
        double DGOptoutDB[1] = {};
        int DGOptoutInt[1] = {};
        bool DGOptoutBool[1] = {};
        std::string DGOptoutStr[1] = {};

        readDataFile(DGOptfileName, DGOptLoc, DGOptkeyWordsDouble, DGOptkeyWordsInt, DGOptkeyWordsBool, DGOptkeyWordsStr, DGOptoutDB, DGOptoutInt, DGOptoutBool, DGOptoutStr, 0, 1, 0, 0);

        systemVar::totalProc = DGOptoutInt[0];
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
            MatLoc = systemVar::wD + "/CASES/" + systemVar::caseName + "/Constant";
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
        std::string DGSchemeskeyWordsDouble[1] = {}, DGSchemeskeyWordsInt[1] = {}, DGSchemeskeyWordsBool[1] = {}, DGSchemeskeyWordsStr[2] = {"diffusionTermScheme","solveTMethod"};
        double DGSchemesoutDB[1] = {};
        int DGSchemesoutInt[1] = {};
        bool DGSchemesoutBool[1] = {};
        std::string DGSchemesoutStr[2] = {};

        readDataFile(fileName, Loc, DGSchemeskeyWordsDouble, DGSchemeskeyWordsInt, DGSchemeskeyWordsBool, DGSchemeskeyWordsStr, DGSchemesoutDB, DGSchemesoutInt, DGSchemesoutBool, DGSchemesoutStr, 0, 0, 0, 2);

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

        if (DGSchemesoutStr[1].compare("implicitly")==0)
        {
            systemVar::solveTImplicitly=true;
        }
        else if (DGSchemesoutStr[1].compare("explicitly")==0)
        {
            systemVar::solveTImplicitly=false;
        }
        else {
            std::string str0("solveTMethod '"+DGSchemesoutStr[1]+"' is not available.");
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
							if ((limitVal::limiterName[ilimiter].compare("massDiffusion") == 0) || (limitVal::limiterName[ilimiter].compare("massdiffusion") == 0)) //PositivityPreserving settings
                            {
                                limitVal::massDiffusion = true;
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
        |5.	slip    		|6. temperatureJump	|2. zeroGradient    |
        |   sigmaU sigmaU   |   sigmaT sigmaT   |                   |
        |   v_wall u v w    |   T_wall T        |                   |
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
                                bcValues::uBCFixed[bcGrp - 1] = 0.0;
                                bcValues::vBCFixed[bcGrp - 1] = 0.0;
                                bcValues::wBCFixed[bcGrp - 1] = 0.0;
							}
							else if ((str0.compare("slip") == 0))  //Type slip
							{
								//Use Maxwell-Smoluchovsky boundary condition
                                bcValues::slipBCFlag=true;
								bcValues::UBcType[bcGrp - 1] = 5;
								std::getline(FileFlux, line);
                                std::istringstream fixedUStream1(line);
                                //Read sigmaU
                                fixedUStream1>>tempStr>>bcValues::sigmaU;
                                if ((tempStr.compare("sigmaU") != 0))
                                {
                                    std::cout<<"ERROR: Cannot find key word 'sigmaU' in file U/group "<<bcGrp<<", sigmaU is set to 1.0.\n";
                                    bcValues::sigmaU=1;
                                }

                                std::getline(FileFlux, line);
                                std::istringstream fixedUStream2(line);
                                //Read UWall
                                fixedUStream2 >> tempStr >> bcValues::uBCFixed[bcGrp - 1] >> bcValues::vBCFixed[bcGrp - 1] >> bcValues::wBCFixed[bcGrp - 1];
                                if ((tempStr.compare("v_wall") != 0))
                                {
                                    std::cout<<"ERROR: Cannot find key word 'v_wall' in file U/group "<<bcGrp<<", wall velocity is set to (0 0 0).\n";
                                    bcValues::uBCFixed[bcGrp - 1]=0;
                                    bcValues::vBCFixed[bcGrp - 1]=0;
                                }
							}
                            else if ((str0.compare("movingWall") == 0))  //Type movingWall
							{
								bcValues::UBcType[bcGrp - 1] = 3;
								std::getline(FileFlux, line);
								std::istringstream fixedUStream(line);
                                fixedUStream >> tempStr >> bcValues::uBCFixed[bcGrp - 1] >> bcValues::vBCFixed[bcGrp - 1] >> bcValues::wBCFixed[bcGrp - 1];
                                std::getline(FileFlux, line);
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
                                fixedUStream >> tempStr >> bcValues::uBCFixed[bcGrp - 1] >> bcValues::vBCFixed[bcGrp - 1] >> bcValues::wBCFixed[bcGrp - 1];
							}
							else if ((str0.compare("outFlow") == 0))  //Type inletOutlet
							{
								bcValues::UBcType[bcGrp - 1] = 4;
								std::getline(FileFlux, line);
								std::istringstream fixedUStream(line);
                                fixedUStream >> tempStr >> bcValues::uBCFixed[bcGrp - 1] >> bcValues::vBCFixed[bcGrp - 1] >> bcValues::wBCFixed[bcGrp - 1];
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
                                bcValues::uBCFixed[bcGrp - 1] = 0.0;
                                bcValues::vBCFixed[bcGrp - 1] = 0.0;
                                bcValues::wBCFixed[bcGrp - 1] = 0.0;
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
                                bcValues::uBCFixed[bcGrp - 1] = 0.0;
                                bcValues::vBCFixed[bcGrp - 1] = 0.0;
                                bcValues::wBCFixed[bcGrp - 1] = 0.0;
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
        |5.	slip    		|6. temperatureJump	|2. zeroGradient    |
        |   sigmaU sigmaU   |   sigmaT sigmaT   |                   |
        |   v_wall u v w    |   T_wall T        |                   |
        +-------------------+-------------------+-------------------+
        U:
        + 3:
        movingWall
        velocity        u v w
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
                            str_val >> iniValues::pIni;
						}
						else if (str1.compare("Type") == 0)  //Bc Type
						{
							std::string str0(ptr[1]);

							if (meshVar::BoundaryType[bcGrp - 1][1] == 1)  //WALL
							{
								if ((str0.compare("zeroGradient") == 0))
								{
									bcValues::pBcType[bcGrp - 1] = 2;
                                    bcValues::pBCFixed[bcGrp - 1] = iniValues::pIni;
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
                                    Stream >> tempStr >> bcValues::pBCFixed[bcGrp - 1];
								}
								else if ((str0.compare("zeroGradient") == 0))
								{
									bcValues::pBcType[bcGrp - 1] = 2;
                                    bcValues::pBCFixed[bcGrp - 1] = iniValues::pIni;
								}
								else if ((str0.compare("outFlow") == 0))  //Type inletOutlet
								{
									bcValues::pBcType[bcGrp - 1] = 4;
									std::getline(FileFlux, line);
									std::istringstream Stream(line);
                                    Stream >> tempStr >> bcValues::pBCFixed[bcGrp - 1];
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
                                    bcValues::pBCFixed[bcGrp - 1] = iniValues::pIni;
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
                                    bcValues::pBCFixed[bcGrp - 1] = iniValues::pIni;
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
                            str_val >> iniValues::TIni;
						}
						else if (str1.compare("Type") == 0)  //Bc Type
						{
							std::string str0(ptr[1]);

							if (meshVar::BoundaryType[bcGrp - 1][1] == 1)  //WALL
							{
								if ((str0.compare("WallAdiabatic") == 0))
								{
									bcValues::TBcType[bcGrp - 1] = 3;
                                    bcValues::TBCFixed[bcGrp - 1] = iniValues::TIni;
								}
								else if ((str0.compare("WallIsothermal") == 0))
								{
									bcValues::TBcType[bcGrp - 1] = 2;
									std::getline(FileFlux, line);
									std::istringstream Stream(line);
                                    Stream >> tempStr >> bcValues::TBCFixed[bcGrp - 1];
								}
                                else if ((str0.compare("temperatureJump") == 0))  //Type temperatureJump
								{
                                    bcValues::temperatureJump=true;
									bcValues::TBcType[bcGrp - 1] = 6;
                                    std::getline(FileFlux, line);
                                    std::istringstream Stream1(line);
                                    //Read sigmaT
                                    Stream1>> tempStr >> bcValues::sigmaT;
                                    if ((tempStr.compare("sigmaT") != 0))
                                    {
                                        std::cout<<"ERROR: Cannot find key word 'sigmaT' in file T/group "<<bcGrp<<", sigmaT is set to 1.0.\n";
                                        bcValues::sigmaT=1;
                                    }

									std::getline(FileFlux, line);
                                    std::istringstream Stream2(line);
                                    //Read TWall
                                    Stream2 >> tempStr >> bcValues::TBCFixed[bcGrp - 1];
                                    if ((tempStr.compare("T_wall") != 0))
                                    {
                                        std::cout<<"ERROR: Cannot find key word 'T_wall' in file T/group "<<bcGrp<<", This is an fatal error and DGSolver will exit.\n";
                                        exit(EXIT_FAILURE);
                                    }
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
                                    Stream >> tempStr >> bcValues::TBCFixed[bcGrp - 1];
								}
								else if ((str0.compare("outFlow") == 0))
								{
									bcValues::TBcType[bcGrp - 1] = 4;
									std::getline(FileFlux, line);
									std::istringstream Stream(line);
                                    Stream >> tempStr >> bcValues::TBCFixed[bcGrp - 1];
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
                                    bcValues::TBCFixed[bcGrp - 1] = iniValues::TIni;
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
                                    bcValues::TBCFixed[bcGrp - 1] = iniValues::TIni;
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
            if (bcValues::TBcType[i]!=6 && (bcValues::UBcType[i]==5))
			{
                message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::SlipBcCompatibleError(i+1));
			}
		}
	}

	void residualOutput(double rhoResGlobal, double rhouResGlobal, double rhovResGlobal, double rhoEResGlobal)
	{
        if (systemVar::iterCount<5)
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

        IO::writeResiduals(systemVar::iterCount,rhoResGlobal,rhouResGlobal,rhovResGlobal,rhoEResGlobal);
	}

    void readDecomposedMeshInfor()
    {
        std::string rankOf2DElemLoc(systemVar::pwd+"/Constant/Mesh/rankOf2DElem.txt"),
                Elem2DlocalIdWithRankLoc(systemVar::pwd+"/Constant/Mesh/Elem2DlocalIdWithRank.txt");

        std::ifstream File1Flux(rankOf2DElemLoc.c_str());
        int i(0);
        if (File1Flux)
        {
            File1Flux>>meshVar::rankOf2DElem[i];
            i++;
        }

        i=0;
        std::ifstream File2Flux(Elem2DlocalIdWithRankLoc.c_str());
        if (File2Flux)
        {
            File2Flux>>meshVar::Elem2DlocalIdWithRank[i];
            i++;
        }
    }

    //Ham write cac field roi rac (discreted) xuong file fileName tai folder Loc
    void writeDiscretedFields(std::string Loc, std::string fileName, double **Var)
    {
        std::string fileLoc(Loc + "/" + fileName);
        std::ofstream fileFlux(fileLoc.c_str());
        fileFlux<<systemVar::headerFile<<"numberOfEntities "<<meshVar::nelem2D<<"\n"<<"{\n";
        for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
        {
            for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
            {
                fileFlux << Var[nelem][iorder] << " ";
            }
            fileFlux << "\n";
        }
        fileFlux<<"}/n";
    }

    //Ham read cac field roi rac (discreted) xuong file fileName tai folder Loc
    void readDiscretedFields(std::string Loc, std::string fileName, double **Var)
    {
        std::string fileLoc(Loc + "/" + fileName), checkStr;
        std::ifstream FileFlux(fileLoc.c_str());

        bool startToRead(false);
        if (FileFlux)
        {
            int nelement(0);
            std::string line;
            while (std::getline(FileFlux, line))
            {
                //auxUlti::addRowTo2DIntArray(meshVar::Elements1D, 3);
                if (startToRead==false)
                {
                    std::istringstream lineflux(line);
                    lineflux >> checkStr;
                    if (checkStr.compare("{") == 0)
                    {
                        startToRead=true;
                    }
                }
                else
                {
                    if (nelement<meshVar::nelem2D)
                    {
                        std::istringstream lineflux(line);
                        for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                        {
                            lineflux >> Var[nelement][iorder];
                        }
                        nelement++;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }
        else
        {
            message::writeLog((systemVar::wD + "/CASES/" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, Loc));
        }
    }

    void saveCase()
	{
		std::string iter_str = std::to_string(systemVar::iterCount);
        std::string fileName, Loc, fileLoc;
        Loc = auxUlti::createTimeStepFolder(systemVar::iterCount,"case");

		/*Conservative variables*/
        IO::writeDiscretedFields(Loc,"rho.txt",rho);
        IO::writeDiscretedFields(Loc,"rhou.txt",rhou);
        IO::writeDiscretedFields(Loc,"rhov.txt",rhov);
        IO::writeDiscretedFields(Loc,"rhoE.txt",rhoE);
		/*end of saving conservative variables*/

        //Save TSurface, USurface
        fileName = "TSurface.txt";
        fileLoc = (Loc + "/" + fileName);
        std::ofstream fileFluxTSurface(fileLoc.c_str());
        for (int iedge = 0; iedge < meshVar::numBCEdges; iedge++)
        {
            double TVar(0.0);
            for (int nG = 0; nG <= mathVar::nGauss; nG++)
            {
                TVar+=SurfaceBCFields::TBc[iedge][nG];
            }
            fileFluxTSurface <<TVar/(mathVar::nGauss+1)<< "\n";
        }
        fileName = "USurface.txt";
        fileLoc = (Loc + "/" + fileName);
        std::ofstream fileFluxUSurface(fileLoc.c_str());
        for (int iedge = 0; iedge < meshVar::numBCEdges; iedge++)
        {
            double uVar(0.0), vVar(0.0);
            for (int nG = 0; nG <= mathVar::nGauss; nG++)
            {
                uVar+=SurfaceBCFields::uBc[iedge][nG];
                vVar+=SurfaceBCFields::vBc[iedge][nG];
            }
            fileFluxUSurface <<uVar/(mathVar::nGauss+1)<< " "<<vVar/(mathVar::nGauss+1)<<"\n";
        }

		/*Residual normalized coeffs*/
		fileName = "ResidualNormCoeffs.txt";
        fileLoc = (Loc + "/" + fileName);
		std::ofstream fileFluxResNorm(fileLoc.c_str());
        fileFluxResNorm << systemVar::rhoResNorm << " " << systemVar::rhouResNorm << " " << systemVar::rhovResNorm << " " << systemVar::rhoEResNorm << "\n";

		/*Time informations*/
        if (systemVar::currentProc==0)
        {
            Loc = systemVar::wD + "/CASES/" + systemVar::caseName;
            fileName = "time.txt";
            fileLoc = (Loc + "/" + fileName);
            std::ofstream fileFluxTime(fileLoc.c_str());
            fileFluxTime << systemVar::iterCount << " "<<runTime<< "\n";
        }
        MPI_Barrier(MPI_COMM_WORLD);
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
            lineflux >> systemVar::iterCount>>runTime;
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

        IO::readDiscretedFields(Loc,"rho.txt",rho);
        IO::readDiscretedFields(Loc,"rhou.txt",rhou);
        IO::readDiscretedFields(Loc,"rhov.txt",rhov);
        IO::readDiscretedFields(Loc,"rhoE.txt",rhoE);

        //Read file TSurface & USurface
        fileName = "TSurface.txt";
        fileLoc = (Loc + "/" + fileName);
        std::ifstream FileFluxTSurface(fileLoc.c_str());
        if (FileFluxTSurface)
        {
            int nEdge(0);
            double TVar(0.0);
            std::string line_T;

            //Gia tri doc vao cua Surface fields la gia tri trung binh tren toan edge
            while (std::getline(FileFluxTSurface, line_T))
            {
                std::istringstream line_Tflux(line_T);
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    line_Tflux>>TVar;
                    for (int nG=0; nG<=mathVar::nGauss; nG++)
                    {
                        SurfaceBCFields::TBc[nEdge][nG]=TVar;
                    }

                }
                nEdge++;
            }
        }
        else
        {
            if (systemVar::currentProc==0)
            {
                std::cout<<"Cannot find file TSurface.txt at folder "<<systemVar::iterCount<<", solver will load surface fields from folder 0\n";
            }
            //process::setIniSurfaceBCValues();
        }

        //Read file TSurface & USurface
        fileName = "USurface.txt";
        fileLoc = (Loc + "/" + fileName);
        std::ifstream FileFluxUSurface(fileLoc.c_str());
        if (FileFluxUSurface)
        {
            int nEdge(0);
            double uVar(0.0), vVar(0.0);
            std::string line_U;

            //Gia tri doc vao cua Surface fields la gia tri trung binh tren toan edge
            while (std::getline(FileFluxTSurface, line_U))
            {
                std::istringstream line_Uflux(line_U);
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    line_Uflux>>uVar>>vVar;
                    for (int nG=0; nG<=mathVar::nGauss; nG++)
                    {
                        SurfaceBCFields::uBc[nEdge][nG]=uVar;
                        SurfaceBCFields::vBc[nEdge][nG]=vVar;
                    }

                }
                nEdge++;
            }
        }
        else
        {
            if (systemVar::currentProc==0)
            {
                std::cout<<"Cannot find file USurface.txt at folder "<<systemVar::iterCount<<", solver will load surface fields from folder 0\n";
            }
            //process::setIniSurfaceBCValues();
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

    void write2DDoubleArrayToFile(std::vector<std::vector<double>> &array, std::string loc, std::string name, int numRow, int numCol)
	{
		std::ofstream Flux(loc.c_str());
		if (Flux)
		{
            Flux<<message::headerFile()<<name<<"\n"<<"NumberOfEntities "<<numRow<<"\n"<<"{\n";
            for (int i = 0; i < numRow; i++)
			{
                for (int j = 0; j < numCol; j++)
				{
					Flux << array[i][j] << " ";
				}
				Flux << "\n";
			}
            Flux<<"}\n";
		}
		else
		{
			std::cout<<"Can not open file "<<loc<<" to write.\n";
			std::cout << "DGSolver will exit after you hit return.\n";
			exit(EXIT_FAILURE);
		}
	}

    void writeResiduals(int iter, double rhoRes, double rhouRes, double rhovRes, double rhoERes)
    {
        IO::openFileToAppend((systemVar::pwd+"/ResidualOutput.txt"),std::to_string(iter)+" "
                             +std::to_string(rhoRes)+" "
                             +std::to_string(rhouRes)+" "
                             +std::to_string(rhovRes)+" "
                             +std::to_string(rhoERes)+"\n");
    }

    void write2DIntArrayToFile(std::vector<std::vector<int>> &array, std::string loc, std::string name, int numRow, int numCol)
	{
		std::ofstream Flux(loc.c_str());
		if (Flux)
		{
            Flux<<message::headerFile()<<name<<"\n"<<"NumberOfEntities "<<numRow<<"\n"<<"{\n";
            for (int i = 0; i < numRow; i++)
			{
                for (int j = 0; j < numCol; j++)
				{
					Flux << array[i][j] << " ";
				}
				Flux << "\n";
			}
            Flux<<"}\n";
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
    //Phan nay k can
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
                //mappSourceToCurrent(fileLoc, rho);
                if (systemVar::currentProc==0)
                {
                    std::cout << "	rho\n";
                }

				fileName = "rhou.txt";
                fileLoc = (sourceLoc + "/" + processor + "/" + std::to_string(time) + "/" + fileName);
                //mappSourceToCurrent(fileLoc, rhou);
                if (systemVar::currentProc==0)
                {
                    std::cout << "	rhou\n";
                }

				fileName = "rhov.txt";
                fileLoc = (sourceLoc + "/" + processor + "/" + std::to_string(time) + "/" + fileName);
                //mappSourceToCurrent(fileLoc, rhov);
                if (systemVar::currentProc==0)
                {
                    std::cout << "	rhov\n";
                }

				fileName = "rhoE.txt";
                fileLoc = (sourceLoc + "/" + processor + "/" + std::to_string(time) + "/" + fileName);
                //mappSourceToCurrent(fileLoc, rhoE);
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

    void write2DDoubleVectorToFile(std::string location, std::string fileName, std::vector<std::vector<double>> &vector)
    {
        std::ofstream FileFlux((location+"/"+fileName).c_str());
        if (FileFlux)
        {
            FileFlux<<message::headerFile()<<fileName<<"\n"<<"NumberOfEntities "<<vector.size()<<"{\n";
            for (int irow=0; irow<vector.size(); irow++)
            {
                for (int icol=0; icol<vector[irow].size(); icol++)
                {
                    FileFlux<<vector[irow][icol]<<" ";
                }
                FileFlux<<"\n";
            }
            FileFlux<<"}\n";
        }
        else
        {
            std::cout<<"Cannot open file at location "<<(location+"/"+fileName)<<" to write.\n";
        }
    }

    void write1DDoubleVectorToFile(std::string location, std::string fileName, double *vector, int length)
    {
        std::ofstream FileFlux((location+"/"+fileName).c_str());
        FileFlux<<message::headerFile()<<fileName<<"\n"<<"NumberOfEntities "<<length<<"{\n";
        if (FileFlux)
        {
            for (int irow=0; irow<length; irow++)
            {
                FileFlux<<vector[irow]<<"\n";
            }
            FileFlux<<"}\n";
        }
        else
        {
            std::cout<<"Cannot open file at location "<<(location+"/"+fileName)<<" to write.\n";
        }
    }

    void write1DIntVectorToFile(std::string location, std::string fileName, int *vector, int length)
    {
        std::ofstream FileFlux((location+"/"+fileName).c_str());
        if (FileFlux)
        {
            FileFlux<<message::headerFile()<<fileName<<"\n"<<"NumberOfEntities "<<length<<"\n"<<"{\n";
            for (int irow=0; irow<length; irow++)
            {
                FileFlux<<vector[irow]<<"\n";
            }
            FileFlux<<"}\n";
        }
        else
        {
            std::cout<<"Cannot open file at location "<<(location+"/"+fileName)<<" to write.\n";
        }
    }
}
