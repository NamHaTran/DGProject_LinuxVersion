#include "DGMessagesLib.h"
#include "VarDeclaration.h"
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>

//include chrono and ctime to get time data
#include <chrono>
#include <ctime>

namespace message
{
	std::string headerFile()
	{
		std::string headerStr(" ");
		headerStr = R"(+------------------------------ CFD LOVES C++ ---------------------------+
|--------------------DISCONTINUOS GALERKIN METHOD SOLVER-----------------|
|                                  Author                                |
|   Nam Ha Tran.                                                         |
|   Ver 1.00                                                             |
+------------------------------------------------------------------------+
|   This program uses Discontinous Galerkin method to solve 2D problems  |
|   on structural and unstructural mesh.                                 |
+------------------------------------------------------------------------+
)";
		std::string Timestr(getTime());
		headerStr += "	Program ran at (d-m-y_h-m-s) " + Timestr + "\n";
		return headerStr;
	}

	std::string undfKeyW(std::string keyW, std::string location)
	{
		std::string str("Cannot find key word <" + keyW + "> at " + location);
		return str;
	}

	std::string opFError(std::string fileName, std::string location)
	{
		std::string str("Cannot open file <" + fileName + "> located at " + location);
		return str;
	}

	void writeLog(std::string location, std::string caseName, std::string str)
	{
		std::string strTime(getTime());
		std::string logFile(location + "\\log_" + caseName + ".txt");
		std::ofstream logFlux(logFile.c_str(), std::ios_base::app | std::ios_base::out);

		//Get time
		std::string crash_time(message::getTime());

		std::cout << "ERROR: " << str << std::endl;
		if (systemVar::wrtLog==true)
		{
			logFlux << "log was created at " << crash_time << std::endl << str << std::endl;
		}
		std::cout << "DGSolver will exit after you hit return.\n";
		exit(EXIT_FAILURE);
	}

	std::string getTime() //This function is referenced from internet
	{
		auto t = std::time(nullptr);
		auto tm = *std::localtime(&t);

		std::ostringstream oss;
		oss << std::put_time(&tm, "%d-%m-%Y_%H-%M-%S");
		auto str = oss.str();

		return str;
	}

	std::string undfBcType(std::string keyW, std::string fileName, std::string bcType)
	{
		std::string str("Boundary condition <" + keyW + "> of file " + fileName + " is not valid for type <" + bcType + ">");
		return str;
	}

	std::string headerpTU(std::string file)
	{
		std::string headerStr(" ");
		if (file.compare("p") == 0)
		{
			headerStr = R"(Pressure conditions (Pa)
initialValue				0
)";
		}
		else if (file.compare("T") == 0)
		{
			headerStr = R"(Temperature conditions (K)
initialValue				0
)";
		}
		else if (file.compare("U") == 0)
		{
			headerStr = R"(Velocity conditions (m/s)
initialValue			0 0 0
)";
		}
		return headerStr;
	}

	std::string SlipBcCompatibleError(int bcGrp)
	{
		std::string Str(" ");
		std::string bcGrpStr(std::to_string(bcGrp));
		Str = "Boundary condition of group " + bcGrpStr + " in file U is slip, but in file T it is not temperatureJump. This is considered as a fatal error.";
		return Str;
	}

	std::string nGaussOrderElemError()
	{
		std::string Str("orderElem should be less or equal to nGauss^2");
		return Str;
	}

	std::string BcCompatibleError(int edgeGrp)
	{
		std::string Str(" ");
		std::string bcGrpStr(std::to_string(edgeGrp));
		Str = "Boundary incompatibility error is detected, please check condition of boundary group " + bcGrpStr + " in folder /0/. This is considered as a fatal error.";
		return Str;
	}

	/*Help functions----------------------------------------------------------------*/

	void UnvReaderHelp()
	{
		std::string Str(" ");
		Str = R"(	UnvReader ultility's help:
To convert unv mesh format to DG2D readable format, do following task step by step:
	- Put mesh file into case folder.
	- Open DG2D solver, enter command <UnvToDG> to console.
	- Input mesh file name to console.
	- Edit boundary conditions, system settings and run simulation.
)";
		std::cout << Str;
	}

	void BCsHelp()
	{
		std::string Str(" ");
		Str = R"(DG Solver supports the following boundary conditions:
		+-------------------+-------------------+-------------------+
		|U					|T					|p					|
		+-------------------+-------------------+-------------------+
		|1. inOutFlow		|1. inOutFlow		|1. inOutFlow		|
		|	Value u v w		|	Value T			|	Value p			|
		+-------------------+-------------------+-------------------+
		|2. noSlip			|2. WallIsothermal	|2. zeroGradient	|
		|					|	Value T			|					|
		+-------------------+-------------------+-------------------+
		|2. noSlip			|3. WallAdiabatic	|2. zeroGradient	|
		+-------------------+-------------------+-------------------+
		|3.	fixedValue		|4. fixedValue		|3. fixedValue		|
		|	Value u v w		|	Value T			|	Value p			|
		+-------------------+-------------------+-------------------+
)";
		std::cout << Str;
	}
}

void exitDG(std::string str)
{
	std::cout << "ERROR: " << str << ". This is considered as fatal error." << std::endl;
	std::cout << "DGSolver will exit after you hit return.\n";
	system("pause");
	exit(EXIT_FAILURE);
}
