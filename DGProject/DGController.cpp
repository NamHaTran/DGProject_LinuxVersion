#include "DGController.h"
#include "DGIOLib.h"
#include "DGMessagesLib.h"
#include "DGMeshReaderLib.h"
#include "DGPostProcessLib.h"
#include "VarDeclaration.h"
#include "VarDeclaration.h"
#include "CommandCheck.h"
#include "DGProcLib.h"
#include "DGAuxUltilitiesLib.h"
#include <iostream>
#include <windows.h>

void Executer(std::string cmd)
{
	if (preProcessKey::checkUnvReader(cmd))
	{
		std::string UnvReaderLoc(systemVar::wD + "\\Ultilities\\MeshReader\\UnvReader.exe");
		auxUlti::openFileEXE(UnvReaderLoc);
	}
	else if (processKey::checkDGRun(cmd))
	{
		PreProcessing();
		Processing();
	}
	else if (postProcessKey::checkExit(cmd))
	{
		systemVar::endKey = true;
	}
	else if (preProcessKey::checkUnvHelper(cmd))
	{
		message::UnvReaderHelp();
	}
	else if (preProcessKey::checkBCsHelper(cmd))
	{
		message::BCsHelp();
	}
	else if (preProcessKey::reSubmit(cmd))
	{
		IO::getCase();
	}
	else
	{
		std::cout << "	ERROR: unknow command <" << cmd << ">!!!\n";
	}
}

void Processing()
{
	/*SET INITIAL VALUES*/
	process::setIniValues();
	//auxUlti::ConserToPri();

	//std::cout << meshVar::BoundaryType;
	std::cout << " \n" << "Simulation is started\n";

	//Calculate initial limiter coefficients
	process::limiter::limiter();

	while (process::checkRunningCond())
	{
		//SOLVE AUXILARY EQUATION
		process::auxEq::solveAuxEquation();

		//SOLVE NSF EQUATION
		process::NSFEq::solveNSFEquation();

		//UPDATE VARIABLES
		process::NSFEq::updateVariables();

		//APPLY LIMITER
		process::limiter::limiter();

		debugTool::checkPointValue(9);
		debugTool::checkPointValue(2338);
		debugTool::checkPointValue(2240);
		debugTool::checkPointValue(2254);
		debugTool::checkPointValue(2221);
		debugTool::checkPointValue(2231);
		debugTool::checkPointValue(2193);
		debugTool::checkPointValue(2212);
		//debugTool::checkPointValue(5538);
	}
}

void PreProcessing()
{
	/*LOAD MESH*/
	IO::loadMesh();

	/*LOAD CONSTANTS*/
	IO::loadConstants();

	/*LOAD p T U*/
	IO::loadpTU();
	//Check subsonic
	refValues::subsonic = auxUlti::checkSubSonic();

	/*PROCESS MESH*/
	MshReader::meshProcess();

	/*CALCULATE JACOBIAN, BASIS FUNCTION AND GAUSSIAN*/
	meshParam::GaussParam();
	meshParam::basisFcParam();
	meshParam::JacobianParam();

	/*CALCULATE CELL METRICS*/
	meshParam::calcCellMetrics();

	/*CALCULATE COORDINATES DERIVATIVES*/
	meshParam::derivCoordinates();

	/*RESIZE ARRAYS*/
	auxUlti::resizeDGArrays();

	auxUlti::mappingEdges();
	//debugTool::checkElemInfor(2029);
}

void PostProcessing()
{

}