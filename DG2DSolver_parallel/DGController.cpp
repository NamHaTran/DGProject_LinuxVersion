#include "DGController.h"
#include "DGIOLib.h"
#include "DGMessagesLib.h"
#include "DGMeshReaderLib.h"
#include "DGPostProcessLib.h"
#include "VarDeclaration.h"
#include "CommandCheck.h"
#include "DGProcLib.h"
#include "DGAuxUltilitiesLib.h"
#include <iostream>
#include "dynamicVarDeclaration.h"
#include "DGLimiterLib.h"

void Executer(std::string cmd)
{
	if (preProcessKey::checkUnvReader(cmd))
	{
        std::string UnvReaderLoc(systemVar::wD + "/Ultilities/MeshReader/UnvToDG");
		auxUlti::openFileEXE(UnvReaderLoc);
	}
	else if (processKey::checkDGRun(cmd))
	{
		if (systemVar::runPreProcess == false)
		{
			PreProcessing();
		}
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
		//PreProcessing();
	}
	else if (preProcessKey::mappResults(cmd))
	{
		PreProcessing();
		meshParam::calcStiffMatrixCoeffs();
		IO::importCase::importResultsFromAnotherCase();
		systemVar::initializedOrNot = true;
	}
    else if (preProcessKey::exportMeshToMetis(cmd))
    {
        IO::loadMesh();
        std::cout<<"Exporting DG mesh to Metis's format\n";
        MshExporter::exportMeshToMetis();
    }
    else if (preProcessKey::testMeshPartitionResult(cmd))
    {
        std::string key("n");
        std::cout<<"Do you want to read mesh? <y/n>: ";
        std::cin>>key;
        if ((key.compare("y") == 0) || (cmd.compare("Y") == 0))
        {
            IO::loadMesh();
        }
        std::cout<<"Creating techplot file\n";
        MshExporter::testMeshPartitionResult();
    }
    else if (preProcessKey::debug::checkElement(cmd))
    {
        int input(-1);
        std::cout << "Input element ID (ID is supplied by SALOME): ";
        std::cin >> input;
        std::cout << " \n";
        debugTool::checkElemInfor(input);
    }
    else if (preProcessKey::decomposeCase(cmd))
	{
        /*LOAD MESH*/
        IO::loadMesh();

        /*PROCESS MESH*/
        MshReader::meshProcess();

        decomposeMesh::decomposingMesh();
	}
	else
	{
		std::cout << "	ERROR: unknow command <" << cmd << ">!!!\n";
	}
}

void Processing()
{
	if (systemVar::loadSavedCase)
	{
		meshParam::calcStiffMatrixCoeffs();
		std::cout << "Loading case...\n" << std::endl;
		IO::loadCase();
	}
	else
	{
		/*SET INITIAL VALUES*/
		if (systemVar::initializedOrNot == false)
		{
			meshParam::calcStiffMatrixCoeffs();
			process::setIniValues();
		}
	}

	std::cout << " \n" << "Simulation is started\n";

	//APPLY LIMITER
	limiter::mathForLimiter::getNeighborElements();
	limitVal::numOfLimitCell = 0;
	limiter::limiter();

	int loadConstCount(0);
	while (process::checkRunningCond())
	{
		systemVar::iterCount++;
        std::cout << "Iteration " << systemVar::iterCount << std::endl
                  << "Total time = "<<runTime<<std::endl;

		//CALCULATE TIME STEP
		process::timeDiscretization::calcGlobalTimeStep();

		//SOLVE TIME MARCHING BY USING TVDRK3
		process::timeDiscretization::TVDRK3();

		//COMPUTE RESIDUALS
		process::timeDiscretization::globalErrorEstimate();
	
		systemVar::savingCout++;
		if (systemVar::savingCout == systemVar::wrtI)
		{
			std::cout << "Saving case...\n" << std::endl;
			IO::saveCase();
			std::cout << "Exporting data to Tecplot...\n" << std::endl;
			DG2Tecplot::exportCellCenteredData(systemVar::iterCount);
			systemVar::savingCout = 0;
		}

		loadConstCount++;
		if (loadConstCount == 10)
		{
			IO::loadConstants();
			loadConstCount = 0;
		}
	}
}

void PreProcessing()
{
	/*LOAD MESH*/
	IO::loadMesh();

	/*LOAD CONSTANTS*/
	IO::loadConstants();
	IO::loadLimiterSettings();

	/*LOAD p T U*/
	IO::loadpTU();

	//Check subsonic
    refValues::subsonic = auxUlti::checkSubSonic();
    //Check case's information
    message::checkCaseInformations();

	/*PROCESS MESH*/
	MshReader::meshProcess();

    /*RESIZE ARRAYS*/
    auxUlti::resizeDGArrays();

	/*CALCULATE JACOBIAN, BASIS FUNCTION AND GAUSSIAN*/
	meshParam::GaussParam();
	meshParam::basisFcParam();
	meshParam::JacobianParam();

	/*CALCULATE CELL METRICS*/
	meshParam::calcCellMetrics();
    meshParam::calcEdgeLength();

	/*CALCULATE COORDINATES DERIVATIVES*/
	meshParam::derivCoordinates();

	auxUlti::mappingEdges();

	systemVar::runPreProcess = true;
}

void PostProcessing()
{
	
}
