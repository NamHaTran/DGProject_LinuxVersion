#include "nonEqmBCs_GenFuncs.h"
#include "nonEqmBCs_Vars.h"
#include "DGAuxUltilitiesLib.h"
#include "VarDeclaration.h"
#include "DGIOLib.h"
#include "dynamicVarDeclaration.h"
#include <iostream>
#include <sstream>

#include "../../bcVariables.h"

//Time-varying boundary conditions
#include "MaxwellSlip/u_MaxwellSlip.h"
#include "SmoluchowskyTJump/T_SmoluchowskyTJump.h"

namespace nonEquilibriumBCs {
    void updateBCs()
    {
        /*
         * Ham update gia tri tren cac field cua surfaceBCFields, chi dung cho cac BC bien thien theo thoi gian.
         * Hien tai, ham su dung cho temperatureJump va slip conditions
        */
        if (auxUlti::checkTimeVaryingBCAvailable() && flowProperties::viscous)
        {
            if (systemVar::currentProc==0)
            {
                std::cout<<"Updating time varying nonequilibrium BCs.\n";
            }

            int globleEdge(0);
            for (int ilocalEdge=0; ilocalEdge<meshVar::numBCEdges; ilocalEdge++)
            {
                globleEdge=auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(ilocalEdge);
                int edgeGrp(auxUlti::getGrpOfEdge(globleEdge));
                int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]);

                for (int nG=0; nG<=mathVar::nGauss; nG++)
                {
                    //Update Temperature
                    if (TType == BCVars::temperatureBCId::SmoluchowskyTJump)
                    {
                        //SmoluchowskyTJump::calcTJump_DGTypeExplicit(globleEdge,edgeGrp,nG);
                        SmoluchowskyTJump::calcTJump_FDMTypeImplicit(globleEdge,edgeGrp,nG);
                    }

                    //Update Velocity
                    if (UType == BCVars::velocityBCId::MaxwellSlip)
                    {
                        //MaxwellSlip::calcUSlip_DGTypeExplicit(globleEdge,edgeGrp,nG);
                        MaxwellSlip::calcUSlip_FDMTypeImplicit(globleEdge,edgeGrp,nG);
                    }
                }
            }
        }
    }

    void resizeSurfaceFields()
    {
        if (auxUlti::checkTimeVaryingBCAvailable())
        {
            nonEqmSurfaceField::TBc = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss+1,iniValues::TIni);
            nonEqmSurfaceField::uBc = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss+1,iniValues::uIni);
            nonEqmSurfaceField::vBc = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss+1,iniValues::vIni);
            nonEqmSurfaceField::pBc = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss+1,iniValues::pIni);

            meshVar::GaussPtsOnBCEdge_x = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss+1,0.0);
            meshVar::GaussPtsOnBCEdge_y = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss+1,0.0);
            meshVar::GaussPtsOnBCEdge_unitVector_x = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss+1,0.0);
            meshVar::GaussPtsOnBCEdge_unitVector_y = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss+1,0.0);
            meshVar::distanceFromGaussPtsToCentroid = auxUlti::resize2DArray(meshVar::numBCEdges,mathVar::nGauss+1,0.0);
        }
    }
}

namespace nonEquilibriumBCs_IO {
    void readSurfaceValues(std::string Loc)
    {
        if (auxUlti::checkTimeVaryingBCAvailable())
        {
            std::string fileLoc(""), fileName("");

            //Read file TSurface & USurface
            fileName = "TSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            nonEqmSurfaceField::TBc=nonEquilibriumBCs_IO::readFromFile(fileLoc,fileName,false);
            fileName = "uSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            nonEqmSurfaceField::uBc=nonEquilibriumBCs_IO::readFromFile(fileLoc,fileName,false);
            fileName = "vSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            nonEqmSurfaceField::vBc=nonEquilibriumBCs_IO::readFromFile(fileLoc,fileName,false);
        }
    }

    void writeSurfaceValues(std::string Loc)
    {
        if (auxUlti::checkTimeVaryingBCAvailable())
        {
            std::string fileLoc(""), fileName("");

            //Write file TSurface & USurface
            fileName = "TSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            nonEquilibriumBCs_IO::writeToFile(fileLoc,fileName,nonEqmSurfaceField::TBc);
            fileName = "uSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            nonEquilibriumBCs_IO::writeToFile(fileLoc,fileName,nonEqmSurfaceField::uBc);
            fileName = "vSurface.txt";
            fileLoc = (Loc + "/" + fileName);
            nonEquilibriumBCs_IO::writeToFile(fileLoc,fileName,nonEqmSurfaceField::vBc);
        }
    }

    void writeToFile(std::string loc, std::string name, double **array)
    {
        std::ofstream Flux(loc.c_str());
        int numRow(meshVar::numBCEdges), numCol(mathVar::nGauss+1);
        if (Flux)
        {
            Flux<<message::headerFile()<<name<<"\n"<<"NumberOfRow "<<numRow<<"\n"<<"NumberOfCol "<<numCol<<"\n"<<"{\n";
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

    double** readFromFile(std::string loc, std::string fileName, bool exitWhenFileNotFound)
    {
        int row(meshVar::numBCEdges), col(1), expectNumCol(mathVar::nGauss+1);
        double**array;
        bool isColAvail(true);
        std::ifstream Flux(loc.c_str());
        if (Flux)
        {
            std::string line(" "), checkStr;
            int iElem(0);
            bool startToRead(false);
            //std::tie(row,std::ignore)=auxUlti::lookForDataOfKeyword(loc,"NumberOfRow");
            std::tie(col,isColAvail)=auxUlti::lookForDataOfKeyword(loc,"NumberOfCol");
            array=auxUlti::resize2DArray(row,expectNumCol,0.0);

            if (!isColAvail)
                col=1;

            //Check whether col # nGauss
            bool calcMeanVal(false);
            if (expectNumCol!=col)
                calcMeanVal=true;

            while (std::getline(Flux, line))
            {
                if (startToRead==false)
                {
                    std::istringstream data(line);
                    data >> checkStr;
                    if (checkStr.compare("{") == 0)
                    {
                        startToRead=true;
                    }
                }
                else
                {
                    if (iElem<row)
                    {
                        std::istringstream data(line);
                        if (!calcMeanVal)
                        {
                            for (int i=0; i<col; i++)
                            {
                                data >> array[iElem][i];
                            }
                        }
                        else
                        {
                            double number(0.0),meanVar(0.0),sum(0.0);
                            for (int i=0; i<col; i++)
                            {
                                data >> number;
                                sum+=number;
                            }
                            meanVar=sum/col;

                            for (int j=0; j<expectNumCol; j++)
                            {
                                array[iElem][j] = meanVar;
                            }
                        }
                        iElem++;
                    }
                    else
                    {
                        break;
                    }
                }
            }
            return array;
        }
        else
        {
            if (exitWhenFileNotFound)
            {
                message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError(fileName, loc));
            }
            return 0;
        }
    }

    void writeGaussPtsCoor(std::string loc, std::string name)
    {
        std::ofstream Flux(loc.c_str());
        int numRow(meshVar::numBCEdges), numCol(2);
        if (Flux)
        {
            Flux<<message::headerFile()<<name<<"\n"<<"NumberOfRow "<<numRow<<"\n"<<"NumberOfCol "<<numCol<<"\n"<<"{\n";
            for (int i = 0; i < numRow; i++)
            {
                Flux << meshVar::GaussPtsOnBCEdge_x[i] << " " << meshVar::GaussPtsOnBCEdge_y[i] << "\n";
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
}
