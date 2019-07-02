#include "DGPostProcessLib.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "DGMath.h"
#include "DGMessagesLib.h"
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>

namespace debugTool
{
	void checkElemSurPt(int ipoin)
	{
		std::cout << "In SALOME, point ";
		std::cout << ipoin << " is surrounded by elements \n";
		ipoin--;
		for (int ii = meshVar::esup2[ipoin]+1; ii <= meshVar::esup2[ipoin+1]; ii++)
		{
			std::cout << meshVar::esup1[ii]+ meshVar::nelem1D + 1 << std::endl;
		}
	}

	void checkPtsSurPt(int ipoin)
	{
		std::cout << "In SALOME, point ";
		std::cout << ipoin << " is surrounded by points \n";
		ipoin--;
		for (int ii = meshVar::psup2[ipoin] + 1; ii <= meshVar::psup2[ipoin + 1]; ii++)
		{
			std::cout << meshVar::psup1[ii] + 1 << std::endl;
		}
	}

	void checkElemsSurElem(int ielem)
	{
		std::cout << "In SALOME, element ";
		std::cout << ielem << " is surrounded by elements \n";
		ielem = ielem - 1 - meshVar::nelem1D;
		for (int i = 0; i < 4; i++)
		{
			if (meshVar::esuel[i][ielem]>=0)
			{
				std::cout << meshVar::esuel[i][ielem] + 1 + meshVar::nelem1D << std::endl;
			}
		}
	}

	void checkElemInfor(int elem)
	{
		std::cout << "In SALOME, element ";
		std::cout << elem << " has the following edges \n";
		elem = elem - 1 - meshVar::nelem1D;
		int edgeName(0), point1(0), point2(0), grp(0), bcType(0), edgeOrder(0), elemType(0);
		double nx(0.0), ny(0.0);
		bool master(true);
		elemType = auxUlti::checkType(elem);
		
		if (elem>0)
		{
			for (int iedge = 0; iedge < elemType; iedge++)
			{
				edgeName = meshVar::inedel[iedge][elem];
				point1 = meshVar::inpoed[0][edgeName] + 1;
				point2 = meshVar::inpoed[1][edgeName] + 1;
				grp = auxUlti::getGrpOfEdge(edgeName);
				bcType = auxUlti::getBCType(edgeName);
				nx = (auxUlti::getNormVectorComp(elem, edgeName, 1));
				ny = (auxUlti::getNormVectorComp(elem, edgeName, 2));
				edgeOrder = auxUlti::findEdgeOrder(elem, edgeName);
				master = auxUlti::checkMaster(elem, edgeName);
				std::cout << "+ Edge " << edgeName << " created by two points: " << point1 << " and " << point2 << std::endl
					<< "	Coordinates of normal vector are: " << nx << ", " << ny << std::endl
					<< "	Edge belongs to group " << grp << " and has bc type is " << bcType << std::endl
					<< "	Order of edge is " << edgeOrder << std::endl
					<< "	Edge is shared by two elements that are " << meshVar::ineled[0][edgeName] + 1 + meshVar::nelem1D << " and " << meshVar::ineled[1][edgeName] + 1 + meshVar::nelem1D << std::endl;
				if (master)
				{
					std::cout << "	Considering element is a master of edge\n";
				}
				else
				{
					std::cout << "	Considering element is not a master of edge\n";
				}
			}
		}
		else
		{
			std::cout << "	Considering element is not a polygonal element!\n";
		}

	}

	void checkPointValue(int element)
	{
		double rhoVal(math::pointValue(element, 0.0, 0.0, 1, 1)),
			uVal(math::pointValue(element, 0.0, 0.0, 2, 1)),
			vVal(math::pointValue(element, 0.0, 0.0, 3, 1)),
			//eVal(math::pointValue(element, 0.0, 0.0, 4, 1)),
			pVal(math::pointValue(element, 0.0, 0.0, 5, 1)),
			TVal(math::pointValue(element, 0.0, 0.0, 6, 1));
			//muVal(math::pointValue(element, 0.0, 0.0, 7, 1));
		std::cout << "- Element: " << element + meshVar::nelem1D + 1 << std::endl
			<< "- rho: " << rhoVal << std::endl
			<< "- u: " << uVal << std::endl
			<< "- v: " << vVal << std::endl
			//<< "- e: " << eVal << std::endl
			<< "- p: " << pVal << std::endl
			<< "- T: " << TVal << std::endl << std::endl;
			//<< "- mu: " << muVal << std::endl << std::endl;
	}

	/*
	void writeMinRho_MinRhoe()
	{
		//forDebugging
		std::string iter_str = std::to_string(systemVar::iterCount);
		std::string fileName_rho("rho_" + iter_str + ".txt"), fileName_rhoe("rhoe_" + iter_str + ".txt"),
            Loc(systemVar::wD + "/CASES/" + systemVar::caseName + "/forDebugging"), code;

        std::string fileLoc_rho(Loc + "/" + fileName_rho), fileLoc_rhoe(Loc + "/" + fileName_rhoe);
		std::ofstream fileFlux_rho(fileLoc_rho.c_str()), fileFlux_rhoe(fileLoc_rhoe.c_str());

		if (fileFlux_rho && fileFlux_rhoe)
		{
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				fileFlux_rho << debug::minRhoArr[i] << std::endl;
				fileFlux_rhoe << debug::minRhoeArr[i] << std::endl;
			}
		}
		else
		{
			std::cout << "Cannot open neither rho_.txt nor rhoe_.txt in folder <forDebugging>\n" << "DGSolver will exit after you hit return.\n";
			system("pause");
			exit(EXIT_FAILURE);
		}
	}
	*/
}

namespace DG2Tecplot
{
	std::vector<double> calcNodeValues(int valType)
	{
		//Only for primary variables
		std::vector<int>elemSurPt;
		std::vector<double>nodeValues(meshVar::npoin,0.0), U(4, 0.0);
		int numElemSurPt(0), elemId(-1), ptAtBCId(-2);
		double aP(0.0), bP(0.0), localValue(0.0);
		for (int ipoint = 0; ipoint < meshVar::npoin; ipoint++)
		{
			ptAtBCId = meshVar::markPointsAtBC[ipoint] - 1;
			if (ptAtBCId < 0)
			{
				localValue = 0.0;
				elemSurPt = auxUlti::postProcess::getElementsSurroundingPoint(ipoint);
				numElemSurPt = elemSurPt.size();
				for (int nelem = 0; nelem < numElemSurPt; nelem++)
				{
					elemId = elemSurPt[nelem];
					std::tie(aP, bP) = auxUlti::postProcess::findPointCoorInStandardSpace(ipoint, elemId);
					localValue += math::pointValue(elemId, aP, bP, valType, 1);
				}
				nodeValues[ipoint] = (localValue / numElemSurPt);
			}
			else
			{
				U = DG2Tecplot::calcNodeValuesAtBC(ptAtBCId);
				if (valType == 1)  //rho
				{
					nodeValues[ipoint] = U[0];
				}
				else if (valType == 2)  //u
				{
					double rhoVal(U[0]), rhouVal(U[1]);
					nodeValues[ipoint] = rhouVal / rhoVal;
				}
				else if (valType == 3)  //v
				{
					double rhoVal(U[0]), rhovVal(U[2]);
					nodeValues[ipoint] = rhovVal / rhoVal;
				}
				else if (valType == 4)  //e
				{
					double rhoVal(U[0]), rhouVal(U[1]), rhovVal(U[2]), rhoEVal(U[3]);
					nodeValues[ipoint] = material::Cv*math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
				}
				else if (valType == 5)  //p
				{
					double rhoVal(U[0]), rhouVal(U[1]), rhovVal(U[2]), rhoEVal(U[3]);
					double TVal(math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal));
					nodeValues[ipoint] = math::CalcP(TVal, rhoVal);
				}
				else if (valType == 6)  //T
				{
					double rhoVal(U[0]), rhouVal(U[1]), rhovVal(U[2]), rhoEVal(U[3]);
					nodeValues[ipoint] = math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
				}
				else if (valType == 7)  //mu
				{
					double rhoVal(U[0]), rhouVal(U[1]), rhovVal(U[2]), rhoEVal(U[3]);
					double TVal(math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal));
					if ((TVal < 0) && fabs(TVal) < 0.001)
					{
						TVal = fabs(TVal);
					}
					nodeValues[ipoint] = math::CalcVisCoef(TVal);
				}
			}
		}
		return nodeValues;
	}

	std::vector<double> calcNodeValuesAtBC(int ptAtBCId)
	{
		std::vector<double> U(4, 0.0);
		int BCEdgeId(0), BCElem(-1);
		double a(0.0), b(0.0);
		for (int iedge = 0; iedge < 2; iedge++)
		{
			BCEdgeId = SurfaceBCFields::BCPointsInfor[ptAtBCId][iedge];
			BCElem = meshVar::MasterElemOfEdge[BCEdgeId];
			int edgeGrp(auxUlti::getGrpOfEdge(BCEdgeId));
			int UType(bcValues::UBcType[edgeGrp - 1]), TType(bcValues::TBcType[edgeGrp - 1]), pType(bcValues::pBcType[edgeGrp - 1]), method(meshVar::BoundaryType[edgeGrp - 1][2]);

			std::tie(a, b) = auxUlti::postProcess::findPointCoorInStandardSpace(BCEdgeId, BCElem);

			if (UType == 1 && TType == 1 && pType == 1)
			{
				U = DG2Tecplot::calcNodeValueAtBCChildFuncs::patch::inFlow(BCElem, edgeGrp, a, b);
			}
			else if (UType == 4 && TType == 4 && pType == 4)
			{
				U = DG2Tecplot::calcNodeValueAtBCChildFuncs::patch::outFlow(BCElem, edgeGrp, a, b);
			}
			else if (UType == 2 && TType == 2 && pType == 2)
			{
				U = DG2Tecplot::calcNodeValueAtBCChildFuncs::wall::noSlipIsoThermal(BCElem, edgeGrp, a, b);
			}
			else if (UType == 2 && TType == 3 && pType == 2)
			{
				U = DG2Tecplot::calcNodeValueAtBCChildFuncs::wall::noSlipAdiabatic(BCElem, a, b);
			}
			else if (UType == 7 && TType == 7 && pType == 7)
			{
				U = DG2Tecplot::calcNodeValueAtBCChildFuncs::Symmetry(BCElem, BCEdgeId, a, b);
			}
			else
			{
				std::string errorStr = message::BcCompatibleError(edgeGrp);
				message::writeLog(systemVar::pwd, systemVar::caseName, errorStr);
			}
		}
		return U;
	}

	std::vector<double> calcCellCenteredValues(int valType)
	{
		std::vector<double>cellCenteredValues(meshVar::nelem2D, 0.0);
		double xC(0.0), yC(0.0), aC(0.0), bC(0.0);
		for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
			std::tie(xC, yC) = auxUlti::getCellCentroid(ielem);
			std::tie(aC, bC) = math::inverseMapping(ielem, xC, yC);
			cellCenteredValues[ielem] = math::pointValue(ielem, aC, bC, valType, 1);
		}
		return cellCenteredValues;
	}

	void exportNodeData(int iter)
	{
		std::vector<double>nodeRho(meshVar::npoin, 0.0), node_u(meshVar::npoin, 0.0), node_v(meshVar::npoin, 0.0), node_p(meshVar::npoin, 0.0), node_T(meshVar::npoin, 0.0), node_uMag(meshVar::npoin, 0.0);
		nodeRho = DG2Tecplot::calcNodeValues(1);
		node_u = DG2Tecplot::calcNodeValues(2);
		node_v = DG2Tecplot::calcNodeValues(3);
		node_p = DG2Tecplot::calcNodeValues(5);
		node_T = DG2Tecplot::calcNodeValues(6);
		for (int i = 0; i < meshVar::npoin; i++)
		{
			node_uMag[i] = sqrt(pow(node_u[i], 2) + pow(node_v[i], 2));
		}

		std::string iter_str = std::to_string(iter);
        std::string fileName(systemVar::caseName + "_node.dat"), Loc(systemVar::wD + "/CASES/" + systemVar::caseName + "/TecplotFile/" + iter_str), code;
        auxUlti::createFolder(Loc);
		
        std::string fileLoc(Loc + "/" + fileName);
		std::ofstream fileFlux(fileLoc.c_str());

		code = R"(
TITLE     = "DG2D to Tecplot"
VARIABLES = "X", "Y", "RHO", "U" ,"V", "VELOCITY_MAG", "P", "T"
ZONE T="ZONE 1"
ZONETYPE=FEQUADRILATERAL
DATAPACKING=BLOCK)";

		if (fileFlux)
		{
			fileFlux << code << std::endl << "NODES=" << std::to_string(meshVar::npoin) << ", " << "ELEMENTS=" << std::to_string(meshVar::nelem2D) << std::endl;
			//X
			int counter(0);
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << meshVar::Points[i][0] << " ";
			}
			fileFlux << std::endl;

			//Y
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << meshVar::Points[i][1] << " ";
			}
			fileFlux << std::endl;

			//RHO
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << nodeRho[i] << " ";
			}
			fileFlux << std::endl;

			//U
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_u[i] << " ";
			}
			fileFlux << std::endl;

			//U
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_v[i] << " ";
			}
			fileFlux << std::endl;

			//VELOCITY_MAG
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_uMag[i] << " ";
			}
			fileFlux << std::endl;

			//P
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_p[i] << " ";
			}
			fileFlux << std::endl;

			//T
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_T[i] << " ";
			}
			fileFlux << std::endl;

			//CONNECTIVITY
			for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
			{
				int elemType(auxUlti::checkType(ielem));
				for (int ipoin = 0; ipoin < 3; ipoin++)
				{
					fileFlux << meshVar::Elements2D[ielem][ipoin] + 1 << " ";
				}
				switch (elemType)
				{
				case 3:
				{
					fileFlux << meshVar::Elements2D[ielem][0] + 1 << std::endl;
					break;
				}
				case 4:
				{
					fileFlux << meshVar::Elements2D[ielem][3] + 1 << std::endl;
					break;
				}
				default:
					break;
				}
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError(systemVar::caseName + ".dat", fileLoc));
		}
	}

	void exportCellCenteredData(int iter)
	{
		std::vector<double>nodeRho(meshVar::nelem2D, 0.0), node_u(meshVar::nelem2D, 0.0), node_v(meshVar::nelem2D, 0.0), node_p(meshVar::nelem2D, 0.0), node_T(meshVar::nelem2D, 0.0), node_uMag(meshVar::nelem2D, 0.0);
		nodeRho = DG2Tecplot::calcCellCenteredValues(1);
		node_u = DG2Tecplot::calcCellCenteredValues(2);
		node_v = DG2Tecplot::calcCellCenteredValues(3);
		node_p = DG2Tecplot::calcCellCenteredValues(5);
		node_T = DG2Tecplot::calcCellCenteredValues(6);
		for (int i = 0; i < meshVar::nelem2D; i++)
		{
			node_uMag[i] = sqrt(pow(node_u[i], 2) + pow(node_v[i], 2));
		}

		std::string iter_str = std::to_string(iter);
        std::string fileName(systemVar::caseName + "cellCentered.dat"), Loc(systemVar::wD + "/CASES/" + systemVar::caseName + "/TecplotFile/" + iter_str), code;
        auxUlti::createFolder(Loc);

        std::string fileLoc(Loc + "/" + fileName);
		std::ofstream fileFlux(fileLoc.c_str());

		code = R"(
TITLE     = "DG2D to Tecplot"
VARIABLES = "X"
"Y"
"RHO"
"U"
"V"
"VELOCITY_MAG"
"P"
"T"
"THETA1"
"THETA2"
ZONE T="ZONE 1"
STRANDID=0
ZONETYPE=FEQuadrilateral
DATAPACKING=BLOCK
VARLOCATION=([3-10]=CELLCENTERED)
DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)
)";

		if (fileFlux)
		{
			fileFlux << code << std::endl << "Nodes=" << std::to_string(meshVar::npoin) << ", " << "Elements=" << std::to_string(meshVar::nelem2D) << std::endl;
			//X
			int counter(0);
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << meshVar::Points[i][0] << " ";
			}
			fileFlux << std::endl;

			//Y
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << meshVar::Points[i][1] << " ";
			}
			fileFlux << std::endl;

			//RHO
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << nodeRho[i] << " ";
			}
			fileFlux << std::endl;

			//U
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_u[i] << " ";
			}
			fileFlux << std::endl;

			//V
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_v[i] << " ";
			}
			fileFlux << std::endl;

			//VELOCITY_MAG
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_uMag[i] << " ";
			}
			fileFlux << std::endl;

			//P
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_p[i] << " ";
			}
			fileFlux << std::endl;

			//T
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_T[i] << " ";
			}
			fileFlux << std::endl;

			//Theta1
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << theta1Arr[i] << " ";
			}
			fileFlux << std::endl;

			//Theta2
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 5)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << theta2Arr[i] << " ";
			}
			fileFlux << std::endl;

			//CONNECTIVITY
			for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
			{
				int elemType(auxUlti::checkType(ielem));
				for (int ipoin = 0; ipoin < 3; ipoin++)
				{
					fileFlux << meshVar::Elements2D[ielem][ipoin] + 1 << " ";
				}
				switch (elemType)
				{
				case 3:
				{
					fileFlux << meshVar::Elements2D[ielem][0] + 1 << std::endl;
					break;
				}
				case 4:
				{
					fileFlux << meshVar::Elements2D[ielem][3] + 1 << std::endl;
					break;
				}
				default:
					break;
				}
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError(systemVar::caseName + ".dat", fileLoc));
		}
	}

	namespace calcNodeValueAtBCChildFuncs
	{
		namespace patch
		{
			std::vector<double> inFlow(int element, int edgeGrp, int a, int b)
			{
				std::vector<double> U(4, 0.0), UPlus(4, 0.0), UMinus(4, 0.0);

				for (int i = 0; i < 4; i++)
				{
					UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
				}

				//Apply weak Riemann infinite value
				UMinus[0] = (bcValues::pBC[edgeGrp - 1] / (material::R*bcValues::TBC[edgeGrp - 1]));
				UMinus[1] = UMinus[0] * bcValues::uBC[edgeGrp - 1];
				UMinus[2] = UMinus[0] * bcValues::vBC[edgeGrp - 1];
				UMinus[3] = UMinus[0] * (bcValues::TBC[edgeGrp - 1] * material::Cv + 0.5*(pow(bcValues::uBC[edgeGrp - 1], 2) + pow(bcValues::vBC[edgeGrp - 1], 2)));

				for (int i = 0; i < 4; i++)
				{
					U[i] = 0.5*(UPlus[i] + UMinus[i]) / 2.0;
				}
				return U;
			}

			std::vector<double> outFlow(int element, int edgeGrp, int a, int b)
			{
				std::vector<double> U(4, 0.0), UPlus(4, 0.0), UMinus(4, 0.0);

				for (int i = 0; i < 4; i++)
				{
					UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
				}

				double uPlus(UPlus[1] / UPlus[0]), vPlus(UPlus[2] / UPlus[0]);

				//Apply PNR (2), R (1)
				int implementation(1);
				switch (implementation)
				{
				case 1: //R
				{
					if (refValues::subsonic)
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
						UMinus[3] = bcValues::pBC[edgeGrp - 1] / (material::gamma - 1) + 0.5*UPlus[0] * (pow(uPlus, 2) + pow(vPlus, 2));
					}
					else
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
						UMinus[3] = UPlus[3];
					}
				}
				break;
				case 2: //PNR
				{
					double TInternal(math::CalcTFromConsvVar(UPlus[0], UPlus[1], UPlus[2], UPlus[3]));
					double pInternal(0);
					pInternal = UPlus[0] * material::R*TInternal;
					if (refValues::subsonic)
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
						UMinus[3] = (2 * bcValues::pBC[edgeGrp - 1] - pInternal) / (material::gamma - 1) + 0.5*UPlus[0] * (pow(uPlus, 2) + pow(vPlus, 2));
					}
					else
					{
						UMinus[0] = UPlus[0];
						UMinus[1] = UPlus[1];
						UMinus[2] = UPlus[2];
						UMinus[3] = UPlus[3];
					}
				}
				break;
				default:
					break;
				}

				for (int i = 0; i < 4; i++)
				{
					U[i] = 0.5*(UPlus[i] + UMinus[i]) / 2.0;
				}
				return U;
			}
		}

		namespace wall
		{
			std::vector <double> noSlipIsoThermal(int element, int edgeGrp, int a, int b)
			{
				std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
				std::vector<double> U(4, 0.0), UMinus(4, 0.0), UPlus(4, 0.0);

				//A1 approach
				for (int i = 0; i < 4; i++)
				{
					UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
				}

				UMinus[0] = UPlus[0];
				UMinus[1] = 0.0;
				UMinus[2] = 0.0;
				UMinus[3] = UPlus[0] * material::Cv*bcValues::TBC[edgeGrp - 1];

				for (int i = 0; i < 4; i++)
				{
					U[i] = 0.5*(UPlus[i] + UMinus[i]) / 2.0;
				}
				return U;
			}

			std::vector <double> noSlipAdiabatic(int element, int a, int b)
			{
				std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
				std::vector<double> U(4, 0.0), UMinus(4, 0.0), UPlus(4, 0.0);

				//A1 approach
				for (int i = 0; i < 4; i++)
				{
					UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
				}

				UMinus[0] = UPlus[0];
				UMinus[1] = 0.0;
				UMinus[2] = 0.0;
				UMinus[3] = UPlus[0] * material::Cv*math::CalcTFromConsvVar(UPlus[0], UPlus[1], UPlus[2], UPlus[3]);

				for (int i = 0; i < 4; i++)
				{
					U[i] = 0.5*(UPlus[i] + UMinus[i]) / 2.0;
				}
				return U;
			}
		}

		std::vector <double> Symmetry(int element, int edge, int a, int b)
		{
			std::vector<double> U(4, 0.0), UPlus(4, 0.0), UMinus(4, 0.0);
			double nx(auxUlti::getNormVectorComp(element, edge, 1)), ny(auxUlti::getNormVectorComp(element, edge, 2));
			for (int i = 0; i < 4; i++)
			{
				UPlus[i] = math::pointValue(element, a, b, i + 1, 2);
			}
			UMinus[0] = UPlus[0];
			UMinus[1] = UPlus[1] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*nx;
			UMinus[2] = UPlus[2] - 2 * (UPlus[1] * nx + UPlus[2] * ny)*ny;
			UMinus[3] = UPlus[3];
			for (int i = 0; i < 4; i++)
			{
				U[i] = 0.5*(UPlus[i] + UMinus[i]) / 2.0;
			}
			return U;
		}
	}
}
