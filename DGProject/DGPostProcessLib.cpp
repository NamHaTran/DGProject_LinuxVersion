#include "DGPostProcessLib.h"
#include "VarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "DGMath.h"
#include <string>
#include <iostream>

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
}