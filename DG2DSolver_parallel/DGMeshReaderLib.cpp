#include "DGMeshReaderLib.h"
#include "DGIOLib.h"
#include "ConstDeclaration.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "DGMath.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <tuple>  //Include this for returning multiple values in function
namespace MshReader
{
	/*Declare local variables of DGMeshReaderLib, these variables are used throughout library*/
	int ninpoed(-1), edgesOfPoint[20][pointsArrSize];

	void meshProcess()
	{
		/*Elements surrounding point*/
		EleSurPt();

		/*Points surrounding point*/
		PtsSurPt();

		/*Elements surrounding element*/
		ElemsSurElem();

		/*Edges's informations*/
		EdgesInfor();
		EdgesOfElem();

		/*Calculate normal vector of each face (edge)*/
		GetNormalVector();

		/*Get points at boundary (for postProcessing)*/
		getBoundaryPoints();

		/*Save mesh data*/
		IO::SaveMeshInfor();
	}

	void EleSurPt()
	{
		//Default element type is quad
		/*CREATE INPOEL MATRIX*/
		//Default element type is quadrature
		for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
            meshVar::inpoel[ielem][0] = meshVar::Elements2D[ielem][0];
            meshVar::inpoel[ielem][1] = meshVar::Elements2D[ielem][1];
            meshVar::inpoel[ielem][2] = meshVar::Elements2D[ielem][2];
            meshVar::inpoel[ielem][3] = meshVar::Elements2D[ielem][3];
			if (meshVar::Elements2D[ielem][3] < 0)
			{
                meshVar::inpoel[ielem][4] = 3; //Type triangle
			}
			else
			{
                meshVar::inpoel[ielem][4] = 4; //Type quadrature
			}
		}

		/*ELEMENTS SURROUNDING POINTS*/
		//Element pass 1
		int ipoi1(0);
		for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
			for (int inode = 0; inode < meshVar::nnode; inode++)
			{
                ipoi1 = meshVar::inpoel[ielem][inode] + 1;
				if (ipoi1 >= 0)
				{
					meshVar::esup2[ipoi1] = meshVar::esup2[ipoi1] + 1;
				}
			}
		}
		//Storage pass 1
		for (int ipoin = 1; ipoin < meshVar::npoin + 1; ipoin++)
		{
			meshVar::esup2[ipoin] = meshVar::esup2[ipoin] + meshVar::esup2[ipoin - 1];
		}
		//Element pass 2
		int ipoin(0), istor(0);
		for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
			for (int inode = 0; inode < meshVar::nnode; inode++)
			{
                ipoin = meshVar::inpoel[ielem][inode];
				if (ipoin >= 0)
				{
					istor = meshVar::esup2[ipoin] + 1;
					meshVar::esup2[ipoin] = istor;
					meshVar::esup1[istor] = ielem;
				}
			}
		}
		//Storage pass 2
		for (int ipoin = meshVar::npoin + 1; ipoin >= 1; ipoin--)
		{
			meshVar::esup2[ipoin] = meshVar::esup2[ipoin - 1];
		}
		meshVar::esup2[0] = 0;
	}

	void PtsSurPt()
	{
		int lpoin[pointsArrSize] = {};
		int istor(0), ielem(0), jpoin(0);

		for (int r = 0; r < pointsArrSize; r++)
		{
			lpoin[r] = -2;  //abitraly non-zero value
		}

		for (int ipoin = 0; ipoin < meshVar::npoin; ipoin++)
		{
			for (int iesup = meshVar::esup2[ipoin] + 1; iesup <= meshVar::esup2[ipoin + 1]; iesup++)
			{
				ielem = meshVar::esup1[iesup];
				for (int inode = 0; inode < meshVar::nnode; inode++)
				{
                    jpoin = meshVar::inpoel[ielem][inode];
					if (jpoin >= 0)
					{
						if ((jpoin != ipoin)&(lpoin[jpoin] != ipoin))
						{
							istor++;
							meshVar::psup1[istor] = jpoin;
							lpoin[jpoin] = ipoin;
						}
					}
				}
			}
			meshVar::psup2[ipoin + 1] = istor;
		}
	}

	void ElemsSurElem()
	{
		int lpoin[pointsArrSize] = {};
		int nfael(4),
			lpofa[2][4] = {},
			lnofa[4] = {};
		int  nnofa(0), ipoin(0);
		int lhelp[2] = {};
		int jelem(0), nnofj(0), jpoin(0);

		int nfaelQuad(4), //A default element has 4 face
			lnofaQuad[4] = {}, //A default element has 4 face, each face has 2 node
			lpofaQuad[2][4] = {};

		//Set initial value for array
		for (int row = 0; row < 4; row++)
		{
			lnofaQuad[row] = 2;
		}

		lpofaQuad[0][0] = 3;
		lpofaQuad[0][1] = 0;
		lpofaQuad[0][2] = 1;
		lpofaQuad[0][3] = 2;

		lpofaQuad[1][0] = 0;
		lpofaQuad[1][1] = 1;
		lpofaQuad[1][2] = 2;
		lpofaQuad[1][3] = 3;
		//----------------------------
		int nfaelTri(3),
			lpofaTri[2][3] = {},
			lnofaTri[3] = {};
		for (int row = 0; row < 3; row++)
		{
			lnofaTri[row] = 2;
		}

		lpofaTri[0][0] = 2;
		lpofaTri[0][1] = 0;
		lpofaTri[0][2] = 1;

		lpofaTri[1][0] = 0;
		lpofaTri[1][1] = 1;
		lpofaTri[1][2] = 2;
		//----------------------------

		//Set initial value of meshVar::esuel
        for (int row = 0; row < elements2DArrSize; row++)
		{
            for (int column = 0; column < 4; column++)
			{
				meshVar::esuel[row][column] = -22;
			}
		}

		for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
            if (meshVar::inpoel[ielem][4] == 4)  //Type Quad
			{
				nfael = 4;
				for (int row = 0; row < 4; row++)
				{
					lnofa[row] = 2;
				}

				lpofa[0][0] = 3;
				lpofa[0][1] = 0;
				lpofa[0][2] = 1;
				lpofa[0][3] = 2;

				lpofa[1][0] = 0;
				lpofa[1][1] = 1;
				lpofa[1][2] = 2;
				lpofa[1][3] = 3;
			}
            else if (meshVar::inpoel[ielem][4] == 3)  //Type Tri
			{
				nfael = 3;
				for (int row = 0; row < 3; row++)
				{
					lnofa[row] = 2;
				}

				lpofa[0][0] = 2;
				lpofa[0][1] = 0;
				lpofa[0][2] = 1;

				lpofa[1][0] = 0;
				lpofa[1][1] = 1;
				lpofa[1][2] = 2;
			}

			for (int ifael = 0; ifael < nfael; ifael++)
			{
				nnofa = lnofa[ifael];
				for (int inofa = 0; inofa < nnofa; inofa++)
				{
                    lhelp[inofa] = meshVar::inpoel[ielem][lpofa[inofa][ifael]];
					lpoin[lhelp[inofa]] = 1;
				}
				ipoin = lhelp[0];
				for (int istor = meshVar::esup2[ipoin]+1; istor <= meshVar::esup2[ipoin+1]; istor++)
				{
					jelem = meshVar::esup1[istor];
                    if (meshVar::inpoel[jelem][4]==4)  //checked element is quad
					{
						if (jelem!=ielem)
						{
							for (int jfael = 0; jfael < nfaelQuad; jfael++)
							{
								nnofj = lnofaQuad[jfael];
								if (nnofj==nnofa)
								{
									int icoun(0);
									for (int jnofa = 0; jnofa < nnofa; jnofa++)
									{
                                        jpoin = meshVar::inpoel[jelem][lpofaQuad[jnofa][jfael]];
										if (jpoin>=0)
										{
											icoun = icoun + lpoin[jpoin];
										}
										if (icoun==nnofa)
										{
                                            meshVar::esuel[ielem][ifael] = jelem;
										}
									}
								}
							}
						}
					}
                    else if (meshVar::inpoel[jelem][4] == 3)  //checked element is quad
					{
						if (jelem!=ielem)
						{
							for (int jfael = 0; jfael < nfaelTri; jfael++)
							{
								nnofj = lnofaTri[jfael];
								if (nnofj==nnofa)
								{
									int icoun(0);
									for (int jnofa = 0; jnofa < nnofa; jnofa++)
									{
                                        jpoin= meshVar::inpoel[jelem][lpofaTri[jnofa][jfael]];
										if (jpoin>=0)
										{
											icoun = icoun + lpoin[jpoin];
										}
										if (icoun==nnofa)
										{
                                            meshVar::esuel[ielem][ifael] = jelem;
										}
									}
								}
							}
						}
					}
				}
				for (int inofa = 0; inofa < nnofa; inofa++)
				{
                    lhelp[inofa] = meshVar::inpoel[ielem][lpofa[inofa][ifael]];
					lpoin[lhelp[inofa]] = 0;
				}
			}
		}
	}

	void EdgesInfor()
	{
        int helpArrIndexI(0), helpArrIndexJ(0), jpoin(0), ipoinIndex(0), jpoinIndex(0);
        //int helpArray[20][pointsArrSize] = {};
        std::vector<std::vector<int>> helpArray(20,std::vector<int>(pointsArrSize,0));
		int iHelpArray[20];
        int flag(0), flag2(0);

		//Set initial values for helpArray array
		for (int row = 0; row < 20; row++)
		{
			for (int col = 0; col < pointsArrSize; col++)
			{
				helpArray[row][col] = -1;
				edgesOfPoint[row][col] = -1;
			}
		}

		for (int ipoin = 0; ipoin < meshVar::npoin; ipoin++)  //scan every ipoin, consider it as base point
		{
			helpArrIndexI = helpArray[19][ipoin];
			for (int row = 0; row < 20; row++)
			{
				iHelpArray[row] = helpArray[row][ipoin];
			}
			for (int isupoin = (meshVar::psup2[ipoin]+1); isupoin <= (meshVar::psup2[ipoin+1]); isupoin++)  //isupoin: i surrounding point
			{
				jpoin = meshVar::psup1[isupoin];  //scan all points which surrounding base point
				flag = CheckConnection(jpoin, iHelpArray, 20);  //check if jpoin has already been connected with ipoin
				helpArrIndexJ = helpArray[19][jpoin];
				if (flag==1)
				{
					flag2 = checkIndividualEdge(ipoin, jpoin);
					//update helpArrIndex
					helpArrIndexI++;
					helpArrIndexJ++;
					helpArray[19][ipoin]++;
					helpArray[19][jpoin]++;

					helpArray[helpArrIndexI][ipoin] = jpoin;
					helpArray[helpArrIndexJ][jpoin] = ipoin;

					if (flag2 == 0)
					{
						if (jpoin<ipoin)
						{
							ninpoed++;
                            meshVar::inpoed[ninpoed][0] = jpoin;
                            meshVar::inpoed[ninpoed][1] = ipoin;

                            ipoinIndex = edgesOfPoint[19][ipoin] + 1;
                            jpoinIndex = edgesOfPoint[19][jpoin] + 1;
							edgesOfPoint[ipoinIndex][ipoin] = ninpoed;
							edgesOfPoint[jpoinIndex][ipoin] = ninpoed;
                            edgesOfPoint[19][ipoin] = ipoinIndex;
                            edgesOfPoint[19][jpoin] = jpoinIndex;

							getBcGrpTp(ipoin, jpoin, ninpoed);
						}
						else
						{
							ninpoed++;
                            meshVar::inpoed[ninpoed][0] = ipoin;
                            meshVar::inpoed[ninpoed][1] = jpoin;

							ipoinIndex = edgesOfPoint[19][ipoin] + 1;
							jpoinIndex = edgesOfPoint[19][jpoin] + 1;
							edgesOfPoint[ipoinIndex][ipoin] = ninpoed;
							edgesOfPoint[jpoinIndex][jpoin] = ninpoed;
							edgesOfPoint[19][ipoin] = ipoinIndex;
							edgesOfPoint[19][jpoin] = jpoinIndex;

							getBcGrpTp(ipoin, jpoin, ninpoed);
						}
					}
				}
			}
		}
		meshVar::inpoedCount = ninpoed + 1;
        auxUlti::clear2DIntVector(helpArray);
	}

	void EdgesOfElem()
	{
		int numEdOfPoin(0), edgeName(0), pointTip(0), edgePointer(0);
		int helpArrayMarker[20] = {}, hlpArrIindex(-1), needReset(0);

		//Set initial values for helpArrayMarker
		for (int i = 0; i < 20; i++)
		{
			helpArrayMarker[i] = -1;
		}

		//Set initial values for ineled
		for (int c = 0; c <= ninpoed; c++)
		{
			for (int r = 0; r < 2; r++)
			{
                meshVar::ineled[c][r] = -meshVar::nelem1D;
			}
            meshVar::ineled[c][2] = -1;
		}
		int helpArray[4 * elements2DArrSize] = {}, index(-1), pointBase(0);

		for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
			//Reset value of helpArray to 0
			for (int ii = 0; ii < 20; ii++)
			{
				needReset = helpArrayMarker[ii];
				if (needReset>=0)
				{
					helpArray[needReset] = 0;
				}
			}

			//Reset value of helpArrayMarker and hlpArrIindex
			hlpArrIindex = -1;
			for (int i = 0; i < 20; i++)
			{
				helpArrayMarker[i] = -1;
			}

			//Reset index
			index = -1;
			for (int ipoin1 = 0; ipoin1 < 4; ipoin1++)  //scan all points of ielem
			{
				pointBase = meshVar::Elements2D[ielem][ipoin1];
				if (pointBase>=0)
				{
					numEdOfPoin = edgesOfPoint[19][pointBase];
					for (int iedge1 = 0; iedge1 <= numEdOfPoin; iedge1++)  //scan all edges content ipoin1
					{
						edgeName = edgesOfPoint[iedge1][pointBase];
						for (int ipoin2 = 0; ipoin2 < 2; ipoin2++)  //scan 2 points of edge
						{
                            pointTip = meshVar::inpoed[edgeName][ipoin2];
							if (pointTip!=pointBase)
							{
								for (int jpoin1 = 0; jpoin1 < 4; jpoin1++)
								{
									if (pointTip==meshVar::Elements2D[ielem][jpoin1])
									{
										if (helpArray[edgeName]!=1)
										{
											index++;
                                            meshVar::inedel[ielem][index] = edgeName;
                                            edgePointer = meshVar::ineled[edgeName][2];
											edgePointer++;
                                            meshVar::ineled[edgeName][edgePointer] = ielem;
                                            meshVar::ineled[edgeName][2] = edgePointer;
											helpArray[edgeName] = 1;

											//Use helpArrayMarker to know where in helpArray needs to be reset value
											hlpArrIindex++;
											helpArrayMarker[hlpArrIindex] = edgeName;
											break;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	void GetNormalVector()
	{
		int elem1(0), elem2(0), masterElem(0);
		int point1(0), point2(0), point1Indice(0), point2Indice(0), inforArr[4] = {};
		double normX(0.0), normY(0.0);
		for (int iedge = 0; iedge <= ninpoed; iedge++)
		{
			auxUlti::addRowTo2DDoubleArray(meshVar::normalVector, 2);

            elem1 = meshVar::ineled[iedge][0];
            elem2 = meshVar::ineled[iedge][1];

			//Find master element of iedge
			if (elem1>elem2)
			{
				masterElem = elem1;
			}
			else if (elem1<elem2)
			{
				masterElem = elem2;
			}
			meshVar::MasterElemOfEdge.push_back(masterElem);

            point1 = meshVar::inpoed[iedge][0];
            point2 = meshVar::inpoed[iedge][1];
			//Determine type of element
			if (meshVar::Elements2D[masterElem][3]>=0)  //Element type quad
			{
				for (int i = 0; i < 4; i++)
				{
					inforArr[i] = meshVar::Elements2D[masterElem][i];
				}
				point1Indice = findIndex(point1, inforArr, 4);
				point2Indice = findIndex(point2, inforArr, 4);
				std::tie(normX, normY) = calcNormVector(point1, point1Indice, point2, point2Indice, 4);  //Use this synstax to get result from tuble function
			}
			else
			{
				for (int i = 0; i < 3; i++)
				{
					inforArr[i] = meshVar::Elements2D[masterElem][i];
				}
				point1Indice = findIndex(point1, inforArr, 4);
				point2Indice = findIndex(point2, inforArr, 4);
				std::tie(normX, normY) = calcNormVector(point1, point1Indice, point2, point2Indice, 3);  //Use this synstax to get result from tuble function
			}
			meshVar::normalVector[iedge][0] = normX;
			meshVar::normalVector[iedge][1] = normY;
		}
	}

	/*Child functions------------------------------------------------*/
	
	/*Function supports for EdgesInfor(), it returns flag=0 if jpoin (input point) has already been connected with ipoin (helpArray of ipoin)*/
	int CheckConnection(int point, int helpArray[], int length)
	{
		int flag(0);
		for (int i = 0; i < length-1; i++)
		{
			if (point==helpArray[i])
			{
				flag = 0;
				break;
			}
			else
			{
				flag = 1;
			}
		}
		return flag;
	}

	/*Function supports for EdgesInfor(), it returns flag=1 if checking edge is individual edge*/
	int checkIndividualEdge(int rootPt, int tipPt)
	{
		int flag(0), index(0), elem(0);
		int elemSuPoint[20] = {};
		int pointArray[4] = {};
		int rootPosition(0), tipPosition(0);
		
		for (int ii = (meshVar::esup2[rootPt]+1); ii <= meshVar::esup2[rootPt+1]; ii++)
		{
			elemSuPoint[index] = meshVar::esup1[ii];
			index++;
		}
		for (int ielem = 0; ielem < index; ielem++)
		{
			elem = elemSuPoint[ielem];
			for (int row = 0; row < 4; row++)
			{
				pointArray[row] = meshVar::Elements2D[elem][row];
			}
			if (pointArray[3]>=0)  //Only quad elements have individual edge
			{
				rootPosition = findIndex(rootPt, pointArray, 4);
				tipPosition = findIndex(tipPt, pointArray, 4);
				if (tipPosition>=0)
				{
					if (((rootPosition==0)&(tipPosition==3||tipPosition==1)) || ((rootPosition == 1)&(tipPosition == 0 || tipPosition == 2)) || ((rootPosition == 2)&(tipPosition == 1 || tipPosition == 3)) || ((rootPosition == 3)&(tipPosition == 2 || tipPosition == 0)))
					{
						flag = 0;
					}
					else
					{
						flag = 1;
						break;
					}
				}
				else
				{
					flag = 0;
				}
			}
			else
			{
				flag = 0;
				break;
			}
		}
		return flag;
	}

	/*Function supports for EdgesInfor(), it returns index of checking number in input iarray*/
	int findIndex(int number, int iarray[], int size)
	{
		int index(0);
		for (int i = 0; i < size; i++)
		{
			if (number==iarray[i])
			{
				index = i;
				break;
			}
			else
			{
				index = -1;
			}
		}
		return index;
	}

	/*Function supports for EdgesInfor(), it gets BcGroup and Bc type of considering edge*/
	void getBcGrpTp(int ipoin, int jpoin, int ninpoed)
	{
		for (int iedge = 0; iedge < meshVar::nelem1D; iedge++)
		{
			if (((ipoin == meshVar::Elements1D[iedge][0])&(jpoin == meshVar::Elements1D[iedge][1])) || ((ipoin == meshVar::Elements1D[iedge][1])&(jpoin == meshVar::Elements1D[iedge][0])))
			{
				int BcGroup(meshVar::Elements1D[iedge][2]);  //Get group which edge belongs to
                meshVar::inpoed[ninpoed][2] = BcGroup;
				if (BcGroup != 0)  //Edge is not belong to internal group
				{
                    meshVar::inpoed[ninpoed][3] = meshVar::BoundaryType[BcGroup - 1][1];  //Get boundary type
					meshVar::adressOfBCVals.push_back(ninpoed);
					meshVar::numBCEdges++;
					break;
				}
			}
		}
	}

	std::tuple<double, double> calcNormVector(int point1, int point1Indice, int point2, int point2Indice, int type)  //Use std::tube<dataType, dataTye, ...> functionName to create function returns multiple variables
	{
		double deltaX(0.0), deltaY(0.0), vectorLength(0.0), normX(0.0), normY(0.0);
		if (std::abs(point1Indice-point2Indice)==1)
		{
			if (point2Indice>point1Indice)
			{
				deltaX = meshVar::Points[point2][0] - meshVar::Points[point1][0];
				deltaY = meshVar::Points[point2][1] - meshVar::Points[point1][1];
			}
			else
			{
				deltaX = meshVar::Points[point1][0] - meshVar::Points[point2][0];
				deltaY = meshVar::Points[point1][1] - meshVar::Points[point2][1];
			}
		}
		else
		{
			if ((point2Indice==0)&(point1Indice==type-1))
			{
				deltaX = meshVar::Points[point2][0] - meshVar::Points[point1][0];
				deltaY = meshVar::Points[point2][1] - meshVar::Points[point1][1];
			}
			else if ((point1Indice == 0)&(point2Indice == type-1))
			{
				deltaX = meshVar::Points[point1][0] - meshVar::Points[point2][0];
				deltaY = meshVar::Points[point1][1] - meshVar::Points[point2][1];
			}
		}
		vectorLength = sqrt(pow(deltaX, 2) + pow(deltaY, 2));
		normX = deltaY / vectorLength;
		normY = -deltaX / vectorLength;
		return std::make_tuple(normX, normY);
	}

	//Note: run this function AFTER mesh processing has been DONE
	void getBoundaryPoints()
	{
		meshVar::markPointsAtBC.resize(meshVar::npoin);
		std::vector<int> helpArr(meshVar::npoin, 0);
		int edgeId(-1), pt1(-1), pt2(-1), BCPtsId(0);
		for (int i = 0; i < meshVar::numBCEdges; i++)
		{
			edgeId = meshVar::adressOfBCVals[i];
            pt1 = meshVar::inpoed[edgeId][0];
            pt2 = meshVar::inpoed[edgeId][1];
			if (helpArr[pt1] == 0)
			{
				auxUlti::addRowTo2DIntArray(SurfaceBCFields::BCPointsInfor, 2);
				//SurfaceBCFields::BCPoints[numBCPts][0] = pt1;
				meshVar::markPointsAtBC[pt1] = BCPtsId + 1;
				SurfaceBCFields::BCPointsInfor[BCPtsId][0] = edgeId;
				helpArr[pt1] = 1;
				BCPtsId++;
			}
			else
			{
				SurfaceBCFields::BCPointsInfor[meshVar::markPointsAtBC[pt1] - 1][1] = edgeId;
			}

			if (helpArr[pt2] == 0)
			{
				auxUlti::addRowTo2DIntArray(SurfaceBCFields::BCPointsInfor, 2);
				//SurfaceBCFields::BCPoints[numBCPts][0] = pt2;
				meshVar::markPointsAtBC[pt2] = BCPtsId + 1;
				SurfaceBCFields::BCPointsInfor[BCPtsId][0] = edgeId;
				helpArr[pt2] = 1;
				BCPtsId++;
			}
			else
			{
				SurfaceBCFields::BCPointsInfor[meshVar::markPointsAtBC[pt2] - 1][1] = edgeId;
			}
		}
		auxUlti::clear1DIntVector(helpArr);
	}

	void sortPointsOfElements()
	{
		std::vector<int> sortedElement; 
		double xCG(0.0), yCG(0.0), xTranslate(0.0), yTranslate(0.0);
		int elemType(0), pointAId(-1);
		for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
			std::vector<double> filterdPointsCoor;
			std::vector<int> vectorPts, filterdPoints;
			elemType = auxUlti::checkType(ielem);
			std::vector<double> xCoor(elemType, 0.0),
				yCoor(elemType, 0.0);
			for (int i = 0; i < elemType; i++)
			{
				std::tie(xCoor[i], yCoor[i]) = auxUlti::getElemCornerCoord(ielem, i);
			}
			std::tie(xCG, yCG) = math::geometricOp::calcGeoCenter(xCoor, yCoor, elemType);

			//Filter 1:
			for (int i = 0; i < elemType; i++)
			{
				xTranslate = xCoor[i] - xCG;
				yTranslate = yCoor[i] - yCG;
				if (yTranslate <= 0.0)
				{
					filterdPoints.push_back(i);
					filterdPointsCoor.push_back(xTranslate);
				}
			}

			//Filter 2:
			if (filterdPoints.size() > 1)
			{
				double minXCoor(*std::min_element(filterdPointsCoor.begin(), filterdPointsCoor.end()));
                for (int id = 0; id < static_cast<int>(filterdPoints.size()); id++)
				{
                    if (fabs(minXCoor - filterdPointsCoor[id]) < 1e-10)
					{
						pointAId = filterdPoints[id];
						goto jumpHere;
					}
				}
			}
			else
			{
				pointAId = filterdPoints[0];
			}

		jumpHere:

			if (pointAId != 0)
			{
				sortedElement.push_back(ielem);
				for (int i = pointAId; i < elemType; i++)
				{
					vectorPts.push_back(meshVar::Elements2D[ielem][i]);
				}
				for (int i = 0; i < pointAId; i++)
				{
					vectorPts.push_back(meshVar::Elements2D[ielem][i]);
				}

				//Put point ids back to Elements2D array
				//std::cout << "Element " << ielem + meshVar::nelem1D + 1 << ". new order: \n";
				for (int i = 0; i < elemType; i++)
				{
					meshVar::Elements2D[ielem][i] = vectorPts[i];
					//std::cout << vectorPts[i] << " ";
				}
				//std::cout << std::endl;
			}
		}

		//Write information for debugging
		std::string  sortedElemLoc = systemVar::pwd + "\\Constant\\Mesh\\sortedElements.txt";
		std::ofstream Flux(sortedElemLoc.c_str());
		if (Flux)
		{
            for (int i = 0; i < static_cast<int>(sortedElement.size()); i++)
			{
				Flux << sortedElement[i] << std::endl;
			}
		}
		else
		{
			std::cout << "Error of writting files\n";
		}
	}
}
