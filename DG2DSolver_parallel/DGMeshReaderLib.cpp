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
#include <sstream>
#include <tuple>  //Include this for returning multiple values in function
#include "DGMessagesLib.h"
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

        if (!systemVar::runDecomposeCaseFnc)
        {
            /*Save mesh data*/
            IO::SaveMeshInfor();
        }
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
                    //Counting number of BC groups
                    if (BcGroup>meshVar::numBCGrp)
                    {
                        meshVar::numBCGrp=BcGroup;
                    }
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

namespace MshExporter
{
    void exportMeshToMetis()
    {
        std::string fileLoc(systemVar::wD + "/CASES/" + systemVar::caseName + "/DGMesh.mesh");
        std::ofstream fileFlux(fileLoc.c_str());
        fileFlux << meshVar::nelem2D<<" 1\n";
        for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
        {
            for (int i = 0; i <= 2; i++)
            {
                fileFlux << meshVar::Elements2D[nelem][i]+1 << " ";
            }
            if (meshVar::Elements2D[nelem][3]>=0)
            {
                fileFlux << meshVar::Elements2D[nelem][3]+1 << " ";
            }
            fileFlux << std::endl;
        }
    }

    void testMeshPartitionResult()
    {
        std::string fileName(systemVar::caseName + "ProcIdValues.dat"), Loc(systemVar::wD + "/CASES/" + systemVar::caseName + "/TecplotFile/PartitionedMesh"), code;
        auxUlti::createFolder(Loc);

        std::string fileLoc(Loc + "/" + fileName);
        std::ofstream fileFlux(fileLoc.c_str());

        code = R"(
TITLE     = "DG2D to Tecplot"
VARIABLES = "X"
"Y"
"PROC_RANK"
ZONE T="ZONE 1"
STRANDID=0
ZONETYPE=FEQuadrilateral
DATAPACKING=BLOCK
VARLOCATION=([3]=CELLCENTERED)
DT=(SINGLE SINGLE SINGLE)
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

            /*Declare loading locations*/
            std::string  procIdLoc = systemVar::pwd + "/DGMesh.mesh.epart";

            /*Load Points*/
            std::ifstream procFlux(procIdLoc.c_str());
            if (procFlux)
            {
                std::string line(" ");
                int nproc;
                meshVar::npoin = 0;
                while (std::getline(procFlux, line))
                {
                    std::istringstream procData(line);
                    procData >> nproc;
                    fileFlux<<nproc<<std::endl;
                }
            }
            else
            {
                std::cout<<"Can not open file DGMesh.mesh.epart.\n";
                std::cout << "DGSolver will exit after you hit return.\n";
                exit(EXIT_FAILURE);
            }

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
            std::cout<<"Can not open file .dat to write data.\n";
            std::cout << "DGSolver will exit after you hit return.\n";
            exit(EXIT_FAILURE);
        }
    }
}

namespace decomposeMesh {
    //NOTE!! ALL FUNCTIONS OF THIS NAMESPACE MUST BE EXCUTED AFTER LOADING MESH

    std::vector<int> loadPartitionedMesh()
    {
        meshVar::rankOf2DElem.resize(meshVar::nelem2D);
        meshVar::Elem2DlocalIdWithRank.resize(meshVar::nelem2D);
        std::vector<int> idMarker(systemVar::totalProc,0);
        /*Declare loading locations*/
        std::string  procIdLoc = systemVar::pwd + "/DGMesh.mesh.epart";

        /*Load process id*/
        std::ifstream procFlux(procIdLoc.c_str());
        if (procFlux)
        {
            std::string line(" ");
            int elemId(0), nproc(0);
            while (std::getline(procFlux, line))
            {
                std::istringstream procData(line);
                procData >> nproc;
                meshVar::rankOf2DElem[elemId]=nproc;
                meshVar::Elem2DlocalIdWithRank[elemId]=idMarker[nproc];
                elemId++;
                idMarker[nproc]++;
            }
        }
        else
        {
            std::cout<<"Can not open file DGMesh.mesh.epart.\n";
            std::cout << "DGSolver will exit after you hit return.\n";
            exit(EXIT_FAILURE);
        }
        return idMarker;
    }

    std::vector<int> findLocalIdOfPts()
    {
        //NOTE: RUN AFTER LOADING PARTITIONED MESH
        int ielem(-1), elemRank(-1);
        //std::vector<int> counter(systemVar::totalProc,0);
        auxUlti::resize2DIntArray(meshVar::PointslocalIdWithRank, meshVar::npoin, systemVar::totalProc);
        std::vector<int> idMarker(systemVar::totalProc,0), helpArray(systemVar::totalProc,0);
        for (int ipoin = 0; ipoin < meshVar::npoin; ipoin++)
        {
            for (int iesup = meshVar::esup2[ipoin] + 1; iesup <= meshVar::esup2[ipoin + 1]; iesup++)
            {
                ielem = meshVar::esup1[iesup];
                elemRank=meshVar::rankOf2DElem[ielem];
                if (helpArray[elemRank]==0)
                {
                    helpArray[elemRank]=1;
                }
            }

            for (int irank=0;irank<systemVar::totalProc;irank++) {
                if (helpArray[irank]==1)
                {
                    meshVar::PointslocalIdWithRank[ipoin][irank]=idMarker[irank];
                    idMarker[irank]++;
                }
                else {
                    meshVar::PointslocalIdWithRank[ipoin][irank]=-1;
                }
                helpArray[irank]=0;
            }
        }
        return idMarker;
    }

    void findEdgeWithBCTypeMatched()
    {
    	//Run this function after loadPartitionedMesh()
        int masterCell, slaveCell, masterCellRank, slaveCellRank;
        for (int iedge=0;iedge<meshVar::inpoedCount;iedge++) {
            std::tie(masterCell,slaveCell)=auxUlti::getMasterServantOfEdge(iedge);
            if (masterCell>=0)
            {
                masterCellRank=meshVar::rankOf2DElem[masterCell];
            }
            if (slaveCell>=0)
            {
                slaveCellRank=meshVar::rankOf2DElem[slaveCell];
            }
            else {
                slaveCellRank=-1;
            }

            if ((masterCellRank!=slaveCellRank)&&(slaveCellRank>=0))
            {
                meshVar::inpoed[iedge][3]=10;
            }
        }
    }

    std::vector<int> getMeshInforOfRanks(std::vector<std::vector<std::vector<double>>>&Points, std::vector<std::vector<std::vector<int>>>&Elem1D, std::vector<std::vector<std::vector<int>>>&Elem2D, std::vector<std::vector<std::vector<int>>>&meshConnection)
    {
        int localptsId;
        for (int ipoint = 0; ipoint < meshVar::npoin; ++ipoint) {
            for (int irank = 0; irank < systemVar::totalProc; ++irank) {
                localptsId=meshVar::PointslocalIdWithRank[ipoint][irank];
                if (localptsId>=0)
                {
                    Points[irank][localptsId][0]=localptsId+1;
                    for (int i = 0; i < 2; ++i) {
                        Points[irank][localptsId][i+1]=meshVar::Points[ipoint][i];
                    }
                    Points[irank][localptsId][3]=0.0;
                }
            }
        }

        int elemRank(0), ptId;
        for (int ielem = 0; ielem < meshVar::nelem2D; ++ielem) {
            elemRank=meshVar::rankOf2DElem[ielem];
            int localId(meshVar::Elem2DlocalIdWithRank[ielem]);
            Elem2D[elemRank][localId][0]=meshVar::Elem2DlocalIdWithRank[ielem]+1;
            for (int i = 0; i < 4; ++i) {
                ptId=meshVar::Elements2D[ielem][i];
                if (ptId>=0)
                {
                    Elem2D[elemRank][localId][i+1]=meshVar::PointslocalIdWithRank[ptId][elemRank]+1;
                }
                else {
                    Elem2D[elemRank][localId][i+1]=-23;
                }
            }
        }

        int pt1, pt2, BCType, BCGrp, masterElem, slaveElem, masterRank, slaveRank, id;
        std::vector<int> idMarker(systemVar::totalProc,0);
        for (int iedge = 0; iedge < meshVar::inpoedCount; ++iedge) {
            pt1=meshVar::inpoed[iedge][0];
            pt2=meshVar::inpoed[iedge][1];
            BCType=auxUlti::getBCType(iedge);
            BCGrp=auxUlti::getGrpOfEdge(iedge);
            std::tie(masterElem,slaveElem)=auxUlti::getMasterServantOfEdge(iedge);
            masterRank=meshVar::rankOf2DElem[masterElem];

            if (slaveElem>=0)
            {
                slaveRank=meshVar::rankOf2DElem[slaveElem];
            }
            else {
                slaveRank=-1;
            }

            if (BCType!=0)
            {
                for (int irank = 0; irank < systemVar::totalProc; ++irank) {
                    if ((meshVar::PointslocalIdWithRank[pt1][irank]>=0) && (meshVar::PointslocalIdWithRank[pt2][irank]>=0))
                    {
                        id=idMarker[irank];
                        Elem1D[irank][id][0]=idMarker[irank]+1;
                        Elem1D[irank][id][1]=meshVar::PointslocalIdWithRank[pt1][irank]+1;
                        Elem1D[irank][id][2]=meshVar::PointslocalIdWithRank[pt2][irank]+1;

                        if (BCType==10)
                        {
                            Elem1D[irank][id][3]=meshVar::numBCGrp+1;
                            //Save information to meshConnection array
                            if (masterRank==irank)
                            {
                                meshConnection[irank][id][0]=slaveRank;
                                meshConnection[irank][id][1]=slaveElem;
                            }
                            else if (slaveRank==irank) {
                                meshConnection[irank][id][0]=masterRank;
                                meshConnection[irank][id][1]=masterElem;
                            }
                        }
                        else {
                            Elem1D[irank][id][3]=BCGrp;
                            meshConnection[irank][id][0]=-1;
                            meshConnection[irank][id][1]=-1;
                        }
                        idMarker[irank]++;
                    }
                }
            }
        }
        return idMarker;
    }

    void decomposingMesh()
    {
    	int maxNumOfElem2D(0), maxNumOfPts(0);
    	std::vector<int> maxElem2DIdOfRanks(systemVar::totalProc,0),
    	maxElem1DIdOfRanks(systemVar::totalProc,0),
    	maxPtsIdOfRanks(systemVar::totalProc,0);

    	//Load partitioned mesh
    	maxElem2DIdOfRanks=decomposeMesh::loadPartitionedMesh();
        maxNumOfElem2D=*std::max_element(maxElem2DIdOfRanks.begin(), maxElem2DIdOfRanks.end());

    	//Find local id of points
    	maxPtsIdOfRanks=decomposeMesh::findLocalIdOfPts();
        maxNumOfPts=*std::max_element(maxPtsIdOfRanks.begin(), maxPtsIdOfRanks.end());

    	//Mark BCType "matched"
    	decomposeMesh::findEdgeWithBCTypeMatched();

    	std::vector<std::vector<std::vector<int>>> Elem2D (systemVar::totalProc,std::vector<std::vector<int>>(maxNumOfElem2D,std::vector <int>(4,0))),
        Elem1D (systemVar::totalProc,std::vector<std::vector<int>>(systemVar::totalProc*meshVar::inpoedCount,std::vector <int>(4,0))),
        meshConnection (systemVar::totalProc,std::vector<std::vector<int>>(systemVar::totalProc*meshVar::inpoedCount,std::vector <int>(2,0)));
    	std::vector<std::vector<std::vector<double>>> Points (systemVar::totalProc,std::vector<std::vector<double>>(maxNumOfPts,std::vector <double>(4,0)));

    	maxElem1DIdOfRanks=decomposeMesh::getMeshInforOfRanks(Points,Elem1D,Elem2D,meshConnection);

    	/*CREATE DECOMPOSED CASE*/
        //Delete previous processor
        std::string key;
        std::cout<<"Do you want to delete previous Processor folders! <y/n>: ";
        std::cin>>key;
        if (key.compare("y") == 0)
        {
            std::cout<<"Deleting Processor folders...\n";
            std::string command("rm -rf "+systemVar::wD + "/CASES/" + systemVar::caseName + "/Processor*");
            const int dir_err = system(command.c_str());
            if (-1 == dir_err)
            {
                printf("Error of deleting Processor folders\n");
                exit(1);
            }
        }

        std::cout<<"Decomposing case...\n";
        for (int irank = 0; irank < systemVar::totalProc; ++irank)
		{
			//Create processor folders
            std::cout<<"Processor "<<irank<<std::endl;
	    	std::string rank = std::to_string(irank);
	        std::string Loc(systemVar::wD + "/CASES/" + systemVar::caseName + "/Processor" + rank);
	        auxUlti::createFolder(Loc);
	        auxUlti::createFolder(Loc+"/Constant");
	        auxUlti::createFolder(Loc+"/Constant/Mesh");

            auxUlti::createFolder(Loc+"/0");
            auxUlti::copyFolder(systemVar::wD + "/CASES/" + systemVar::caseName+"/0",systemVar::wD + "/CASES/" + systemVar::caseName+ "/Processor" + rank + "/0");

            std::string  elems1DLoc = systemVar::pwd + "/Processor" + rank + "/Constant/Mesh/Elements1D.txt",
            elems2DLoc = systemVar::pwd + "/Processor" + rank + "/Constant/Mesh/Elements2D.txt",
            pointsLoc = systemVar::pwd + "/Processor" + rank + "/Constant/Mesh/Points.txt",
            meshConnectionLoc = systemVar::pwd + "/Processor" + rank + "/Constant/Mesh/meshConnection.txt";

			//Elements2D
            IO::write2DIntArrayToFile(Elem2D[irank], elems2DLoc, maxElem2DIdOfRanks[irank], 4);

			//Points
            IO::write2DDoubleArrayToFile(Points[irank], pointsLoc, maxPtsIdOfRanks[irank], 4);

            //Elements1D
            IO::write2DIntArrayToFile(Elem1D[irank], elems1DLoc, maxElem1DIdOfRanks[irank], 4);

            //Elements1D
            IO::write2DIntArrayToFile(meshConnection[irank], meshConnectionLoc, maxElem1DIdOfRanks[irank], 2);

            decomposeMesh::exportPartitionedMesh(irank,maxPtsIdOfRanks[irank],maxElem2DIdOfRanks[irank],Points[irank],Elem2D[irank]);

            decomposeMesh::decomposingTime0(Loc);

            //Create file boundaryPatch and Material
            auxUlti::copyFile(systemVar::wD + "/CASES/" + systemVar::caseName + "/Constant/boundaryPatch.txt",Loc+"/Constant");
            auxUlti::copyFile(systemVar::wD + "/CASES/" + systemVar::caseName + "/Constant/Material.txt",Loc+"/Constant");
            std::string content1 = R"(matchedBoundary
{
        Group )", content2 = R"(
        Type				matched
}
    )";
            IO::openFileToAppend(Loc+"/Constant/boundaryPatch.txt",content1+std::to_string(meshVar::numBCGrp+1)+content2);
		}
        std::cout<<"DONE!\n";
    }

    void decomposingTime0(std::string Loc)
    {
        std::string  pLoc = Loc + "/0/p.txt",
        TLoc = Loc + "/0/T.txt", ULoc = Loc + "/0/U.txt";

        std::string content1 = R"(matchedBoundary
{
        Group )", content2 = R"(
        Type				matched
}
)";
        IO::openFileToAppend(pLoc,content1+std::to_string(meshVar::numBCGrp+1)+content2);
        IO::openFileToAppend(TLoc,content1+std::to_string(meshVar::numBCGrp+1)+content2);
        IO::openFileToAppend(ULoc,content1+std::to_string(meshVar::numBCGrp+1)+content2);
    }

    void exportPartitionedMesh(int rank, int npoin, int nelem2D, std::vector<std::vector<double>>&Points, std::vector<std::vector<int>>&Elements2D)
    {
        //DIDN'T CHECK FOR QUADRILATERAL ELEMENTS

        std::string fileName(systemVar::caseName + "proc"+std::to_string(rank)+".dat"), Loc(systemVar::wD + "/CASES/" + systemVar::caseName + "/TecplotFile/PartitionedMesh"), code;
        auxUlti::createFolder(Loc);

        std::string fileLoc(Loc + "/" + fileName);
        std::ofstream fileFlux(fileLoc.c_str());

        code = R"(
TITLE     = "DG2D to Tecplot"
VARIABLES = "X"
"Y"
"PROC_RANK"
ZONE T="ZONE 1"
STRANDID=0
ZONETYPE=FEQuadrilateral
DATAPACKING=BLOCK
VARLOCATION=([3]=CELLCENTERED)
DT=(SINGLE SINGLE SINGLE)
)";

        if (fileFlux)
        {
            fileFlux << code << std::endl << "Nodes=" << std::to_string(npoin) << ", " << "Elements=" << std::to_string(nelem2D) << std::endl;
            //X
            int counter(0);
            for (int i = 0; i < npoin; i++)
            {
                counter++;
                if (counter == 5)
                {
                    fileFlux << std::endl;
                    counter = 0;
                }
                fileFlux << Points[i][1] << " ";
            }
            fileFlux << std::endl;

            //Y
            counter = 0;
            for (int i = 0; i < npoin; i++)
            {
                counter++;
                if (counter == 5)
                {
                    fileFlux << std::endl;
                    counter = 0;
                }
                fileFlux << Points[i][2] << " ";
            }
            fileFlux << std::endl;

            //RHO
            counter = 0;
            for (int i = 0; i < nelem2D; i++)
            {
                counter++;
                if (counter == 5)
                {
                    fileFlux << std::endl;
                    counter = 0;
                }
                fileFlux << rank << " ";
            }
            fileFlux << std::endl;

            //CONNECTIVITY
            for (int ielem = 0; ielem < nelem2D; ielem++)
            {
                int elemType(3);
                for (int ipoin = 1; ipoin < 4; ipoin++)
                {
                    fileFlux << Elements2D[ielem][ipoin] << " ";
                }
                switch (elemType)
                {
                case 3:
                {
                    fileFlux << Elements2D[ielem][1] << std::endl;
                    break;
                }
                case 4:
                {
                    fileFlux << Elements2D[ielem][4] << std::endl;
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
}
