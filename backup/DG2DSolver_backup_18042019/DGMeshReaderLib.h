#ifndef DGMESHREADERLIB_H_INCLUDED
#define DGMESHREADERLIB_H_INCLUDED
#include <tuple>  //Include this for returning multiple values in function
namespace MshReader
{
	/*Function creates matrixes which save mesh data get from files, by using EleSurPtProcess()*/
	void meshProcess();

	/*Function gets informations of elements surrounding points*/
	void EleSurPt();

	/*Function gets informations of points surrounding points*/
	void PtsSurPt();

	/*Function gets informations of elements surrounding element*/
	void ElemsSurElem();

	/*Function gets informations of edges*/
	void EdgesInfor();
	void EdgesOfElem();

	/*Calculate normal vector of each face (edge)*/
	void GetNormalVector();

	/*Function supports for EdgesInfor(), it returns index of checking number in input iarray*/
	int findIndex(int number, int iarray[], int size);

	//Note: run this function AFTER mesh processing is DONE
	void getBoundaryPoints();

	/*Child functions*/
	int CheckConnection(int point, int helpArray[], int length);
	int checkIndividualEdge(int rootPt, int tipPt);
	int findIndex(int number, int iarray[], int size);
	void getBcGrpTp(int ipoin, int jpoin, int ninpoed);
	std::tuple<double, double> calcNormVector(int point1, int point1Indice, int point2, int point2Indice, int type);

	void sortPointsOfElements();
}
#endif // DGMESHREADERLIB_H_INCLUDED
