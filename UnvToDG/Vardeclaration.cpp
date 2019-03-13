#include <QString>
#include <vector>
#include "VarDeclaration.h"

/*VARIABLE DECLARATION-------------------------------------------------------------------*/
int const meshArrSize(1500000);
int const pointsArrSize(500000);
int const elements2DArrSize(1000000);

QString caseName(" "), wD(" "), pwd(" "), meshFileName(" ");
std::vector<std::vector<double>> Points(pointsArrSize, std::vector<double>(4, 0));

std::vector<std::vector<int>> location(3, std::vector<int>(2, 0)),
Elements1D(pointsArrSize, std::vector<int>(4, 0)),
Elements2D(elements2DArrSize, std::vector<int>(5, 0)),
boundaries(pointsArrSize, std::vector<int>(2, 0));

std::vector<int> boundLocation;  //variable for GetBoundaries
std::vector<QString> boundName, Mesh(meshArrSize);
bool savingFlag(false);  //Flag of saving data
/*Declare size of Points, Elements1D, Elements2D*/
int nodeNumber(0), n1D(0), n2D(0), numOfBoundEdge(1), numOfBound(0), noLines(0);
