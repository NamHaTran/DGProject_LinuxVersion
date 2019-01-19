#ifndef VARDECLARATION_H
#define VARDECLARATION_H
#include <QString>
#include <vector>
/*VARIABLE DECLARATION-------------------------------------------------------------------*/
extern int const meshArrSize;
extern int const pointsArrSize;
extern int const elements2DArrSize;

extern QString caseName, wD, pwd, meshFileName;
extern std::vector<std::vector<double>> Points;

extern std::vector<std::vector<int>> location, Elements1D, Elements2D, boundaries;

extern std::vector<int> boundLocation;  //variable for GetBoundaries
extern std::vector<QString> boundName, Mesh;
extern bool savingFlag;  //Flag of saving data
/*Declare size of Points, Elements1D, Elements2D*/
extern int nodeNumber, n1D, n2D, numOfBoundEdge, numOfBound, noLines;
#endif // VARDECLARATION_H
