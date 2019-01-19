#ifndef DYNAMICVARDECLARATION_H_INCLUDED
#define DYNAMICVARDECLARATION_H_INCLUDED
#include "VarDeclaration.h"
#include <vector>

namespace meshVar
{
	/*Gauss points on edges*/
	extern std::vector<std::vector<double>> edgeGaussPoints_a, edgeGaussPoints_b;
}

/*Conservative variables declaration
extern double rho[elements2DArrSize][maxOrder],
rhou[elements2DArrSize][maxOrder],
rhov[elements2DArrSize][maxOrder],
rhoE[elements2DArrSize][maxOrder];*/
extern std::vector<std::vector<double>> rho, rhou, rhov, rhoE;
extern std::vector<std::vector<double>> rhoN, rhouN, rhovN, rhoEN;

/*Primary variables declaration
extern double u[elements2DArrSize][maxOrder],
v[elements2DArrSize][maxOrder],
e[elements2DArrSize][maxOrder],
p[elements2DArrSize][maxOrder],
T[elements2DArrSize][maxOrder],
mu[elements2DArrSize][maxOrder];*/
extern std::vector<std::vector<double>>u, v, e, p, T, mu;

/*Auxilary variables
//X direction
extern double rhoX[elements2DArrSize][maxOrder],
rhouX[elements2DArrSize][maxOrder],
rhovX[elements2DArrSize][maxOrder],
rhoEX[elements2DArrSize][maxOrder];*/
extern std::vector<std::vector<double>> rhoX, rhouX, rhovX, rhoEX;

/*Y direction
extern double rhoY[elements2DArrSize][maxOrder],
rhouY[elements2DArrSize][maxOrder],
rhovY[elements2DArrSize][maxOrder],
rhoEY[elements2DArrSize][maxOrder];*/
extern std::vector<std::vector<double>> rhoY, rhouY, rhovY, rhoEY;

//time step
extern double dt, runTime;

//Limiting coefficients
extern std::vector<double>
theta1Arr,
theta2Arr;

//Mean values
extern std::vector<std::vector<double>> meanVals;
#endif // DYNAMICVARDECLARATION_H_INCLUDED
