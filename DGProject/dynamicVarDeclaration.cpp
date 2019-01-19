#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include <vector>

namespace meshVar
{
	/*Gauss points on edges*/
	std::vector<std::vector<double>>
		edgeGaussPoints_a(1, std::vector<double>(1, 0.0)),
		edgeGaussPoints_b(1, std::vector<double>(1, 0.0));
}

/*Conservative variables declaration
double rho[meshVar::nelem2D][maxOrder] = {},
rhou[meshVar::nelem2D][maxOrder] = {},
rhov[meshVar::nelem2D][maxOrder] = {},
rhoE[meshVar::nelem2D][maxOrder] = {};*/
std::vector<std::vector<double>>
rho(1, std::vector<double>(1, 0.0)),
rhou(1, std::vector<double>(1, 0.0)),
rhov(1, std::vector<double>(1, 0.0)),
rhoE(1, std::vector<double>(1, 0.0));
std::vector<std::vector<double>>
rhoN(1, std::vector<double>(1, 0.0)),
rhouN(1, std::vector<double>(1, 0.0)),
rhovN(1, std::vector<double>(1, 0.0)),
rhoEN(1, std::vector<double>(1, 0.0));

/*Primary variables declaration
double u[meshVar::nelem2D][maxOrder] = {},
v[meshVar::nelem2D][maxOrder] = {},
e[meshVar::nelem2D][maxOrder] = {},
p[meshVar::nelem2D][maxOrder] = {},
T[meshVar::nelem2D][maxOrder] = {},
mu[meshVar::nelem2D][maxOrder] = {};*/
std::vector<std::vector<double>>
u(1, std::vector<double>(1, 0.0)),
v(1, std::vector<double>(1, 0.0)),
e(1, std::vector<double>(1, 0.0)),
p(1, std::vector<double>(1, 0.0)),
T(1, std::vector<double>(1, 0.0)),
mu(1, std::vector<double>(1, 0.0));

/*Auxilary variables
//X direction
double rhoX[meshVar::nelem2D][maxOrder] = {},
rhouX[meshVar::nelem2D][maxOrder] = {},
rhovX[meshVar::nelem2D][maxOrder] = {},
rhoEX[meshVar::nelem2D][maxOrder] = {};*/
std::vector<std::vector<double>>
rhoX(1, std::vector<double>(1, 0.0)),
rhouX(1, std::vector<double>(1, 0.0)),
rhovX(1, std::vector<double>(1, 0.0)),
rhoEX(1, std::vector<double>(1, 0.0));

/*Y direction
double rhoY[meshVar::nelem2D][maxOrder] = {},
rhouY[meshVar::nelem2D][maxOrder] = {},
rhovY[meshVar::nelem2D][maxOrder] = {},
rhoEY[meshVar::nelem2D][maxOrder] = {};*/
std::vector<std::vector<double>>
rhoY(1, std::vector<double>(1, 0.0)),
rhouY(1, std::vector<double>(1, 0.0)),
rhovY(1, std::vector<double>(1, 0.0)),
rhoEY(1, std::vector<double>(1, 0.0));

//Limiting coefficients
std::vector<double>
theta1Arr(1, 1.0),
theta2Arr(1, 1.0);

/*Mean values
row1: mean rho
row2: mean rhou
row3: mean rhov
row4: mean rhoE*/
std::vector<std::vector<double>> meanVals(1, std::vector<double>(4, 0.0));
