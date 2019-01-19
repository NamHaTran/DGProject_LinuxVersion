#ifndef DGPROCLIB_H_INCLUDED
#define DGPROCLIB_H_INCLUDED
#include <vector>

namespace meshParam
{
	/*Function calculate Gaussian constants*/
	void GaussParam();

	/*Function calculate Jacobian*/
	void JacobianParam();

	/*Function calculate basis function*/
	void basisFcParam();

	/*Function saves coordinates derivatives to array*/
	void derivCoordinates();

	/*Function calculates centroid and cell size of elements*/
	void calcCellMetrics();
}

namespace process
{
	/*Function sets initial value to all elements*/
	void setIniValues();

	/*Function distributes initial value to each order of accuracy of element*/
	std::vector<double> calcIniValues(double iniVal, int element);

	/*Function computes RHS of initial condition equation*/
	std::vector<double> calcIniValuesRHS(int element, double iniVal);

	namespace auxEq
	{
		/*Function solves auxilary equation at all elements for auxilary variables*/
		void solveAuxEquation();

		/*Function calculates right hand side term of auxilary equations at all order*/
		void CalcRHSTerm(int element, std::vector<double> &rhoRHSOx, std::vector<double> &rhoRHSOy, std::vector<double> &rhouRHSOx, std::vector<double> &rhouRHSOy, std::vector<double> &rhovRHSOx, std::vector<double> &rhovRHSOy, std::vector<double> &rhoERHSOx, std::vector<double> &rhoERHSOy);

		/*Function calculates volume integral terms of element*/
		void calcVolumeIntegralTerms(int element, std::vector<double> &rhoRHSTerm, std::vector<double> &rhouRHSTerm, std::vector<double> &rhovRHSTerm, std::vector<double> &rhoERHSTerm, int dir);

		/*Function calculates surface integral terms of element*/
		void calcSurfaceIntegralTerms(int element, std::vector<double> &rhoSurfInt, std::vector<double> &rhouSurfInt, std::vector<double> &rhovSurfInt, std::vector<double> &rhoESurfInt, int dir);

		/*Function returns matrix content Gauss values of <valType> conservative variable at all Gauss points in element volume,
		use only for auxilary equation*/
		std::vector<std::vector<double>> getGaussMatrixOfConserVar(int element, int valType);

		/*Function returns vector of flux of all conservative variables at all faces of element*/
		void getGaussVectorOfConserVar(int element, std::vector<std::vector<double>> &rhoFlux, std::vector<std::vector<double>> &rhouFlux, std::vector<std::vector<double>> &rhovFlux, std::vector<std::vector<double>> &rhoEFlux, int dir);

		/*Function returns vector content Gauss values of flux of <valType> conservative variable at all Gauss points on the edge,
		use only for auxilary equation*/
		std::vector<double> getVectorOfConserVarFluxesAtInternal(int edge, int element, int nG, double nVectorComp);
	}

	namespace NSFEq
	{
		/*Function solves NSEF equation at all elements for conservative variables*/
		void solveNSFEquation();

		/*Function calculates right hand side terms of all conservative variables at ONLY one order*/
		void CalcRHSTerm(int element, std::vector<double> &term1RHS, std::vector<double> &term2RHS, std::vector<double> &term3RHS, std::vector<double> &term4RHS);

		/*Function calculates Inviscid terms at Gauss point (a, b) and returns 2D matrix
		InviscidTerm 2D array has 4 row 2 column:
		- column 1: Ox direction
		- column 2: Oy direction*/
		std::vector<std::vector<double>> calcGaussInviscidTerm(int element, double a, double b);

		/*Function calculates Viscous terms at Gauss point (a, b) and returns 2D matrix
		ViscousTerm 2D array has 4 row 2 column:
		- column 1: Ox direction
		- column 2: Oy direction*/
		std::vector<std::vector<double>> calcGaussViscousTerm(int element, double a, double b);

		/*Function calculates volume integral terms in NSF equation at ONLY ONE ORDER*/
		void calcVolumeIntegralTerms(int element, std::vector<double> &VolIntTerm1, std::vector<double> &VolIntTerm2, std::vector<double> &VolIntTerm3, std::vector<double> &VolIntTerm4);

		/*Function calculates surface integral terms in NSF equation at ONLY ONE ORDER*/
		void calcSurfaceIntegralTerms(int element, std::vector<double> &SurfIntTerm1, std::vector<double> &SurfIntTerm2, std::vector<double> &SurfIntTerm3, std::vector<double> &SurfIntTerm4);

		/*Function calculates flux at nGauss point of all conservative variables at internal egde*/
		std::vector<std::vector<double>> getGaussVectorOfConserVarFluxesAtInternal(int edgeName, int element, int nGauss);

		/*Function applies input ddtSchemme to solve time marching*/
		std::vector<double> solveTimeMarching(std::vector<double> &ddtArr, std::vector<double> &UnArr);
		
		/*Function updates values in conservative variable array*/
		void updateVariables();
	}

	namespace Euler
	{
		/*Function computes local time step*/
		double localTimeStep(int element);

		/*Function computes time derivative at centroid of standard cell*/
		double localErrorEstimate(int element, std::vector<double> &ddtArr);

		/*Function computes global errors*/
		std::tuple <double, double, double, double> globalErrorEstimate(std::vector<double> &RhoError, std::vector<double> &RhouError, std::vector<double> &RhovError, std::vector<double> &RhoEError);
	}

	namespace limiter
	{
		//Function applies limiter to input element
		void limiter();

		namespace Pp
		{
			//Function computes coefficients theta1 theta2 for positivity preserving limiter
			std::tuple<double, double> calcPpLimiterCoef(int element);
		}
	}

	/*Function calculates volume integral terms of auxilary equations of ONLY one order
	User's guide:
	elem: element index
	Ui: array of Gauss values
	direction: index of direction
	1: x direction
	2: y direction*/
	double volumeInte(int elem, std::vector< std::vector<double> > &Ui, int order, int direction);

	/*Function calculates surface integral at 1 surface of auxilary equations of ONLY one order
	User's guide:
	direction: index of direction
	1: x direction
	2: y direction*/
	double surfaceInte(int elem, int edge, std::vector<double> &FluxVector, int order);

	/*Function calculates Auxilary Stiff matrix*/
	std::vector<std::vector<double>> calculateStiffMatrix(int element);

	/*Function calculates value of element of Auxilary Stiff matrix*/
	double calculateStiffMatrixElement(int element, int order1, int order2);

	/*Function return false if running contion is wrong*/
	bool checkRunningCond();
}

#endif // DGPROCLIB_H_INCLUDED