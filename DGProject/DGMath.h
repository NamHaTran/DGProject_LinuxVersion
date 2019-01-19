#ifndef DGMATH_H_INCLUDED
#define DGMATH_H_INCLUDED
#include <tuple>
#include <vector>

namespace math
{
	/*Function calculates Gauss values*/
	void Gauss(int nGauss);

	/*Function calculates Gauss-Lobatto values*/
	void GaussLobatto(int nGauss);

	/*Function calculates basis function*/
	void basisFc(double a, double b);

	/*Function calculates derivatives of basis function respect to a, b*/
	void dBasisFc(double a, double b);

	/*Function calculates 2D Jacobian*/
	double J2DCal(int elem, double a, double b);

	/*Function calculates 1D Jacobian for an edge, it returns master Jacobi and servant Jacobi*/
	std::tuple<double, double> J1DCal(int edge);

	/*Function calculates integral which used to set initial values for all elements of mesh*/
	double iniIntegral(int order);

	/*Function calculates volume integral from Gauss values*/
	double volumeInte(std::vector< std::vector<double> > &Fvalue, int elem);

	/*Function calculates 2D jacobian for quad elements*/
	double jacobianQuad(double xA, double xB, double xC, double xD, double yA, double yB, double yC, double yD, double a, double b);

	/*Function calculates 2D jacobian for tri elements*/
	double jacobianTri(double xA, double xB, double xC, double yA, double yB, double yC, double a, double b);

	/*Function calculates 1D jacobian for 1 edge of quad elements*/
	double jacobian1DQuad(int edgeIndex, double xA, double xB, double xC, double xD, double yA, double yB, double yC, double yD);

	/*Function calculates 1D jacobian for 1 edge of tri elements*/
	double jacobian1DTri(int edgeIndex, double xA, double xB, double xC, double yA, double yB, double yC);

	/*Function solves system of linear equations by using Gauss-Seidel algorithm*/
	std::vector<double> SolveSysEqs(std::vector< std::vector<double> > &a, std::vector<double> &b);

	/*Function calculates error for Gauss-Seidel algorithm*/
	double errorGS(std::vector< std::vector<double> > &a, std::vector<double> &b, std::vector<double> &No);

	/*Function calculates dynamic viscosity coefficient by using Sutherland formular*/
	double CalcVisCoef(double T);

	/*Function calculates T at arbitrary point of arbitrary element from rho, rhou, rhov, rhoE (ver 1)*/
	double CalcT(int elem, double a, double b);

	/*Function calculates T from rho, rhou, rhov, rhoE (ver 2)*/
	double CalcTFromConsvVar(double rho, double rhou, double rhov, double rhoE);

	/*Function calculates p from T and rho*/
	double CalcP(double T, double rho);

	/*Function calculates dxa, dxb, dya, dyb at inputed Gauss point*/
	std::tuple<double, double, double, double> Calc_dxydab(int elem, double a, double b);

	/*Function calculates dBx, dBy at Gauss point (na,nb) of element at order*/
	double Calc_dBxdBy(int elem, int order, int na, int nb, int opt);

	/*Function calculates surface integral (for 2D case, surface integral is equal to line integral)*/
	double surfaceInte(std::vector<double> &Fvalue, int edge, int elem);

	/*Function calculates values of primary and conservative variables at arbitrary point of arbitrary element.
	ValKind 1: primary variable
		valType 1: rho
		2: u
		3: v
		4: e
		5: p
		6: T
		7: mu
	ValKind 2: conservative variable
		valType 1: rho
		2: rhou
		3: rhov
		4: rhoE*/
	double pointValue(int element, double a, double b, int valType, int valKind);
	
	/*Function calculates values of primary and conservative variables at arbitrary point of arbitrary element without applying limiter.
	ValKind 1: primary variable
		valType 1: rho
		2: u
		3: v
		4: e
		5: p
		6: T
		7: mu
	ValKind 2: conservative variable
		valType 1: rho
		2: rhou
		3: rhov
		4: rhoE*/
	double pointValueNoLimiter(int element, double a, double b, int valType, int valKind);

	/*Function calculates dot product of 2 vectors*/
	double vectorDotProduct(std::vector<double> &a, std::vector<double> &b);

	/*Function calculates values of primary and conservative variables on plus side and minus side of edge at arbitrary Gauss point
	ValKind 1: primary variable
	valType 1: rho
	2: u
	3: v
	4: e
	5: p
	6: T
	7: mu
	ValKind 2: conservative variable
	valType 1: rho
	2: rhou
	3: rhov
	4: rhoE*/
	std::tuple <double, double> internalSurfaceValue(int edge, int element, int nG, int valType, int valKind);

	/*Function calculates values of auxilary variables (derivative of conservative variables) on plus side and minus side of edge at arbitrary Gauss point
	valType 1: drho
	2: drhou
	3: drhov
	4: drhoE
	dir 1: Ox
	2: Oy*/
	std::tuple <double, double> internalSurfaceDerivativeValue(int edge, int element, int nG, int valType, int dir);

	/*Function calculates values of primary and conservative variables on plus side of edge at arbitrary Gauss point
	ValKind 1: primary variable
	valType 1: rho
	2: u
	3: v
	4: e
	5: p
	6: T
	7: mu
	ValKind 2: conservative variable
	valType 1: rho
	2: rhou
	3: rhov
	4: rhoE
	NOTE: SHOULD BE USED FOR BOUNDARY EDGEs*/
	double SurfaceValueFromMaster(int edge, int nG, int valType, int valKind);

	/*Function calculates sum of 2 vectors*/
	std::vector<double> vectorSum(std::vector<double> &a, std::vector<double> &b);

	/*Function calculates speed of sound*/
	double CalcSpeedOfSound(double T);

	/*Function calculates derivatives of conservative variables rho, rhou, rhov, rhoE from auxilary variables
	NOTE! direction of derivatives based on input values, for example:
		d(u)/dx = [d(rhou)/dx - d(rho)/dx.(rhou)/rho]/(rho)
		so input values are
		d(rho)/dx	----> dRhoVal
		d(rhou)/dx	----> dRhoUVal
		rhou		----> rhoUVal
		rho			----> rhoVal*/
	double calcRhouvEDeriv(double dRhoUVal, double dRhoVal, double rhoUVal, double rhoVal);

	/*Function calculates derivatives of T from auxilary variables
	NOTE! direction of derivatives based on input values, for example:
	d(T)/dx = f(dE/dx, du/dx, dv/dx, u, v)
	so input values are
	dE/dx	----> dE
	du/dx	----> du
	dv/dx	----> dv*/
	double calcTDeriv(double dE, double du, double dv, double u, double v);

	/*Function calculates values of auxilary variables at arbitrary point of arbitrary element.
	valType 1: drho
	2: drhou
	3: drhov
	4: drhoE
	dir 1: Ox
	dir 2: Oy*/
	double pointAuxValue(int element, double a, double b, int valType, int dir);

	/*Function calculates thermal conductivity*/
	double calcThermalConductivity(double muVal);

	//Function solves quadratic equation
	std::tuple<bool, double, double> solvQuadraticEq(double A, double B, double C);

	//Function computes value of primary and conservative variables at center of cell: WARNING THIS BC CAUSES BLOW UP!!!
	double centerValue(int element, int valType, int valKind);

	//Function computes value of auxilary variables at center of cell: WARNING THIS BC CAUSES BLOW UP!!!
	double centerAuxValue(int element, int valType, int dir);

	//Function maps point coordinates from standard element to real element
	std::tuple<double, double> mappingStdToReal(int element, double aCoor, double bCoor);

	//Function maps point coordinates from real element to standard element
	std::tuple<double, double> mappingRealToStd(int edge, int element, double xCoor, double yCoor);

	/*
	//Function supports for math::mappingRealToStd
	double solve_abQuad(int option, double A, double B, double D, double C, double inVar);

	//Function supports for math::mappingRealToStd
	double solve_abTri(int option, double A, double B, double C, double inVar);*/

	namespace numericalFluxes
	{
		/*Function calculates auxilary flux at Gauss point*/
		double auxFlux(double MinusVal, double PlusVar, double vectorComp);

		/*Function calculates advective flux*/
		double advectiveFlux(double FPlus, double FMinus, double UPlus, double UMinus, double C, double vectorComp);

		/*Function calculates diffusive flux*/
		double diffusiveFlux(double MinusVal, double PlusVar, double vectorComp);

		/*Function calculates fluxes of all inviscid terms of NSF equation from conservative variables*/
		std::vector<std::vector<double>> NSFEqFluxFromConserVars(std::vector<double> &UPlus, std::vector<double> &UMinus, std::vector<double> &dUXPlus, std::vector<double> &dUXMinus, std::vector<double> &dUYPlus, std::vector<double> &dUYMinus, std::vector<double> &normVector);

		/*Function calculates constant C for advective flux*/
		double constantC(double uMagP, double uMagM, double aP, double aM);
	}

	namespace inviscidTerms
	{
		/*Function calculates inviscid terms of NSF equation from primary variables*/
		std::tuple<double, double, double, double> calcInvisTermsFromPriVars(double rhoVal, double uVal, double vVal, double totalE, double pVal, int dir);
	}

	namespace viscousTerms
	{
		/*Function calulates stress
		Formular of stress
		tau_xx = -mu((4/3)*du/dx - (2/3)*dv/dy)
		tau_yy = -mu((4/3)*dv/dy - (2/3)*du/dx)
		tau_xy = -mu(du/dy + dv/dx)

		input variables:
		index			fstDeriv		sndDeriv
		11 (tau_xx)		du/dx			dv/dy
		22 (tau_yy)		dv/dy			du/dx
		12 (tau_xy)		du/dy			dv/dx
		*/
		double calcStressComponent(int index, double mu, double fstDeriv, double sndDeriv);

		/*Function calculates stress tensor from viscosity, conservative variables and derivative of conservative variables, stress tensor has form:
		[tau_xx	    tau_xy		Qx]
		[tau_yx	    tau_yy		Qy]
		*/
		std::vector<std::vector<double>> calcStressTensorAndHeatFlux(double muVal, std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy);

		/*Function calculates heat flux terms Qx, Qy*/
		std::tuple<double, double> calcHeatFluxTerms(double dTx, double dTy, double k);

		/*Function calculates viscous terms of NSF equation from Stress and Heat flux matrix returned from calcStressTensorAndHeatFlux function*/
		std::tuple<double, double, double, double> calcViscousTermsFromStressHeatFluxMatrix(std::vector< std::vector<double> > &StressHeatFlux, double uVal, double vVal, int dir);
	}

	namespace limiter
	{
		/*Function calculates mean value of input valType of quad element
		Available of valType
		valType 1: rho
		2: rhou
		3: rhov
		4: rhoE
		*/
		double calcMeanConsvVarQuad(int element, int valType);

		/*Function calculates mean value of input valType of tri element
		Available of valType
		valType 1: rho
		2: rhou
		3: rhov
		4: rhoE
		*/
		double calcMeanConsvVarTri(int element, int valType);

		//Function calculates minimum value of rho of quad element
		double calcMinRhoQuad(int element);

		//Function calculates minimum value of rho of tri element
		double calcMinRhoTri(int element);

		//Function calculates minimum value of p of tri element
		double calcMinPTri(int element);

		//Function calculates modified value of Rho at abitrary point (for calculating theta2)
		double calcRhoModified(int element, double a, double b, double theta1, double rhoMean);

		/*Function computes value of conservative variables at abitrary point with applying limiter
		valType:
		1: rho
		2: rhou
		3: rhov
		4: rhoE*/
		double calcConsvVarWthLimiter(int element, double a, double b, int valType);

		//Function returns true if element is needed to limit, and value of rho which applied first time limiter
		std::tuple<bool, double> checkLimiterForQuad(int element, double a, double b);

		//Function computes theta1 coefficient and omega for Pp limiter
		std::tuple<double, double> calcTheta1Coeff(double meanRho, double minRho, double meanP);

		//Function computes theta2 at 1 Gauss point in input direction
		double calcTheta2Coeff(int element, int na, int nb, double theta1, double omega, double meanRho, double meanRhou, double meanRhov, double meanRhoE, int dir);
	}

	namespace geometricOp
	{
		//Function computes geometric center of polygon, inputs are coordinates of vertexs of polygon
		std::tuple<double, double> calcGeoCenter(std::vector<double> &xCoor, std::vector<double> &yCoor, int type);

		//Function computes area of polygon
		double calcPolygonArea(std::vector<double> &xCoor, std::vector<double> &yCoor, int type);

		//Function computes centroid of quad elements, for tri elements, centroid is coincident with geometric center
		std::tuple<double, double> calcQuadCentroid(int element, double xCG, double yCG, double area);
	}

	namespace residualManipulation
	{
		//Function calculates normalized coefficients for residuals
		void calcNormResidual(double rhoRes, double rhouRes, double rhovRes, double rhoERes);
	}
}

#endif // DGMATH_H_INCLUDED
