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
	void basisFc(double a, double b, int elemType);

	/*Function calculates derivatives of basis function respect to a, b*/
	void dBasisFc(double a, double b, int elemType);

	/*Function calculates 2D Jacobian*/
	double J2DCal(int elem, double a, double b);

	/*Function calculates 1D Jacobian for an edge, it returns master Jacobi and servant Jacobi*/
	std::tuple<double, double> J1DCal(int edge);

	/*Function calculates integral which used to set initial values for all elements of mesh*/
	double iniIntegral(int elemType, int order);

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

    /*Function calculates T from rho, rhou, rhov, rhoE, rhoX, rhoY in case of mass diffusion*/
    double CalcTFromConsvVar_massDiff(double rho, double rhou, double rhov, double rhoE, double rhox, double rhoy);

	/*Function calculates p from T and rho*/
	double CalcP(double T, double rho);

	/*Function calculates dxa, dxb, dya, dyb at inputed Gauss point*/
	std::tuple<double, double, double, double> Calc_dxydab(int elem, double a, double b);

	/*Function calculates dBx, dBy at Gauss point (na,nb) of element at order*/
	double Calc_dBxdBy(int elem, int order, int na, int nb, int opt);

	/*Function calculates surface integral (for 2D case, surface integral is equal to line integral)*/
    double surfaceInte(std::vector<double> &Fvalue, int edge);

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
	double pointValueNoLimiter(int element, double a, double b, int valType);

    double pointDerivRho(int element, double a, double b, int dir);

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
	std::tuple<double, double>directMapping(int element, double aCoor, double bCoor);

	//Function maps point coordinates from real element to standard element
	std::tuple<double, double> mappingRealToStd(int edge, int element, double xCoor, double yCoor);

	//Function maps point coordinates from real element to standard element
	std::tuple<double, double> inverseMapping(int element, double xCoor, double yCoor);

	/*
	//Function supports for math::mappingRealToStd
	double solve_abQuad(int option, double A, double B, double D, double C, double inVar);

	//Function supports for math::mappingRealToStd
	double solve_abTri(int option, double A, double B, double C, double inVar);*/

	/*Function computes value of conservative variables at abitrary point with applying limiter
			valType:
			1: rho
			2: rhou
			3: rhov
			4: rhoE*/
	double calcConsvVarWthLimiter(int element, double a, double b, int valType);

	//Function returns id of input number in input array
	int findIndex(int number, std::vector<int> InArray);

	//Function finds product of tensor and vector
	std::vector<double> tensorVectorDotProduct(std::vector<std::vector<double>> tensor, std::vector<double> vector);

	//Function finds norm of vector
	double vectorNorm(std::vector<double> vector);

	double calcResidualFromResidualOfOrder(int element, double a, double b, int valType);

    namespace solvePolynomialsEq {
        double NewtonRaphson(std::vector<double> &power, std::vector<double> &coefs, double initialValue);

        double Bisection(std::vector<double> &power, std::vector<double> &coefs, double initialValue);

        double subValToPolynomial(std::vector<double> &power, std::vector<double> &coefs, double Value);
    }

	namespace numericalFluxes
	{
		/*Function calculates auxilary flux at Gauss point*/
		double auxFlux(double MinusVal, double PlusVar, double vectorComp);

		/*Function calculates advective flux*/
		double advectiveFlux(double FPlus, double FMinus, double UPlus, double UMinus, double C, double vectorComp);

		/*Function calculates diffusive flux*/
		double diffusiveFlux(double FPlus, double FMinus, double UPlus, double UMinus, double Beta, double vectorComp);

		/*Function calculates constant C for advective flux*/
		double constantC(double uMagP, double uMagM, double aP, double aM);

		/*Function calculates constant Beta for diffusive flux*/
        double constantBeta(double uMagP, double uMagM, double rhoP, double rhoM, double eP, double eM, double pP, double pM, std::vector<std::vector<double>> stressHeatFluxP, std::vector<std::vector<double>> stressHeatFluxM, std::vector<double> nP);
	}

	namespace inviscidTerms
	{
		/*Function calculates inviscid terms of NSF equation from primary variables*/
        std::tuple<double, double, double, double> calcInvisTermsFromPriVars(double rhoVal, double uVal, double umVal, double vVal, double vmVal, double totalE, double pVal, int dir);
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
		double calcStressComponent(int index, double fstDeriv, double sndDeriv);

		/*Function calculates stress tensor from viscosity, conservative variables and derivative of conservative variables, stress tensor has form:
		[tau_xx	    tau_xy		Qx]
		[tau_yx	    tau_yy		Qy]
		*/
        std::vector<std::vector<double>> calcStressTensorAndHeatFlux(std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy, double TVal);

		/*Function calculates heat flux terms Qx, Qy*/
		std::tuple<double, double> calcHeatFluxTerms(double dTx, double dTy, double k);

		/*Function calculates viscous terms of NSF equation from Stress and Heat flux matrix returned from calcStressTensorAndHeatFlux function*/
        std::tuple<double, double, double, double> calcViscousTermsFromStressHeatFluxMatrix(std::vector< std::vector<double> > &StressHeatFlux, double uVal, double vVal, double dRho, int dir);
	}

	namespace geometricOp
	{
		//Function computes geometric center of polygon, inputs are coordinates of vertexs of polygon
		std::tuple<double, double> calcGeoCenter(std::vector<double> &xCoor, std::vector<double> &yCoor, int type);

		//Function computes area of polygon
		double calcPolygonArea(std::vector<double> &xCoor, std::vector<double> &yCoor, int type);

		//Function computes centroid of quad elements, for tri elements, centroid is coincident with geometric center
		std::tuple<double, double> calcQuadCentroid(int element, double xCG, double yCG, double area);

		//Function calculates edge metrics
		std::tuple<double, double> calEdgeMetric(int edgeId, int elementId);

		//
		std::tuple<double, double> calDifferenceOfElementsCellCenters(int elem1, int elem2);

		//Function computes local cell size
		double calLocalCellSize(int element, double elementArea);

		double calDistBetween2Points(double xPt1, double yPt1, double xPt2, double yPt2);

		double calPolygonPerimeter(std::vector<double> &xCoor, std::vector<double> &yCoor, int numOfEdge);

        std::tuple<double, double> calROfInscribedCircleOfTriElement(std::vector<double> &xCoor, std::vector<double> &yCoor);

        std::tuple<double, double> calSizeOfQuadElement(std::vector<double> &xCoor, std::vector<double> &yCoor);
	}

	namespace residualManipulation
	{
		//Function calculates normalized coefficients for residuals
		//void calcNormResidual(double rhoRes, double rhouRes, double rhovRes, double rhoERes);
	}

	double calcMaxT(int element);

    namespace massDiffusionFncs {
    double calcTotalVelocity(double rho, double advecV, double mudRho);
    }
}

#endif // DGMATH_H_INCLUDED
