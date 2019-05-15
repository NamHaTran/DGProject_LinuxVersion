#ifndef DGLIMITERLIB_H_INCLUDED
#define DGLIMITERLIB_H_INCLUDED
#include <tuple>
#include <vector>
namespace limiter
{
	//Function applies limiter to input element
	void limiter();

	namespace Pp
	{
		namespace triangleCell
		{
			//Function computes coefficients theta1 theta2 for positivity preserving limiter
			std::tuple<double, double> calcPpLimiterCoef(int element);
		}

		namespace quadratureCell
		{
			//Function computes coefficients theta1 theta2 for positivity preserving limiter
			std::tuple<double, double> calcPpLimiterCoef(int element);

			namespace simplifiedVersion
			{
				std::tuple<double, double> calcPpLimiterCoef(int element);
			}
		}
	}

	namespace pAdaptive
	{
		void pAdaptiveLimiter(int element, int valType);

		std::vector<double> pAdaptiveChildFunction_Quad(int element, int valType, double IPlus, double IMinus, double JPlus, double JMinus);

		std::vector<double> pAdaptiveChildFunction_Tri(int element, int valType, double IPlus, double IMinus, double JMinus);
	}

	namespace mathForLimiter
	{
		namespace triangleCell
		{
			//Function calculates minimum value of rho of tri element
			double calcMinRho(int element);

			//Function computes theta1 coefficient and omega for Pp limiter
			std::tuple<double, double> calcTheta1Coeff(double meanRho, double minRho, double meanP);

			//Function calculates modified value of Rho at abitrary point (for calculating theta2)
			double calcRhoModified(int element, double a, double b, double theta1);

			//Function checks condition of running limiter
			bool checkLimiter(int element, double theta1, double omega);

			//Function supports for computing limiter
			std::tuple<double, double> calcXYBySigma(double sigma, double xi, double yi, double xC, double yC);
			//Function supports for computing limiter
			double calcP(int element, double x, double y, double theta1);

			//Function computes theta2 for Pp limiter
			double calcTheta2Coeff(int element, double theta1, double omega);
		}

		namespace quadratureCell
		{
			//Function calculates minimum value of rho of quad element
			double calcMinRhoQuad(int element);

			//Function calculates modified value of Rho at abitrary point (for calculating theta2)
			double calcRhoModified(int element, double a, double b, double theta1);

			//Function returns true if element is needed to limit, and value of rho which applied first time limiter
			std::tuple<bool, double> checkLimiterForQuad(int element, double a, double b);

			//Function computes theta1 coefficient and omega for Pp limiter
			std::tuple<double, double> calcTheta1Coeff(double meanRho, double minRho, double meanP);

			//Function computes theta2 at 1 Gauss point in input direction
			double calcTheta2Coeff(int element, double aG, double bG, double theta1, double omega, double meanRho, double meanRhou, double meanRhov, double meanRhoE);

			namespace simplifiedVersion
			{
				double calcMinRhoeQuad(int element, double theta1);

				double calRhoeFromConserVars(double rho, double rhou, double rhov, double rhoE);

				double calcTheta1Coeff(double minRho, double meanRho);

				double calcTheta2Coeff(double meanRhoe, double minRhoe, double meanRho);
			}
		}

		//Function gets sign of input double number (return -1, 0, 1 if input number is negative, zero or positive, respectively)
		int getSignOfDouble(double input);

		//Minmod function
		double minmod(std::vector<double> vectorArguments);

		//Modified minmon function
		double modifiedMinmod(std::vector<double> vectorArguments, double condition);

		double calM(int element, int valType);

		void getNeighborElements();
	}
}

namespace troubledCellDetection
{
	bool checkTroubledCell(std::vector<double> InputVector, double condition);
}
#endif // DGLIMITERLIB_H_INCLUDED
