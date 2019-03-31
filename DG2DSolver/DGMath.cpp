#include "DGMath.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "DGMessagesLib.h"
#include <math.h>
#include <tuple>  //Include this for returning multiple values in function
#include <algorithm>
#include <iostream>

namespace math
{
	void Gauss(int nGauss)
	{
		if (nGauss==0)
		{
			mathVar::xGauss[nGauss] = 0.0;
			mathVar::wGauss[nGauss] = 2.0;
		}
		else if (nGauss==1)
		{
			mathVar::xGauss[nGauss - 1] = -0.577350269189625764509148780502;
			mathVar::xGauss[nGauss] = 0.577350269189625764509148780502;

			mathVar::wGauss[nGauss - 1] = 1.0;
			mathVar::wGauss[nGauss] = 1.0;
		}
		else if (nGauss==2)
		{
			mathVar::xGauss[nGauss - 2] = -0.774596669241483377035853079956;
			mathVar::xGauss[nGauss - 1] = 0.0;
			mathVar::xGauss[nGauss] = 0.774596669241483377035853079956;

			mathVar::wGauss[nGauss - 2] = 0.555555555555555555555;
			mathVar::wGauss[nGauss - 1] = 0.888888888888888888888;
			mathVar::wGauss[nGauss] = 0.55555555555555555555555;
		}
        else if (nGauss==3)
        {
            mathVar::xGauss[nGauss - 3] = -0.8611363115940526;
            mathVar::xGauss[nGauss - 2] = -0.3399810435848563;
            mathVar::xGauss[nGauss - 1] = 0.3399810435848563;
            mathVar::xGauss[nGauss] = 0.8611363115940526;

            mathVar::wGauss[nGauss - 3] = 0.3478548451374538;
            mathVar::wGauss[nGauss - 2] = 0.6521451548625461;
            mathVar::wGauss[nGauss - 1] = 0.6521451548625461;
            mathVar::wGauss[nGauss] = 0.3478548451374538	;
        }
        else if (nGauss==4)
        {
            mathVar::xGauss[nGauss - 4] = -0.9061798459386640;
            mathVar::xGauss[nGauss - 3] = -0.5384693101056831;
            mathVar::xGauss[nGauss - 2] = 0.0000000000000000;
            mathVar::xGauss[nGauss - 1] = 0.5384693101056831;
            mathVar::xGauss[nGauss] = 0.9061798459386640;

            mathVar::wGauss[nGauss - 4] = 0.2369268850561891;
            mathVar::wGauss[nGauss - 3] = 0.4786286704993665;
            mathVar::wGauss[nGauss - 2] = 0.5688888888888889;
            mathVar::wGauss[nGauss - 1] = 0.4786286704993665;
            mathVar::wGauss[nGauss] = 0.2369268850561891;
        }
        else if (nGauss==5)
        {
            mathVar::xGauss[nGauss - 5] = -0.9324695142031521;
            mathVar::xGauss[nGauss - 4] = -0.6612093864662645;
            mathVar::xGauss[nGauss - 3] = -0.2386191860831969;
            mathVar::xGauss[nGauss - 2] = 0.2386191860831969;
            mathVar::xGauss[nGauss - 1] = 0.6612093864662645;
            mathVar::xGauss[nGauss] = 0.9324695142031521;

            mathVar::wGauss[nGauss - 5] = 0.1713244923791704;
            mathVar::wGauss[nGauss - 4] = 0.3607615730481386;
            mathVar::wGauss[nGauss - 3] = 0.4679139345726910;
            mathVar::wGauss[nGauss - 2] = 0.4679139345726910;
            mathVar::wGauss[nGauss - 1] = 0.3607615730481386;
            mathVar::wGauss[nGauss] = 0.1713244923791704;
        }
	}

	void GaussLobatto(int nGauss)
	{
		if (nGauss == 0)
		{
			mathVar::xGaussLobatto[nGauss] = 0.0;
			mathVar::wGaussLobatto[nGauss] = 2.0;
		}
		else if (nGauss == 1)
		{
			mathVar::xGaussLobatto[nGauss - 1] = -1.0;
			mathVar::xGaussLobatto[nGauss] = 1.0;

			mathVar::wGaussLobatto[nGauss - 1] = 1.0;
			mathVar::wGaussLobatto[nGauss] = 1.0;
		}
		else if (nGauss == 2)
		{
			mathVar::xGaussLobatto[nGauss - 2] = -1.0;
			mathVar::xGaussLobatto[nGauss - 1] = 0.0;
			mathVar::xGaussLobatto[nGauss] = 1.0;

			mathVar::wGaussLobatto[nGauss - 2] = 1.0/3.0;
			mathVar::wGaussLobatto[nGauss - 1] = 4.0/3.0;
			mathVar::wGaussLobatto[nGauss] = 1.0/3.0;
		}
	}

	void basisFc(double a, double b, int elemType)
	{
		switch (elemType)
		{
		case 3:
		{
			for (int i = 0; i <= mathVar::orderElem; i++)
			{
				if (i == 0)
				{
					mathVar::B[i] = 1.0;
				}
				else if (i == 1)
				{
					mathVar::B[i] = (3.0*b + 1.0) / 2.0;
				}
				else if (i == 2)
				{
					mathVar::B[i] = a * (1 - b);
				}
				else if (i == 3)
				{
					mathVar::B[i] = -0.5 + b + 5.0*pow(b, 2) / 2.0;
				}
                else if (i == 4)
                {
                    mathVar::B[i] = a*(1-b)*(1.5+2.5*b);
                }
                else if (i == 5)
                {
                    mathVar::B[i] = (1.5*a*a-0.5)*(1-b)*(1-b);
                }
                else if (i == 6)
                {
                    mathVar::B[i] = -(3.0/8.0)-(15.0/8.0)*b+(15.0/8.0)*b*b+(35.0/8.0)*b*b*b;
                }
                else if (i == 7)
                {
                    mathVar::B[i] = a*(1-b)*(0.25+4.5*b+(21.0/4.0)*b*b);
                }
                else if (i == 8)
                {
                    mathVar::B[i] = (1.5*a*a-0.5*(1-b)*(1-b))*(2.5+3.5*b);
                }
                else if (i == 9)
                {
                    mathVar::B[i] = (2.5*a*a*a-1.5)*pow((1-b),3.0);
                }
			}
		}
		break;
		case 4:
		{
			for (int i = 0; i <= mathVar::orderElem; i++)
			{
				if (i == 0)
				{
					mathVar::B[i] = 1.0;
				}
				else if (i == 1)
				{
					mathVar::B[i] = a;
				}
				else if (i == 2)
				{
					mathVar::B[i] = b;
				}
				else if (i == 3)
				{
					mathVar::B[i] = a * b;
				}
			}
		}
		break;
		default:
			break;
		}
	}

	void dBasisFc(double a, double b, int elemType)
	{
		switch (elemType)
		{
		case 3:
		{
			for (int i = 0; i <= mathVar::orderElem; i++)
			{
				if (i == 0)
				{
					mathVar::dBa[i] = 0;
					mathVar::dBb[i] = 0;
				}
				else if (i == 1)
				{
					mathVar::dBa[i] = 0;
					mathVar::dBb[i] = 3.0 / 2.0;
				}
				else if (i == 2)
				{
					mathVar::dBa[i] = 1 - b;
					mathVar::dBb[i] = -a;
				}
                else if (i == 3)
                {
                    mathVar::dBa[i] = 0;
                    mathVar::dBb[i] = 1+5*b;
                }
                else if (i == 4)
                {
                    mathVar::dBa[i] = 1.5+b-2.5*b*b;
                    mathVar::dBb[i] = a*(1-5*b);
                }
                else if (i == 5)
                {
                    mathVar::dBa[i] = 3*a*(1-b)*(1-b);
                    mathVar::dBb[i] = 2*(b-1)*(1.5*a*a-0.5);
                }
                else if (i == 6)
                {
                    mathVar::dBa[i] = 0;
                    mathVar::dBb[i] = (-15.0/8.0)+(15.0/4.0)*b+(35.0/24.0)*b*b;
                }
                else if (i == 7)
                {
                    mathVar::dBa[i] = (1-b)*(0.25+4.5*b+(21.0/4.0)*b*b);
                    mathVar::dBb[i] = a*(17.0/4.0+3*b/8.0+7*b*b/4.0);
                }
                else if (i == 8)
                {
                    mathVar::dBa[i] = 3*a*(2.5+3.5*b);
                    mathVar::dBb[i] = (21.0/4.0)*(a*a-b*b)+9.0*b/2.0+3.0/4.0;
                }
                else if (i == 9)
                {
                    mathVar::dBa[i] = (15.0/2.0)*a*a*pow(1-b,3.0);
                    mathVar::dBb[i] = -3*pow(1-b,2.0)*(5*a*a*a/2.0-1.5);
                }
			}
		}
		break;
		case 4:
		{
			for (int i = 0; i <= mathVar::orderElem; i++)
			{
				if (i == 0)
				{
					mathVar::dBa[i] = 0.0;
					mathVar::dBb[i] = 0.0;
				}
				else if (i == 1)
				{
					mathVar::dBa[i] = 1.0;
					mathVar::dBb[i] = 0.0;
				}
				else if (i == 2)
				{
					mathVar::dBa[i] = 0.0;
					mathVar::dBb[i] = 1.0;
				}
				else if (i == 3)
				{
					mathVar::dBa[i] = b;
					mathVar::dBb[i] = a;
				}
			}
		}
		break;
		default:
			break;
		}
	}

	double J2DCal(int elem, double a, double b)
	{
		double Jacobi(0.0), xA(0.0), xB(0.0), xC(0.0), xD(0.0), yA(0.0), yB(0.0), yC(0.0), yD(0.0);
		int elemType(auxUlti::checkType(elem));
		if (elemType == 4)  //Quad
		{
			std::tie(xA, yA) = auxUlti::getElemCornerCoord(elem, 0);
			std::tie(xB, yB) = auxUlti::getElemCornerCoord(elem, 1);
			std::tie(xC, yC) = auxUlti::getElemCornerCoord(elem, 2);
			std::tie(xD, yD) = auxUlti::getElemCornerCoord(elem, 3);

			Jacobi = math::jacobianQuad(xA, xB, xC, xD, yA, yB, yC, yD, a, b);
		}
		else if (elemType == 3)  //Tri
		{
			std::tie(xA, yA) = auxUlti::getElemCornerCoord(elem, 0);
			std::tie(xB, yB) = auxUlti::getElemCornerCoord(elem, 1);
			std::tie(xC, yC) = auxUlti::getElemCornerCoord(elem, 2);

			Jacobi = math::jacobianTri(xA, xB, xC, yA, yB, yC, a, b);
		}
		return Jacobi;
	}

	std::tuple<double, double> J1DCal(int edge)
	{
		int master(0), servant(0), masterIndex(0), servantIndex(0), masterType(0), servantType(0);
		double JMaster(0.0), JServant(0.0), xA(0.0), xB(0.0), xC(0.0), xD(0.0), yA(0.0), yB(0.0), yC(0.0), yD(0.0);

		if (meshVar::ineled[0][edge]>meshVar::ineled[1][edge])
		{
			master = meshVar::ineled[0][edge];
			servant = meshVar::ineled[1][edge];
		}
		else if (meshVar::ineled[0][edge]<meshVar::ineled[1][edge])
		{
			master = meshVar::ineled[1][edge];
			servant = meshVar::ineled[0][edge];
		}

		masterIndex = auxUlti::findEdgeOrder(master, edge);
		masterType = auxUlti::checkType(master);

		if (masterType==4)
		{
			std::tie(xA, yA) = auxUlti::getElemCornerCoord(master, 0);
			std::tie(xB, yB) = auxUlti::getElemCornerCoord(master, 1);
			std::tie(xC, yC) = auxUlti::getElemCornerCoord(master, 2);
			std::tie(xD, yD) = auxUlti::getElemCornerCoord(master, 3);

			JMaster = math::jacobian1DQuad(masterIndex, xA, xB, xC, xD, yA, yB, yC, yD);
		}
		else if (masterType==3)
		{
			std::tie(xA, yA) = auxUlti::getElemCornerCoord(master, 0);
			std::tie(xB, yB) = auxUlti::getElemCornerCoord(master, 1);
			std::tie(xC, yC) = auxUlti::getElemCornerCoord(master, 2);

			JMaster = math::jacobian1DTri(masterIndex, xA, xB, xC, yA, yB, yC);
		}

		if (servant >= 0)
		{
			servantIndex = auxUlti::findEdgeOrder(servant, edge);
			servantType = auxUlti::checkType(servant);
			if (servantType == 4)
			{
				std::tie(xA, yA) = auxUlti::getElemCornerCoord(servant, 0);
				std::tie(xB, yB) = auxUlti::getElemCornerCoord(servant, 1);
				std::tie(xC, yC) = auxUlti::getElemCornerCoord(servant, 2);
				std::tie(xD, yD) = auxUlti::getElemCornerCoord(servant, 3);

				JServant = math::jacobian1DQuad(servantIndex, xA, xB, xC, xD, yA, yB, yC, yD);
			}
			else if (servantType == 3)
			{
				std::tie(xA, yA) = auxUlti::getElemCornerCoord(servant, 0);
				std::tie(xB, yB) = auxUlti::getElemCornerCoord(servant, 1);
				std::tie(xC, yC) = auxUlti::getElemCornerCoord(servant, 2);

				JServant = math::jacobian1DTri(servantIndex, xA, xB, xC, yA, yB, yC);
			}
		}
		else
		{
			JServant = 0.0;
		}

		return std::make_tuple(JMaster, JServant);
	}

	double iniIntegral(int elemType, int order)
	{
		double w1(0.0), w2(0.0), integral(0.0);
		switch (elemType)
		{
		case 3:
		{
			for (int na = 0; na <= mathVar::nGauss; na++)
			{
				for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					w1 = mathVar::wGaussPts[na][nb][0];
					w2 = mathVar::wGaussPts[na][nb][1];
					integral += w1 * w2* mathVar::BPts_Tri[order][na][nb];
				}
			}
		}
		break;
		case 4:
		{
			for (int na = 0; na <= mathVar::nGauss; na++)
			{
				for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					w1 = mathVar::wGaussPts[na][nb][0];
					w2 = mathVar::wGaussPts[na][nb][1];
					integral += w1 * w2* mathVar::BPts_Quad[order][na][nb];
				}
			}
		}
		break;
		default:
			break;
		}
		return integral;
	}

	double volumeInte(std::vector< std::vector<double> > &Fvalue, int elem)
	{
        double J2D(0.0), w1(0.0), w2(0.0), integral(0.0);
		for (int na = 0; na <= mathVar::nGauss; na++)
		{
			for (int nb = 0; nb <= mathVar::nGauss; nb++)
			{
				w1 = mathVar::wGaussPts[na][nb][0];
				w2 = mathVar::wGaussPts[na][nb][1];
                J2D = meshVar::J2D[elem][na][nb];
                integral += w1 * w2*(J2D)*Fvalue[na][nb];
			}
		}
		return integral;
	}

	double jacobianQuad(double xA, double xB, double xC, double xD, double yA, double yB, double yC, double yD, double a, double b)
	{
		double dxa(0.0), dxb(0.0), dya(0.0), dyb(0.0);
		double jQuad(0.0);

		dxa = (1.0 / 4.0)*(-xA + xB - xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*b;
		dxb = (1.0 / 4.0)*(-xA - xB + xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*a;
		dya = (1.0 / 4.0)*(-yA + yB - yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*b;
		dyb = (1.0 / 4.0)*(-yA - yB + yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*a;
		jQuad = dxa * dyb - dxb * dya;
		return jQuad;
	}

	double jacobianTri(double xA, double xB, double xC, double yA, double yB, double yC, double a, double b)
	{
		double dxa(0.0), dxb(0.0), dya(0.0), dyb(0.0);
		double jTri(0.0);

		dxa = (1 - b)*(xB - xA) / 4.0;
		dxb = a * (xA - xB) / 4.0 + (-xA - xB + 2 * xC) / 4.0;
		dya = (1 - b)*(yB - yA) / 4.0;
		dyb = a * (yA - yB) / 4.0 + (-yA - yB + 2 * yC) / 4.0;
		jTri = dxa * dyb - dxb * dya;
		return jTri;
	}

	double jacobian1DQuad(int edgeIndex, double xA, double xB, double xC, double xD, double yA, double yB, double yC, double yD)
	{
		double dx(0.0), dy(0.0), C(0.0);
		double jacobi(0.0);
		if ((edgeIndex == 0) || (edgeIndex == 2))  //AB or CD
		{
			if (edgeIndex == 0)
			{
				C = -1.0;
			}
			else if (edgeIndex == 2)
			{
				C = 1.0;
			}
			dx = (1.0 / 4.0)*(-xA + xB - xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*C;
			dy = (1.0 / 4.0)*(-yA + yB - yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*C;
		}
		else if ((edgeIndex == 1) || (edgeIndex == 3))  //BC or DA
		{
			if (edgeIndex == 1)
			{
				C = 1.0;
			}
			else if (edgeIndex == 3)
			{
				C = -1.0;
			}
			dx = (1.0 / 4.0)*(-xA - xB + xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*C;
			dy = (1.0 / 4.0)*(-yA - yB + yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*C;
		}

		jacobi = std::sqrt(pow(dx, 2) + pow(dy, 2));
		return jacobi;
	}

	double jacobian1DTri(int edgeIndex, double xA, double xB, double xC, double yA, double yB, double yC)
	{
		double dx(0.0), dy(0.0), C(0.0);
		double jacobi(0.0);
		if ((edgeIndex == 1) || (edgeIndex == 2))  //BC or CA
		{
			if (edgeIndex == 1)
			{
				C = 1.0;
			}
			else if (edgeIndex == 2)
			{
				C = -1.0;
			}
			dx = C * (xA - xB) / 4.0 + (-xA - xB + 2 * xC) / 4.0;
			dy = C * (yA - yB) / 4.0 + (-yA - yB + 2 * yC) / 4.0;
		}
        else if (edgeIndex == 0)  //AB
		{
			C = -1.0;
			dx = (1 - C)*(xB - xA) / 4.0;
			dy = (1 - C)*(yB - yA) / 4.0;
		}

		jacobi = std::sqrt(pow(dx, 2) + pow(dy, 2));
		return jacobi;
	}

	std::vector<double> SolveSysEqs(std::vector< std::vector<double> > &a, std::vector<double> &b)
	{
		int n = static_cast<int>(b.size());
		double eMax(1e-9), e(1.0), sum(0.0), xi(0.0);
		std::vector<double> results(n, 1.0);
		int counter(0);
		//eMax = fabs(*std::min_element(b.begin(), b.end())) / 1e3;

		while (e>eMax && counter<=15)
		{
			for (int i = 0; i < n; i++)
			{
				sum = b[i];
				for (int j = 0; j < n; j++)
				{
					if (i != j)
					{
						sum -= a[i][j] * results[j];
					}
				}
				xi = sum / a[i][i];
				results[i] = xi;
			}
			e = math::errorGS(a, b, results);
			counter++;
		}
		return results;
	}

	double errorGS(std::vector< std::vector<double> > &a, std::vector<double> &b, std::vector<double> &No)
	{
		/* Calculate a*No-b */
		int n = static_cast<int>(b.size());
        double error(1.0);
		std::vector<double> R(n, 0.0);

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				R[i] += a[i][j] * No[j];
			}
			R[i] = fabs(R[i] - b[i]);
		}
		error = *std::max_element(R.begin(), R.end());  //find max value of vector
		return error;
	}

	double CalcVisCoef(double T)
	{
		double mu(0.0);
		mu = material::As*pow(T, 1.5) / (T + material::Ts);
		return mu;
	}

	double CalcT(int elem, double a, double b)
	{
		double T(0.0), rhoGs(0.0), rhouGs(0.0), rhovGs(0.0), rhoEGs(0.0);
		rhoGs = math::pointValue(elem, a, b, 1, 2);
		rhouGs = math::pointValue(elem, a, b, 2, 2);
		rhovGs = math::pointValue(elem, a, b, 3, 2);
		rhoEGs = math::pointValue(elem, a, b, 4, 2);
		T = (material::gamma - 1)*(rhoEGs - 0.5*(pow(rhouGs, 2) + pow(rhovGs, 2)) / rhoGs) / (material::R*rhoGs);
		return T;
	}

	double CalcTFromConsvVar(double rho, double rhou, double rhov, double rhoE)
	{
		double T((material::gamma - 1)*(rhoE - 0.5*(pow(rhou, 2) + pow(rhov, 2)) / rho) / (material::R*rho));
		if ((T <= 0) && (fabs(T) < 0.001))
		{
			//std::cout << "Warning!!! limiting T " << T <<std::endl;
			//T = limitVal::TDwn;
			//limitVal::limitTOrNot = true;
			//system("pause");
			T = fabs(T);
		}
		return T;
	}

	double CalcP(double T, double rho)
	{
		double p = T * (material::R*rho);
		return p;
	}

	std::tuple<double, double, double, double> Calc_dxydab(int elem, double a, double b)
	{
		double dxa(0.0), dxb(0.0), dya(0.0), dyb(0.0);
		double xA(0.0), xB(0.0), xC(0.0), xD(0.0), yA(0.0), yB(0.0), yC(0.0), yD(0.0);
		int elemType(0);

		elemType = auxUlti::checkType(elem);
		std::tie(xA, yA) = auxUlti::getElemCornerCoord(elem, 0);
		std::tie(xB, yB) = auxUlti::getElemCornerCoord(elem, 1);
		std::tie(xC, yC) = auxUlti::getElemCornerCoord(elem, 2);
		
		if (elemType==4)  //Quad
		{
			std::tie(xD, yD) = auxUlti::getElemCornerCoord(elem, 3);
            dxa = (1.0 / 4.0)*(-xA + xB - xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*b;
            dxb = (1.0 / 4.0)*(-xA - xB + xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*a;
            dya = (1.0 / 4.0)*(-yA + yB - yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*b;
            dyb = (1.0 / 4.0)*(-yA - yB + yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*a;
		}
		else if (elemType == 3)  //Tri
		{
            dxa = (1 - b)*(xB - xA) / 4.0;
            dxb = a * (xA - xB) / 4.0 + (-xA - xB + 2 * xC) / 4.0;
            dya = (1 - b)*(yB - yA) / 4.0;
            dyb = a * (yA - yB) / 4.0 + (-yA - yB + 2 * yC) / 4.0;
		}
		
		return std::make_tuple(dxa, dxb, dya, dyb);
	}

	double Calc_dBxdBy(int elem, int order, int na, int nb, int opt)
	{
		double dB(0.0);
		int elemType(auxUlti::checkType(elem));
		switch (elemType)
		{
		case 3:
		{
			if (opt == 1)  //x direction
			{
				dB = (1 / meshVar::J2D[elem][na][nb]) * (mathVar::dBaPts_Tri[order][na][nb] * meshVar::dyb[elem][na][nb] - mathVar::dBbPts_Tri[order][na][nb] * meshVar::dya[elem][na][nb]);
			}
			else if (opt == 2)  //y direction
			{
				dB = (1 / meshVar::J2D[elem][na][nb]) * (mathVar::dBbPts_Tri[order][na][nb] * meshVar::dxa[elem][na][nb] - mathVar::dBaPts_Tri[order][na][nb] * meshVar::dxb[elem][na][nb]);
			}
		}
		break;
		case 4:
		{
			if (opt == 1)  //x direction
			{
				dB = (1 / meshVar::J2D[elem][na][nb]) * (mathVar::dBaPts_Quad[order][na][nb] * meshVar::dyb[elem][na][nb] - mathVar::dBbPts_Quad[order][na][nb] * meshVar::dya[elem][na][nb]);
			}
			else if (opt == 2)  //y direction
			{
				dB = (1 / meshVar::J2D[elem][na][nb]) * (mathVar::dBbPts_Quad[order][na][nb] * meshVar::dxa[elem][na][nb] - mathVar::dBaPts_Quad[order][na][nb] * meshVar::dxb[elem][na][nb]);
			}
		}
		break;
		default:
			break;
		}
		return dB;
	}

	double surfaceInte(std::vector<double> &Fvalue, int edge, int elem)
	{
        double inte(0.0), J(meshVar::J1D[edge]), w(0.0);
		for (int nG = 0; nG <= mathVar::nGauss; nG++)
		{
			w = mathVar::wGauss[nG];
			inte += w * Fvalue[nG] * J;
		}
		return inte;
	}

	double pointValue(int element, double a, double b, int valType, int valKind)
	{
		//Compute primary variables from conservative variables
		double out(0.0);
		if (valKind==1)  //primary variables
		{
			if (valType == 1)  //rho
			{
				out= math::calcConsvVarWthLimiter(element, a, b, valType);
			}
			else if (valType == 2)  //u
			{
				double rhoVal(math::calcConsvVarWthLimiter(element, a, b, 1)),
					rhouVal(math::calcConsvVarWthLimiter(element, a, b, 2));
				out = rhouVal / rhoVal;
			}
			else if (valType == 3)  //v
			{
				double rhoVal(math::calcConsvVarWthLimiter(element, a, b, 1)),
					rhovVal(math::calcConsvVarWthLimiter(element, a, b, 3));
				out = rhovVal / rhoVal;
			}
			else if (valType == 4)  //e
			{
				double rhoVal(math::calcConsvVarWthLimiter(element, a, b, 1)),
					rhouVal(math::calcConsvVarWthLimiter(element, a, b, 2)),
					rhovVal(math::calcConsvVarWthLimiter(element, a, b, 3)),
					rhoEVal(math::calcConsvVarWthLimiter(element, a, b, 4));
				out = material::Cv*math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
			}
			else if (valType == 5)  //p
			{
				double rhoVal(math::calcConsvVarWthLimiter(element, a, b, 1)),
					rhouVal(math::calcConsvVarWthLimiter(element, a, b, 2)),
					rhovVal(math::calcConsvVarWthLimiter(element, a, b, 3)),
					rhoEVal(math::calcConsvVarWthLimiter(element, a, b, 4));
				double TVal(math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal));
				out = math::CalcP(TVal, rhoVal);
			}
			else if (valType == 6)  //T
			{
				double rhoVal(math::calcConsvVarWthLimiter(element, a, b, 1)),
					rhouVal(math::calcConsvVarWthLimiter(element, a, b, 2)),
					rhovVal(math::calcConsvVarWthLimiter(element, a, b, 3)),
					rhoEVal(math::calcConsvVarWthLimiter(element, a, b, 4));
				out = math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
				if ((out < 0) && (fabs(out) < 0.001))
				{
					std::cout << "Negative T" << out << " at cell " << element + meshVar::nelem1D + 1 << std::endl;
					//system("pause");
					std::cout << theta1Arr[element] << std::endl;
					std::cout << theta2Arr[element] << std::endl;
					std::cout << a << ", " << b << std::endl;
					out = fabs(out);
				}
			}
			else if (valType == 7)  //mu
			{
				double rhoVal(math::calcConsvVarWthLimiter(element, a, b, 1)),
					rhouVal(math::calcConsvVarWthLimiter(element, a, b, 2)),
					rhovVal(math::calcConsvVarWthLimiter(element, a, b, 3)),
					rhoEVal(math::calcConsvVarWthLimiter(element, a, b, 4));
				double TVal(math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal));
				if ((TVal<0) && fabs(TVal) < 0.001)
				{
					TVal = fabs(TVal);
				}
				out = math::CalcVisCoef(TVal);
				if (out < 0 || out != out)
				{
					std::cout << "unphysical mu at cell " << element + meshVar::nelem1D + 1 << std::endl;
					std::cout << "TVal " << TVal << std::endl;
					std::cout << "rhoVal " << rhoVal << std::endl;
					std::cout << "rhouVal " << rhouVal << std::endl;
					std::cout << "rhovVal " << rhovVal << std::endl;
					std::cout << "rhoEVal " << rhoEVal << std::endl;
					std::cout << theta1Arr[element] << std::endl;
					std::cout << theta2Arr[element] << std::endl;
					std::cout << a << ", " << b << std::endl;
					system("pause");
				}
			}
		}
		else if (valKind==2)  //conservative variables
		{	
			out = math::calcConsvVarWthLimiter(element, a, b, valType);
		}
		return out;
	}

	double pointValueNoLimiter(int element, double a, double b, int valType)
	{
		double out(0.0);
		std::vector<double> Value(mathVar::orderElem + 1, 0.0);
		Value = auxUlti::getElementConserValuesOfOrder(element, valType);

		math::basisFc(a, b, auxUlti::checkType(element));
		for (int order = 0; order <= mathVar::orderElem; order++)
		{
			out += Value[order] * mathVar::B[order];
		}

		return out;
	}

	double vectorDotProduct(std::vector<double> &a, std::vector<double> &b)
	{
		double out(0.0);
		out = a[0] * b[0] + a[1] * b[1];
		return out;
	}

	std::tuple <double, double> internalSurfaceValue(int edge, int element, int nG, int valType, int valKind)
	{
		int masterElem(0), servantElem(0);
		std::tie(masterElem, servantElem) = auxUlti::getMasterServantOfEdge(edge);
		double valPlus(0.0), valMinus(0.0), aMaster(0.0), bMaster(0.0), aServant(0.0), bServant(0.0);

		std::tie(aMaster, bMaster) = auxUlti::getGaussSurfCoor(edge, masterElem, nG);
		std::tie(aServant, bServant) = auxUlti::getGaussSurfCoor(edge, servantElem, nG);

		if (masterElem==element)  //considering element is master
		{
			valPlus = math::pointValue(masterElem, aMaster, bMaster, valType, valKind);
			valMinus = math::pointValue(servantElem, aServant, bServant, valType, valKind);
		}
		else
		{
			valMinus = math::pointValue(masterElem, aMaster, bMaster, valType, valKind);
			valPlus = math::pointValue(servantElem, aServant, bServant, valType, valKind);
		}

		return std::make_tuple(valPlus, valMinus);
	}

	std::tuple <double, double> internalSurfaceDerivativeValue(int edge, int element, int nG, int valType, int dir)
	{
		int masterElem(0), servantElem(0);
		std::tie(masterElem, servantElem) = auxUlti::getMasterServantOfEdge(edge);
		double valPlus(0.0), valMinus(0.0), aMaster(0.0), bMaster(0.0), aServant(0.0), bServant(0.0);

		std::tie(aMaster, bMaster) = auxUlti::getGaussSurfCoor(edge, masterElem, nG);
		std::tie(aServant, bServant) = auxUlti::getGaussSurfCoor(edge, servantElem, nG);

		if (masterElem == element)  //considering element is master
		{
			valPlus = math::pointAuxValue(masterElem, aMaster, bMaster, valType, dir);
			valMinus = math::pointAuxValue(servantElem, aServant, bServant, valType, dir);
		}
		else
		{
			valMinus = math::pointAuxValue(masterElem, aMaster, bMaster, valType, dir);
			valPlus = math::pointAuxValue(servantElem, aServant, bServant, valType, dir);
		}

		return std::make_tuple(valPlus, valMinus);
	}

	double SurfaceValueFromMaster(int edge, int nG, int valType, int valKind)
	{
		int masterElem(0), servantElem(0);
		std::tie(masterElem, servantElem) = auxUlti::getMasterServantOfEdge(edge);
		double valPlus(0.0), aMaster(0.0), bMaster(0.0);
		std::tie(aMaster, bMaster) = auxUlti::getGaussSurfCoor(edge, masterElem, nG);
		valPlus = math::pointValue(masterElem, aMaster, bMaster, valType, valKind);
		return valPlus;
	}

	double CalcSpeedOfSound(double T)
	{
		double C(sqrt(material::gamma*material::R*T));
		return C;
	}

	std::vector<double> vectorSum(std::vector<double> &a, std::vector<double> &b)
	{
		static int size(a.size());
		std::vector<double> out(size, 0.0);
		for (int i = 0; i < size; i++)
		{
			out[i] = a[i] + b[i];
		}
		return out;
	}

	double calcRhouvEDeriv(double dRhoUVal, double dRhoVal, double rhoUVal, double rhoVal)
	{
		/*NOTE! direction of derivatives based on input values, for example:
		d(u)/dx = [d(rhou)/dx - d(rho)/dx.(rhou)/rho]/(rho)
		so input values are
		d(rho)/dx	----> dRhoVal
		d(rhou)/dx	----> dRhoUVal
		rhou		----> rhoUVal
		rho			----> rhoVal*/
		double outVal(0.0);
		outVal = (dRhoUVal - dRhoVal * rhoUVal / rhoVal) / rhoVal;
		return outVal;
	}

	double calcTDeriv(double dE, double du, double dv, double u, double v)
	{
		/*NOTE! direction of derivatives based on input values, for example:
		d(T)/dx = f(dE/dx, du/dx, dv/dx, u, v)
		so input values are
		dE/dx	----> dE
		du/dx	----> du
		dv/dx	----> dv*/
		double outVal(0.0);
		outVal = (material::gamma - 1)*(dE - (u*du + v * dv)) / material::R;
		return outVal;
	}

	double pointAuxValue(int element, double a, double b, int valType, int dir)
	{
		double out(0.0);
		std::vector<double> Value(mathVar::orderElem + 1, 0.0);

		Value = auxUlti::getElementAuxValuesOfOrder(element, valType, dir);

		math::basisFc(a, b, auxUlti::checkType(element));
		for (int order = 0; order <= mathVar::orderElem; order++)
		{
			out += Value[order] * mathVar::B[order];
		}
		//out = out / muVal;
		return out;
	}

	double calcThermalConductivity(double muVal)
	{
		double k(0.0);
		k = material::Cp*muVal / material::Pr;
		return k;
	}

	std::tuple<bool, double, double> solvQuadraticEq(double A, double B, double C)
	{
        double delta(B*B - 4 * A*C), root1(0.0), root2(0.0);
		bool realRoot(true);

        if (A != 0.0)
		{
			if (delta>0)
			{
				root1 = ((-B + sqrt(delta)) / (2 * A));
				root2 = ((-B - sqrt(delta)) / (2 * A));
			}
            else if (delta == 0.0)
			{
				root1 = (-B / (2 * A));
				root2 = root1;
			}
			else
			{
				realRoot = false;
			}
		}
        else if (A == 0.0 && B !=0.0)
		{
			root1 = -C / B;
			root2 = root1;
		}
		else
		{
			realRoot = false;
		}
		return std::make_tuple(realRoot, root1, root2);
	}

	double centerValue(int element, int valType, int valKind)
	{
		double xC(-1.0 / 3.0), yC(1.0 / 3.0), output(0.0);

		if (auxUlti::checkType(element) == 4)
		{
			xC = 0.0;
			yC = 0.0;
		}
		output = math::pointValue(element, xC, yC, valType, valKind);
		return output;
	}

	double centerAuxValue(int element, int valType, int dir)
	{
		double xC(-1.0 / 3.0), yC(1.0 / 3.0), output(0.0);

		if (auxUlti::checkType(element) == 4)
		{
			xC = 0.0;
			yC = 0.0;
		}
		output = math::pointAuxValue(element, xC, yC, valType, dir);
		return output;
	}

	std::tuple<double, double> directMapping(int element, double aCoor, double bCoor)
	{
		int elemType(auxUlti::checkType(element));
		double xCoor(0.0), yCoor(0.0), xA(0.0), xB(0.0), xC(0.0),
			yA(0.0), yB(0.0), yC(0.0);

		std::tie(xA, yA) = auxUlti::getElemCornerCoord(element, 0);
		std::tie(xB, yB) = auxUlti::getElemCornerCoord(element, 1);
		std::tie(xC, yC) = auxUlti::getElemCornerCoord(element, 2);

		if (elemType == 3) //Tri element
		{
			xCoor = 0.25*(1 - aCoor)*(1 - bCoor)*xA + 0.25*(1 + aCoor)*(1 - bCoor)*xB + 0.5*(1 + bCoor)*xC;
			yCoor = 0.25*(1 - aCoor)*(1 - bCoor)*yA + 0.25*(1 + aCoor)*(1 - bCoor)*yB + 0.5*(1 + bCoor)*yC;
		}
		else if (elemType == 4) //Quad element
		{
			double xD(0.0), yD(0.0);
			std::tie(xD, yD) = auxUlti::getElemCornerCoord(element, 3);

			xCoor = 0.25*(1 - aCoor)*(1 - bCoor)*xA + 0.25*(1 + aCoor)*(1 - bCoor)*xB + 0.25*(1 - aCoor)*(1 + bCoor)*xD + 0.25*(1 + aCoor)*(1 + bCoor)*xC;
			yCoor = 0.25*(1 - aCoor)*(1 - bCoor)*yA + 0.25*(1 + aCoor)*(1 - bCoor)*yB + 0.25*(1 - aCoor)*(1 + bCoor)*yD + 0.25*(1 + aCoor)*(1 + bCoor)*yC;
		}
		return std::make_tuple(xCoor, yCoor);
	}

	std::tuple<double, double> mappingRealToStd(int edge, int element, double xCoor, double yCoor)
	{
        double xA(0.0), xB(0.0), xC(0.0),
			yA(0.0), yB(0.0), yC(0.0),
			aCoor(0.0), bCoor(0.0),
			xCal(0.0), yCal(0.0), eX(100.0), eY(100.0);
		int elemType(auxUlti::checkType(element));
		std::vector<std::vector<double>> vectorGaussPoints(mathVar::nGauss + 1, std::vector<double>(2, 0.0));

		std::tie(xA, yA) = auxUlti::getElemCornerCoord(element, 0);
		std::tie(xB, yB) = auxUlti::getElemCornerCoord(element, 1);
		std::tie(xC, yC) = auxUlti::getElemCornerCoord(element, 2);

		vectorGaussPoints = auxUlti::getVectorGaussSurfCoor(edge, element);

		if (elemType==3) //tri
		{
			for (int nG = 0; nG <= mathVar::nGauss; nG++)
			{
				aCoor = vectorGaussPoints[nG][0];
				bCoor = vectorGaussPoints[nG][1];
				xCal = 0.25*(1 - aCoor)*(1 - bCoor)*xA + 0.25*(1 + aCoor)*(1 - bCoor)*xB + 0.5*(1 + bCoor)*xC;
				yCal = 0.25*(1 - aCoor)*(1 - bCoor)*yA + 0.25*(1 + aCoor)*(1 - bCoor)*yB + 0.5*(1 + bCoor)*yC;
				eX = fabs(xCal - xCoor) * 100 / xCoor;
				eY = fabs(yCal - yCoor) * 100 / yCoor;
				if ((eX < 0.005) && (eY < 0.005))
				{
					break;
				}
			}
		}
		else if (elemType==4) //quad
		{
			double xD(0.0), yD(0.0);
			std::tie(xD, yD) = auxUlti::getElemCornerCoord(element, 3);

			for (int nG = 0; nG <= mathVar::nGauss; nG++)
			{
				aCoor = vectorGaussPoints[nG][0];
				bCoor = vectorGaussPoints[nG][1];
				xCal = 0.25*(1 - aCoor)*(1 - bCoor)*xA + 0.25*(1 + aCoor)*(1 - bCoor)*xB + 0.25*(1 - aCoor)*(1 + bCoor)*xD + 0.25*(1 + aCoor)*(1 + bCoor)*xC;
				yCal = 0.25*(1 - aCoor)*(1 - bCoor)*yA + 0.25*(1 + aCoor)*(1 - bCoor)*yB + 0.25*(1 - aCoor)*(1 + bCoor)*yD + 0.25*(1 + aCoor)*(1 + bCoor)*yC;
				eX = fabs((xCal - xCoor) * 100 / xCoor);
				eY = fabs((yCal - yCoor) * 100 / yCoor);
				if ((eX < 0.05) && (eY < 0.05))
				{
					break;
				}
			}
		}

		return std::make_tuple(aCoor, bCoor);
	}

	std::tuple<double, double> inverseMapping(int element, double xCoor, double yCoor)
	{
		int elemType(auxUlti::checkType(element));
		double aCoor(0.0), bCoor(0.0), xA(0.0), xB(0.0), xC(0.0),
			yA(0.0), yB(0.0), yC(0.0),
			Ax(0.0), Bx(0.0), Cx(0.0), Dx(0.0),
			Ay(0.0), By(0.0), Cy(0.0), Dy(0.0), AA(0.0), BB(0.0), CC(0.0);
		bool realRoot(true);

		std::tie(xA, yA) = auxUlti::getElemCornerCoord(element, 0);
		std::tie(xB, yB) = auxUlti::getElemCornerCoord(element, 1);
		std::tie(xC, yC) = auxUlti::getElemCornerCoord(element, 2);

		if (elemType == 3) //Tri element
		{
			Ax = (xA - xB) / 4.0;
			Bx = (xB - xA) / 4.0;
			Cx = (2 * xC - xA - xB) / 4.0;
			Dx = (2 * xC + xA + xB) / 4.0 - xCoor;
			Ay = (yA - yB) / 4.0;
			By = (yB - yA) / 4.0;
			Cy = (2 * yC - yA - yB) / 4.0;
			Dy = (2 * yC + yA + yB) / 4.0 - yCoor;
		}
		else if (elemType == 4) //Quad element
		{
			double xD(0.0), yD(0.0);
			std::tie(xD, yD) = auxUlti::getElemCornerCoord(element, 3);

			Ax = (xA - xB - xD + xC) / 4.0;
			Bx = (- xA + xB - xD + xC) / 4.0;
			Cx = (- xA - xB + xD + xC) / 4.0;
			Dx = (xA + xB + xD + xC) / 4.0 - xCoor;
			Ay = (yA - yB - yD + yC) / 4.0;
			By = (-yA + yB - yD + yC) / 4.0;
			Cy = (-yA - yB + yD + yC) / 4.0;
			Dy = (yA + yB + yD + yC) / 4.0 - yCoor;
		}

		//solve a
		AA = -Ax * By + Ay * Bx;
		BB = -Ax * Dy + Bx * Cy - Cx * By + Ay * Dx;
		CC = Dx * Cy - Cx * Dy;
		if (fabs(AA) < 1e-12 && fabs(BB) < 1e-12 && fabs(CC) < 1e-12)
		{
			aCoor = -1;
			bCoor = 1;
		}
		else
		{
			std::vector<double> vectorRoot_a(2, 0.0);
			std::tie(realRoot, vectorRoot_a[0], vectorRoot_a[1]) = math::solvQuadraticEq(AA, BB, CC);
			if (realRoot)
			{
				for (int i = 0; i < 2; i++)
				{
					double error(fabs(fabs(vectorRoot_a[i]) - 1) * 100);
					if ((fabs(vectorRoot_a[i]) <= 1) || (error <= 0.00001))
					{
						aCoor = vectorRoot_a[i];
						break;
					}
				}
			}
			else
			{
				/*
				std::cout << "Element " << element << std::endl;
				std::cout << "A: " << AA << " B: " << BB << " C: " << CC << std::endl;
				std::cout << "x: " << xCoor << " y: " << yCoor << std::endl;
				*/
				std::string str("mapping error occured_a");
				exitDG(str);
			}
			
			if (fabs(aCoor*Ay + Cy) <= 1e-12)
			{
				//std::cout << "aCoor: " << aCoor << " criteria: " << fabs((aCoor*By + Dy) - (aCoor*Ay + Cy)) << std::endl;
				//bCoor = 2;
				AA = -Ay * Bx + Ax * By;
				BB = -Ay * Dx + By * Cx - Cy * Bx + Ax * Dy;
				CC = Dy * Cx - Cy * Dx;
				if (fabs(AA) < 1e-12 && fabs(BB) < 1e-12 && fabs(CC) < 1e-12)
				{
					aCoor = -1;
					bCoor = 1;
				}
				else
				{
					std::vector<double> vectorRoot_a(2, 0.0);
					std::tie(realRoot, vectorRoot_a[0], vectorRoot_a[1]) = math::solvQuadraticEq(AA, BB, CC);
					if (realRoot)
					{
						for (int i = 0; i < 2; i++)
						{
							double error(fabs(fabs(vectorRoot_a[i]) - 1) * 100);
							if ((fabs(vectorRoot_a[i]) <= 1) || (error <= 0.00001))
							{
								aCoor = vectorRoot_a[i];
								break;
							}
						}
					}
					else
					{
						/*
						std::cout << "Element " << element << std::endl;
						std::cout << "A: " << AA << " B: " << BB << " C: " << CC << std::endl;
						std::cout << "x: " << xCoor << " y: " << yCoor << std::endl;
						*/
						std::string str("mapping error occured");
						exitDG(str);
					}

					if (fabs(aCoor*Ax + Cx) == 0)
					{
						std::string str("mapping error occured");
						exitDG(str);
					}
					else
					{
						bCoor = -(aCoor*Bx + Dx) / (aCoor*Ax + Cx);
					}
				}
			}
			else
			{
				bCoor = -(aCoor*By + Dy) / (aCoor*Ay + Cy);
			}
		}
		
		return std::make_tuple(aCoor, bCoor);
	}

	std::vector<double> tensorVectorDotProduct(std::vector<std::vector<double>> tensor, std::vector<double> vector)
	{
		int size(vector.size());
		std::vector<double> product(size, 0.0);
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				product[i] += tensor[i][j] * vector[j];
			}
		}
		return product;
	}

	double vectorNorm(std::vector<double> vector)
	{
		int size(vector.size());
		double norm(0.0);
		for (int i = 0; i < size; i++)
		{
			norm += pow(vector[i], 2);
		}
		norm = pow(norm, 0.5);
		return norm;
	}

	namespace numericalFluxes
	{
		double auxFlux(double MinusVal, double PlusVar, double vectorComp)
		{
			/*use central numerical flux*/
			double flux(0.5*(MinusVal + PlusVar)*vectorComp);
			return flux;
		}

		double advectiveFlux(double FPlus, double FMinus, double UPlus, double UMinus, double C, double vectorComp)
		{
			/*use Lax - Friedrich numerical flux*/
			/*Plus is inside, minus is outside*/
            double flux(0.5*((FPlus + FMinus)*vectorComp - 0.5 * C * (UMinus - UPlus)));
			return flux;
		}

		double diffusiveFlux(double FPlus, double FMinus, double UPlus, double UMinus, double Beta, double vectorComp)
		{
			/*use central numerical flux*/
			Beta = 0.0;
			double flux(0.5*((FPlus + FMinus)*vectorComp + Beta * (UMinus - UPlus)));
			return flux;
		}

		std::vector<std::vector<double>> NSFEqAdvDiffFluxFromConserVars(int edge, std::vector<double> &UPlus, std::vector<double> &UMinus, std::vector<double> &dUXPlus, std::vector<double> &dUXMinus, std::vector<double> &dUYPlus, std::vector<double> &dUYMinus, std::vector<double> &normVector)
		{
			/*Fluxes array has the following form:
			- column 0: advective fluxes
			- column 1: diffusive fluxes*/
			std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));

			/*StressHeat matrix has form:
			[tauXx		tauXy		Qx]
			[tauYx		tauYy		Qy]
			*/
			std::vector<std::vector<double>> StressHeatP(2, std::vector<double>(3, 0.0));
			std::vector<std::vector<double>> StressHeatM(2, std::vector<double>(3, 0.0));

			double rhoPlus(UPlus[0]), rhouPlus(UPlus[1]), rhovPlus(UPlus[2]), rhoEPlus(UPlus[3]),
				rhoMinus(UMinus[0]), rhouMinus(UMinus[1]), rhovMinus(UMinus[2]), rhoEMinus(UMinus[3]),
				nx(normVector[0]), ny(normVector[1]);

			double
				uPlus(rhouPlus / rhoPlus),
				uMinus(rhouMinus / rhoMinus),

				vPlus(rhovPlus / rhoPlus),
				vMinus(rhovMinus / rhoMinus),

				totalEPlus(rhoEPlus / rhoPlus),
				totalEMinus(rhoEMinus / rhoMinus),

				TPlus(0.0),
				TMinus(0.0),

				pPlus(0.0),
				pMinus(0.0);

				//muPlus(0.0),
				//muMinus(0.0);

			double
				termX1P(0.0), termX1M(0.0),  //(rho*u)					or 0
				termX2P(0.0), termX2M(0.0),  //(rho*u^2 + p)			or tauxx
				termX3P(0.0), termX3M(0.0),  //(rho*u*v)				or tauxy
				termX4P(0.0), termX4M(0.0),  //(rho*totalE + p)*u		or tauxx*u + tauxy*v + Qx

				termY1P(0.0), termY1M(0.0),  //(rho*v)					or 0
				termY2P(0.0), termY2M(0.0),  //(rho*u*v)				or tauxy
				termY3P(0.0), termY3M(0.0),  //(rho*v^2 + p)			or tauyy
				termY4P(0.0), termY4M(0.0);  //(rho*totalE + p)*v		or tauxy*u + tauyy*v + Qy

			double C(0.0), Beta(0.0),
				uMagP(0.0),
				uMagM(0.0),
				aP(0.0),
				aM(0.0);

			/*INVISCID TERMS*/
			/*calculate velocity magnitude*/
			uMagP = sqrt(pow(uPlus, 2) + pow(vPlus, 2));
			uMagM = sqrt(pow(uMinus, 2) + pow(vMinus, 2));

			/*calculate T and P*/
			TPlus = math::CalcTFromConsvVar(rhoPlus, rhouPlus, rhovPlus, rhoEPlus);
			TMinus = math::CalcTFromConsvVar(rhoMinus, rhouMinus, rhovMinus, rhoEMinus);
			pPlus = math::CalcP(TPlus, rhoPlus);
			pMinus = math::CalcP(TMinus, rhoMinus);
			//muPlus = math::CalcVisCoef(TPlus);
			//muMinus = math::CalcVisCoef(TMinus);

			/*calculate speed of sound*/
			aP = math::CalcSpeedOfSound(TPlus);
			aM = math::CalcSpeedOfSound(TMinus);

			if (auxUlti::getBCType(edge) == 0)
			{
				C = LxFConst[edge];
			}
			else
			{
				C = math::numericalFluxes::constantC(uMagP, uMagM, aP, aM);
			}

			/*calculate inviscid terms*/
			std::tie(termX1P, termX2P, termX3P, termX4P) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoPlus, uPlus, vPlus, totalEPlus, pPlus, 1);
			std::tie(termY1P, termY2P, termY3P, termY4P) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoPlus, uPlus, vPlus, totalEPlus, pPlus, 2);

			std::tie(termX1M, termX2M, termX3M, termX4M) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoMinus, uMinus, vMinus, totalEMinus, pMinus, 1);
			std::tie(termY1M, termY2M, termY3M, termY4M) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoMinus, uMinus, vMinus, totalEMinus, pMinus, 2);

			/*Calculate fluxes*/
			Fluxes[0][0] = math::numericalFluxes::advectiveFlux(termX1P, termX1M, rhoPlus, rhoMinus, C, nx) + math::numericalFluxes::advectiveFlux(termY1P, termY1M, rhoPlus, rhoMinus, C, ny);
			Fluxes[1][0] = math::numericalFluxes::advectiveFlux(termX2P, termX2M, rhouPlus, rhouMinus, C, nx) + math::numericalFluxes::advectiveFlux(termY2P, termY2M, rhouPlus, rhouMinus, C, ny);
			Fluxes[2][0] = math::numericalFluxes::advectiveFlux(termX3P, termX3M, rhovPlus, rhovMinus, C, nx) + math::numericalFluxes::advectiveFlux(termY3P, termY3M, rhovPlus, rhovMinus, C, ny);
			Fluxes[3][0] = math::numericalFluxes::advectiveFlux(termX4P, termX4M, rhoEPlus, rhoEMinus, C, nx) + math::numericalFluxes::advectiveFlux(termY4P, termY4M, rhoEPlus, rhoEMinus, C, ny);

			/*VISCOUS TERMS*/
			/*calculate inviscid terms*/
			StressHeatP = math::viscousTerms::calcStressTensorAndHeatFlux(UPlus, dUXPlus, dUYPlus);
			StressHeatM = math::viscousTerms::calcStressTensorAndHeatFlux(UMinus, dUXMinus, dUYMinus);
			/*
			if (auxUlti::getBCType(edge) == 0)
			{
				Beta=DiffusiveFluxConst[edge];
			}
			else
			{
				double eP(material::Cv*TPlus), eM(material::Cv*TMinus);
				std::vector<double> nP(2, 0.0);
				nP[0] = nx;
				nP[0] = ny;
				Beta = math::numericalFluxes::constantBeta(uMagP, uMagM, UPlus[0], UMinus[0], eP, eM, pPlus, pMinus, StressHeatP, StressHeatM, nP);
			}
			*/
			std::tie(termX1P, termX2P, termX3P, termX4P) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uPlus, vPlus, 1);
			std::tie(termY1P, termY2P, termY3P, termY4P) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uPlus, vPlus, 2);
			std::tie(termX1M, termX2M, termX3M, termX4M) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uMinus, vMinus, 1);
			std::tie(termY1M, termY2M, termY3M, termY4M) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uMinus, vMinus, 2);

			/*Calculate fluxes*/
			Fluxes[0][1] = math::numericalFluxes::diffusiveFlux(termX1M, termX1P, rhoPlus, rhoMinus, Beta, nx) + math::numericalFluxes::diffusiveFlux(termY1M, termY1P, rhoPlus, rhoMinus, Beta, ny);
			Fluxes[1][1] = math::numericalFluxes::diffusiveFlux(termX2M, termX2P, rhouPlus, rhouMinus, Beta, nx) + math::numericalFluxes::diffusiveFlux(termY2M, termY2P, rhouPlus, rhouMinus, Beta, ny);
			Fluxes[2][1] = math::numericalFluxes::diffusiveFlux(termX3M, termX3P, rhovPlus, rhovMinus, Beta, nx) + math::numericalFluxes::diffusiveFlux(termY3M, termY3P, rhovPlus, rhovMinus, Beta, ny);
			Fluxes[3][1] = math::numericalFluxes::diffusiveFlux(termX4M, termX4P, rhoEPlus, rhoEMinus, Beta, nx) + math::numericalFluxes::diffusiveFlux(termY4M, termY4P, rhoEPlus, rhoEMinus, Beta, ny);

			return Fluxes;
		}

		double constantC(double uMagP, double uMagM, double aP, double aM)
		{
			std::vector<double> CArray(2, 0.0);
			double C(0.0);
			CArray[0] = uMagP + aP;
			CArray[1] = uMagM + aM;
			C = *std::max_element(CArray.begin(), CArray.end());  //find max value of vector
			return C;
		}

		double constantBeta(double uMagP, double uMagM, double rhoP, double rhoM, double eP, double eM, double pP, double pM, std::vector<std::vector<double>> stressHeatFluxP, std::vector<std::vector<double>> stressHeatFluxM, std::vector<double> nP)
		{
			std::vector<std::vector<double>> stressP(2, std::vector<double>(2, 0.0)), stressM(2, std::vector<double>(2, 0.0));
			std::vector<double> heatP(2, 0.0), heatM(2, 0.0),
				productStressNP(2, 0.0), productStressNM(2, 0.0), nM(2, 0.0), pMultiN_P(2, 0.0) , pMultiN_M(2, 0.0);
			double productHeatNP(0.0), productStressNP2(0.0), productHeatNM(0.0), productStressNM2(0.0), BM(0.0), BP(0.0), BOut(0.0);
			//update variables
			heatP[0] = stressHeatFluxP[0][2];
			heatP[1] = stressHeatFluxP[1][2];
			heatM[0] = stressHeatFluxM[0][2];
			heatM[1] = stressHeatFluxM[1][2];
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					stressP[i][j] = stressHeatFluxP[i][j];
					stressM[i][j] = stressHeatFluxM[i][j];
				}
			}
			nM[0] = -nP[0];
			nM[1] = -nP[1];

			pMultiN_P[0] = -pP * nP[0];
			pMultiN_P[1] = -pP * nP[1];

			pMultiN_M[0] = -pM * nM[0];
			pMultiN_M[1] = -pM * nM[1];

			//calculating
			productHeatNP = math::vectorDotProduct(heatP, nP);
			productHeatNM = math::vectorDotProduct(heatM, nM);

			productStressNP = math::tensorVectorDotProduct(stressP, nP);
			//productStressNP = math::vectorSum(productStressNP, pMultiN_P);

			productStressNM = math::tensorVectorDotProduct(stressM, nM);
			//productStressNM = math::vectorSum(productStressNM, pMultiN_M);

			productStressNP2 = math::vectorNorm(productStressNP);
			productStressNM2 = math::vectorNorm(productStressNM);

			BP = (pow(pow(rhoP*productHeatNP, 2) + 2 * rhoP*rhoP*eP*pow(productStressNP2, 2), 0.5) + rhoP * fabs(productHeatNP)) / (2 * rhoP*rhoP*eP);
			BM = (pow(pow(rhoP*productHeatNM, 2) + 2 * rhoM*rhoM*eM*pow(productStressNM2, 2), 0.5) + rhoM * fabs(productHeatNM)) / (2 * rhoM*rhoM*eM);
			if (BP > BM)
			{
				BOut = BP;
			}
			else
			{
				BOut = BM;
			}
			return BOut;
		}
	}//end of namespace numericalFluxes

	namespace inviscidTerms
	{
		std::tuple<double, double, double, double> calcInvisTermsFromPriVars(double rhoVal, double uVal, double vVal, double totalE, double pVal, int dir)
		{
			double term1(0.0), term2(0.0), term3(0.0), term4(0.0);

			if (dir==1)  //Ox direction
			{ 
				term1 = rhoVal * uVal;
				term2 = rhoVal * pow(uVal, 2) + pVal;
				term3 = rhoVal * uVal*vVal;
				term4 = (rhoVal*totalE + pVal)*uVal;
			}
			else if (dir==2)  //Oy direction
			{
				term1 = rhoVal * vVal;
				term2 = rhoVal * uVal*vVal;
				term3 = rhoVal * pow(vVal, 2) + pVal;
				term4 = (rhoVal*totalE + pVal)*vVal;
			}
			return std::make_tuple(term1, term2, term3, term4);
		}
	}//end of namespace invicidTerms

	namespace viscousTerms
	{
		std::vector<std::vector<double>> calcStressTensorAndHeatFlux(std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy)
		{
			/*Output matrix has form:
			[tauXx		tauXy		Qx]
			[tauYx		tauYy		Qy]
			*/
			std::vector<std::vector<double>> OutputMatrix(2, std::vector<double>(3, 0.0));

			double dux(0.0), duy(0.0), dvx(0.0), dvy(0.0), dEx(0.0), dEy(0.0), rhouVal(0.0), rhovVal, rhoVal(0.0), rhoEVal(0.0);
			double drhox(0.0), drhoy(0.0),
				drhoux(0.0), drhouy(0.0),
				drhovx(0.0), drhovy(0.0),
				drhoEx(0.0), drhoEy(0.0),
				dTx(0.0), dTy(0.0);
			double uVal(0.0), vVal(0.0), Qx(0.0), Qy(0.0), k(0.0);
			int index(0);

			rhoVal = U[0];
			rhouVal = U[1];
			rhovVal = U[2];
			rhoEVal = U[3];

			uVal = rhouVal / rhoVal;
			vVal = rhovVal / rhoVal;

			drhox = dUx[0];
			drhoy = dUy[0];

			drhoux = dUx[1];
			drhouy = dUy[1];

			drhovx = dUx[2];
			drhovy = dUy[2];

			drhoEx = dUx[3];
			drhoEy = dUy[3];

			dux = math::calcRhouvEDeriv(drhoux, drhox, rhouVal, rhoVal);
			duy = math::calcRhouvEDeriv(drhouy, drhoy, rhouVal, rhoVal);

			dvx = math::calcRhouvEDeriv(drhovx, drhox, rhovVal, rhoVal);
			dvy = math::calcRhouvEDeriv(drhovy, drhoy, rhovVal, rhoVal);

			dEx = math::calcRhouvEDeriv(drhoEx, drhox, rhoEVal, rhoVal);
			dEy = math::calcRhouvEDeriv(drhoEy, drhoy, rhoEVal, rhoVal);

			/*calculate stresses*/
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					index = 10 * (i + 1) + (j + 1);
					if (index == 11)  //x normal stress (tau_xx)
					{
						OutputMatrix[i][j] = math::viscousTerms::calcStressComponent(index, dux, dvy);
					}
					else if (index == 22)  //y normal stress (tau_yy)
					{
						OutputMatrix[i][j] = math::viscousTerms::calcStressComponent(index, dvy, dux);
					}
					else  //shear stress (tau_xy)
					{
						OutputMatrix[i][j] = math::viscousTerms::calcStressComponent(index, duy, dvx);
					}
				}
			}

			/*calculate heat flux*/
			dTx = math::calcTDeriv(dEx, dux, dvx, uVal, vVal);
			dTy = math::calcTDeriv(dEy, duy, dvy, uVal, vVal);
			//k = math::calcThermalConductivity(muVal);
			k = material::Cp / material::Pr;
			std::tie(Qx, Qy) = math::viscousTerms::calcHeatFluxTerms(dTx, dTy, k);

			OutputMatrix[0][2] = Qx;
			OutputMatrix[1][2] = Qy;
			return OutputMatrix;
		}

		double calcStressComponent(int index, double fstDeriv, double sndDeriv)
		{
			/*Formular of stress
			tau_xx = -mu((4/3)*du/dx - (2/3)*dv/dy)
			tau_yy = -mu((4/3)*dv/dy - (2/3)*du/dx)
			tau_xy = -mu(du/dy + dv/dx)

			input variables:
			index			fstDeriv		sndDeriv
			11 (tau_xx)		du/dx			dv/dy
			22 (tau_yy)		dv/dy			du/dx
			12 (tau_xy)		du/dy			dv/dx
			*/
			double tau(0.0);
			if (index == 11 || index == 22)
			{
				tau = -((4.0 / 3.0)*fstDeriv - (2.0 / 3.0)*sndDeriv);
			}
			else if (index == 12 || index == 21)
			{
				tau = -(fstDeriv + sndDeriv);
			}
			return tau;
		}

		std::tuple<double, double> calcHeatFluxTerms(double dTx, double dTy, double k)
		{
			double Qx(0.0), Qy(0.0);
			Qx = -k * dTx;
			Qy = -k * dTy;
			return std::make_tuple(Qx, Qy);
		}

		std::tuple<double, double, double, double> calcViscousTermsFromStressHeatFluxMatrix(std::vector< std::vector<double> > &StressHeatFlux, double uVal, double vVal, int dir)
		{
			/*StressHeatFlux is 2D array contents stress and heat flux component, has the following form:
			[tauXx		tauXy		Qx]
			[tauYx		tauYy		Qy]*/

			double tauXy(0.0), tauXx(0.0), tauYy(0.0), Qx(0.0), Qy(0.0);
			double viscTerm1(0.0), viscTerm2(0.0), viscTerm3(0.0), viscTerm4(0.0);

			tauXx = StressHeatFlux[0][0];
			tauYy = StressHeatFlux[1][1];
			tauXy = StressHeatFlux[0][1];
			Qx = StressHeatFlux[0][2];
			Qy = StressHeatFlux[1][2];

			if (dir==1)
			{
				/*1. Ox direction*/
				viscTerm1 = 0.0;
				viscTerm2 = tauXx;
				viscTerm3 = tauXy;
				viscTerm4 = tauXx * uVal + tauXy * vVal + Qx;
			}
			else if (dir==2)
			{
				/*2. Oy direction*/
				viscTerm1 = 0.0;
				viscTerm2 = tauXy;
				viscTerm3 = tauYy;
				viscTerm4 = tauXy * uVal + tauYy * vVal + Qy;
			}
			return std::make_tuple(viscTerm1, viscTerm2, viscTerm3, viscTerm4);
		}
	}//end of namespace viscousTerms

	/*Function computes value of conservative variables at abitrary point with applying limiter
			valType:
			1: rho
			2: rhou
			3: rhov
			4: rhoE*/
	double calcConsvVarWthLimiter(int element, double a, double b, int valType)
	{
		double out(0.0);
		std::vector<double> Value(mathVar::orderElem + 1, 0.0);
		Value = auxUlti::getElementConserValuesOfOrder(element, valType);

		//Compute value at point (a, b) without limiter
		math::basisFc(a, b, auxUlti::checkType(element));

		if (valType == 1)
		{
			for (int order = 1; order <= mathVar::orderElem; order++)
			{
				out += Value[order] * mathVar::B[order] * theta2Arr[element] * theta1Arr[element];
			}
		}
		else
		{
			for (int order = 1; order <= mathVar::orderElem; order++)
			{
				out += Value[order] * mathVar::B[order] * theta2Arr[element];
			}
		}
		out += Value[0];

		return out;
	}

	double calcResidualFromResidualOfOrder(int element, double a, double b, int valType)
	{
		double out(0.0);
		std::vector<double> Value(mathVar::orderElem + 1, 0.0);
		Value = auxUlti::getResidualValuesOfOrder(element, valType);

		//Compute value at point (a, b) without limiter
		math::basisFc(a, b, auxUlti::checkType(element));

		for (int order = 1; order <= mathVar::orderElem; order++)
		{
			out += Value[order] * mathVar::B[order];
		}
		out += Value[0];

		return fabs(out);
	}

	namespace geometricOp
	{
		std::tuple<double, double> calcGeoCenter(std::vector<double> &xCoor, std::vector<double> &yCoor, int type)
		{
			double x(0.0), y(0.0), xCG(0.0), yCG(0.0);
			for (int i = 0; i < type; i++)
			{
				x = xCoor[i];
				y = yCoor[i];
				xCG += x;
				yCG += y;
			}

			xCG = xCG / type;
			yCG = yCG / type;
			return std::make_tuple(xCG, yCG);
		}

		double calcPolygonArea(std::vector<double> &xCoor, std::vector<double> &yCoor, int type)
		{
			double x1(xCoor[0]), x2(xCoor[1]), x3(xCoor[2]), x4(0.0), y1(yCoor[0]), y2(yCoor[1]), y3(yCoor[2]), y4(0.0),Area(0.0);
			if (type == 3)
			{
				Area = fabs(0.5*(x1 * y2 - x2 * y1 + x2 * y3 - x3 * y2 + x3 * y1 - x1 * y3));
			}
			else if (type == 4)
			{
				x4 = xCoor[3];
				y4 = yCoor[3];

				Area = fabs(0.5*(x1 * y2 - x2 * y1 + x2 * y3 - x3 * y2 + x3 * y4 - x4 * y3 + x4 * y1 - x1 * y4));
			}
			return Area;
		}

		std::tuple<double, double> calcQuadCentroid(int element, double xCG, double yCG, double area)
		{
			std::vector<double> xSubTriCoor(3, 0.0), ySubTriCoor(3, 0.0);
			std::vector<double> xCGSubTri(4, 0.0), yCGSubTri(4, 0.0), subTriArea(4,0.0);
			double xC(0.0), yC(0.0);
			//1. point 0, 1
			std::tie(xSubTriCoor[0], ySubTriCoor[0]) = auxUlti::getElemCornerCoord(element, 0);
			std::tie(xSubTriCoor[1], ySubTriCoor[1]) = auxUlti::getElemCornerCoord(element, 1);
			xSubTriCoor[2] = xCG;
			ySubTriCoor[2] = yCG;
			std::tie(xCGSubTri[0], yCGSubTri[0]) = math::geometricOp::calcGeoCenter(xSubTriCoor, ySubTriCoor, 3);
			subTriArea[0] = math::geometricOp::calcPolygonArea(xSubTriCoor, ySubTriCoor, 3);

			//2. point 1, 2
			std::tie(xSubTriCoor[0], ySubTriCoor[0]) = auxUlti::getElemCornerCoord(element, 1);
			std::tie(xSubTriCoor[1], ySubTriCoor[1]) = auxUlti::getElemCornerCoord(element, 2);
			xSubTriCoor[2] = xCG;
			ySubTriCoor[2] = yCG;
			std::tie(xCGSubTri[1], yCGSubTri[1]) = math::geometricOp::calcGeoCenter(xSubTriCoor, ySubTriCoor, 3);
			subTriArea[1] = math::geometricOp::calcPolygonArea(xSubTriCoor, ySubTriCoor, 3);

			//3. point 2, 3
			std::tie(xSubTriCoor[0], ySubTriCoor[0]) = auxUlti::getElemCornerCoord(element, 2);
			std::tie(xSubTriCoor[1], ySubTriCoor[1]) = auxUlti::getElemCornerCoord(element, 3);
			xSubTriCoor[2] = xCG;
			ySubTriCoor[2] = yCG;
			std::tie(xCGSubTri[2], yCGSubTri[2]) = math::geometricOp::calcGeoCenter(xSubTriCoor, ySubTriCoor, 3);
			subTriArea[2] = math::geometricOp::calcPolygonArea(xSubTriCoor, ySubTriCoor, 3);

			//4. point 3, 0
			std::tie(xSubTriCoor[0], ySubTriCoor[0]) = auxUlti::getElemCornerCoord(element, 3);
			std::tie(xSubTriCoor[1], ySubTriCoor[1]) = auxUlti::getElemCornerCoord(element, 0);
			xSubTriCoor[2] = xCG;
			ySubTriCoor[2] = yCG;
			std::tie(xCGSubTri[3], yCGSubTri[3]) = math::geometricOp::calcGeoCenter(xSubTriCoor, ySubTriCoor, 3);
			subTriArea[3] = math::geometricOp::calcPolygonArea(xSubTriCoor, ySubTriCoor, 3);

			for (int i = 0; i < 4; i++)
			{
				xC += subTriArea[i] * xCGSubTri[i];
				yC += subTriArea[i] * yCGSubTri[i];
			}
			xC = xC / area;
			yC = yC / area;
			return std::make_tuple(xC, yC);
		}

		double calLocalCellSize(int element, double elementArea)
		{
			int elemType(auxUlti::checkType(element)), edgeId(0), neighborElemId(0);
			double deltaXe(0.0), deltaYe(0.0), Lx(0.0), Ly(0.0), Lxy(0.0), deltaXc(0.0), deltaYc(0.0);

			for (int e = 0; e < elemType; e++)
			{
				edgeId = meshVar::inedel[e][element];
				neighborElemId = auxUlti::getNeighborElement(element, edgeId);
				std::tie(deltaXe, deltaYe) = math::geometricOp::calEdgeMetric(edgeId, element);
				if (neighborElemId >= 0)
				{
					std::tie(deltaXc, deltaYc) = math::geometricOp::calDifferenceOfElementsCellCenters(element, neighborElemId);
				}
				else
				{
					deltaXc = 0.0;
					deltaYc = 0.0;
				}

				Lx += deltaXc * fabs(deltaXc)*deltaYe;
				Ly += deltaYc * fabs(deltaYc)*deltaXe;
			}
			Lx *= (1.0 / elementArea);
			Ly *= (-1.0 / elementArea);
			Lxy = sqrt(Lx*Lx + Ly*Ly);
			return Lxy;
		}

		std::tuple<double, double> calEdgeMetric(int edgeId, int elementId)
		{
			int point1(meshVar::inpoed[0][edgeId]), point2(meshVar::inpoed[1][edgeId]);
			double deltaX(0.0), deltaY(0.0);
			if (auxUlti::findVertexOrder(point1, elementId) > auxUlti::findVertexOrder(point2, elementId))
			{
				deltaX = meshVar::Points[point1][0] - meshVar::Points[point2][0];
				deltaY = meshVar::Points[point1][1] - meshVar::Points[point2][1];
			}
			else
			{
				deltaX = meshVar::Points[point2][0] - meshVar::Points[point1][0];
				deltaY = meshVar::Points[point2][1] - meshVar::Points[point1][1];
			}
			return std::make_tuple(deltaX, deltaY);
		}

		std::tuple<double, double> calDifferenceOfElementsCellCenters(int elem1, int elem2)
		{
			double deltaXc(0.0), deltaYc(0.0);
			//elem2 is neighbor
			deltaXc = (meshVar::geoCenter[elem2][0] - meshVar::geoCenter[elem1][0]);
			deltaYc = (meshVar::geoCenter[elem2][1] - meshVar::geoCenter[elem1][1]);
			return std::make_tuple(deltaXc, deltaYc);
		}

		double calDistBetween2Points(double xPt1, double yPt1, double xPt2, double yPt2) {
			double length(sqrt(pow(xPt1 - xPt2, 2) + pow(yPt1 - yPt2, 2)));
			return length;
		}

		double calPolygonPerimeter(std::vector<double> &xCoor, std::vector<double> &yCoor, int numOfEdge) {
			double perimeter(0.0), dx(0.0), dy(0.0);
			for (int iedge = 0; iedge < numOfEdge - 1; iedge++)
			{
				perimeter += math::geometricOp::calDistBetween2Points(xCoor[iedge], yCoor[iedge], xCoor[iedge + 1], yCoor[iedge + 1]);
			}
			perimeter += math::geometricOp::calDistBetween2Points(xCoor[numOfEdge -1], yCoor[numOfEdge - 1], xCoor[0], yCoor[0]);
			return perimeter;
		}

        std::tuple<double, double> calROfInscribedCircleOfTriElement(std::vector<double> &xCoor, std::vector<double> &yCoor) {
			double rIn(0.0), cellArea(0.0);
            cellArea = math::geometricOp::calcPolygonArea(xCoor, yCoor, 3);
            rIn = 2 * cellArea / math::geometricOp::calPolygonPerimeter(xCoor, yCoor, 3);
			return std::make_tuple(rIn, cellArea);
		}

        std::tuple<double, double> calSizeOfQuadElement(std::vector<double> &xCoor, std::vector<double> &yCoor) {
            int elemType(4);
            double size(0.0), size1(0.0), size2(0.0), cellArea(0.0);
            cellArea = math::geometricOp::calcPolygonArea(xCoor, yCoor, elemType);
            size1 = math::geometricOp::calDistBetween2Points(xCoor[0], yCoor[0], xCoor[2], yCoor[2]);
            size2 = math::geometricOp::calDistBetween2Points(xCoor[1], yCoor[1], xCoor[3], yCoor[3]);
            if (size1>size2)
            {
                size=size2;
            }
            else
            {
                size=size1;
            }
            return std::make_tuple(size, cellArea);
        }
	}

	namespace residualManipulation
	{
		/*
		void calcNormResidual(double rhoRes, double rhouRes, double rhovRes, double rhoERes)
		{
			systemVar::rhoResNormVector[systemVar::iterCount - 1] = rhoRes;
			systemVar::rhouResNormVector[systemVar::iterCount - 1] = rhouRes;
			systemVar::rhovResNormVector[systemVar::iterCount - 1] = rhovRes;
			systemVar::rhoEResNormVector[systemVar::iterCount - 1] = rhoERes;

			systemVar::rhoResNorm = *std::max_element(systemVar::rhoResNormVector.begin(), systemVar::rhoResNormVector.end());  //find min value of vector
			systemVar::rhouResNorm = *std::max_element(systemVar::rhouResNormVector.begin(), systemVar::rhouResNormVector.end());
			systemVar::rhovResNorm = *std::max_element(systemVar::rhovResNormVector.begin(), systemVar::rhovResNormVector.end());
			systemVar::rhoEResNorm = *std::max_element(systemVar::rhoEResNormVector.begin(), systemVar::rhoEResNormVector.end());
		}
		*/
	}

	double calcMaxT(int element)
	{
		std::vector<double> vectorT(2 * (mathVar::nGauss + 1) * (mathVar::nGauss + 1), 0.0);
		double aG(0.0), bG(0.0), aGL(0.0), bGL(0.0), min(0.0);
		int index(0);
		for (int na = 0; na <= mathVar::nGauss; na++)
		{
			for (int nb = 0; nb <= mathVar::nGauss; nb++)
			{
				aG = mathVar::GaussPts[na][nb][0];
				bG = mathVar::GaussPts[na][nb][1];

				aGL = mathVar::GaussLobattoPts[na][nb][0];
				bGL = mathVar::GaussLobattoPts[na][nb][1];

				vectorT[index] = math::pointValue(element, aG, bGL, 6, 1);
				index++;
				vectorT[index] = math::pointValue(element, aGL, bG, 6, 1);
				index++;
			}
		}
		min = *std::max_element(vectorT.begin(), vectorT.end());  //find min value of vector
		return min;
	}

	int findIndex(int number, std::vector<int> InArray)
	{
		int index(0), size(InArray.size());
		for (int i = 0; i < size; i++)
		{
			if (number == InArray[i])
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
}//end of namespace math
