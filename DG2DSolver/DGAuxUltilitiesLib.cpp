#include "DGAuxUltilitiesLib.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGMath.h"
#include <vector>
#include <math.h>
#include <tuple>
#include <QString>
#include <QDir>
#include <QProcess>

namespace auxUlti
{
    int findEdgeOrder(int element, int edge)
	{
		int order(0);
        int pt1(meshVar::inpoed[edge][0]), pt2(meshVar::inpoed[edge][1]);
		int typeElem(checkType(element));
		int ABpt1(0), ABpt2(0), BCpt1(0), BCpt2(0), CDpt1(0), CDpt2(0), DApt1(0), DApt2(0), CApt1(0), CApt2(0);

		if (typeElem==4)  //Quad element
		{
			ABpt1 = meshVar::Elements2D[element][0];
			ABpt2 = meshVar::Elements2D[element][1];

            BCpt1 = ABpt2;
			BCpt2 = meshVar::Elements2D[element][2];

			CDpt1 = BCpt2;
			CDpt2 = meshVar::Elements2D[element][3];

			DApt1 = CDpt2;
			DApt2 = ABpt1;

			if ((pt1 == ABpt1 && pt2 == ABpt2) || (pt1 == ABpt2 && pt2 == ABpt1))
			{
				order = 0;
			}
			else if ((pt1 == BCpt1 && pt2 == BCpt2) || (pt1 == BCpt2 && pt2 == BCpt1))
			{
				order = 1;
			}
			else if ((pt1 == CDpt1 && pt2 == CDpt2) || (pt1 == CDpt2 && pt2 == CDpt1))
			{
				order = 2;
			}
			else if ((pt1 == DApt1 && pt2 == DApt2) || (pt1 == DApt2 && pt2 == DApt1))
			{
				order = 3;
			}
		}
		else if (typeElem == 3)  //Tri element
		{
			ABpt1 = meshVar::Elements2D[element][0];
			ABpt2 = meshVar::Elements2D[element][1];

			BCpt1 = ABpt2;
			BCpt2 = meshVar::Elements2D[element][2];

			CApt1 = BCpt2;
			CApt2 = ABpt1;

			if ((pt1 == ABpt1 && pt2 == ABpt2) || (pt1 == ABpt2 && pt2 == ABpt1))
			{
				order = 0;
			}
			else if ((pt1 == BCpt1 && pt2 == BCpt2) || (pt1 == BCpt2 && pt2 == BCpt1))
			{
				order = 1;
			}
			else if ((pt1 == CApt1 && pt2 == CApt2) || (pt1 == CApt2 && pt2 == CApt1))
			{
				order = 2;
			}
		}
		return order;
	}

	int checkType(int element)
	{
		int typeElem(0);
		int typeFlag(meshVar::Elements2D[element][3]);
		if (typeFlag<0)
		{
			typeElem = 3;
		}
		else
		{
			typeElem = 4;
		}
		return typeElem;
	}

    std::tuple<double, double> getElemCornerCoord(int elem, int index)
	{
		/*Note: index starts from 0*/
		int pt(meshVar::Elements2D[elem][index]);
		double x(meshVar::Points[pt][0]), y(meshVar::Points[pt][1]);
		return std::make_tuple(x, y);
	}
	
	std::string workingdir()
	{
        QString workingDirectory(QDir::currentPath());
        //Convert QString to std::string
        std::string out = workingDirectory.toUtf8().constData();
        return out;
	}

	bool checkMaster(int elem, int edge)
	{
		bool master(true);
        if (meshVar::MasterElemOfEdge[edge]==elem)
		{
			master = true;
		}
		else
		{
			master = false;
		}
		return master;
	}

	std::tuple<double, double> getGaussSurfCoor(int edge, int elem, int nG)
	{
		double a(0.0), b(0.0);
		bool isMaster(auxUlti::checkMaster(elem, edge));

		if (isMaster)
		{
			a = meshVar::edgeGaussPoints_a[edge][nG];
			b = meshVar::edgeGaussPoints_b[edge][nG];
		}
		else
		{
			a = meshVar::edgeGaussPoints_a[edge][nG + mathVar::nGauss + 1];
			b = meshVar::edgeGaussPoints_b[edge][nG + mathVar::nGauss + 1];
		}

		return std::make_tuple(a, b);
	}

	std::tuple<double, double> getGaussSurfCoorMaster(int edge, int elem, int nG)
	{
		double a(0.0), b(0.0);
		int edgeOrder(auxUlti::findEdgeOrder(elem, edge));
		int elemType(auxUlti::checkType(elem));

		if (elemType == 4)  //quad element
		{
			if (edgeOrder == 0)
			{
				a = mathVar::xGauss[nG];
				b = -1.0;
			}
			else if (edgeOrder == 1)
			{
				a = 1.0;
				b = mathVar::xGauss[nG];
			}
			else if (edgeOrder == 2)
			{
				a = mathVar::xGauss[nG];
				b = 1.0;
			}
			else if (edgeOrder == 3)
			{
				a = -1.0;
				b = mathVar::xGauss[nG];
			}
		}
		else if (elemType == 3)  //tri element
		{
			if (edgeOrder == 0)
			{
				a = mathVar::xGauss[nG];
				b = -1.0;
			}
			else if (edgeOrder == 1)
			{
				a = 1.0;
				b = mathVar::xGauss[nG];
			}
			else if (edgeOrder == 2)
			{
				a = -1.0;
				b = mathVar::xGauss[nG];
			}
		}
		return std::make_tuple(a, b);
	}

	//This function supports for inverse coodinates mapping
	std::vector<std::vector<double>> getVectorGaussSurfCoor(int edge, int elem)
	{
		std::vector<std::vector<double>> vectorGaussPoints(mathVar::nGauss + 1, std::vector<double>(2, 0.0));
        for (int nG = 0; nG <= mathVar::nGauss; nG++)
		{
			std::tie(vectorGaussPoints[nG][0], vectorGaussPoints[nG][1]) = auxUlti::getGaussSurfCoorMaster(edge, elem, nG);
		}
		return vectorGaussPoints;
	}

	double getNormVectorComp(int elem, int edge, int dir)
	{
		double n(0.0);
		bool master(auxUlti::checkMaster(elem, edge));
		if (dir == 1)  //x direction
		{
			n = meshVar::normalVector[edge][0];
		}
		else if (dir == 2)  //y direction
		{
			n = meshVar::normalVector[edge][1];
		}

		if (master == false)
		{
			n = -1.0 * n;
		}
		return n;
	}

	void openFileEXE(std::string location)
	{
        QProcess *myProcess = new QProcess();
        myProcess->start(QString::fromUtf8(location.c_str()));
	}

	std::vector<double> getElementConserValuesOfOrder(int element, int type)
	{
		std::vector<double> Out(mathVar::orderElem + 1, 0.0);
		if (type == 1)  //rho
		{
            for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
			{
				Out[iorder] = rho[element][iorder];
			}
		}
		else if (type == 2)  //rhou
		{
            for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
			{
				Out[iorder] = rhou[element][iorder];
			}
		}
		else if (type == 3)  //rhov
		{
            for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
			{
				Out[iorder] = rhov[element][iorder];
			}
		}
		else if (type == 4)  //rhoE
		{
            for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
			{
				Out[iorder] = rhoE[element][iorder];
			}
		}
		
		return Out;
	}

	std::vector<double> getResidualValuesOfOrder(int element, int type)
	{
		std::vector<double> Out(mathVar::orderElem + 1, 0.0);
		if (type == 1)  //rho
		{
            for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
			{
				Out[iorder] = rhoResArr[element][iorder];
			}
		}
		else if (type == 2)  //rhou
		{
            for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
			{
				Out[iorder] = rhouResArr[element][iorder];
			}
		}
		else if (type == 3)  //rhov
		{
            for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
			{
				Out[iorder] = rhovResArr[element][iorder];
			}
		}
		else if (type == 4)  //rhoE
		{
            for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
			{
				Out[iorder] = rhoEResArr[element][iorder];
			}
		}

		return Out;
	}

	std::vector<double> getElementAuxValuesOfOrder(int element, int type, int dir)
	{
		std::vector<double> Out(mathVar::orderElem + 1, 0.0);
		if (dir==1)  //Ox direction
		{
			if (type == 1)  //d(rho)x
			{
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhoX[element][iorder];
				}
			}
			else if (type == 2)  //d(rhou)x
			{
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhouX[element][iorder];
				}
			}
			else if (type == 3)  //d(rhov)x
			{
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhovX[element][iorder];
				}
			}
			else if (type == 4)  //d(rhoE)x
			{
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhoEX[element][iorder];
				}
			}
		}
		else if (dir==2)  //Oy direction
		{
			if (type == 1)  //d(rho)y
			{
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhoY[element][iorder];
				}
			}
			else if (type == 2)  //d(rhou)y
			{
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhouY[element][iorder];
				}
			}
			else if (type == 3)  //d(rhov)y
			{
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhovY[element][iorder];
				}
			}
			else if (type == 4)  //d(rhoE)y
			{
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhoEY[element][iorder];
				}
			}
		}

		return Out;
	}

	std::tuple<double, double> getGaussCoor(int na, int nb)
	{
		double a(0.0), b(0.0);
		a = mathVar::GaussPts[na][nb][0];
		b= mathVar::GaussPts[na][nb][1];
		return std::make_tuple(a, b);
	}

	int getGrpOfEdge(int edge)
	{
        int grp(meshVar::inpoed[edge][2]);
		return grp;
	}

	int getBCType(int edge)
	{
        int bcType(meshVar::inpoed[edge][3]);
		return bcType;
	}

	bool checkSubSonic()
	{
		double uInf(0.0), vInf(0.0), TInf(0.0), SpeedOfSound(0.0), Mach(0.0), Velocity(0.0);
		bool Out(true);
        for (int i = 0; i < meshVar::nBc; i++)
		{
			if (bcValues::UBcType[i]==1 || bcValues::UBcType[i] == 4)
			{
				TInf = bcValues::TBC[i];
				uInf = bcValues::uBC[i];
				vInf = bcValues::vBC[i];
				SpeedOfSound = (sqrt(material::gamma*material::R*TInf));
				Velocity = (sqrt(uInf*uInf + vInf * vInf));
				Mach = Velocity / SpeedOfSound;
				if (Mach >= 1.0)
				{
					Out = false;
					break;
				}
			}
		}
		return Out;
	}

	bool checkSubSonicLocally(double TVal, double uVal, double vVal)
	{
		bool Out(true);
		double SpeedOfSound(sqrt(material::gamma*material::R*TVal)),
			Velocity(sqrt(uVal*uVal + vVal * vVal));

		double Mach(Velocity / SpeedOfSound);
		if (Mach >= 1.0)
		{
			Out = false;
		}
		return Out;
	}

	std::tuple<int, int> getMasterServantOfEdge(int edge)
	{
        int master(meshVar::MasterElemOfEdge[edge]), servant(0), elem1(meshVar::ineled[edge][0]), elem2(meshVar::ineled[edge][1]);
		if (master==elem1)
		{
			servant = elem2;
		}
		else if (master==elem2)
		{
			servant = elem1;
		}
		return std::make_tuple(master, servant);
	}

    void resize2DArray(std::vector<std::vector<double>> &Array, int row, int column)
	{
		Array.resize(row);
        for (int i = 0; i < row; ++i)
		{
			Array[i].resize(column);
		}
	}

    void resize2DIntArray(std::vector<std::vector<int>> &Array, int row, int column)
	{
		Array.resize(row);
        for (int i = 0; i < row; ++i)
		{
			Array[i].resize(column);
		}
	}

    void resize3DArray(std::vector<std::vector<std::vector<double>>> &Array, int direct1, int direct2, int direct3)
	{
		Array.resize(direct1);
        for (int i = 0; i < direct1; ++i)
		{
			Array[i].resize(direct2);
            for (int j = 0; j < direct2; j++)
			{
				Array[i][j].resize(direct3);
			}
		}
	}

    void addRowTo2DIntArray(std::vector<std::vector<int>> &Array, int numCol)
	{
        int length(Array.size());
		Array.push_back(std::vector<int>());
        for (int icol = 0; icol < numCol; icol++)
		{
			Array[length].push_back(-1);
		}
	}

    void addRowTo2DDoubleArray(std::vector<std::vector<double>> &Array, int numCol)
	{
		int length(Array.size());
		Array.push_back(std::vector<double>());
        for (int icol = 0; icol < numCol; icol++)
		{
			Array[length].push_back(0.0);
		}
	}

	//Function returns cell centroid coordinates and size (cell area)
	std::tuple<double, double, double> getCellMetrics(int element)
	{
		double xC(meshVar::geoCenter[element][0]), yC(meshVar::geoCenter[element][1]), size(meshVar::cellSize[element]);
		return std::make_tuple(xC, yC, size);
	}

	void mappingEdges()
	{
        int masterElem(0), servantElem(0), bcType(0);
		double aMaster(0.0), bMaster(0.0), aServant(0.0), bServant(0.0),
			xMaster(0.0), yMaster(0.0);

        for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
		{
 			std::tie(masterElem, servantElem) = auxUlti::getMasterServantOfEdge(iedge);
			bcType = auxUlti::getBCType(iedge);
            for (int nG = 0; nG <= mathVar::nGauss; nG++)
			{
				std::tie(aMaster, bMaster) = auxUlti::getGaussSurfCoorMaster(iedge, masterElem, nG);
				std::tie(xMaster, yMaster) = math::directMapping(masterElem, aMaster, bMaster);
				meshVar::edgeGaussPoints_a[iedge][nG] = aMaster;
				meshVar::edgeGaussPoints_b[iedge][nG] = bMaster;
				if (bcType != 0)
				{
					aServant = 0.0;
					bServant = 0.0;
				}
				else
				{
					std::tie(aServant, bServant) = math::inverseMapping(servantElem, xMaster, yMaster);
				}
				meshVar::edgeGaussPoints_a[iedge][nG + mathVar::nGauss + 1] = aServant;
				meshVar::edgeGaussPoints_b[iedge][nG + mathVar::nGauss + 1] = bServant;
			}
		}
	}

	void resizeDGArrays()
	{
        //Resize mesh arrays
        auxUlti::resize2DArray(meshVar::geoCenter,meshVar::nelem2D,2);
        meshVar::cellArea.resize(meshVar::nelem2D);
        meshVar::cellSize.resize(meshVar::nelem2D);
        meshVar::localCellSize.resize(meshVar::nelem2D);
        auxUlti::resize3DArray(meshVar::dxa, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(meshVar::dya, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(meshVar::dxb, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(meshVar::dyb, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(meshVar::J2D, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);
        meshVar::J1D.resize(meshVar::inpoedCount);
        //auxUlti::resize2DArray(meshVar::J1D, meshVar::inpoedCount, 2);

        auxUlti::resize2DArray(meshVar::edgeGaussPoints_a, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(meshVar::edgeGaussPoints_b, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));

        auxUlti::resize2DArray(rho, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhou, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhov, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhoE, meshVar::nelem2D, mathVar::orderElem + 1);

        auxUlti::resize2DArray(rho0, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhou0, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhov0, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhoE0, meshVar::nelem2D, mathVar::orderElem + 1);

        auxUlti::resize2DArray(rhoResArr, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhouResArr, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhovResArr, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhoEResArr, meshVar::nelem2D, mathVar::orderElem + 1);

        auxUlti::resize2DArray(SurfaceBCFields::rhoBc, mathVar::nGauss + 1, meshVar::numBCEdges);
        auxUlti::resize2DArray(SurfaceBCFields::rhouBc, mathVar::nGauss + 1, meshVar::numBCEdges);
        auxUlti::resize2DArray(SurfaceBCFields::rhovBc, mathVar::nGauss + 1, meshVar::numBCEdges);
        auxUlti::resize2DArray(SurfaceBCFields::rhoEBc, mathVar::nGauss + 1, meshVar::numBCEdges);

        auxUlti::resize2DArray(rhoN, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhouN, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhovN, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhoEN, meshVar::nelem2D, mathVar::orderElem + 1);

        auxUlti::resize2DArray(aux_interface_rho, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(aux_interface_rhou, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(aux_interface_rhov, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(aux_interface_rhoE, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));

        auxUlti::resize2DArray(interface_rho, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(interface_rhou, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(interface_rhov, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(interface_rhoE, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));

        auxUlti::resize2DArray(rhoX, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhouX, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhovX, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhoEX, meshVar::nelem2D, mathVar::orderElem + 1);

        auxUlti::resize2DArray(rhoY, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhouY, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhovY, meshVar::nelem2D, mathVar::orderElem + 1);
        auxUlti::resize2DArray(rhoEY, meshVar::nelem2D, mathVar::orderElem + 1);

        auxUlti::resize2DArray(invis_interface_rhoX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(invis_interface_rhouX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(invis_interface_rhovX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(invis_interface_rhoEX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));

        auxUlti::resize2DArray(invis_interface_rhoY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(invis_interface_rhouY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(invis_interface_rhovY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(invis_interface_rhoEY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));

        auxUlti::resize2DArray(Vis_interface_rhoX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(Vis_interface_rhouX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(Vis_interface_rhovX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(Vis_interface_rhoEX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));

        auxUlti::resize2DArray(Vis_interface_rhoY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(Vis_interface_rhouY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(Vis_interface_rhovY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(Vis_interface_rhoEY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));

        theta1Arr.resize(meshVar::nelem2D);
        theta2Arr.resize(meshVar::nelem2D);
        //debug::minRhoArr.resize(meshVar::nelem2D);
        //debug::minRhoeArr.resize(meshVar::nelem2D);

        LxFConst.resize(meshVar::inpoedCount);
        DiffusiveFluxConst.resize(meshVar::inpoedCount);

        auxUlti::resize2DArray(stiffMatrixCoeffs, meshVar::nelem2D, mathVar::orderElem + 1);

        auxUlti::resize3DArray(rhoVolGauss, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(rhouVolGauss, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(rhovVolGauss, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(rhoEVolGauss, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);

        //meshVar::adressOfBCVals.resize(meshVar::numBCEdges);
        auxUlti::resize2DIntArray(meshVar::neighboringElements, meshVar::nelem2D, 4);

        //Resize mathVar array
        mathVar::wGauss.resize(mathVar::nGauss+1);
        mathVar::xGauss.resize(mathVar::nGauss+1);
        mathVar::wGaussLobatto.resize(mathVar::nGauss+1);
        mathVar::xGaussLobatto.resize(mathVar::nGauss+1);
        mathVar::B.resize(mathVar::nGauss+1);
        mathVar::dBa.resize(mathVar::nGauss+1);
        mathVar::dBb.resize(mathVar::nGauss+1);
        auxUlti::resize3DArray(mathVar::BPts_Quad, mathVar::orderElem+1, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(mathVar::dBaPts_Quad, mathVar::orderElem+1, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(mathVar::dBbPts_Quad, mathVar::orderElem+1, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(mathVar::BPts_Tri, mathVar::orderElem+1, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(mathVar::dBaPts_Tri, mathVar::orderElem+1, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(mathVar::dBbPts_Tri, mathVar::orderElem+1, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(mathVar::GaussPts, mathVar::nGauss + 1, mathVar::nGauss + 1, 2);
        auxUlti::resize3DArray(mathVar::wGaussPts, mathVar::nGauss + 1, mathVar::nGauss + 1, 2);
        auxUlti::resize3DArray(mathVar::GaussLobattoPts, mathVar::nGauss + 1, mathVar::nGauss + 1, 2);
        auxUlti::resize3DArray(mathVar::wGaussLobattoPts, mathVar::nGauss + 1, mathVar::nGauss + 1, 2);
	}

	int getAdressOfBCEdgesOnBCValsArray(int edge)
	{
		int locate(0);
        for (int i = 0; i < meshVar::numBCEdges; i++)
		{
			if (edge == meshVar::adressOfBCVals[i])
			{
				locate = i;
				break;
			}
		}
		return locate;
	}

	std::tuple<double, double> getCellCentroid(int element)
	{
		double xC(meshVar::geoCenter[element][0]), yC(meshVar::geoCenter[element][1]);
		return std::make_tuple(xC, yC);
	}

	void clear1DIntVector(std::vector<int>&vector)
	{
		int vectorLenth(vector.size());
		/*Use erase function to clear vector*/
		vector.erase(vector.begin(), vector.begin() + vectorLenth);
		/*Shrink to fit*/
		vector.shrink_to_fit();
	}

    void clear2DIntVector(std::vector<std::vector<int>>&vector)
    {
        int numRow(static_cast<int>(vector.size()));
        /*Use erase function to clear vector*/
        for (int row = 0; row < numRow; ++row) {
            int length(static_cast<int>(vector[row].size()));
            vector[row].erase(vector[row].begin(), vector[row].begin() + length);
            vector[row].shrink_to_fit();
        }

        /*Shrink to fit*/
        vector.shrink_to_fit();
    }

	int findVertexOrder(int point, int element)
	{
		int elemType(auxUlti::checkType(element)), order(-1);
        for (int i = 0; i < elemType; i++)
		{
			if (point == meshVar::Elements2D[element][i])
			{
				order = i;
				break;
			}
		}
		return order;
	}

	namespace postProcess
	{
		std::vector<int> getElementsSurroundingPoint(int point)
		{
			std::vector<int>ElSurPt;
            for (int iesup = meshVar::esup2[point] + 1; iesup <= meshVar::esup2[point + 1]; iesup++)
			{
				ElSurPt.push_back(meshVar::esup1[iesup]);
			}
			return ElSurPt;
		}

		std::tuple<double, double> findPointCoorInStandardSpace(int point, int element)
		{
			std::vector<int> iarray;
			int index(0), elemType(auxUlti::checkType(element));
			double a(0.0), b(0.0);
            for (int ipoin = 0; ipoin < elemType; ipoin++)
			{
				iarray.push_back(meshVar::Elements2D[element][ipoin]);
			}

            for (int i = 0; i < elemType; i++)
			{
				if (point == iarray[i])
				{
					index = i;
					break;
				}
				else
				{
					index = -1;
				}
			}

			switch (index)
			{
			case 0:
			{
				a = -1.0;
				b = -1.0;
				break;
			}
			case 1:
			{
				a = 1.0;
				b = -1.0;
				break;
			}
			case 2:
			{
				if (elemType==3)
				{
					a = -1.0;
					b = 1.0;
				}
				else
				{
					a = 1.0;
					b = 1.0;
				}
				break;
			}
			case 3:
			{
				a = -1.0;
				b = 1.0;
				break;
			}
			default:
				break;
			}
			return std::make_tuple(a, b);
		}
	}

	int getNeighborElement(int element, int edge)
	{
		int neighbor(0);
		if (auxUlti::getBCType(edge) == 0)
		{
            if (meshVar::ineled[edge][0] == element)
			{
                neighbor = meshVar::ineled[edge][1];
			}
			else
			{
                neighbor = meshVar::ineled[edge][0];
			}
		}
		else
		{
			neighbor = -1;
		}
		return neighbor;
	}

	int getEdgeHasInputOrderOfElement(int element, int inputEdgeOrder)
	{
		int elemType(auxUlti::checkType(element)), edgeId(0), edgeOrder(0), outputEdgeId(0);
        for (int i = 0; i < elemType; i++)
		{
            edgeId = meshVar::inedel[element][i];
			edgeOrder = auxUlti::findEdgeOrder(element, edgeId);
			if (edgeOrder == inputEdgeOrder)
			{
				outputEdgeId = edgeId;
				break;
			}
		}
		return outputEdgeId;
	}

    void createFolder(std::string location)
    {
        std::string command("mkdir -p "+location);
        const int dir_err = system(command.c_str());
        if (-1 == dir_err)
        {
            printf("Error creating directory!n");
            exit(1);
        }
    }
}
