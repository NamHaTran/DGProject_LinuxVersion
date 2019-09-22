#include "DGAuxUltilitiesLib.h"
#include "DGMessagesLib.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGMath.h"
#include <vector>
#include <math.h>
#include <tuple>
#include <QString>
#include <QDir>
#include <QProcess>
#include <mpi.h>

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
                    Out[iorder] = BR1Vars::rhoX[element][iorder];
                }
            }
            else if (type == 2)  //d(rhou)x
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Out[iorder] = BR1Vars::rhouX[element][iorder];
                }
            }
            else if (type == 3)  //d(rhov)x
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Out[iorder] = BR1Vars::rhovX[element][iorder];
                }
            }
            else if (type == 4)  //d(rhoE)x
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Out[iorder] = BR1Vars::rhoEX[element][iorder];
                }
            }
        }
        else if (dir==2)  //Oy direction
        {
            if (type == 1)  //d(rho)y
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Out[iorder] = BR1Vars::rhoY[element][iorder];
                }
            }
            else if (type == 2)  //d(rhou)y
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Out[iorder] = BR1Vars::rhouY[element][iorder];
                }
            }
            else if (type == 3)  //d(rhov)y
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Out[iorder] = BR1Vars::rhovY[element][iorder];
                }
            }
            else if (type == 4)  //d(rhoE)y
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Out[iorder] = BR1Vars::rhoEY[element][iorder];
                }
            }
        }

		return Out;
	}

    std::vector<double> getElementAuxValuesOfOrder_BR2_vol(int element, int type, int dir)
    {
        std::vector<double> Out(mathVar::orderElem + 1, 0.0);
        if (dir==1) //Ox
        {
            switch (type) {
            case 1:
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Out[iorder] = BR2Vars::rhoXVol[element][iorder];
                }
            }
                break;
            case 2:
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Out[iorder] = BR2Vars::rhouXVol[element][iorder];
                }
            }
                break;
            case 3:
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Out[iorder] = BR2Vars::rhovXVol[element][iorder];
                }
            }
                break;
            case 4:
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Out[iorder] = BR2Vars::rhoEXVol[element][iorder];
                }
            }
                break;
            default:
                break;
            }
        }
        else { //Oy
            switch (type) {
            case 1:
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Out[iorder] = BR2Vars::rhoYVol[element][iorder];
                }
            }
                break;
            case 2:
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Out[iorder] = BR2Vars::rhouYVol[element][iorder];
                }
            }
                break;
            case 3:
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Out[iorder] = BR2Vars::rhovYVol[element][iorder];
                }
            }
                break;
            case 4:
            {
                for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                {
                    Out[iorder] = BR2Vars::rhoEYVol[element][iorder];
                }
            }
                break;
            default:
                break;
            }
        }

        return Out;
    }

     std::vector<double> getElementAuxValuesOfOrder_BR2_sur(int edge, int element, int type, int dir)
     {
         bool isMaster(auxUlti::checkMaster(element,edge));
         std::vector<double> Out(mathVar::orderElem + 1, 0.0);
         if (isMaster)
         {
             if (dir==1) //Ox
             {
                 switch (type) {
                 case 1:
                 {
                     for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                     {
                         Out[iorder] = BR2Vars::rhoXSurMaster[edge][iorder];
                     }
                 }
                     break;
                 case 2:
                 {
                     for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                     {
                         Out[iorder] = BR2Vars::rhouXSurMaster[edge][iorder];
                     }
                 }
                     break;
                 case 3:
                 {
                     for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                     {
                         Out[iorder] = BR2Vars::rhovXSurMaster[edge][iorder];
                     }
                 }
                     break;
                 case 4:
                 {
                     for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                     {
                         Out[iorder] = BR2Vars::rhoEXSurMaster[edge][iorder];
                     }
                 }
                     break;
                 default:
                     break;
                 }
             }
             else { //Oy
                 switch (type) {
                 case 1:
                 {
                     for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                     {
                         Out[iorder] = BR2Vars::rhoYSurMaster[edge][iorder];
                     }
                 }
                     break;
                 case 2:
                 {
                     for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                     {
                         Out[iorder] = BR2Vars::rhouYSurMaster[edge][iorder];
                     }
                 }
                     break;
                 case 3:
                 {
                     for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                     {
                         Out[iorder] = BR2Vars::rhovYSurMaster[edge][iorder];
                     }
                 }
                     break;
                 case 4:
                 {
                     for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                     {
                         Out[iorder] = BR2Vars::rhoEYSurMaster[edge][iorder];
                     }
                 }
                     break;
                 default:
                     break;
                 }
             }
         }
         else {
             if (dir==1) //Ox
              {
                  switch (type) {
                  case 1:
                  {
                      for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                      {
                          Out[iorder] = BR2Vars::rhoXSurSlave[edge][iorder];
                      }
                  }
                      break;
                  case 2:
                  {
                      for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                      {
                          Out[iorder] = BR2Vars::rhouXSurSlave[edge][iorder];
                      }
                  }
                      break;
                  case 3:
                  {
                      for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                      {
                          Out[iorder] = BR2Vars::rhovXSurSlave[edge][iorder];
                      }
                  }
                      break;
                  case 4:
                  {
                      for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                      {
                          Out[iorder] = BR2Vars::rhoEXSurSlave[edge][iorder];
                      }
                  }
                      break;
                  default:
                      break;
                  }
              }
              else { //Oy
                  switch (type) {
                  case 1:
                  {
                      for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                      {
                          Out[iorder] = BR2Vars::rhoYSurSlave[edge][iorder];
                      }
                  }
                      break;
                  case 2:
                  {
                      for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                      {
                          Out[iorder] = BR2Vars::rhouYSurSlave[edge][iorder];
                      }
                  }
                      break;
                  case 3:
                  {
                      for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                      {
                          Out[iorder] = BR2Vars::rhovYSurSlave[edge][iorder];
                      }
                  }
                      break;
                  case 4:
                  {
                      for (int iorder = 0; iorder <= mathVar::orderElem; iorder++)
                      {
                          Out[iorder] = BR2Vars::rhoEYSurSlave[edge][iorder];
                      }
                  }
                      break;
                  default:
                      break;
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
                    if (systemVar::parallelMode && bcType==10)
                    {
                        int localEdgeIdOnElem1DArray=auxUlti::getAdressOfBCEdgesOnBCValsArray(iedge);
                        std::tie(parallelBuffer::aCoor[localEdgeIdOnElem1DArray][nG], parallelBuffer::bCoor[localEdgeIdOnElem1DArray][nG]) = math::inverseMapping_ForParallel(iedge, xMaster, yMaster);
                    }
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

        auxUlti::resize2DArray(surfaceFields::aux_rho, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::aux_rhou, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::aux_rhov, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::aux_rhoE, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));

        auxUlti::resize2DArray(surfaceFields::rho, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::rhou, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::rhov, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::rhoE, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));

        if (systemVar::auxVariables==1)
        {
            auxUlti::resize2DArray(BR1Vars::rhoX, meshVar::nelem2D, mathVar::orderElem + 1);
            auxUlti::resize2DArray(BR1Vars::rhouX, meshVar::nelem2D, mathVar::orderElem + 1);
            auxUlti::resize2DArray(BR1Vars::rhovX, meshVar::nelem2D, mathVar::orderElem + 1);
            auxUlti::resize2DArray(BR1Vars::rhoEX, meshVar::nelem2D, mathVar::orderElem + 1);

            auxUlti::resize2DArray(BR1Vars::rhoY, meshVar::nelem2D, mathVar::orderElem + 1);
            auxUlti::resize2DArray(BR1Vars::rhouY, meshVar::nelem2D, mathVar::orderElem + 1);
            auxUlti::resize2DArray(BR1Vars::rhovY, meshVar::nelem2D, mathVar::orderElem + 1);
            auxUlti::resize2DArray(BR1Vars::rhoEY, meshVar::nelem2D, mathVar::orderElem + 1);

        }
        else if (systemVar::auxVariables==2)
        {
            auxUlti::resize2DArray(BR2Vars::rhoXVol, meshVar::nelem2D, mathVar::orderElem + 1);
            auxUlti::resize2DArray(BR2Vars::rhouXVol, meshVar::nelem2D, mathVar::orderElem + 1);
            auxUlti::resize2DArray(BR2Vars::rhovXVol, meshVar::nelem2D, mathVar::orderElem + 1);
            auxUlti::resize2DArray(BR2Vars::rhoEXVol, meshVar::nelem2D, mathVar::orderElem + 1);

            auxUlti::resize2DArray(BR2Vars::rhoYVol, meshVar::nelem2D, mathVar::orderElem + 1);
            auxUlti::resize2DArray(BR2Vars::rhouYVol, meshVar::nelem2D, mathVar::orderElem + 1);
            auxUlti::resize2DArray(BR2Vars::rhovYVol, meshVar::nelem2D, mathVar::orderElem + 1);
            auxUlti::resize2DArray(BR2Vars::rhoEYVol, meshVar::nelem2D, mathVar::orderElem + 1);

            auxUlti::resize2DArray(BR2Vars::rhoXSurMaster,meshVar::inpoedCount,(mathVar::orderElem + 1));
            auxUlti::resize2DArray(BR2Vars::rhouXSurMaster,meshVar::inpoedCount,(mathVar::orderElem + 1));
            auxUlti::resize2DArray(BR2Vars::rhovXSurMaster,meshVar::inpoedCount,(mathVar::orderElem + 1));
            auxUlti::resize2DArray(BR2Vars::rhoEXSurMaster,meshVar::inpoedCount,(mathVar::orderElem + 1));

            auxUlti::resize2DArray(BR2Vars::rhoYSurMaster,meshVar::inpoedCount,(mathVar::orderElem + 1));
            auxUlti::resize2DArray(BR2Vars::rhouYSurMaster,meshVar::inpoedCount,(mathVar::orderElem + 1));
            auxUlti::resize2DArray(BR2Vars::rhovYSurMaster,meshVar::inpoedCount,(mathVar::orderElem + 1));
            auxUlti::resize2DArray(BR2Vars::rhoEYSurMaster,meshVar::inpoedCount,(mathVar::orderElem + 1));

            auxUlti::resize2DArray(BR2Vars::rhoXSurSlave,meshVar::inpoedCount,(mathVar::orderElem + 1));
            auxUlti::resize2DArray(BR2Vars::rhouXSurSlave,meshVar::inpoedCount,(mathVar::orderElem + 1));
            auxUlti::resize2DArray(BR2Vars::rhovXSurSlave,meshVar::inpoedCount,(mathVar::orderElem + 1));
            auxUlti::resize2DArray(BR2Vars::rhoEXSurSlave,meshVar::inpoedCount,(mathVar::orderElem + 1));

            auxUlti::resize2DArray(BR2Vars::rhoYSurSlave,meshVar::inpoedCount,(mathVar::orderElem + 1));
            auxUlti::resize2DArray(BR2Vars::rhouYSurSlave,meshVar::inpoedCount,(mathVar::orderElem + 1));
            auxUlti::resize2DArray(BR2Vars::rhovYSurSlave,meshVar::inpoedCount,(mathVar::orderElem + 1));
            auxUlti::resize2DArray(BR2Vars::rhoEYSurSlave,meshVar::inpoedCount,(mathVar::orderElem + 1));
        }

        auxUlti::resize2DArray(surfaceFields::invis_rhoX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::invis_rhouX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::invis_rhovX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::invis_rhoEX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));

        auxUlti::resize2DArray(surfaceFields::invis_rhoY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::invis_rhouY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::invis_rhovY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::invis_rhoEY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));

        auxUlti::resize2DArray(surfaceFields::Vis_rhoX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::Vis_rhouX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::Vis_rhovX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::Vis_rhoEX, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));

        auxUlti::resize2DArray(surfaceFields::Vis_rhoY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::Vis_rhouY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::Vis_rhovY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::Vis_rhoEY, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));
        auxUlti::resize2DArray(surfaceFields::T, meshVar::inpoedCount, 2 * (mathVar::nGauss + 1));

        theta1Arr.resize(meshVar::nelem2D);
        theta2Arr.resize(meshVar::nelem2D);
        //debug::minRhoArr.resize(meshVar::nelem2D);
        //debug::minRhoeArr.resize(meshVar::nelem2D);

        LxFConst.resize(meshVar::inpoedCount);
        DiffusiveFluxConst.resize(meshVar::inpoedCount);

        auxUlti::resize2DArray(stiffMatrixCoeffs, meshVar::nelem2D, mathVar::orderElem + 1);

        auxUlti::resize3DArray(volumeFields::rhoVolGauss, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(volumeFields::rhouVolGauss, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(volumeFields::rhovVolGauss, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(volumeFields::rhoEVolGauss, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(volumeFields::drhoXVolGauss, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(volumeFields::drhoYVolGauss, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);
        auxUlti::resize3DArray(volumeFields::T, meshVar::nelem2D, mathVar::nGauss + 1, mathVar::nGauss + 1);

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

        //Buffer (parallel computing)
        //Conservative variables
        auxUlti::resize2DArray(parallelBuffer::rho, meshVar::numBCEdges, mathVar::orderElem + 1);
        auxUlti::resize2DArray(parallelBuffer::rhou, meshVar::numBCEdges, mathVar::orderElem + 1);
        auxUlti::resize2DArray(parallelBuffer::rhov, meshVar::numBCEdges, mathVar::orderElem + 1);
        auxUlti::resize2DArray(parallelBuffer::rhoE, meshVar::numBCEdges, mathVar::orderElem + 1);

        //Auxilary variables
        auxUlti::resize2DArray(parallelBuffer::drhoX, meshVar::numBCEdges, mathVar::orderElem + 1);
        auxUlti::resize2DArray(parallelBuffer::drhouX, meshVar::numBCEdges, mathVar::orderElem + 1);
        auxUlti::resize2DArray(parallelBuffer::drhovX, meshVar::numBCEdges, mathVar::orderElem + 1);
        auxUlti::resize2DArray(parallelBuffer::drhoEX, meshVar::numBCEdges, mathVar::orderElem + 1);
        auxUlti::resize2DArray(parallelBuffer::drhoY, meshVar::numBCEdges, mathVar::orderElem + 1);
        auxUlti::resize2DArray(parallelBuffer::drhouY, meshVar::numBCEdges, mathVar::orderElem + 1);
        auxUlti::resize2DArray(parallelBuffer::drhovY, meshVar::numBCEdges, mathVar::orderElem + 1);
        auxUlti::resize2DArray(parallelBuffer::drhoEY, meshVar::numBCEdges, mathVar::orderElem + 1);

        //Surface Gauss point coordinates
        auxUlti::resize2DArray(parallelBuffer::aCoor, meshVar::numBCEdges, mathVar::nGauss + 1);
        auxUlti::resize2DArray(parallelBuffer::bCoor, meshVar::numBCEdges, mathVar::nGauss + 1);
        //Neighbor cell vertexes
        auxUlti::resize2DArray(parallelBuffer::xCoor, meshVar::numBCEdges, 4);
        auxUlti::resize2DArray(parallelBuffer::yCoor, meshVar::numBCEdges, 4);

        parallelBuffer::theta1.resize(meshVar::numBCEdges);
        parallelBuffer::theta2.resize(meshVar::numBCEdges);
        parallelBuffer::elemType.resize(meshVar::numBCEdges);

        //Mesh connection array
        //auxUlti::resize2DIntArray(meshVar::meshConnection,meshVar::numBCEdges,2);
	}

    int getAdressOfBCEdgesOnBCValsArray(int globalEdgeId)
	{
        return meshVar::inpoed[globalEdgeId][4];
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

    void createFolder(std::string location, bool passExit)
    {
        std::string command("mkdir -p "+location);
        const int dir_err = system(command.c_str());
        if (-1 == dir_err)
        {
            printf("Error of creating directory\n");
            if (!passExit)
            {
                exit(1);
            }
        }
    }

    std::string createTimeStepFolder(int iter, std::string option)
    {
        std::string Loc;
        if (systemVar::currentProc==0)
        {
            std::string iter_str = std::to_string(iter), tecplotName("/");
            if (option.compare("tecplot")==0)
            {
                tecplotName="/TecplotFile/";
            }

            if (systemVar::parallelMode)
            {
                Loc = systemVar::wD + "/CASES/" + systemVar::caseName + "/Processor0" + tecplotName + iter_str;
                auxUlti::createFolder(Loc,true);
                for (int iproc = 1; iproc < systemVar::totalProc; iproc++)
                {
                    Loc = systemVar::wD + "/CASES/" + systemVar::caseName + "/Processor" + std::to_string(iproc) + tecplotName + iter_str;
                    auxUlti::createFolder(Loc,true);
                    auxUlti::functionsOfParallelComputing::sendString(Loc,iproc,1);
                }
            }
            else {
                Loc = systemVar::wD + "/CASES/" + systemVar::caseName + tecplotName + iter_str;
                auxUlti::createFolder(Loc,true);
            }
        }
        else
        {
            Loc=auxUlti::functionsOfParallelComputing::receiveString(0,1);
        }
        return Loc;
    }

    void copyFolder(std::string source, std::string destination)
    {
        std::string command("cp -a "+source+"/. "+destination+"/");
        const int dir_err = system(command.c_str());
        if (-1 == dir_err)
        {
            printf("Error of copying directory\n");
            exit(1);
        }
    }

    void copyFile(std::string source, std::string destination)
    {
        std::string command("cp -r "+source + " " +destination+"/");
        const int dir_err = system(command.c_str());
        if (-1 == dir_err)
        {
            printf("Error of copying file\n");
            exit(1);
        }
    }

    std::tuple<double, double> getUAtInterfaces(int edge, int element, int nG, int valType)
    {
        bool isMaster(auxUlti::checkMaster(element, edge));
        int locationPlus(-1), locationMinus(-1);
        double valPlus(0.0), valMinus(0.0);
        if (isMaster)
        {
            locationPlus = nG;
            locationMinus = nG + mathVar::nGauss + 1;
        }
        else
        {
            locationPlus = nG + mathVar::nGauss + 1;
            locationMinus = nG;
        }

        switch (valType)
        {
        case 1: //rho
        {
            valPlus = surfaceFields::rho[edge][locationPlus];
            valMinus = surfaceFields::rho[edge][locationMinus];
        }
        break;
        case 2: //rhou
        {
            valPlus = surfaceFields::rhou[edge][locationPlus];
            valMinus = surfaceFields::rhou[edge][locationMinus];
        }
        break;
        case 3: //rhov
        {
            valPlus = surfaceFields::rhov[edge][locationPlus];
            valMinus = surfaceFields::rhov[edge][locationMinus];
        }
        break;
        case 4: //rhoE
        {
            valPlus = surfaceFields::rhoE[edge][locationPlus];
            valMinus = surfaceFields::rhoE[edge][locationMinus];
        }
        break;
        default:
            break;
        }
        return std::make_tuple(valPlus, valMinus);
    }

    double getUPlusAtBC(int edge, int nG, int valType)
    {
        double val(0.0);
        switch (valType)
        {
        case 1: //rho
        {
            val = surfaceFields::rho[edge][nG];
        }
        break;
        case 2: //rhou
        {
            val = surfaceFields::rhou[edge][nG];
        }
        break;
        case 3: //rhov
        {
            val = surfaceFields::rhov[edge][nG];
        }
        break;
        case 4: //rhoE
        {
            val = surfaceFields::rhoE[edge][nG];
        }
        break;
        default:
            break;
        }
        return val;
    }


    std::tuple<double, double> getTAtInterfaces(int edge, int element, int nG)
    {
        bool isMaster(auxUlti::checkMaster(element, edge));
        int locationPlus(-1), locationMinus(-1);
        double TPlus(0.0), TMinus(0.0);
        if (isMaster)
        {
            locationPlus = nG;
            locationMinus = nG + mathVar::nGauss + 1;
        }
        else
        {
            locationPlus = nG + mathVar::nGauss + 1;
            locationMinus = nG;
        }

        TPlus = surfaceFields::T[edge][locationPlus];
        TMinus = surfaceFields::T[edge][locationMinus];
        return std::make_tuple(TPlus, TMinus);
    }

    double getTPlusAtBC(int edge, int nG)
    {
        return surfaceFields::T[edge][nG];
    }

    void saveUAtBCToSurfaceFields(int edge, int nG, std::vector<double>&UPlus, std::vector<double>&UMinus)
    {
        surfaceFields::rho[edge][nG]=UPlus[0];
        surfaceFields::rhou[edge][nG]=UPlus[1];
        surfaceFields::rhov[edge][nG]=UPlus[2];
        surfaceFields::rhoE[edge][nG]=UPlus[3];
        surfaceFields::rho[edge][nG+mathVar::nGauss+1]=UMinus[0];
        surfaceFields::rhou[edge][nG+mathVar::nGauss+1]=UMinus[1];
        surfaceFields::rhov[edge][nG+mathVar::nGauss+1]=UMinus[2];
        surfaceFields::rhoE[edge][nG+mathVar::nGauss+1]=UMinus[3];
    }

    namespace functionsOfParallelComputing {
    std::tuple<double,double> getGaussPointCoorsOfNeighborCell(int loc, int nG)
    {
        double a(parallelBuffer::aCoor[loc][nG]), b(parallelBuffer::bCoor[loc][nG]);
        return std::make_tuple(a,b);
    }

    void prepareParallelCase()
    {
        MPI_Init(NULL, NULL);
        int maxNode;
        MPI_Comm_size(MPI_COMM_WORLD, &maxNode);
        if (maxNode<systemVar::totalProc)
        {
            std::cout<<"Number of processors ("<<systemVar::totalProc<<") exceeds maximum number of available processors of system ("<<maxNode<<").\n";
            std::cout << "DGSolver will exit after you hit return.\n";
            exit(EXIT_FAILURE);
        }
        MPI_Comm_rank(MPI_COMM_WORLD, &systemVar::currentProc);
    }

    void sendReceiveU()
    {
        /*
         * List of tag:
         * - 100: rho -> 10i tuong duong iorder
         * - 200: rhou
         * - 300: rhov
         * - 400: rhoE
         * - 500: drhoX
         * - 600: drhouX
         * - 700: drhovX
         * - 800: drhoEX
         * - 900: drhoY
         * - 110: drhouY
         * - 120: drhovY
         * - 130: drhoEY
         * - 140: theta1
         * - 141: theta2
        */
        MPI_Request request_send, request_recv;
        MPI_Barrier(MPI_COMM_WORLD);
        int coef(1);
        int destination, source, cellId, neighborCellId, tag_sent, tag_recv;

        //Send
        for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge) {
            destination=meshVar::meshConnection[iBCedge][1];
            cellId=meshVar::meshConnection[iBCedge][0];
            neighborCellId=meshVar::meshConnection[iBCedge][2];

            //tag = neighborCellId*10 + destination
            if (destination>=0)
            {
                if (systemVar::totalProc<10)
                {
                    coef=10;
                }
                else if ((systemVar::totalProc>=10)&&(systemVar::totalProc<100))
                {
                    coef=100;
                }

                //define tags
                tag_sent = neighborCellId*coef + destination;

                for (int iorder = 0; iorder <= mathVar::orderElem; ++iorder) {
                    //rho
                    MPI_Isend(&rho[cellId][iorder], 1, MPI_DOUBLE, destination, tag_sent*1000+100+iorder, MPI_COMM_WORLD, &request_send);
                    //rhou
                    MPI_Isend(&rhou[cellId][iorder], 1, MPI_DOUBLE, destination, tag_sent*1000+200+iorder, MPI_COMM_WORLD, &request_send);
                    //rhov
                    MPI_Isend(&rhov[cellId][iorder], 1, MPI_DOUBLE, destination, tag_sent*1000+300+iorder, MPI_COMM_WORLD, &request_send);
                    //rhoE
                    MPI_Isend(&rhoE[cellId][iorder], 1, MPI_DOUBLE, destination, tag_sent*1000+400+iorder, MPI_COMM_WORLD, &request_send);
                }
                //theta1
                MPI_Isend(&theta1Arr[cellId], 1, MPI_DOUBLE, destination, tag_sent*1000+140, MPI_COMM_WORLD, &request_send);
                //theta2
                MPI_Isend(&theta2Arr[cellId], 1, MPI_DOUBLE, destination, tag_sent*1000+141, MPI_COMM_WORLD, &request_send);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        /*Chi receive khi nao toan bo cac processor send data xong*/

        //Receive
        for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge) {
            source=meshVar::meshConnection[iBCedge][1];
            cellId=meshVar::meshConnection[iBCedge][0];
            neighborCellId=meshVar::meshConnection[iBCedge][2];

            //tag = neighborCellId*10 + destination
            if (source>=0)
            {
                if (systemVar::totalProc<10)
                {
                    coef=10;
                }
                else if ((systemVar::totalProc>=10)&&(systemVar::totalProc<100))
                {
                    coef=100;
                }

                //define tags
                tag_recv = cellId*coef + systemVar::currentProc;

                for (int iorder = 0; iorder <= mathVar::orderElem; ++iorder) {
                    //rho
                    MPI_Irecv(&parallelBuffer::rho[iBCedge][iorder], 1, MPI_DOUBLE, source, tag_recv*1000+100+iorder, MPI_COMM_WORLD, &request_recv);
                    //rhou
                    MPI_Irecv(&parallelBuffer::rhou[iBCedge][iorder], 1, MPI_DOUBLE, source, tag_recv*1000+200+iorder, MPI_COMM_WORLD, &request_recv);
                    //rhov
                    MPI_Irecv(&parallelBuffer::rhov[iBCedge][iorder], 1, MPI_DOUBLE, source, tag_recv*1000+300+iorder, MPI_COMM_WORLD, &request_recv);
                    //rhoE
                    MPI_Irecv(&parallelBuffer::rhoE[iBCedge][iorder], 1, MPI_DOUBLE, source, tag_recv*1000+400+iorder, MPI_COMM_WORLD, &request_recv);
                }
                //theta1
                MPI_Irecv(&parallelBuffer::theta1[iBCedge], 1, MPI_DOUBLE, source, tag_recv*1000+140, MPI_COMM_WORLD, &request_recv);
                //theta2
                MPI_Irecv(&parallelBuffer::theta1[iBCedge], 1, MPI_DOUBLE, source, tag_recv*1000+141, MPI_COMM_WORLD, &request_recv);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    void sendReceivedU()
    {
        /*
         * List of tag:
         * - 100: rho -> 10i tuong duong iorder
         * - 200: rhou
         * - 300: rhov
         * - 400: rhoE
         * - 500: drhoX
         * - 600: drhouX
         * - 700: drhovX
         * - 800: drhoEX
         * - 900: drhoY
         * - 110: drhouY
         * - 120: drhovY
         * - 130: drhoEY
        */
        MPI_Request request_send, request_recv;
        MPI_Barrier(MPI_COMM_WORLD);
        int coef(1);
        int destination, source, cellId, neighborCellId, tag_sent, tag_recv;

        //Send
        for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge) {
            destination=meshVar::meshConnection[iBCedge][1];
            cellId=meshVar::meshConnection[iBCedge][0];
            neighborCellId=meshVar::meshConnection[iBCedge][2];

            //tag = neighborCellId*10 + destination
            if (destination>=0)
            {
                if (systemVar::totalProc<10)
                {
                    coef=10;
                }
                else if ((systemVar::totalProc>=10)&&(systemVar::totalProc<100))
                {
                    coef=100;
                }

                //define tags
                tag_sent = neighborCellId*coef + destination;

                for (int iorder = 0; iorder <= mathVar::orderElem; ++iorder) {
                    if (!flowProperties::massDiffusion)
                    {
                        if (systemVar::auxVariables==1)
                        {
                            //drhoX
                            MPI_Isend(&BR1Vars::rhoX[cellId][iorder], 1, MPI_DOUBLE, destination, tag_sent*1000+500+iorder, MPI_COMM_WORLD, &request_send);

                            //drhoY
                            MPI_Isend(&BR1Vars::rhoY[cellId][iorder], 1, MPI_DOUBLE, destination, tag_sent*1000+900+iorder, MPI_COMM_WORLD, &request_send);
                        }
                        else if (systemVar::auxVariables==2)
                        {
                            //Chua lam cho method BR2
                        }
                    }

                    if (systemVar::auxVariables==1)
                    {
                        //drhouX
                        MPI_Isend(&BR1Vars::rhouX[cellId][iorder], 1, MPI_DOUBLE, destination, tag_sent*1000+600+iorder, MPI_COMM_WORLD, &request_send);
                        //drhovX
                        MPI_Isend(&BR1Vars::rhovX[cellId][iorder], 1, MPI_DOUBLE, destination, tag_sent*1000+700+iorder, MPI_COMM_WORLD, &request_send);
                        //drhoEX
                        MPI_Isend(&BR1Vars::rhoEX[cellId][iorder], 1, MPI_DOUBLE, destination, tag_sent*1000+800+iorder, MPI_COMM_WORLD, &request_send);

                        //drhouY
                        MPI_Isend(&BR1Vars::rhouY[cellId][iorder], 1, MPI_DOUBLE, destination, tag_sent*1000+110+iorder, MPI_COMM_WORLD, &request_send);
                        //drhovY
                        MPI_Isend(&BR1Vars::rhovY[cellId][iorder], 1, MPI_DOUBLE, destination, tag_sent*1000+120+iorder, MPI_COMM_WORLD, &request_send);
                        //drhoEY
                        MPI_Isend(&BR1Vars::rhoEY[cellId][iorder], 1, MPI_DOUBLE, destination, tag_sent*1000+130+iorder, MPI_COMM_WORLD, &request_send);
                    }
                    else if (systemVar::auxVariables==2)
                    {
                        //Chua lam cho method BR2
                    }
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        /*Chi receive khi nao toan bo cac processor send data xong*/

        //Receive
        for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge) {
            source=meshVar::meshConnection[iBCedge][1];
            cellId=meshVar::meshConnection[iBCedge][0];
            neighborCellId=meshVar::meshConnection[iBCedge][2];

            //tag = neighborCellId*10 + destination
            if (source>=0)
            {
                if (systemVar::totalProc<10)
                {
                    coef=10;
                }
                else if ((systemVar::totalProc>=10)&&(systemVar::totalProc<100))
                {
                    coef=100;
                }

                //define tags
                tag_recv = cellId*coef + systemVar::currentProc;

                for (int iorder = 0; iorder <= mathVar::orderElem; ++iorder) {
                    if (!flowProperties::massDiffusion)
                    {
                        if (systemVar::auxVariables==1)
                        {
                            //drhoX
                            MPI_Irecv(&parallelBuffer::drhoX[iBCedge][iorder], 1, MPI_DOUBLE, source, tag_recv*1000+500+iorder, MPI_COMM_WORLD, &request_recv);

                            //drhoY
                            MPI_Irecv(&parallelBuffer::drhoY[iBCedge][iorder], 1, MPI_DOUBLE, source, tag_recv*1000+900+iorder, MPI_COMM_WORLD, &request_recv);
                        }
                        else if (systemVar::auxVariables==2)
                        {
                            //Chua lam cho method BR2
                        }
                    }

                    if (systemVar::auxVariables==1)
                    {
                        //drhouX
                        MPI_Irecv(&parallelBuffer::drhouX[iBCedge][iorder], 1, MPI_DOUBLE, source, tag_recv*1000+600+iorder, MPI_COMM_WORLD, &request_recv);
                        //drhovX
                        MPI_Irecv(&parallelBuffer::drhovX[iBCedge][iorder], 1, MPI_DOUBLE, source, tag_recv*1000+700+iorder, MPI_COMM_WORLD, &request_recv);
                        //drhoEX
                        MPI_Irecv(&parallelBuffer::drhoEX[iBCedge][iorder], 1, MPI_DOUBLE, source, tag_recv*1000+800+iorder, MPI_COMM_WORLD, &request_recv);

                        //drhouY
                        MPI_Irecv(&parallelBuffer::drhouY[iBCedge][iorder], 1, MPI_DOUBLE, source, tag_recv*1000+110+iorder, MPI_COMM_WORLD, &request_recv);
                        //drhovY
                        MPI_Irecv(&parallelBuffer::drhovY[iBCedge][iorder], 1, MPI_DOUBLE, source, tag_recv*1000+120+iorder, MPI_COMM_WORLD, &request_recv);
                        //drhoEY
                        MPI_Irecv(&parallelBuffer::drhoEY[iBCedge][iorder], 1, MPI_DOUBLE, source, tag_recv*1000+130+iorder, MPI_COMM_WORLD, &request_recv);
                    }
                    else if (systemVar::auxVariables==2)
                    {
                        //Chua lam cho method BR2
                    }
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    void sendReceivedRho()
    {
        /*
         * List of tag:
         * - 100: rho -> 10i tuong duong iorder
         * - 500: drhoX
         * - 900: drhoY
        */
        MPI_Request request_send, request_recv;
        MPI_Barrier(MPI_COMM_WORLD);
        int coef(1);
        int destination, source, cellId, neighborCellId, tag_sent, tag_recv;

        //Send
        for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge) {
            destination=meshVar::meshConnection[iBCedge][1];
            cellId=meshVar::meshConnection[iBCedge][0];
            neighborCellId=meshVar::meshConnection[iBCedge][2];

            //tag = neighborCellId*10 + destination
            if (destination>=0)
            {
                if (systemVar::totalProc<10)
                {
                    coef=10;
                }
                else if ((systemVar::totalProc>=10)&&(systemVar::totalProc<100))
                {
                    coef=100;
                }

                //define tags
                tag_sent = neighborCellId*coef + destination;

                for (int iorder = 0; iorder <= mathVar::orderElem; ++iorder) {
                    if (systemVar::auxVariables==1)
                    {
                        //drhoX
                        MPI_Isend(&BR1Vars::rhoX[iBCedge][iorder], 1, MPI_DOUBLE, destination, tag_sent*1000+500+iorder, MPI_COMM_WORLD, &request_send);

                        //drhoY
                        MPI_Isend(&BR1Vars::rhoY[iBCedge][iorder], 1, MPI_DOUBLE, destination, tag_sent*1000+900+iorder, MPI_COMM_WORLD, &request_send);
                    }
                    else if (systemVar::auxVariables==2)
                    {
                        //Chua lam cho method BR2
                    }
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        /*Chi receive khi nao toan bo cac processor send data xong*/

        //Receive
        for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge) {
            source=meshVar::meshConnection[iBCedge][1];
            cellId=meshVar::meshConnection[iBCedge][0];
            neighborCellId=meshVar::meshConnection[iBCedge][2];

            //tag = neighborCellId*10 + destination
            if (source>=0)
            {
                if (systemVar::totalProc<10)
                {
                    coef=10;
                }
                else if ((systemVar::totalProc>=10)&&(systemVar::totalProc<100))
                {
                    coef=100;
                }

                //define tags
                tag_recv = cellId*coef + systemVar::currentProc;

                for (int iorder = 0; iorder <= mathVar::orderElem; ++iorder) {
                    if (systemVar::auxVariables==1)
                    {
                        //drhoX
                        MPI_Irecv(&parallelBuffer::drhoX[iBCedge][iorder], 1, MPI_DOUBLE, source, tag_recv*1000+500+iorder, MPI_COMM_WORLD, &request_recv);

                        //drhoY
                        MPI_Irecv(&parallelBuffer::drhoY[iBCedge][iorder], 1, MPI_DOUBLE, source, tag_recv*1000+900+iorder, MPI_COMM_WORLD, &request_recv);
                    }
                    else if (systemVar::auxVariables==2)
                    {
                        //Chua lam cho method BR2
                    }
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    void sendReceiveMeshData()
    {
        MPI_Request request_send, request_recv;
        MPI_Barrier(MPI_COMM_WORLD);
        int coef(1);
        int destination, source, cellId, neighborCellId, tag_sent, tag_recv;

        //Send
        for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge) {
            destination=meshVar::meshConnection[iBCedge][1];
            cellId=meshVar::meshConnection[iBCedge][0];
            neighborCellId=meshVar::meshConnection[iBCedge][2];

            //tag = neighborCellId*10 + destination
            if (destination>=0)
            {
                if (systemVar::totalProc<10)
                {
                    coef=10;
                }
                else if ((systemVar::totalProc>=10)&&(systemVar::totalProc<100))
                {
                    coef=100;
                }

                //define tags
                tag_sent = neighborCellId*coef + destination;

                int elemType(auxUlti::checkType(cellId)), ptId(0);

                //send toa do 3 dinh cua cell
                for (int vertex = 0; vertex < elemType; ++vertex) {
                    ptId=meshVar::Elements2D[cellId][vertex];
                    //tag*10 + 1: gui toa do x
                    //tag*10 + 2: gui toa do y
                    MPI_Isend(&meshVar::Points[ptId][0], 1, MPI_DOUBLE, destination, tag_sent*10+1, MPI_COMM_WORLD, &request_send);
                    MPI_Isend(&meshVar::Points[ptId][1], 1, MPI_DOUBLE, destination, tag_sent*10+2, MPI_COMM_WORLD, &request_send);
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        /*Chi receive khi nao toan bo cac processor send data xong*/

        //Receive
        for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge) {
            source=meshVar::meshConnection[iBCedge][1];
            cellId=meshVar::meshConnection[iBCedge][0];
            neighborCellId=meshVar::meshConnection[iBCedge][2];

            //tag = neighborCellId*10 + destination
            if (source>=0)
            {
                if (systemVar::totalProc<10)
                {
                    coef=10;
                }
                else if ((systemVar::totalProc>=10)&&(systemVar::totalProc<100))
                {
                    coef=100;
                }

                //define tags
                tag_recv = cellId*coef + systemVar::currentProc;
                int elemType(auxUlti::checkType(cellId));

                //receive toa do 3 dinh cua cell
                for (int vertex = 0; vertex < elemType; ++vertex) {
                    //tag*10 + 1: gui toa do x
                    //tag*10 + 2: gui toa do y
                    MPI_Irecv(&parallelBuffer::xCoor[iBCedge][vertex], 1, MPI_DOUBLE, source, tag_recv*10+1, MPI_COMM_WORLD, &request_recv);
                    MPI_Irecv(&parallelBuffer::yCoor[iBCedge][vertex], 1, MPI_DOUBLE, source, tag_recv*10+2, MPI_COMM_WORLD, &request_recv);
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    void sendString(std::string content, int destination, int tag)
    {
        MPI_Send(content.c_str(), content.size(), MPI_CHAR, destination, tag, MPI_COMM_WORLD);
    }

    std::string receiveString(int source, int tag)
    {
        //From internet
        MPI::Status status;
        MPI::COMM_WORLD.Probe(source, tag, status);
        int l = status.Get_count(MPI_CHAR);
        char *buf = new char[l];
        MPI::COMM_WORLD.Recv(buf, l, MPI_CHAR, source, tag, status);
        std::string str(buf, l);
        delete[] buf;
        return str;
    }
    }

    void checkInforBeforeRunning()
    {
        std::string runOrNot;
        if (systemVar::currentProc==0)
        {
            //Check subsonic
            refValues::subsonic = auxUlti::checkSubSonic();
            //Check case's information
            message::showCaseInformations();
            std::cout<<"Do you want to continue? <y/n> ";
            std::cin>>runOrNot;
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                auxUlti::functionsOfParallelComputing::sendString(runOrNot,irank,2);
            }
        }
        else {
            runOrNot=auxUlti::functionsOfParallelComputing::receiveString(0,2);
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        if (runOrNot.compare("n") == 0){
            std::cout << "DGSolver is exitting.\n";
            exit(EXIT_SUCCESS);
        }
    }

    void getCommand()
    {
        if (systemVar::currentProc==0)
        {
            std::cout << ">> ";
            std::cin >> systemVar::cmd;

            //Send command to all processors
            for (int irank=1;irank<systemVar::totalProc;irank++) {
                auxUlti::functionsOfParallelComputing::sendString(systemVar::cmd,irank,1);
            }
        }
        else {
            systemVar::cmd=auxUlti::functionsOfParallelComputing::receiveString(0,1);
        }
    }
}
