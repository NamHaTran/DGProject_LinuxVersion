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
#include <fstream>
#include <sstream>
#include <mpi.h>

#include <iostream>

int calcArrId(int id1, int id2, int length)
{
    return (id1*length+id2);
}

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

    int lookForDataOfKeyword(std::string fileLoc, std::string inputKeyWord)
    {
        std::ifstream FileFlux(fileLoc.c_str());
        std::string line(" "), keyWord;
        int length(0);
        while (std::getline(FileFlux, line))
        {
            std::istringstream line2str(line);
            std::vector<std::string> ptr;
            //Split <line2str> stringstream into array of words
            while ((line2str >> keyWord))
            {
                ptr.push_back(keyWord);
            }

            int numWrd = static_cast<int>(ptr.size());

            if (numWrd >= 2)
            {
                std::istringstream strdata(ptr[1]);
                if (ptr[0].compare(inputKeyWord) == 0)
                {
                    strdata >> length;
                    break;
                }
            }
        }
        return length;
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
        int nanb(calcArrId(na,nb,mathVar::nGauss+1));
        a = mathVar::GaussPts[nanb][0];
        b= mathVar::GaussPts[nanb][1];
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
                TInf = bcValues::TBCFixed[i];
                uInf = bcValues::uBCFixed[i];
                vInf = bcValues::vBCFixed[i];
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

    /*
    void resize1DArray(double*Array, int row, double initialValue)
    {
        Array = new double [row];
        for (int i=0; i<row; i++)
        {
            Array[i]=initialValue;
        }
    }

    void resize1DIntArray(int*Array, int row, int initialValue)
    {
        Array = new int [row];
        for (int i=0; i<row; i++)
        {
            Array[i]=initialValue;
        }
    }
    */
    void initialize1DArray(double*Array, int row, double initialValue)
    {
        for (int i=0; i<row; i++)
        {
         Array[i]=initialValue;
        }
    }

    void initialize1DIntArray(int*Array, int row, int initialValue)
    {
        for (int i=0; i<row; i++)
        {
         Array[i]=initialValue;
        }
    }

    double** resize2DArray(int row, int column, double initialValue)
	{
        /*Ham tao mang 2D gom cac o nho lien tuc (contiguous memory location)*/
        /*
        double**Array = new double*[row];
        Array[0]= new double[row*column];
        for(int j= 1; j < row; j++) {
          Array[j]=Array[j-1] + column;
        }
        */
        double** Array = new double*[row];
        for(int i = 0; i < row; ++i)
            Array[i] = new double[column];

        //Dat gia tri ban dau cua array la 0.0
        for(int j= 0; j < row; j++) {
            for(int k= 0; k < column; k++) {
              Array[j][k]= initialValue;
            }
        }
        return Array;
	}

    int** resize2DIntArray(int row, int column, int initialValue)
	{
        /*
        int**Array = new int*[row];
        Array[0]= new int[row*column];
        for(int j= 1; j < row; j++) {
          //Array[j]= &Array[0][j*column];
            Array[j]=Array[j-1] + column;
        }*/
        int** Array = new int*[row];
        for(int i = 0; i < row; ++i)
            Array[i] = new int[column];

        //Dat gia tri ban dau cua array la 0
        for(int j= 0; j < row; j++) {
            for(int k= 0; k < column; k++) {
              Array[j][k]= initialValue;
            }
        }
        return Array;
	}

    void resize3DArray(std::vector<std::vector<std::vector<double>>> &Array, int direct1, int direct2, int direct3)
	{
		Array.resize(direct1);
        for (int i = 0; i < direct1; ++i)
		{
            Array[i].resize(direct2);
            for (int j = 0; j < direct2; j++)
			{
                Array[i][j].resize(direct3, 0.0);
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
                    if (bcType==4) //boundary type 4 at file boundaryPatch: matched
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
        meshVar::geoCenter=auxUlti::resize2DArray(meshVar::nelem2D,2,0.0);

        meshVar::cellArea = new double [meshVar::nelem2D];
        meshVar::cellSize = new double [meshVar::nelem2D];
        meshVar::localCellSize = new double [meshVar::nelem2D];

        meshVar::dxa = auxUlti::resize2DArray(meshVar::nelem2D, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);
        meshVar::dya = auxUlti::resize2DArray(meshVar::nelem2D, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);
        meshVar::dxb = auxUlti::resize2DArray(meshVar::nelem2D, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);
        meshVar::dyb = auxUlti::resize2DArray(meshVar::nelem2D, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);

        meshVar::J2D = auxUlti::resize2DArray(meshVar::nelem2D, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);

        //meshVar::J1D = new double [meshVar::inpoedCount];
        meshVar::J1D = new double [meshVar::inpoedCount];
        auxUlti::initialize1DArray(meshVar::J1D, meshVar::inpoedCount, 0.0);
        //auxUlti::resize2DArray(meshVar::J1D, meshVar::inpoedCount, 2);

        meshVar::edgeGaussPoints_a = auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);
        meshVar::edgeGaussPoints_b = auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);

        rho = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
        rhou = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
        rhov = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
        rhoE = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);

        rho0 = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
        rhou0 = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
        rhov0 = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
        rhoE0 = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);

        rhoResArr = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
        rhouResArr = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
        rhovResArr = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
        rhoEResArr = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);

        SurfaceBCFields::uBc = new double [meshVar::numBCEdges];
        auxUlti::initialize1DArray(SurfaceBCFields::uBc, meshVar::numBCEdges, 0.0);
        SurfaceBCFields::vBc = new double [meshVar::numBCEdges];
        auxUlti::initialize1DArray(SurfaceBCFields::vBc, meshVar::numBCEdges, 0.0);
        SurfaceBCFields::TBc = new double [meshVar::numBCEdges];
        auxUlti::initialize1DArray(SurfaceBCFields::TBc, meshVar::numBCEdges, 0.0);
        meshVar::distanceFromCentroidToBCEdge = new double [meshVar::numBCEdges];

        auxUlti::initialize1DArray(meshVar::distanceFromCentroidToBCEdge, meshVar::numBCEdges, 0.0);

        rhoN = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
        rhouN = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
        rhovN = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
        rhoEN = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);

        /*Khong giai phuong trinh phu khi dong inviscid*/
        if (flowProperties::viscous)
        {
            surfaceFields::aux_rho = auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);
            surfaceFields::aux_rhou = auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);
            surfaceFields::aux_rhov = auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);
            surfaceFields::aux_rhoE = auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);
        }

        surfaceFields::rho= auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);
        surfaceFields::rhou= auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);
        surfaceFields::rhov= auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);
        surfaceFields::rhoE= auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);

        if (flowProperties::viscous)
        {
            if (systemVar::auxVariables==1)
            {
                BR1Vars::rhoX = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
                BR1Vars::rhouX = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
                BR1Vars::rhovX = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
                BR1Vars::rhoEX = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);

                BR1Vars::rhoY = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
                BR1Vars::rhouY = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
                BR1Vars::rhovY = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
                BR1Vars::rhoEY = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);
            }
            else if (systemVar::auxVariables==2)
            {
                //Bo BR2
            }
        }

        surfaceFields::invis_rho = auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);
        surfaceFields::invis_rhou = auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);
        surfaceFields::invis_rhov = auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);
        surfaceFields::invis_rhoE = auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);

        if (flowProperties::viscous)
        {
            surfaceFields::Vis_rho = auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);
            surfaceFields::Vis_rhou = auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);
            surfaceFields::Vis_rhov = auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);
            surfaceFields::Vis_rhoE = auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);
        }

        surfaceFields::T = auxUlti::resize2DArray(meshVar::inpoedCount, 2 * (mathVar::nGauss + 1),0.0);

        //theta1Arr=new double [meshVar::nelem2D];
        //theta2Arr=new double [meshVar::nelem2D];
        //debug::minRhoArr.resize(meshVar::nelem2D);
        //debug::minRhoeArr.resize(meshVar::nelem2D);
        limitVal::troubleCellsMarker=new bool [meshVar::nelem2D];
        theta1Arr=new double[meshVar::nelem2D];
        auxUlti::initialize1DArray(theta1Arr, meshVar::nelem2D, 1.0);
        theta2Arr=new double[meshVar::nelem2D];
        auxUlti::initialize1DArray(theta2Arr, meshVar::nelem2D, 1.0);

        LxFConst = new double [meshVar::inpoedCount];
        DiffusiveFluxConst = new double [meshVar::inpoedCount];

        stiffMatrixCoeffs = auxUlti::resize2DArray(meshVar::nelem2D, mathVar::orderElem + 1,0.0);

        volumeFields::rhoVolGauss = auxUlti::resize2DArray(meshVar::nelem2D, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);
        volumeFields::rhouVolGauss = auxUlti::resize2DArray(meshVar::nelem2D, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);
        volumeFields::rhovVolGauss = auxUlti::resize2DArray(meshVar::nelem2D, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);
        volumeFields::rhoEVolGauss = auxUlti::resize2DArray(meshVar::nelem2D, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);
        volumeFields::drhoXVolGauss = auxUlti::resize2DArray(meshVar::nelem2D, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);
        volumeFields::drhoYVolGauss = auxUlti::resize2DArray(meshVar::nelem2D, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);
        volumeFields::T = auxUlti::resize2DArray(meshVar::nelem2D, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);

        //meshVar::adressOfBCVals.resize(meshVar::numBCEdges);
        meshVar::neighboringElements = auxUlti::resize2DIntArray(meshVar::nelem2D, 4,0.0);

        //Resize mathVar array
        mathVar::wGauss=new double [mathVar::nGauss+1];
        mathVar::xGauss=new double [mathVar::nGauss+1];
        mathVar::wGaussLobatto=new double [mathVar::nGauss+1];
        mathVar::xGaussLobatto=new double [mathVar::nGauss+1];
        mathVar::B=new double [mathVar::orderElem+1];
        mathVar::dBa=new double [mathVar::orderElem+1];
        mathVar::dBb=new double [mathVar::orderElem+1];

        mathVar::BPts_Quad = auxUlti::resize2DArray(mathVar::orderElem+1, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);
        mathVar::dBaPts_Quad = auxUlti::resize2DArray(mathVar::orderElem+1, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);
        mathVar::dBbPts_Quad = auxUlti::resize2DArray(mathVar::orderElem+1, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);
        mathVar::BPts_Tri = auxUlti::resize2DArray(mathVar::orderElem+1, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);
        mathVar::dBaPts_Tri = auxUlti::resize2DArray(mathVar::orderElem+1, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);
        mathVar::dBbPts_Tri = auxUlti::resize2DArray(mathVar::orderElem+1, (mathVar::nGauss + 1)*(mathVar::nGauss + 1),0.0);
        mathVar::GaussPts = auxUlti::resize2DArray((mathVar::nGauss + 1)*(mathVar::nGauss + 1),2,0.0);
        mathVar::wGaussPts = auxUlti::resize2DArray((mathVar::nGauss + 1)*(mathVar::nGauss + 1),2,0.0);
        mathVar::GaussLobattoPts = auxUlti::resize2DArray((mathVar::nGauss + 1)*(mathVar::nGauss + 1),2,0.0);
        mathVar::wGaussLobattoPts = auxUlti::resize2DArray((mathVar::nGauss + 1)*(mathVar::nGauss + 1),2,0.0);

        //Buffer (parallel computing)
        //Conservative variables
        parallelBuffer::rho= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::orderElem + 1,0.0);
        parallelBuffer::rhou= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::orderElem + 1,0.0);
        parallelBuffer::rhov= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::orderElem + 1,0.0);
        parallelBuffer::rhoE= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::orderElem + 1,0.0);

        //Auxilary variables
        if (flowProperties::viscous)
        {
            parallelBuffer::drhoX= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::orderElem + 1,0.0);
            parallelBuffer::drhouX= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::orderElem + 1,0.0);
            parallelBuffer::drhovX= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::orderElem + 1,0.0);
            parallelBuffer::drhoEX= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::orderElem + 1,0.0);
            parallelBuffer::drhoY= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::orderElem + 1,0.0);
            parallelBuffer::drhouY= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::orderElem + 1,0.0);
            parallelBuffer::drhovY= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::orderElem + 1,0.0);
            parallelBuffer::drhoEY= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::orderElem + 1,0.0);
        }

        //Surface Gauss point coordinates
        parallelBuffer::aCoor= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,0.0);
        parallelBuffer::bCoor= auxUlti::resize2DArray(meshVar::numBCEdges, mathVar::nGauss + 1,0.0);
        //Neighbor cell vertexes
        parallelBuffer::xCoor= auxUlti::resize2DArray(meshVar::numBCEdges, 4,0.0);
        parallelBuffer::yCoor= auxUlti::resize2DArray(meshVar::numBCEdges, 4,0.0);

        //parallelBuffer::theta1 = new double [meshVar::numBCEdges];
        //parallelBuffer::theta2 = new double [meshVar::numBCEdges];
        //parallelBuffer::elemType = new int [meshVar::numBCEdges];
        parallelBuffer::theta1 = new double [meshVar::numBCEdges];
        auxUlti::initialize1DArray(parallelBuffer::theta1, meshVar::numBCEdges, 1.0);
        parallelBuffer::theta2 = new double [meshVar::numBCEdges];
        auxUlti::initialize1DArray(parallelBuffer::theta2, meshVar::numBCEdges, 1.0);
        parallelBuffer::elemType = new int [meshVar::numBCEdges];
        auxUlti::initialize1DIntArray(parallelBuffer::elemType,meshVar::numBCEdges,3);

        //Mesh connection array
        //auxUlti::resize2DIntArray(meshVar::meshConnection,meshVar::numBCEdges,2);

        //For Maxwell-Smoluchowsky BC
        meshVar::normProjectionOfCenterToBCEdge_realSysCoor = auxUlti::resize2DArray(meshVar::numBCEdges, 2,0.0);
        meshVar::normProjectionOfCenterToBCEdge_standardSysCoor = auxUlti::resize2DArray(meshVar::numBCEdges, 2,0.0);
	}

    int getAdressOfBCEdgesOnBCValsArray(int globalEdgeId)
	{
        return meshVar::inpoed[globalEdgeId][4];
	}

    int getGlobalEdgeIdFromLocalBCEdgeId(int localBCEdgeId)
    {
        return SurfaceBCFields::localGlobalBCEdgesMatching[localBCEdgeId];
    }

	std::tuple<double, double> getCellCentroid(int element)
	{
		double xC(meshVar::geoCenter[element][0]), yC(meshVar::geoCenter[element][1]);
		return std::make_tuple(xC, yC);
	}

	void clear1DIntVector(std::vector<int>&vector)
	{
		//Use erase function to clear vector
		vector.erase(vector.begin(), vector.end());
		//Shrink to fit
		vector.shrink_to_fit();
	}

    void clear2DIntVector(std::vector<std::vector<int>>&vector)
    {
        int numRow(static_cast<int>(vector.size()));
        //Use erase function to clear vector
        for (int row = 0; row < numRow; row++) {
            int length(static_cast<int>(vector[row].size()));
            vector[row].erase(vector[row].begin(), vector[row].begin() + length);
            vector[row].shrink_to_fit();
        }

        //Shrink to fit
        vector.erase(vector.begin(), vector.end());
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
                for (int iproc = 1; iproc < systemVar::totalProc; iproc++)
                {
                    Loc = systemVar::wD + "/CASES/" + systemVar::caseName + "/Processor" + std::to_string(iproc) + tecplotName + iter_str;
                    auxUlti::createFolder(Loc,true);
                    auxUlti::functionsOfParallelComputing::sendString(Loc,iproc,1);
                }
                Loc = systemVar::wD + "/CASES/" + systemVar::caseName + "/Processor0" + tecplotName + iter_str;
                auxUlti::createFolder(Loc,true);
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

    void getVectorUMeanOfCell(int element, std::vector<double> &U)
    {
        U[0]=rho[element][0];
        U[1]=rhou[element][0];
        U[2]=rhov[element][0];
        U[3]=rhoE[element][0];
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

    /*void sendRecvDiscretedVar(double**Var,double**Buffer,int order)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        //Declare dynamic array
        MPI_Request *request_send = new MPI_Request[meshVar::numBCEdges*(mathVar::orderElem+1)],
                *request_recv = new MPI_Request[meshVar::numBCEdges*(mathVar::orderElem+1)];

        int coef(1);
        int destination, source, cellId, neighborCellId, tag_sent, tag_recv, sentCount(0), recvCount(0);

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

                MPI_Isend(&Var[cellId][order], 1, MPI_DOUBLE, destination, tag_sent*1000+100+order, MPI_COMM_WORLD, &request_send[sentCount]);
                sentCount++;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        //Chi receive khi nao toan bo cac processor send data xong

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

                MPI_Irecv(&Buffer[iBCedge][order], 1, MPI_DOUBLE, source, tag_recv*1000+100+order, MPI_COMM_WORLD, &request_recv[recvCount]);
                recvCount++;
            }
        }

        //Wait for all request fullfilled
        MPI_Waitall(sentCount,request_send,MPI_STATUSES_IGNORE);
        MPI_Waitall(recvCount,request_recv,MPI_STATUSES_IGNORE);
        //MPI_Barrier(MPI_COMM_WORLD);

        delete [] request_send;
        delete [] request_recv;
    }*/

    void sendRecvDiscretedVar(int sendingProc, int receivingProc, double**Var,double**Buffer,int order)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        //Declare dynamic array
        MPI_Request *request = new MPI_Request[meshVar::numBCEdges];

        int destination, source, cellId, neighborCellId, tag_sent, tag_recv, requestCount(0);

        if (systemVar::currentProc==sendingProc) //current proc send, neighbor receive
        {
            //Send
            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge)
            {
                destination=meshVar::meshConnection[iBCedge][1];
                cellId=meshVar::meshConnection[iBCedge][0];
                neighborCellId=meshVar::meshConnection[iBCedge][2];

                //tag = neighborCellId*10 + destination
                if (destination==receivingProc)
                {
                    //define tags
                    tag_sent = neighborCellId;

                    MPI_Isend(&Var[cellId][order], 1, MPI_DOUBLE, destination, tag_sent*1000+100+order, MPI_COMM_WORLD, &request[requestCount]);
                    requestCount++;
                }
            }
            //Wait for all request fullfilled
            MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
        }
        else if (systemVar::currentProc==receivingProc)
        {
            //Receive
            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge) {
                source=meshVar::meshConnection[iBCedge][1];
                cellId=meshVar::meshConnection[iBCedge][0];
                neighborCellId=meshVar::meshConnection[iBCedge][2];

                //tag = neighborCellId*10 + destination
                if (source==sendingProc)
                {
                    //define tags
                    tag_recv = cellId;

                    MPI_Irecv(&Buffer[iBCedge][order], 1, MPI_DOUBLE, source, tag_recv*1000+100+order, MPI_COMM_WORLD, &request[requestCount]);
                    requestCount++;
                }
            }
            //Wait for all request fullfilled
            MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        delete [] request;
    }

    /*void sendRecvTheta(double*thetaArray,double*Buffer)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        int coef(1);
        int destination, source, cellId, neighborCellId, tag_sent, tag_recv, sentCount(0), recvCount(0);

        //Declare dynamic array
        MPI_Request *request_send = new MPI_Request[meshVar::numBCEdges],
                *request_recv = new MPI_Request[meshVar::numBCEdges];

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

                MPI_Isend(&thetaArray[cellId], 1, MPI_DOUBLE, destination, tag_sent*1000+140, MPI_COMM_WORLD, &request_send[sentCount]);
                sentCount++;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        //Chi receive khi nao toan bo cac processor send data xong

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

                MPI_Irecv(&Buffer[iBCedge], 1, MPI_DOUBLE, source, tag_recv*1000+140, MPI_COMM_WORLD, &request_recv[recvCount]);
                recvCount++;
            }
        }

        //Wait for all request fullfilled
        MPI_Waitall(sentCount,request_send,MPI_STATUSES_IGNORE);
        MPI_Waitall(recvCount,request_recv,MPI_STATUSES_IGNORE);
        //MPI_Barrier(MPI_COMM_WORLD);

        delete [] request_send;
        delete [] request_recv;
    }*/

    void sendRecvTheta(int sendingProc, int receivingProc, double*thetaArray, double*Buffer)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        //Declare dynamic array
        MPI_Request *request = new MPI_Request[meshVar::numBCEdges];

        int destination, source, cellId, neighborCellId, tag_sent, tag_recv, requestCount(0);

        if (systemVar::currentProc==sendingProc) //current proc send, neighbor receive
        {
            //Send
            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge)
            {
                destination=meshVar::meshConnection[iBCedge][1];
                cellId=meshVar::meshConnection[iBCedge][0];
                neighborCellId=meshVar::meshConnection[iBCedge][2];

                //tag = neighborCellId*10 + destination
                if (destination==receivingProc)
                {
                    //define tags
                    tag_sent = neighborCellId;

                    MPI_Isend(&thetaArray[cellId], 1, MPI_DOUBLE, destination, tag_sent*1000+140, MPI_COMM_WORLD, &request[requestCount]);
                    requestCount++;
                }
            }
            //Wait for all request fullfilled
            MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
        }
        else if (systemVar::currentProc==receivingProc)
        {
            //Receive
            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge) {
                source=meshVar::meshConnection[iBCedge][1];
                cellId=meshVar::meshConnection[iBCedge][0];
                neighborCellId=meshVar::meshConnection[iBCedge][2];

                //tag = neighborCellId*10 + destination
                if (source==sendingProc)
                {
                    //define tags
                    tag_recv = cellId;

                    MPI_Irecv(&Buffer[iBCedge], 1, MPI_DOUBLE, source, tag_recv*1000+140, MPI_COMM_WORLD, &request[requestCount]);
                    requestCount++;
                }
            }
            //Wait for all request fullfilled
            MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        delete [] request;
    }

    void sendReceiveU()
    {
        int sendingProc, receivingProc;
        for (int i=0; i<systemVar::sendRecvOrder_length; i++)
        {
            sendingProc=systemVar::sendRecvOrder[i][0];
            receivingProc=systemVar::sendRecvOrder[i][1];
            for (int iorder = 0; iorder <= mathVar::orderElem; ++iorder)
            {
                auxUlti::functionsOfParallelComputing::sendRecvDiscretedVar(sendingProc, receivingProc, rho,parallelBuffer::rho,iorder);
                auxUlti::functionsOfParallelComputing::sendRecvDiscretedVar(sendingProc, receivingProc, rhou,parallelBuffer::rhou,iorder);
                auxUlti::functionsOfParallelComputing::sendRecvDiscretedVar(sendingProc, receivingProc, rhov,parallelBuffer::rhov,iorder);
                auxUlti::functionsOfParallelComputing::sendRecvDiscretedVar(sendingProc, receivingProc, rhoE,parallelBuffer::rhoE,iorder);
            }
        }
    }

    void sendReceivedU()
    {
        if (systemVar::auxVariables==1)
        {
            int sendingProc, receivingProc;
            for (int i=0; i<systemVar::sendRecvOrder_length; i++)
            {
                sendingProc=systemVar::sendRecvOrder[i][0];
                receivingProc=systemVar::sendRecvOrder[i][1];
                for (int iorder = 0; iorder <= mathVar::orderElem; ++iorder)
                {
                    auxUlti::functionsOfParallelComputing::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhoX,parallelBuffer::drhoX,iorder);
                    auxUlti::functionsOfParallelComputing::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhoY,parallelBuffer::drhoY,iorder);

                    auxUlti::functionsOfParallelComputing::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhouX,parallelBuffer::drhouX,iorder);
                    auxUlti::functionsOfParallelComputing::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhouY,parallelBuffer::drhouY,iorder);

                    auxUlti::functionsOfParallelComputing::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhovX,parallelBuffer::drhovX,iorder);
                    auxUlti::functionsOfParallelComputing::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhovY,parallelBuffer::drhovY,iorder);

                    auxUlti::functionsOfParallelComputing::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhoEX,parallelBuffer::drhoEX,iorder);
                    auxUlti::functionsOfParallelComputing::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhoEY,parallelBuffer::drhoEY,iorder);
                }
            }
        }
        else if (systemVar::auxVariables==2)
        {
            //Chua lam cho method BR2
        }
    }

    void sendReceivedRho()
    {
        if (systemVar::auxVariables==1)
        {
            int sendingProc, receivingProc;
            for (int i=0; i<systemVar::sendRecvOrder_length; i++)
            {
                sendingProc=systemVar::sendRecvOrder[i][0];
                receivingProc=systemVar::sendRecvOrder[i][1];
                for (int iorder = 0; iorder <= mathVar::orderElem; ++iorder)
                {
                    auxUlti::functionsOfParallelComputing::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhoX,parallelBuffer::drhoX,iorder);
                    auxUlti::functionsOfParallelComputing::sendRecvDiscretedVar(sendingProc, receivingProc, BR1Vars::rhoY,parallelBuffer::drhoY,iorder);
                }
            }
        }
        else if (systemVar::auxVariables==2)
        {
            //Chua lam cho method BR2
        }
    }

    /*void sendReceiveMeshData(int vertex, int dir, double**Buffer)
    {
        //dir = 1: xCoor
        //dir = 2: yCoor

        //MPI_Request request_send, request_recv;
        MPI_Request *request_send = new MPI_Request[meshVar::numBCEdges*4], //Phai nhan 4 vi phan tu nhieu node nhat la quad
                *request_recv = new MPI_Request[meshVar::numBCEdges*4];

        MPI_Barrier(MPI_COMM_WORLD);
        int coef(1);
        int destination, source, cellId, neighborCellId, tag_sent, tag_recv, sentCount(0), recvCount(0);

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

                //send toa do dinh cua cell
                if (vertex<elemType)
                {
                    ptId=meshVar::Elements2D[cellId][vertex];
                    MPI_Isend(&meshVar::Points[ptId][dir-1], 1, MPI_DOUBLE, destination, tag_sent*100+10+vertex, MPI_COMM_WORLD, &request_send[sentCount]);
                    sentCount++;
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        //Chi receive khi nao toan bo cac processor send data xong

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

                //elemType nay la kieu phan tu cua cell neighbor chu k phai cell dang xet (cellId) ---> BUG!!!!!
                int elemType(auxUlti::checkType(cellId));

                //receive toa do dinh cua cell
                if (vertex<elemType)
                {
                    MPI_Irecv(&Buffer[iBCedge][vertex], 1, MPI_DOUBLE, source, tag_recv*100+10+vertex, MPI_COMM_WORLD, &request_recv[recvCount]);
                    recvCount++;
                }
            }
        }

        MPI_Waitall(sentCount,request_send,MPI_STATUSES_IGNORE);
        MPI_Waitall(recvCount,request_recv,MPI_STATUSES_IGNORE);

        MPI_Barrier(MPI_COMM_WORLD);
        delete [] request_send;
        delete [] request_recv;
    }*/

    void sendReceiveMeshData(int vertex, int dir, double**Buffer)
    {
        //dir = 1: xCoor
        //dir = 2: yCoor

        int sendingProc, receivingProc;
        for (int i=0; i<systemVar::sendRecvOrder_length; i++)
        {
            sendingProc=systemVar::sendRecvOrder[i][0];
            receivingProc=systemVar::sendRecvOrder[i][1];

            MPI_Barrier(MPI_COMM_WORLD);
            //Declare dynamic array
            MPI_Request *request = new MPI_Request[meshVar::numBCEdges];

            int destination, source, cellId, neighborCellId, tag_sent, tag_recv, requestCount(0);
            if (systemVar::currentProc==sendingProc) //current proc send, neighbor receive
            {
                //Send
                for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge)
                {
                    destination=meshVar::meshConnection[iBCedge][1];
                    cellId=meshVar::meshConnection[iBCedge][0];
                    neighborCellId=meshVar::meshConnection[iBCedge][2];

                    int elemType(auxUlti::checkType(cellId)), ptId(0);

                    if (destination==receivingProc)
                    {
                        tag_sent=neighborCellId;
                        //send toa do dinh cua cell
                        if (vertex<elemType)
                        {
                            ptId=meshVar::Elements2D[cellId][vertex];
                            MPI_Isend(&meshVar::Points[ptId][dir-1], 1, MPI_DOUBLE, destination, tag_sent*100+10+vertex, MPI_COMM_WORLD, &request[requestCount]);
                            requestCount++;
                        }
                    }
                }
                //Wait for all request fullfilled
                MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
            }
            else if (systemVar::currentProc==receivingProc)
            {
                //Receive
                for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; ++iBCedge) {
                    source=meshVar::meshConnection[iBCedge][1];
                    cellId=meshVar::meshConnection[iBCedge][0];
                    neighborCellId=meshVar::meshConnection[iBCedge][2];

                    //tag = neighborCellId*10 + destination
                    if (source==sendingProc)
                    {
                        int elemType(auxUlti::checkType(cellId));

                        tag_recv=cellId;
                        //receive toa do dinh cua cell
                        if (vertex<elemType)
                        {
                            MPI_Irecv(&Buffer[iBCedge][vertex], 1, MPI_DOUBLE, source, tag_recv*100+10+vertex, MPI_COMM_WORLD, &request[requestCount]);
                            requestCount++;
                        }
                    }
                }
                //Wait for all request fullfilled
                MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
            }

            MPI_Barrier(MPI_COMM_WORLD);
            delete [] request;
        }
    }

    void sendString(std::string content, int destination, int tag)
    {
        MPI_Send(content.c_str(), content.size(), MPI_CHAR, destination, tag, MPI_COMM_WORLD);
    }

    std::string receiveString(int source, int tag)
    {
        //From internet (giu nguyen khong sua)
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
        //Check subsonic
        refValues::subsonic = auxUlti::checkSubSonic();
        if (systemVar::currentProc==0)
        {
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

    bool checkTimeVaryingBCAvailable()
    {
        if (bcValues::slipBCFlag || bcValues::temperatureJump)
        {
            return true;
        }
        else
        {
            return false;
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

    void shrink2DIntVector(std::vector<std::vector<int>>&vector, int numRow)
    {
        int maxNumRow(static_cast<int>(vector.size()));
        /*Use erase function to clear vector*/
        //numRow+2 de chac an la khong xoa nham phan chua data
        for (int row = numRow+2; row < maxNumRow; row++) {
            vector[row].erase(vector[row].begin(), vector[row].end());
            vector[row].shrink_to_fit();
        }
        vector.erase(vector.begin()+numRow+2, vector.end());
        vector.shrink_to_fit();
    }

    void releaseMemory()
    {
        /*
        Ham release bo nho sau khi doc va xu luoi
        */
        //meshVar::esup1.clear();
        //meshVar::esup2.clear();
        //meshVar::psup1.clear();
        //meshVar::psup2.clear();

        //shrink vector to fit its size
        //auxUlti::shrink2DIntVector(meshVar::inpoel,meshVar::nelem2D);
        //auxUlti::shrink2DIntVector(meshVar::esuel,meshVar::nelem2D);
        //auxUlti::shrink2DIntVector(meshVar::inpoed,meshVar::inpoedCount);
        //auxUlti::shrink2DIntVector(meshVar::inedel,meshVar::nelem2D);
        //auxUlti::shrink2DIntVector(meshVar::ineled,meshVar::inpoedCount);
    }

    void resizeTemporaryArrays()
    {
        //Resize array
        meshVar::esup1.resize(6*meshVar::npoin);
        meshVar::esup2.resize(meshVar::npoin+1);

        meshVar::inpoel = auxUlti::resize2DIntArray(meshVar::nelem2D,5,0);

        meshVar::psup1.resize(meshVar::npoin*6);
        meshVar::psup2.resize(meshVar::npoin+1);

        meshVar::esuel = auxUlti::resize2DIntArray(meshVar::nelem2D,4,0);

        meshVar::inpoed = auxUlti::resize2DIntArray(meshVar::nelem2D*4,5,0);

        meshVar::inedel = auxUlti::resize2DIntArray(meshVar::nelem2D,4,0);

        meshVar::ineled = auxUlti::resize2DIntArray(meshVar::nelem2D*4,5,0);

        SurfaceBCFields::localGlobalBCEdgesMatching = new int [meshVar::nelem1D];
    }
}
