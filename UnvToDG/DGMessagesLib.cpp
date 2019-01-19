#include "DGMessagesLib.h"
#include <iostream>
#include <QString>
#include <QTextStream>
#include <QFile>
#include <iomanip>

namespace message
{
    QString headerFile()
	{
        QString headerStr(" ");
		headerStr = R"(+------------------------------ CFD LOVES C++ ---------------------------+
|--------------------DISCONTINUOS GALERKIN METHOD SOLVER-----------------|
|                                  Author                                |
|   Nam Ha Tran.                                                         |
|   Ver 1.00                                                             |
+------------------------------------------------------------------------+
|   This program uses Discontinous Galerkin method to solve 2D problems  |
|   on structural and unstructural mesh.                                 |
+------------------------------------------------------------------------+
)";
        QString Timestr(getTime());
		headerStr += "	Program ran at (d-m-y_h-m-s) " + Timestr + "\n";
		return headerStr;
	}

    QString undfKeyW(QString keyW, QString location)
	{
        QString str("Cannot find key word <" + keyW + "> at " + location);
		return str;
	}

    QString opFError(QString fileName, QString location)
	{
        QString str("Cannot open file <" + fileName + "> located at " + location);
		return str;
	}

    void writeLog (QString location, QString caseName, QString str)
	{
        QString logFile(location + "/log_" + caseName + ".txt");

        //Do not write log file, waste of time
        std::cout << "ERROR: " << str.toUtf8().constData() << std::endl;
        std::cout << "DGSolver will exit after you hit return.\n";
        //system("pause");
		exit(EXIT_FAILURE);
	}

    QString getTime() //This function is referenced from internet
	{
		/*
		auto t = std::time(nullptr);
		auto tm = *std::localtime(&t);

		std::ostringstream oss;
		oss << std::put_time(&tm, "%d-%m-%Y_%H-%M-%S");
		auto str = oss.str();
		*/
        QString str("time");
		return str;
	}

    QString headerpTU(QString file)
	{
        QString headerStr(" ");
		if (file.compare("p") == 0)
		{
			headerStr = R"(Pressure conditions (Pa)
initialValue				0
)";
		}
		else if (file.compare("T") == 0)
		{
			headerStr = R"(Temperature conditions (K)
initialValue				0
)";
		}
		else if (file.compare("U") == 0)
		{
			headerStr = R"(Velocity conditions (m/s)
initialValue			0 0 0
)";
		}
		return headerStr;
	}
}
