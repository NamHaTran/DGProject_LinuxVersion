#include "Unv2DReaderLib.h"
#include "DGMessagesLib.h"
#include "DGAuxUltilitiesLib.h"
#include "varDeclaration.h"
#include <QTextStream>
#include <iostream>
#include <QString>
#include <vector>
#include <QFile>

extern QString caseName, wD, pwd;

void dispLogo()
{
	std::cout << "+------------------------------ CFD LOVES C++ ---------------------------+" << std::endl
		<< "|-------------------UNV MESH FORMAT READER (2D version)------------------|" << std::endl
		<< "|                                 Author(s)                              |" << std::endl
		<< "|   Nam Ha Tran.                                                         |" << std::endl
		<< "|   Ver 1.00                                                             |" << std::endl
		<< "+------------------------------------------------------------------------+" << std::endl
		<< "|   This tool is developed as a part of project 'Discontinous Galerkin   |" << std::endl
		<< "|   Method for 2D problems'.                                             |" << std::endl
		<< "+------------------------------------------------------------------------+" << std::endl;
}

void getCase()
{
	std::cout << "***Getting case's information***\n";
	wD = auxUlti::workingdir();

	/*Get caseName from SubmitCase*/
    QString submitLoc(wD);  //submitLoc contents location of submitingCase
	//submitLoc.erase(submitLoc.end() - 22, submitLoc.end());
    submitLoc += "/CASES/SubmitCase.txt";

    QFile submitCaseFlux(submitLoc);
    if (submitCaseFlux.open(QIODevice::ReadOnly))
	{
        QTextStream line(&submitCaseFlux);
        QString keyWord;
        int keyWFlag(0);
        while (!line.atEnd())  //Read submitCaseFlux line by line
		{
            QString line2str = line.readLine();
            QStringList ptr = line2str.split(" ", QString::SkipEmptyParts);  //Split line into list of strings by spliter "space"

			int numWrd = ptr.size();
			if (numWrd == 2)
			{
                QString str1(ptr[0]), str2("caseName");
                if (!str1.compare(str2))
				{
					caseName = ptr[1];
                    keyWFlag = 1;
                    std::cout << "	Case " << caseName.toStdString() << " has been submitted\n";
				}
			}
		}
		if (keyWFlag == 0)
		{
            std::cout << message::undfKeyW("caseName", submitLoc).toStdString() << std::endl;
		}
	}
	else
	{
        message::writeLog((wD + "/CASES/"), "", message::opFError("SubmitCase.txt", submitLoc));
	}
    pwd = wD + "/CASES/" + caseName;
}

void GetMesh(QString meshName)
{
	int iline(0);

    QFile meshFlux(meshName);
    if (meshFlux.open(QIODevice::ReadOnly))  //open file to read
	{
        QTextStream line(&meshFlux);
        QString keyWord;
        while (!line.atEnd())  //Read meshFlux line by line
        {
            Mesh[iline] = line.readLine();  //Read line and convert it into Qstring type and pass it to Mesh array
            QStringList ptr = Mesh[iline].split(" ", QString::SkipEmptyParts);  //Split line into list of strings by spliter "space"
            iline++;

            int sizePtr(static_cast<int>(ptr.size()));
            if (sizePtr == 1)
            {
                int elemType = ptr[0].toInt();
                if (elemType == 2411)  //Nodes
                {
                    location[0][0] = iline;
                    location[0][1] = elemType;
                }
                else if (elemType == 2412)  //Elements
                {
                    location[1][0] = iline;
                    location[1][1] = elemType;
                }
                else if (elemType == 2467)  //Boundaries
                {
                    location[2][0] = iline;
                    location[2][1] = elemType;
                }
            }
            noLines += 1;
        }
	}
	else
	{
        QString logFile(wD + "/CASES/" + caseName);
		message::writeLog(pwd, caseName, message::opFError(meshName, logFile));
	}
}

void GetNodes()
{
	int firstNodeLine(location[0][0]), lastNodeLine(location[1][0] - 3);
	nodeNumber = ((lastNodeLine - firstNodeLine) / 2);
	std::cout << "**Reading nodes coordinates**" << std::endl
		<< "-----Mesh has " << nodeNumber << " nodes." << std::endl;
	int index(0);
	for (int nline = firstNodeLine + 1; nline < lastNodeLine + 1; nline += 2)
	{
		Points[index][0] = index + 1;
        QStringList nodeDataStr = Mesh[nline].split(" ", QString::SkipEmptyParts);
        for (int i = 1; i < 4; ++i) {
            Points[index][i] = nodeDataStr[i - 1].toDouble();
        }
		index++;
	}
	std::cout << "DONE!" << std::endl << " " << std::endl;
}

void GetElements()
{
	std::cout << "**Reading information of elements**" << std::endl;

    int firstElemLine(location[1][0]), lastElemLine(location[2][0] - 3);
	int elemLocation[4]; //From row 1 to 4: starting and ending lines of 1D and 2D elements
	elemLocation[0] = firstElemLine + 1;
	for (int nline = firstElemLine; nline < lastElemLine; nline += 3)
	{
        QStringList elemProperties(Mesh[nline].split(" ", QString::SkipEmptyParts));
        if (elemProperties[1].toInt() != 11)
		{
			elemLocation[1] = nline;
			break;
		}
	}
	elemLocation[2] = elemLocation[1] + 1;
	elemLocation[3] = lastElemLine;
	n1D = (elemLocation[1] - elemLocation[0] + 1) / 3;
	n2D = (elemLocation[3] - elemLocation[2] + 1) / 2;
	std::cout << "-----Mesh has " << n1D << " 1D elements." << std::endl
		<< "-----Mesh has " << n2D << " 2D elements." << std::endl;

	//Set default value for Elements2D array
	for (int nRow = 0; nRow < n2D; nRow++)
	{
		Elements2D[nRow][4] = -22;
	}

	//1D ELEMENTS
	int index(0);
	for (int nline = elemLocation[0] - 1; nline < elemLocation[1] - 1; nline += 3)
	{
        QStringList elemProperties(Mesh[nline].split(" ", QString::SkipEmptyParts));
        Elements1D[index][0] = elemProperties[0].toInt();

        elemProperties = Mesh[nline + 2].split(" ", QString::SkipEmptyParts);
        Elements1D[index][1] = elemProperties[0].toInt();
        Elements1D[index][2] = elemProperties[1].toInt();
		index++;
	}

	//2D ELEMENTS
	index = 0;
	for (int nline = elemLocation[2] - 1; nline < elemLocation[3] - 1; nline += 2)
	{
		int typeElem(0);
        QStringList elemProperties (Mesh[nline].split(" ", QString::SkipEmptyParts));
        Elements2D[index][0] = elemProperties[0].toInt();
        typeElem = elemProperties[1].toInt();

        elemProperties = Mesh[nline + 1].split(" ", QString::SkipEmptyParts);
		if (typeElem == 44)  //quad element
		{
            for (int i = 0; i < 4; ++i) {
                Elements2D[index][i+1]=elemProperties[i].toInt();
            }
		}
		else if (typeElem == 41)  //tri element
		{
            for (int i = 0; i < 3; ++i) {
                Elements2D[index][i+1]=elemProperties[i].toInt();
            }
		}
		index++;
	}
	std::cout << "DONE!" << std::endl << " " << std::endl;
}

void GetBoundaries()
{
	std::cout << "**Reading information of boundaries**" << std::endl;
	int firstBoundLine(location[2][0] + 1), lastBoundLine(noLines - 2);

    QString BCName;
	for (int nline = firstBoundLine; nline < lastBoundLine; nline++)
	{
        QStringList lineProperties(Mesh[nline].split(" ", QString::SkipEmptyParts));
        if (!(lineProperties[0].toInt()))  //Check input is string or number
		{
			/*NOTE: when line lineProperties>>checkNumber is excuted, a pointer of stream lineProperties
			is changed (it's not located at the first data of stream), if we use stream lineProperties again,
			(ex: lineProperties>>BCname), lineProperties will transfer data at second location to BCname*/

            boundName.push_back(Mesh[nline]);
			boundLocation.push_back(nline);
			numOfBound += 1;
		}
	}
	numOfBound = static_cast<int>(boundLocation.size());
	if (numOfBound != 0)
	{
		savingFlag = true;
		std::cout << "-----Mesh has " << numOfBound << " boundaries:" << std::endl;
		for (int i = 0; i < numOfBound; i++)
		{
            std::cout << "     - " << boundName[i].toUtf8().constData() << std::endl;
		}

		//Set all values of boundaries to -22
		for (int nRow = 0; nRow < n1D; nRow++)
		{
			for (int nColumn = 0; nColumn < 2; nColumn++)
			{
				boundaries[nRow][nColumn] = -22;
			}
		}

		int index(0), group(0);
		for (int nline = firstBoundLine; nline < lastBoundLine + 1; nline++)
		{
			std::vector<int> linePropertiesNum;
            QStringList lineProperties(Mesh[nline].split(" ", QString::SkipEmptyParts));

            for (int i = 0; i < static_cast<int>(lineProperties.size()); ++i) {
                linePropertiesNum.push_back(lineProperties[i].toInt());
            }

            if (nline<lastBoundLine && (static_cast<int>(linePropertiesNum.size())== 0))
			{
                QStringList linePropertiesStrNext(Mesh[nline + 1].split(" ", QString::SkipEmptyParts));
				std::vector<int> linePropertiesNextNum;

                for (int i = 0; i < static_cast<int>(linePropertiesStrNext.size()); ++i) {
                    linePropertiesNextNum.push_back(linePropertiesStrNext[i].toInt());
                }

				if (static_cast<int>(linePropertiesNextNum.size()) != 0 && linePropertiesNextNum[0] == 8)
				{
					group += 1;
				}
			}

			if (static_cast<int>(linePropertiesNum.size()) == 8 && linePropertiesNum[0] == 8 && linePropertiesNum[4] == 8)
			{
				//Line contents boundary edges
				boundaries[index][0] = linePropertiesNum[1];
				boundaries[index + 1][0] = linePropertiesNum[5];
				boundaries[index][1] = group;
				boundaries[index + 1][1] = group;
				index += 2;
				numOfBoundEdge += 2;
			}
			else if (static_cast<int>(linePropertiesNum.size()) == 4 && linePropertiesNum[0] == 8)
			{
				boundaries[index][0] = linePropertiesNum[1];
				boundaries[index][1] = group;
				index += 1;
				numOfBoundEdge += 1;
			}
		}
		std::cout << "DONE!" << std::endl << " " << std::endl;
	}
	else
	{
		savingFlag = false;
        QString errorStr(" ");
		errorStr = R"(-----WARNING: Mesh has no boundary!!!
     Boundaries must be created before simulation.
     Reading mesh has finished with no boundary information.
)";
		message::writeLog(pwd, caseName, errorStr);
	}
}

void MarkBoundary()
{
	for (int i = 0; i < numOfBoundEdge; i++)
	{
		for (int j = 0; j < n1D; j++)
		{
			if (Elements1D[j][3] == 0)
			{
				if (Elements1D[j][0] == boundaries[i][0])
				{
					Elements1D[j][3] = boundaries[i][1];
				}
			}
		}
	}
}

void SaveMeshData()
{
	std::cout << "**Saving data**\n";
    extern std::vector<QString> boundName;

	//Declare saving locations
    QString bcPatchPath(wD + "/CASES/" + caseName + "/Constant/boundaryPatch.txt");
    QString PointsPatchPath(wD + "/CASES/" + caseName + "/Constant/Mesh/Points.txt");
    QString Elem1DPatchPath(wD + "/CASES/" + caseName + "/Constant/Mesh/Elements1D.txt");
    QString Elem2DPatchPath(wD + "/CASES/" + caseName + "/Constant/Mesh/Elements2D.txt");

    QFile FluxBC(bcPatchPath);
    if (FluxBC.open(QIODevice::ReadWrite | QIODevice::Truncate | QIODevice::Text))
	{
        QString header(message::headerFile());
        QTextStream FluxBCstream(&FluxBC);
        FluxBCstream << header
            << "  \n"
			<< "Boundary condition definitions:\n";

		for (int i = 0; i < numOfBound; i++)
		{
            FluxBCstream << boundName[i] << endl
				<< "{\n"
                << "	Group				" << i + 1 << endl
				<< "	Type				wall\n"
				<< "	Method				weakRiemann\n"
				<< "}\n" << " \n";
		}

	}
	else
	{
		message::writeLog(pwd, caseName, message::opFError("boundaryPatch.txt", bcPatchPath));
	}

	//Save points data
    QFile FluxPt(PointsPatchPath);
    if (FluxPt.open(QIODevice::ReadWrite | QIODevice::Truncate | QIODevice::Text))
	{
        QTextStream FluxPtstream(&FluxPt);
		for (int i = 0; i < nodeNumber; i++)
		{
			for (int j = 0; j < 4; j++)
			{
                FluxPtstream << Points[i][j] << " ";
			}
            FluxPtstream << endl;
		}

	}
	else
	{
		message::writeLog(pwd, caseName, message::opFError("Points.txt", PointsPatchPath));
	}

	//Save elements 1D data
    QFile FluxElem1D(Elem1DPatchPath);
    if (FluxElem1D.open(QIODevice::ReadWrite | QIODevice::Truncate | QIODevice::Text))
	{
        QTextStream FluxElem1Dstream(&FluxElem1D);
		for (int i = 0; i < n1D; i++)
		{
			for (int j = 0; j < 4; j++)
			{
                FluxElem1Dstream << Elements1D[i][j] << " ";
			}
            FluxElem1Dstream << endl;
		}

	}
	else
	{
		message::writeLog(pwd, caseName, message::opFError("Elements1D.txt", Elem1DPatchPath));
	}

	//Save element 2D data
    QFile FluxElem2D(Elem2DPatchPath);
    if (FluxElem2D.open(QIODevice::ReadWrite | QIODevice::Truncate | QIODevice::Text))
	{
        QTextStream FluxElem2Dstream(&FluxElem2D);
		for (int i = 0; i < n2D; i++)
		{
			for (int j = 0; j < 5; j++)
			{
                FluxElem2Dstream << Elements2D[i][j] << " ";
			}
            FluxElem2Dstream << endl;
		}

	}
	else
	{
		message::writeLog(pwd, caseName, message::opFError("Elements2D.txt", Elem2DPatchPath));
	}

	std::cout << "DONE!" << std::endl << " " << std::endl;
}

void createFolder0()
{
    extern std::vector<QString> boundName;

	int numOfBC(static_cast<int>(boundName.size()));
    QString TLoc(wD + "/CASES/" + caseName + "/0/T.txt");
    QString ULoc(wD + "/CASES/" + caseName + "/0/U.txt");
    QString PLoc(wD + "/CASES/" + caseName + "/0/p.txt");
    QString header(message::headerFile());
	/*create T.txt*/
    QFile TFlux(TLoc);
    if (TFlux.open(QIODevice::ReadWrite | QIODevice::Truncate | QIODevice::Text))
	{
        QTextStream TFluxstream(&TFlux);
        TFluxstream << header << message::headerpTU("T");

		for (int i = 0; i < numOfBC; i++)
		{
            TFluxstream << boundName[i] << endl
				<< "{\n"
                << "	Group				" << i + 1 << endl
				<< "	Type				WallAdiabatic\n" << "}\n" << " \n";
		}

	}
	else
	{
        QString logFile(wD + "/CASES/" + caseName);
		message::writeLog(pwd, caseName, message::opFError("T.txt", logFile));
	}

	/*Create P.txt*/
    QFile PFlux(PLoc);
    if (PFlux.open(QIODevice::ReadWrite | QIODevice::Truncate | QIODevice::Text))
	{
        QTextStream PFluxstream(&PFlux);
        PFluxstream << header << message::headerpTU("p");

		for (int i = 0; i < numOfBC; i++)
		{
            PFluxstream << boundName[i] << endl
				<< "{\n"
                << "	Group " << i + 1 << endl
				<< "	Type				zeroGradient\n" << "}\n" << " \n";
		}

	}
	else
	{
        QString logFile(wD + "/CASES/" + caseName);
		message::writeLog(pwd, caseName, message::opFError("p.txt", logFile));
	}

	/*Create U.txt*/
    QFile UFlux(ULoc);
    if (UFlux.open(QIODevice::ReadWrite | QIODevice::Truncate | QIODevice::Text))
	{
        QTextStream UFluxstream(&UFlux);
        UFluxstream << header << message::headerpTU("U");

		for (int i = 0; i < numOfBC; i++)
		{
            UFluxstream << boundName[i] << endl
				<< "{\n"
                << "	Group " << i + 1 << endl
				<< "	Type				noSlip\n" << "}\n" << " \n";
		}

	}
	else
	{
        QString logFile(wD + "/CASES/" + caseName);
		message::writeLog(pwd, caseName, message::opFError("U.txt", logFile));
	}
}

void createTemplate()
{
	std::cout << "**Creating template files**\n";

    QString header(message::headerFile());
    QString DGOpLoc(wD + "/CASES/" + caseName + "/System/DGOptions.txt");
    QString MatLoc(wD + "/CASES/" + caseName + "/Constant/Material.txt");

	/*Create file DGOptions*/
    QFile DGOpFlux(DGOpLoc);
    if (DGOpFlux.open(QIODevice::ReadWrite | QIODevice::Truncate | QIODevice::Text))
	{
        QTextStream DGOpFluxstream(&DGOpFlux);
        DGOpFluxstream << header
			<< "	Note: all parameter are in SI unit.\n"
			<< "DGoptions\n"
			<< "{\n"
			<< "	numberOfGaussPoints			2\n"
			<< "	orderOfAccuracy				2\n"
			<< "	CourantNumber				0.5\n"
			<< "	totalTime(s)				10\n"
			<< "	writeInterval				100\n"
			<< "	writeLog					true\n"
			<< "	ddtScheme					Euler\n"
			<< "}";
	}
	else
	{
        QString logFile(wD + "/CASES/" + caseName);
		message::writeLog(pwd, caseName, message::opFError("DGOptions.txt", logFile));
	}

	/*Create file Material*/
    QFile MatFlux(MatLoc);
    if (MatFlux.open(QIODevice::ReadWrite | QIODevice::Truncate | QIODevice::Text))
	{
        QTextStream MatFluxstream(&MatFlux);
        MatFluxstream << header
			<< "	Note: all parameter are in SI unit.\n"
			<< "Material Properties\n"
			<< "{\n"
			<< "	gammaRatio					1.4\n"
			<< "	gasConstant					287\n"
			<< "	PrandtlNumber				0.72\n"
			<< "	SutherlandAs				1.46e-6\n"
			<< "	SutherlandTs				110.4\n"
			<< "}";
	}
	else
	{
        QString logFile(wD + "/CASES/" + caseName);
		message::writeLog(pwd, caseName, message::opFError("Material.txt", logFile));
	}

    /*Create file LimiterSettings*/
    QString limiterSettingLoc(wD + "/CASES/" + caseName + "/System/LimiterSettings.txt");
    QFile limiterSettingFlux(limiterSettingLoc);
    if (limiterSettingFlux.open(QIODevice::ReadWrite | QIODevice::Truncate | QIODevice::Text))
    {
        QTextStream limiterSettingFluxstream(&limiterSettingFlux);
        limiterSettingFluxstream << header
            << "LimiterSettings\n \n"
                << "limiter		PositivityPreserving	PAdaptive\n"
                << "PositivityPreserving\n"
                << "{\n"
                << "	version	simplified\n"
                << "}\n"
                << "PAdaptive\n"
                << "{\n"
                << "	/nothing here\n"
                << "}\n";
        }
        else
        {
            QString logFile(wD + "/CASES/" + caseName);
            message::writeLog(pwd, caseName, message::opFError("LimiterSettings.txt", logFile));
        }

	std::cout << "DONE!\n";
}

void clearVar()
{
	int numOfBC(static_cast<int>(boundName.size()));
	int numOfBCLoc(static_cast<int>(boundLocation.size()));

	/*Use erase function to clear vector*/
	boundName.erase(boundName.begin(), boundName.begin() + numOfBC);
	boundLocation.erase(boundLocation.begin(), boundLocation.begin() + numOfBCLoc);
	/*Shrink to fit
	Requests the container to reduce its capacity to fit its size.*/
	boundName.shrink_to_fit();
	boundLocation.shrink_to_fit();

	/*Reset all variable*/
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			location[i][j] = 0;
		}
	}

	for (size_t i = 0; i < pointsArrSize; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			Elements1D[i][j] = 0;
		}
	}

	for (size_t i = 0; i < elements2DArrSize; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			Elements2D[i][j] = 0;
		}
	}

	for (size_t i = 0; i < pointsArrSize; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			boundaries[i][j] = 0;
		}
	}

	numOfBound = 0;
	noLines = 0;
	savingFlag = false;  //Flag of saving data
	nodeNumber = 0;
	n1D = 0;
	n2D = 0;
	numOfBoundEdge = 1;
}
