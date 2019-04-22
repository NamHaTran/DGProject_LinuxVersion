#include <iostream>
#include <QString>
#include <QTextStream>
#include "Unv2DReaderLib.h"
#include "VarDeclaration.h"
using namespace std;

int main()
{
    //DISPLAY LOGO--------------------------------------------------------------------------
    dispLogo();

    //INPUT FILE NAME-----------------------------------------------------------------------
    getCase();

    cout << "Enter file name: ";
    QTextStream qtin(stdin);
    qtin >> meshFileName;

    QString fullName = pwd + "/" + meshFileName;
    cout << " " << endl
        << "**Unv2DReader is running**" << endl;

    //GET MESH INFORMATION------------------------------------------------------------------
    GetMesh(fullName);

    //GET NODES INFORMATION-----------------------------------------------------------------
    GetNodes();

    //GET ELEMENTS INFORMATION--------------------------------------------------------------
    GetElements();

    //GET BOUNDARIES INFORMATION------------------------------------------------------------
    GetBoundaries();

    //MARK BOUNDARY EDGE--------------------------------------------------------------------
    MarkBoundary();

    //SAVE DATA-----------------------------------------------------------------------------
    SaveMeshData();

    QString cmd("y");
    cout << "Do you want to create template setting files? <y/n>: ";
    qtin >> cmd;
    cout << "\n";
        if (cmd.compare("y") == 0)
        {
            //CREATE FOLDER 0-----------------------------------------------------------------------
            createFolder0();

            //CREATE TEMPLATE FILES-----------------------------------------------------------------
            createTemplate();
        }

    //system("pause");
    return 0;
}
