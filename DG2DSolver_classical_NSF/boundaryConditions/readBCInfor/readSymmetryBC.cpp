#include "readSymmetryBC.h"
#include "VarDeclaration.h"
#include <fstream>
#include <string>
#include <sstream>
#include "supportReadingBCFuncs.h"
#include "../fixedValue.h"
#include "../bcVariables.h"

void readSymmetryBC(std::ifstream &FileFlux, int bcGrp)
{
    readUAppMeth(FileFlux,bcGrp,BCVars::DirichletAppMethGeneralBCStrong);
    readGradUAppMeth(FileFlux,bcGrp,BCVars::NewmannAppMethGradGeneralBCStrong);
}
