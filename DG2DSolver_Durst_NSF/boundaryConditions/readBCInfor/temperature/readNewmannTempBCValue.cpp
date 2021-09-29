#include "readNewmannTempBCValue.h"
#include "VarDeclaration.h"
#include <fstream>
#include <string>
#include <sstream>
#include "../supportReadingBCFuncs.h"
#include "../../fixedValue.h"
#include "../../bcVariables.h"

void readZeroGradT(std::ifstream &FileFlux, int bcGrp)
{
    //readUAppMeth(FileFlux,bcGrp,BCVars::DirichletAppMethUStrong);
    readGradUAppMeth(FileFlux,bcGrp,BCVars::NewmannAppMethGradTStrong,"T");
}
