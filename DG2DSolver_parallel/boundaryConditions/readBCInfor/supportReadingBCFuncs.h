#ifndef SUPPORTREADINGBCFUNCS_H
#define SUPPORTREADINGBCFUNCS_H
#include <sstream>
#include <vector>

void readUAppMeth(std::ifstream &FileFlux, int bcGrp, std::vector<bool> &UMethod);
void readGradUAppMeth(std::ifstream &FileFlux, int bcGrp, std::vector<bool> &gradUMethod);

#endif // SUPPORTREADINGBCFUNCS_H
