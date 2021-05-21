#ifndef READDIRICHLETVELBCVALUE_H
#define READDIRICHLETVELBCVALUE_H
#include <sstream>
#include <vector>
void readNoSlipU(std::ifstream &FileFlux, int bcGrp);
void readFixedValueU(std::ifstream &FileFlux, int bcGrp);
void readMaxwellSlipU(std::ifstream &FileFlux, int bcGrp);
#endif // READDIRICHLETVELBCVALUE_H
