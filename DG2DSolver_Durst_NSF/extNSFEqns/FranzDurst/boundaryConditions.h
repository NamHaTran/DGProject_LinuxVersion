#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H
#include <vector>

namespace bcForExtNSF_Durst {
extern double
//! Variable for correcting boundary condition in case of Durst model. tempdTx is temporary dTx.
tempdTx,
//! Variable for correcting boundary condition in case of Durst model. tempdTy is temporary dTy.
tempdTy;

void dropNormDiffVel(int edge, int edgeGrp, double &dRhoXM, double &dRhoYM, double dRhoXP, double dRhoYP, double rhoBC, double dTBCx, double dTBCy, const std::vector<double> &n);

void checkConditionToAddDiffTerms(int edge);

void resetNeedToRemoveDiffTermFlag();
}
#endif // BOUNDARYCONDITIONS_H
