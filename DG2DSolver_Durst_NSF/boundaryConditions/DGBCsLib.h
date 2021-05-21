#ifndef DGBCSLIB_H_INCLUDED
#define DGBCSLIB_H_INCLUDED
#include <vector>
/*NOTES:
- All boundary functions must return numerical fluxes (rho, rhou, rhov, rhoE fluxes) at surface Gauss points!
*/

/*Function applies advective and diffusive boundary conditions to edges in NSF equation, it returns values of fluxes at Gauss point*/
std::vector<double> NSFEqBCsImplement(int element, int edge, int nG);
/*Cac ham con cua NSFEqBCsImplement*/
void NSFEqBCsForNormalBoundary(std::vector<double> &Fluxes, int element, int edge, int edgeGrp, int nG);

/*Function applies auxilary boundary conditions to edges in auxilary equation, it returns values of auxiraly flux at Gauss point*/
std::vector<std::vector<double>> auxEqBCsImplement(int element, int edge, int nG);
/*Cac ham con cua auxEqBCsImplement*/
void auxEqBCsForNormalBoundary(std::vector<std::vector<double>> &Fluxes, int element, int edge, int edgeGrp, int nG);

//Implement bondary condition of Rho (use when massDiffusion is on)
std::tuple<double, double> rhoBCsImplement(int edge, int nG);

#endif // DGBCSLIB_H_INCLUDED
