#ifndef DURSTMODEL_H
#define DURSTMODEL_H
#include <tuple>
#include <vector>

namespace extNSF_Durst {
    /*Variables*/
    extern bool
    enable,
    diffusionAtWall,
    needToRemoveDiffTerm;
    extern double Dm;


    void correctViscousTerms(std::vector<std::vector<double>> &diffTerms, std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy);

    std::tuple<double, double, double, double> calcSelfDiffusionTerms(std::vector< std::vector<double> > &selfDiffusionTensor, double uVal, double vVal, int dir);

    std::vector<std::vector<double>> calcSelfDiffusionTensor(std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy);

    double calcDiffusionStressComponent(int index, double fstTerm, double sndTerm);

    void correctEnergyEqnVolIntTerm(int element, std::vector<double> &VolIntTerm4);
}

#endif // DURSTMODEL_H
