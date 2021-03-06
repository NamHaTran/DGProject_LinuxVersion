#ifndef FRANZDURSTMODEL_H
#define FRANZDURSTMODEL_H
#include <tuple>
#include <vector>

namespace extNSF_Durst {
    /*Variables*/
    extern bool enable;


    void correctViscousTerms(std::vector<std::vector<double>> &diffTerms, std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy);

    std::tuple<double, double, double, double> calcSelfDiffusionTerms(std::vector< std::vector<double> > &selfDiffusionTensor, double uVal, double vVal, int dir);

    std::vector<std::vector<double>> calcSelfDiffusionTensor(std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy);

    double calcDiffusionStressComponent(int index, double fstTerm, double sndTerm);
}

#endif // FRANZDURSTMODEL_H
