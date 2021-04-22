#ifndef BCSUPPORTFNCS_H
#define BCSUPPORTFNCS_H
#include <vector>
namespace BCSupportFncs
{
    /*Function returns true if at considering face (edge), flow is going into computational domain, returns false if flow is going out of computational domain.
    - nx, ny are surface unit normal vector components.*/
    bool checkInflow(double u, double v, double nx, double ny);

    void decompseU(std::vector<double> &priVars, const std::vector<double> &U, double rhoX, double rhoY, double iniT, bool calcT);

    void decompsedU(std::vector<double> &priVars, const std::vector<double> &U, const std::vector<double> &dU, double T);

    void correctPriVars(int edge, int edgeGrp, std::vector<double> &priVarsM, const std::vector<double> &priVarsP, const std::vector<double> &priVarsMean, const std::vector<double> &n, bool inflow);

    void correctPriVarsGrad(int edgeGrp, std::vector<double> &dpriVarsXM, std::vector<double> &dpriVarsYM, const std::vector<double> &dpriVarsXP, const std::vector<double> &dpriVarsYP, const std::vector<double> &UM, double TM, const std::vector<double> &n, bool inflow);

    void correctDensity(std::vector<double> &priVarsM, int edgeGrp, const std::vector<double> &priVarsP);

    std::tuple <double, double> correctVelocity(int edge, int edgeGrp, double uP, double uMean, double vP, double vMean, const std::vector<double> &n, bool inflow);

    std::tuple <double, double> correctVelocityGrad(int edgeGrp, double duXP, double duYP, bool inflow, const std::vector<double> &n);

    double correctTemperature(int edge, int edgeGrp, double TP, double TMean, bool inflow);

    std::tuple<double, double> correctTemperatureGrad(int edgeGrp, double dTXP, double dTYP, bool inflow, const std::vector<double> &n);

    double correctPressure(int edge, int edgeGrp, double pP, double pMean, bool inflow);

    std::tuple<double, double> correctPressureGrad(int edgeGrp, double dpXP, double dpYP, bool inflow, const std::vector<double> &n);

    void reconstructConvectiveU(std::vector<double> &U, const std::vector<double> &priVars);

    void reconstructdU(std::vector<double> &dU, const std::vector<double> &priVars, const std::vector<double> &U);

    namespace auxilaryBCs {
        void getUPlus(int edge, int nG, std::vector<double> &UPlus);

        void calcUMean(int element, std::vector<double> &UMean);

        /*Ham apply dieu kien back flow tai vi tri bi back (reversed) flow*/
        void calcUReversedFlow(int edge, int nG, std::vector<double> &U);
    }

    namespace NSFEqBCs {
        void findMaxLxFConstantOnEdge(int element, int edge);

        void getdUPlus(int edge, int nG, std::vector<double> &dUXPlus, std::vector<double> &dUYPlus);

        void calcdUMean(int edge, int element, std::vector<double> &dUXMean, std::vector<double> &dUYMean);

        double calcTotalEnergy(const std::vector<double> &UMinus, double T, double drhoX, double drhoY);

        /*Function calculates fluxes of all advective and diffusive terms of NSF equation from conservative variables*/
        void NSFEqFluxes(int edgeId, std::vector<double> &Fluxes, double TPlus, double TMinus, std::vector<double> &UPlus, std::vector<double> &UMinus, std::vector<double> &dUXPlus, std::vector<double> &dUXMinus, std::vector<double> &dUYPlus, std::vector<double> &dUYMinus, std::vector<double> &normVector);
    }

    namespace parallel {
    void getUMinus_total(int edge, int nG, std::vector<double> &UMinus);

    void getUMinus_convective(int edge, int nG, std::vector<double> &UMinus);

    void getdUMinus(int edge, int nG, std::vector<double> &dUMinus, int dir);
    }
}
#endif // BCSUPPORTFNCS_H
