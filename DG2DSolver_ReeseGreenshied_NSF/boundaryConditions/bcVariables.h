#ifndef BCVARIABLES_H
#define BCVARIABLES_H
#include <vector>
namespace BCVars {

    extern std::vector<bool> DirichletAppMethUStrong, NewmannAppMethGradUStrong,
    DirichletAppMethTStrong, NewmannAppMethGradTStrong,
    DirichletAppMethPStrong, NewmannAppMethGradPStrong,
    DirichletAppMethGeneralBCStrong, NewmannAppMethGradGeneralBCStrong;

    extern bool correctRho;

    namespace generalBCId {
        static constexpr int symmetry = 7;
        static constexpr int matched = 10;
    }

    namespace pressureBCId {
        static constexpr int interpFrmDensity = 11;
        static constexpr int zeroGrad = 2;
        static constexpr int fixedValue = 1;
        static constexpr int inletOutlet = 4;
    }

    namespace temperatureBCId {
        static constexpr int temperatureJump = 11;
        static constexpr int zeroGrad = 2;
        static constexpr int fixedValue = 1;
        static constexpr int inletOutlet = 4;
    }

    namespace velocityBCId {
        static constexpr int fixedValue = 1;
        static constexpr int noSlip = 2;
        static constexpr int movingWall = 3;
        static constexpr int slipWall = 5;
        static constexpr int inletOutlet = 4;
        static constexpr int zeroGrad = 6;
    }
}
#endif // BCVARIABLES_H
