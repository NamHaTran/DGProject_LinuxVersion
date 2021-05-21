#ifndef NONEQUILIBRIUMBCSLIB_H
#define NONEQUILIBRIUMBCSLIB_H

/*Function calculates and updates values of surfaceBCFields*/
void updateTimeVaryingBCs();

namespace nonequilibriumBCs
{
    namespace MaxwellSmoluchowskiBC
    {
        void explicitMethod_FVMType(int edge, int edgeGrp);

        void implicit2ndOrderMethod(int edge, int edgeGrp);

        void explicitMethod_DGType(int edge, int edgeGrp);
    }
}
#endif // NONEQUILIBRIUMBCSLIB_H
