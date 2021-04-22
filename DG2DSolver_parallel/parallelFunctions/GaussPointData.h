#ifndef GAUSSPOINTDATA_H
#define GAUSSPOINTDATA_H
#include <vector>
#include <string>
namespace parallelFuncs_GaussPt
{
    void sendRecvMatchedBCGaussPtValues(int sendingProc, int receivingProc, double**Var,double**Buffer,int nG);

    void sendRecvMatchedBCGaussPtCoors(int sendingProc, int receivingProc, double**Var,double**Buffer,int nG);

    /**
     * @brief Function send/recv value at Gauss points on "matched" BC.
     */
    void sendReceiveSurfGaussArray(double**Var,double**Buffer);

    /**
     * @brief Function send/recv U values at Gauss points on "matched" BC.
     */
    void synchSurfGaussU();

    void sendRecvSurfaceBCGaussPtArray(double**Var,double**Buffer);

    void synchGaussPtCoors();
}
#endif // GAUSSPOINTDATA_H
