#ifndef GENERALPARALLELFUNCS_H
#define GENERALPARALLELFUNCS_H
#include <vector>
#include <string>
namespace parallelFuncs_Gen
{
    //std::tuple<double,double> getGaussPointCoorsOfNeighborCell(int loc, int nG);
    void prepareParallelCase();
    void sendString(std::string content, int destination, int tag);
    std::string receiveString(int source, int tag);
    void sendReceiveMeshData(int vertex, int dir, double**Buffer);
    int getMatchedEdgeIdOfCell(int element);

    void resizeMeshParallelBuffers();

    //void resizeGaussPtParallelBuffers();
}
#endif // GENERALPARALLELFUNCS_H
