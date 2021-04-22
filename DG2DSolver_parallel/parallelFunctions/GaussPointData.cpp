#include "GaussPointData.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <iostream>
#include <mpi.h>
#include "parallelVariables.h"
#include "generalParallelFuncs.h"

namespace parallelFuncs_GaussPt
{
    void sendRecvMatchedBCGaussPtCoors(int sendingProc, int receivingProc, double**Var,double**Buffer,int nG)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        //Declare dynamic array
        MPI_Request *request = new MPI_Request[meshVar::numBCEdges];

        int destination, source, cellId, neighborCellId, tag_sent, tag_recv, requestCount(0), edgeId;

        if (systemVar::currentProc==sendingProc) //current proc send, neighbor receive
        {
            //Send
            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; iBCedge++)
            {
                destination=meshVar::meshConnection[iBCedge][1];
                //cellId=meshVar::meshConnection[iBCedge][0];
                neighborCellId=meshVar::meshConnection[iBCedge][2];

                //tag = neighborCellId*10 + destination
                if (destination==receivingProc)
                {
                    edgeId = auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(iBCedge);

                    //define tags
                    tag_sent = neighborCellId;

                    MPI_Isend(&Var[edgeId][mathVar::nGauss + nG + 1], 1, MPI_DOUBLE, destination, tag_sent, MPI_COMM_WORLD, &request[requestCount]);
                    requestCount++;
                }
            }
            //Wait for all request fullfilled
            MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
        }
        else if (systemVar::currentProc==receivingProc)
        {
            //Receive
            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; iBCedge++) {
                source=meshVar::meshConnection[iBCedge][1];
                cellId=meshVar::meshConnection[iBCedge][0];
                //neighborCellId=meshVar::meshConnection[iBCedge][2];

                //tag = neighborCellId*10 + destination
                if (source==sendingProc)
                {
                    edgeId = auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(iBCedge);

                    //define tags
                    tag_recv = cellId;

                    MPI_Irecv(&Buffer[edgeId][nG], 1, MPI_DOUBLE, source, tag_recv, MPI_COMM_WORLD, &request[requestCount]);
                    requestCount++;
                }
            }
            //Wait for all request fullfilled
            MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        delete [] request;
    }

    void sendRecvMatchedBCGaussPtValues(int sendingProc, int receivingProc, double**Var,double**Buffer,int nG)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        //Declare dynamic array
        MPI_Request *request = new MPI_Request[meshVar::numBCEdges];

        int destination, source, cellId, neighborCellId, tag_sent, tag_recv, requestCount(0), edgeId;

        if (systemVar::currentProc==sendingProc) //current proc send, neighbor receive
        {
            //Send
            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; iBCedge++)
            {
                destination=meshVar::meshConnection[iBCedge][1];
                //cellId=meshVar::meshConnection[iBCedge][0];
                neighborCellId=meshVar::meshConnection[iBCedge][2];

                //tag = neighborCellId*10 + destination
                if (destination==receivingProc)
                {
                    edgeId = auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(iBCedge);

                    //define tags
                    tag_sent = neighborCellId;

                    MPI_Isend(&Var[edgeId][nG], 1, MPI_DOUBLE, destination, tag_sent, MPI_COMM_WORLD, &request[requestCount]);
                    requestCount++;
                }
            }
            //Wait for all request fullfilled
            MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
        }
        else if (systemVar::currentProc==receivingProc)
        {
            //Receive
            for (int iBCedge = 0; iBCedge < meshVar::numBCEdges; iBCedge++) {
                source=meshVar::meshConnection[iBCedge][1];
                cellId=meshVar::meshConnection[iBCedge][0];
                //neighborCellId=meshVar::meshConnection[iBCedge][2];

                //tag = neighborCellId*10 + destination
                if (source==sendingProc)
                {
                    edgeId = auxUlti::getGlobalEdgeIdFromLocalBCEdgeId(iBCedge);

                    //define tags
                    tag_recv = cellId;

                    MPI_Irecv(&Buffer[edgeId][mathVar::nGauss + nG + 1], 1, MPI_DOUBLE, source, tag_recv, MPI_COMM_WORLD, &request[requestCount]);
                    requestCount++;
                }
            }
            //Wait for all request fullfilled
            MPI_Waitall(requestCount,request,MPI_STATUSES_IGNORE);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        delete [] request;
    }

    /**
     * @brief Function send/recv value at Gauss points on "matched" BC.
     */
    void sendReceiveSurfGaussArray(double**Var,double**Buffer)
    {
        int sendingProc, receivingProc;
        for (int i=0; i<systemVar::sendRecvOrder_length; i++)
        {
            sendingProc=systemVar::sendRecvOrder[i][0];
            receivingProc=systemVar::sendRecvOrder[i][1];
            for (int nG = 0; nG <= mathVar::nGauss; nG++)
            {
                sendRecvMatchedBCGaussPtValues(sendingProc, receivingProc,Var,Buffer,nG);
            }
        }
    }

    /**
     * @brief Function send/recv U values at Gauss points on "matched" BC.
     */
    void synchSurfGaussU()
    {
        int sendingProc, receivingProc;
        for (int i=0; i<systemVar::sendRecvOrder_length; i++)
        {
            sendingProc=systemVar::sendRecvOrder[i][0];
            receivingProc=systemVar::sendRecvOrder[i][1];

            if (!flowProperties::massDiffusion)
            {
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    sendRecvMatchedBCGaussPtValues(sendingProc, receivingProc, surfaceFields::rho,surfaceFields::rho,nG);
                }
            }
            for (int nG = 0; nG <= mathVar::nGauss; nG++)
            {
                sendRecvMatchedBCGaussPtValues(sendingProc, receivingProc, surfaceFields::rhou,surfaceFields::rhou,nG);
                sendRecvMatchedBCGaussPtValues(sendingProc, receivingProc, surfaceFields::rhov,surfaceFields::rhov,nG);
                sendRecvMatchedBCGaussPtValues(sendingProc, receivingProc, surfaceFields::rhoE,surfaceFields::rhoE,nG);
            }
        }
    }

    void sendRecvSurfaceBCGaussPtArray(double**Var,double**Buffer)
    {
        int sendingProc, receivingProc;
        for (int i=0; i<systemVar::sendRecvOrder_length; i++)
        {
            sendingProc=systemVar::sendRecvOrder[i][0];
            receivingProc=systemVar::sendRecvOrder[i][1];
            for (int nG = 0; nG <= mathVar::nGauss; nG++)
            {
                sendRecvMatchedBCGaussPtValues(sendingProc, receivingProc,Var,Buffer,nG);
            }
        }
    }

    void synchGaussPtCoors()
    {
        int sendingProc, receivingProc;
        for (int i=0; i<systemVar::sendRecvOrder_length; i++)
        {
            sendingProc=systemVar::sendRecvOrder[i][0];
            receivingProc=systemVar::sendRecvOrder[i][1];

            if (sendingProc>receivingProc)
            {
                for (int nG = 0; nG <= mathVar::nGauss; nG++)
                {
                    //sendRecvMatchedBCGaussPtCoor(sendingProc,receivingProc,parallelBuffer::aCoor,meshVar::edgeGaussPoints_a,nG);
                    //sendRecvMatchedBCGaussPtCoor(sendingProc,receivingProc,parallelBuffer::bCoor,meshVar::edgeGaussPoints_b,nG);
                    sendRecvMatchedBCGaussPtCoors(sendingProc,receivingProc,meshVar::edgeGaussPoints_a,meshVar::edgeGaussPoints_a,nG);
                    sendRecvMatchedBCGaussPtCoors(sendingProc,receivingProc,meshVar::edgeGaussPoints_b,meshVar::edgeGaussPoints_b,nG);
                }
            }
        }
    }
}
