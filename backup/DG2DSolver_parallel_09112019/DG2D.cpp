#include "DGIOLib.h"
#include "DGController.h"
#include "VarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <iostream>
#include <mpi.h>

using namespace std;

int main()
{
    auxUlti::functionsOfParallelComputing::prepareParallelCase();

	/*DISPLAY LOGO*/
    if (systemVar::currentProc==0)
    {
        IO::dispLogo();
        std::cout << "                     Welcome to DG2D solver console!!\n";
    }
    IO::getCase();

    //Phai doc so cores truoc tien
    IO::readNumberOfCores();

	while (systemVar::endKey==false)
	{
        auxUlti::getCommand();

        checkCommandLine(systemVar::cmd);

        MPI_Barrier(MPI_COMM_WORLD);
        Executer();
	}
    MPI_Finalize();
	return 0;
}
