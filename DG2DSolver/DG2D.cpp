#include "DGIOLib.h"
#include "DGController.h"
#include "VarDeclaration.h"
#include <iostream>
using namespace std;

int main()
{
	/*DISPLAY LOGO*/
	IO::dispLogo();
	std::cout << "                     Welcome to DG2D solver console!!\n";
	IO::getCase();
	//PreProcessing();
	while (systemVar::endKey==false)
	{
		std::cout << ">> ";
		std::cin >> systemVar::cmd;
		Executer(systemVar::cmd);
	}
	return 0;
}
