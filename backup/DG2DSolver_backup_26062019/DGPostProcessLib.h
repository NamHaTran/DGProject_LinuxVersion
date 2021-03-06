#ifndef DGPOSTPROCESSLIB_H_INCLUDED
#define DGPOSTPROCESSLIB_H_INCLUDED
#include <vector>
namespace debugTool
{
	void checkElemSurPt(int ipoin);
	void checkPtsSurPt(int ipoin);
	void checkElemsSurElem(int ielem);
	void checkElemInfor(int elem);
	void checkPointValue(int element);
	//void writeMinRho_MinRhoe();
}

namespace DG2Tecplot
{
	std::vector<double> calcNodeValues(int valType);

	std::vector<double> calcNodeValuesAtBC(int ptAtBCId);

	std::vector<double> calcCellCenteredValues(int valType);

	void exportNodeData(int iter);

	void exportCellCenteredData(int iter);

	namespace calcNodeValueAtBCChildFuncs
	{
		namespace patch
		{
			std::vector<double> inFlow(int element, int edgeGrp, int a, int b);

			std::vector<double> outFlow(int element, int edgeGrp, int a, int b);
		}

		namespace wall
		{
			std::vector <double> noSlipIsoThermal(int element, int edgeGrp, int a, int b);

			std::vector <double> noSlipAdiabatic(int element, int a, int b);
		}

		std::vector <double> Symmetry(int element, int edge, int a, int b);
	}
}
#endif // DGPOSTPROCESSLIB_H_INCLUDED