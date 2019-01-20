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

	std::vector<double> calcCellCenteredValues(int valType);

	void exportNodeData(int iter);

	void exportCellCenteredData(int iter);
}
#endif // DGPOSTPROCESSLIB_H_INCLUDED