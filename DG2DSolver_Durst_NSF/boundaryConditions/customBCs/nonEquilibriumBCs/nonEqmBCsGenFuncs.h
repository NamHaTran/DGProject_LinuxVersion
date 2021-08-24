#ifndef NONEQMBCSGENFUNCS_H
#define NONEQMBCSGENFUNCS_H
#include <string>

namespace nonEquilibriumBCs {
    void readSurfaceValues(std::string Loc);

    void writeSurfaceValues(std::string Loc);

    void updateBCs();
}
#endif // NONEQMBCSGENFUNCS_H
