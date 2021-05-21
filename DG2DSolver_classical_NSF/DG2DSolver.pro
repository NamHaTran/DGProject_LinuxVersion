QT -= gui

CONFIG += c++11 console
CONFIG -= app_bundle

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    CommandCheck.cpp \
    ConstDeclaration.cpp \
    DG2D.cpp \
    DGAuxUltilitiesLib.cpp \
    DGController.cpp \
    DGIOLib.cpp \
    DGMath.cpp \
    DGMeshReaderLib.cpp \
    DGMessagesLib.cpp \
    DGPostProcessLib.cpp \
    DGProcLib.cpp \
    NonEquilibriumBCsLib.cpp \
    boundaryConditions/BCSupportFncs.cpp \
    boundaryConditions/DGBCsLib.cpp \
    boundaryConditions/bcVariables.cpp \
    boundaryConditions/fixedValue.cpp \
    boundaryConditions/matched.cpp \
    boundaryConditions/readBCInfor/pressure/readDirichletPresBCValue.cpp \
    boundaryConditions/readBCInfor/pressure/readMixedPresBCValue.cpp \
    boundaryConditions/readBCInfor/pressure/readNewmannPresBCValue.cpp \
    boundaryConditions/readBCInfor/readSymmetryBC.cpp \
    boundaryConditions/readBCInfor/supportReadingBCFuncs.cpp \
    boundaryConditions/readBCInfor/temperature/readDirichletTempBCValue.cpp \
    boundaryConditions/readBCInfor/temperature/readMixedTempBCValue.cpp \
    boundaryConditions/readBCInfor/temperature/readNewmannTempBCValue.cpp \
    boundaryConditions/readBCInfor/velocity/readDirichletVelBCValue.cpp \
    boundaryConditions/readBCInfor/velocity/readMixedVelBCValue.cpp \
    boundaryConditions/readBCInfor/velocity/readNewmannVelBCValue.cpp \
    boundaryConditions/symmetry.cpp \
    boundaryConditions/zeroGradient.cpp \
    debuggingFuncs.cpp \
    dynamicVarDeclaration.cpp \
    VarDeclaration.cpp \
    limiters/detectTroubleCells.cpp \
    limiters/limiterController.cpp \
    limiters/massDiffusion/massDiffusion.cpp \
    limiters/mathFunctions.cpp \
    limiters/pAdaptive/pAdaptive.cpp \
    limiters/parallelFuncs.cpp \
    limiters/positivityPreserving/positivityPreserving.cpp \
    limiters/varsOfLimiters.cpp \
    parallelFunctions/GaussPointData.cpp \
    parallelFunctions/cellData.cpp \
    parallelFunctions/generalParallelFuncs.cpp \
    parallelFunctions/parallelVariables.cpp

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
    CommandCheck.h \
    ConstDeclaration.h \
    DGAuxUltilitiesLib.h \
    DGController.h \
    DGIOLib.h \
    DGMath.h \
    DGMeshReaderLib.h \
    DGMessagesLib.h \
    DGPostProcessLib.h \
    DGProcLib.h \
    NonEquilibriumBCsLib.h \
    boundaryConditions/BCSupportFncs.h \
    boundaryConditions/DGBCsLib.h \
    boundaryConditions/bcVariables.h \
    boundaryConditions/fixedValue.h \
    boundaryConditions/matched.h \
    boundaryConditions/readBCInfor/pressure/readDirichletPresBCValue.h \
    boundaryConditions/readBCInfor/pressure/readMixedPresBCValue.h \
    boundaryConditions/readBCInfor/pressure/readNewmannPresBCValue.h \
    boundaryConditions/readBCInfor/readSymmetryBC.h \
    boundaryConditions/readBCInfor/supportReadingBCFuncs.h \
    boundaryConditions/readBCInfor/temperature/readDirichletTempBCValue.h \
    boundaryConditions/readBCInfor/temperature/readMixedTempBCValue.h \
    boundaryConditions/readBCInfor/temperature/readNewmannTempBCValue.h \
    boundaryConditions/readBCInfor/velocity/readDirichletVelBCValue.h \
    boundaryConditions/readBCInfor/velocity/readMixedVelBCValue.h \
    boundaryConditions/readBCInfor/velocity/readNewmannVelBCValue.h \
    boundaryConditions/symmetry.h \
    boundaryConditions/zeroGradient.h \
    debuggingFuncs.h \
    dynamicVarDeclaration.h \
    VarDeclaration.h \
    limiters/detectTroubleCell.h \
    limiters/limiterController.h \
    limiters/massDiffusion/massDiffusion.h \
    limiters/mathFunctions.h \
    limiters/pAdaptive/pAdaptive.h \
    limiters/parallelFuncs.h \
    limiters/positivityPreserving/positivityPreserving.h \
    limiters/varsOfLimiters.h \
    parallelFunctions/GaussPointData.h \
    parallelFunctions/cellData.h \
    parallelFunctions/generalParallelFuncs.h \
    parallelFunctions/parallelVariables.h

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
