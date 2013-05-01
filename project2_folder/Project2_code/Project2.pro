TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    atom.cpp \
    cell.cpp \
    cellcontainer.cpp \
    cellsolver.cpp \
    normal.cpp \
    test.cpp \
    makematrix.cpp \
    simulation.cpp

HEADERS += \
    atom.h \
    cell.h \
    cellcontainer.h \
    cellsolver.h \
    normal.hpp \
    test.h \
    makematrix.h \
    simulation.h

LIBS +=-fopenmp
LIBS += -lconfig++

COMMON_CXXFLAGS +=  -fopenmp
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS
