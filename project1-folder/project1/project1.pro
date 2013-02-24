TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt



SOURCES += main.cpp \
    lattice.cpp \
    atom.cpp \
    facecube.cpp \
    normal.cpp \
    verlet_solver.cpp \
    cellsolver.cpp \
    cell.cpp \
    cellcontainer.cpp \
    atomnode.cpp \
    energytest.cpp \
    timetest.cpp \
    thermostattest.cpp

HEADERS += \
    lattice.h \
    atom.h \
    facecube.h\
    normal.hpp \
    verlet_solver.h \
    cellsolver.h \
    cell.h \
    cellcontainer.h \
    atomnode.h \
    energytest.h \
    timetest.h \
    thermostattest.h

