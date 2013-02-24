#ifndef ENERGYTEST_H
#define ENERGYTEST_H
#include<armadillo>
#include"cellsolver.h"
#include<sstream>
class EnergyTest
{
public:
    EnergyTest();
    void runFor(double dt, int numOfTimeSteps, string trial);
};

#endif // ENERGYTEST_H
