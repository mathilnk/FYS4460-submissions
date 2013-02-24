#ifndef THERMOSTATTEST_H
#define THERMOSTATTEST_H
#include"cellsolver.h"
class ThermostatTest
{
public:
    ThermostatTest();
    void runBandA_Always();
    void runBUpDown(double first, double second, double first_timestep, double second_timestep);
    void runWithoutThermo();
    void varyTmeasureP(double T_start, double T_stop, int numOfPlots);
};
#endif // THERMOSTATTEST_H
