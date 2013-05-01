#ifndef THERMOSTATTEST_H
#define THERMOSTATTEST_H
#include"cellsolver.h"
class Test
{
public:
    Test();
    void runBandA_Always();
    void runBUpDown(double first, double second, double first_timestep, double second_timestep);
    void runWithoutThermo();
    void varyTmeasureP(double T_start, double T_stop, int numOfPlots);
    void parallel_efficiency_1();
    void parallel_efficiency_2();
    void parallel_efficiency_3();
    void testLoad();
    double run_serial_code();
    double run_parallel(int num_pro, double N);
};
#endif // TEST_H
