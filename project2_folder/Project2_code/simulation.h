#ifndef SIMULATION_H
#define SIMULATION_H
#include<makematrix.h>
using namespace std;
class Simulation

{
public:
    Simulation();
    void simulate_first_porous_system();
    void simulate_first_porous_system_fromFile();
    void prepareSystem(string filename);
    void prepareSystem();
    void makePores(int num_pores, double r_min, double r_max);
    void solve(int timesteps);
    void thermalize(double T);
    void try_config();
    void test_ting();
    void thermalize_stuff_with_forces();
};

#endif // SIMULATION_H
