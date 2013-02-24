#ifndef VERLET_SOLVER_H
#define VERLET_SOLVER_H
#include"lattice.h"
#include<math.h>
#include<time.h>

class Verlet_solver
{
public:
    Verlet_solver(Lattice &molecule);
    Lattice molecule;
    Lattice solve_one_time_step(double t, double dt, string filename);
    vec force(vec r);
    void solve(double t_start,int timesteps, double dt, string filename);
private:
    int current_time_step;
    void writeVMDfile_Verlet(char*filename, string comment);
    vec force_between(vec r_1, vec r_2);
    vec force_on(vec r, int index);
    void findForces();
    void cleanForces();
};

#endif // VERLET_SOLVER_H
