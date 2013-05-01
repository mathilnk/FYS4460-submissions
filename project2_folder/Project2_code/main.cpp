#include <iostream>
#include<cmath>

#include"atom.h"

#include <armadillo>
#include"normal.hpp"

#include<time.h>
#include"cellsolver.h"


#include"test.h"
#include"makematrix.h"
#include"simulation.h"


using namespace std;
using namespace arma;


int main(int argc, const char* argv[])
{
    string experiment_name = argv[1];
    string experiment_folder = argv[2];
    string config_file = argv[3];
    CellSolver* s = new CellSolver(experiment_name, experiment_folder, config_file);
    s->solve();
    //new Simulation();
    return 0;
};


