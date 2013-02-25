#include <iostream>
#include<cmath>
#include"lattice.h"
#include"atom.h"

#include <armadillo>
#include"normal.hpp"
#include "verlet_solver.h"
#include<time.h>
#include"cellsolver.h"
#include"atomnode.h"
#include"energytest.h"
#include"timetest.h"
#include"thermostattest.h";


//using namespace std;
//using namespace arma;

int main()
{
    new ThermostatTest();
    //new TimeTest();
    return 0;
};


