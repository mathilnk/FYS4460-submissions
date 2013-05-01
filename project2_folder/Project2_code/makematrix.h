#ifndef MAKEMATRIX_H
#define MAKEMATRIX_H
#include"cellsolver.h"
class MakeMatrix
{
public:
    MakeMatrix();
    CellSolver * mySolver;
    CellSolver* makeAndThermalize(int N,double b, double T);
    void cutCyllinder(vec centrum, double R);
    CellSolver* makeFromFile(string filename, int N, double b, int max_time);
    CellSolver* runCyllinder();
    void makePorousMatrix(double r_max, double r_min, int max_time);
};

#endif // MAKEMATRIX_H
