#ifndef FACECUBE_H
#define FACECUBE_H
#include <armadillo>
#include "atom.h"

using namespace std;
using namespace arma;

class FaceCube
{
public:
    FaceCube(double b, string element);
    double b;
    string element;
    mat mat_cube;

private:
    void makeCube();


};

#endif // FACECUBE_H
