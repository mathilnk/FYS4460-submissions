#ifndef ATOM_H
#define ATOM_H
#include<armadillo>
#include<iostream>
using namespace arma;
using namespace std;
class Atom
{
public:
    Atom();
    Atom(vec position, vec velocity, string element);
    vec  position;
    vec  real_position;
    vec initial_position;
    vec velocity;
    vec force;
    string element;
    int myindexInAllAtoms;
    bool canMove;
    //double myDisplacement;
};

#endif // ATOM_H
