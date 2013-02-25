#ifndef LATTICE_H
#define LATTICE_H
#include<armadillo>
#include<facecube.h>
#include<vector>
#include<iostream>
#include<fstream>
#include"atom.h"
#include"normal.hpp"
//#include<random>
using namespace std;
using namespace arma;
class Lattice
{
public:
    Lattice();
    Lattice(int Nx, int Ny, int Nz, string element, double b, double T, double mass);
    vector<Atom*> allAtoms;
    void writeVMDfile(string Filename, string comment);
    void findPosAndMakeAtoms(vec posBase);
    double T;
    int numberOfAtoms;
    double mass;
    double Lx,Ly,Lz;
    double k;

private:
    int Nx,Ny,Nz,N;
    string element;
    double b;
    void makeLattice();

    double gauss(double s, double mean);
    double s,mean;
    double pi;
    int seed;



};

#endif // LATTICE_H
