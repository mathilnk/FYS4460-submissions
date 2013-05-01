#ifndef CELL_H
#define CELL_H
#include<armadillo>
#include<vector>
#include"atom.h"


class Cell
{
public:
    Cell();
    //vector<Atom*> myAtoms;
    void addAtom(Atom *a);
    bool removeAtom(int index);

    vector<Atom*> myAtoms;

    int myCellNumber;

    bool isEmpty();
    int numberOfAtoms;





};

#endif // CELL_H
