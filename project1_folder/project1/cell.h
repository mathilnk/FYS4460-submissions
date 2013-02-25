#ifndef CELL_H
#define CELL_H
#include<armadillo>
#include<vector>
#include"atom.h"
#include"atomnode.h"

class Cell
{
public:
    Cell();
    //vector<Atom*> myAtoms;
    void addAtom(Atom *a);
    bool removeAtom(int index);
    AtomNode* first;
    vector<Atom*> myAtoms;

    int myCellNumber;

    bool isEmpty();
    int numberOfAtoms;




};

#endif // CELL_H
