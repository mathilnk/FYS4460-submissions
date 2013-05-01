#ifndef CELLCONTAINER_H
#define CELLCONTAINER_H
#include"cell.h"

class CellContainer
{
public:
    CellContainer();
    CellContainer(int Nx, int Ny, int Nz, double r_cut);
    int Nx,Ny,Nz,N, numberOfAtoms, numOfMovableAtoms, numOfNotMovAtoms;
    double r_cut;
    vector<Cell> myCells;
    vector<Atom*> allAtoms;

    int findCellNumberForAtom(Atom * a);
    int findCellNumberForCell(int x, int y, int z);
    void addAtom(Atom * a);
    void removeAtom(int atomIndexInAllAtoms);
    void writeVMDfile(string filename, string comment, string element);
    void writeOnlyMovAtomsVMD(string filename, string comment, string element);
    void writeOnlyNOTMovAtomsVMD(string filename, string comment, string element);
    Col<int> findCellPosition(int cellNumber);
    vector<int> findMyNeighbors(int cellNumber);
    bool first;
    int index;
    void updateCells();
private:
    void makeCells();

};

#endif // CELLCONTAINER_H
