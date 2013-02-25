#include "cellcontainer.h"

CellContainer::CellContainer(int Nx, int Ny, int Nz, double r_cut)
{
    this->Nx = Nx;
    this->Ny = Ny;
    this->Nz = Nz;
    this->N = Nx*Ny*Nz;
    this->numberOfAtoms = 0;
    this->r_cut = r_cut;
    first = true;
    makeCells();

}

void CellContainer::makeCells(){
    /*
      Makes the cells. Should only be used in the constructor
      */
    ////cout<<Nz<<endl;
    for(int i=0;i<N;i++){

        myCells.push_back(Cell());
    }
}


void CellContainer::addAtom(Atom * a){
    int cellNumber = findCellNumberForAtom(a);
    ////cout<<myCells.size()<<endl;
    myCells[cellNumber].addAtom(a);
    if(first){
        allAtoms.push_back(a);
    }
    ////cout<<a->force<<endl;

    numberOfAtoms++;
}

void CellContainer::removeAtom(int cellNumber, int atomNumber){
    Cell thisCell = myCells[cellNumber];
    thisCell.removeAtom(atomNumber);

}

int CellContainer::findCellNumberForAtom(Atom* a){
    /*
      This method finds the correct cellNumber for a given atom. It finds the correct cellCoordinates, and then uses the method findCellNumberForCell.
      */
    Col<int> cellPosition = zeros<Col<int> >(3);

    for(int i=0;i<3;i++){
        cellPosition[i] = a->position[i]/r_cut;
    }
    //cout<<a->position[0]<<endl;//atomet vi sender inn har ingen posisjon
    int cellNumber = findCellNumberForCell(cellPosition[0], cellPosition[1], cellPosition[2]);

    return cellNumber;
}


int CellContainer::findCellNumberForCell(int x, int y, int z){
    /*
      x,y,z are position of the cell in the matrix made of cells.NOT position of the particle.
      */
    return z*Nx*Ny + y*Nx + x;
}

Col<int> CellContainer::findCellPosition(int cellNumber){
    Col<int> cellPosition = zeros<Col<int> >(3);
    cellPosition[0] = cellNumber%Nx;
    cellPosition[1] = ((cellNumber -cellPosition[0])/Nx)%Ny;
    cellPosition[2] = (cellNumber - cellPosition[0] - Nx*cellPosition[1])/(Nx*Ny);
    return cellPosition;
}



vector<int>CellContainer::findMyNeighbors(int cellNumber){
    Col<int> myPosition = findCellPosition(cellNumber);
    vector<int>neighborNumbers;
    int neiPosX, neiPosY, neiPosZ;
    if(myCells.size()<=1){
        return neighborNumbers;
    }
    for(int z=1;z>=0;z--){
        for(int y=1;y>=-z;y--){
            for(int x=-1;x<2;x++){
                if(z==0 && y==0){
                    //legg på en på x
                    x=1;
                    neiPosX = (myPosition[0] + x + Nx)%Nx;
                    neiPosY = (myPosition[1] + y + Ny)%Ny;
                    neiPosZ = (myPosition[2] + z + Nz)%Nz;
                    neighborNumbers.push_back(findCellNumberForCell(neiPosX, neiPosY, neiPosZ));
                    break;

                }
                neiPosX = (myPosition[0] + x + Nx)%Nx;
                neiPosY = (myPosition[1] + y + Ny)%Ny;
                neiPosZ = (myPosition[2] + z + Nz)%Nz;
                neighborNumbers.push_back(findCellNumberForCell(neiPosX, neiPosY, neiPosZ));
            }
        }
    }
    return neighborNumbers;
}

//vector<int> CellContainer::findMyNeighbors(int cellNumber){
//    Col<int> myPosition = findCellPosition(cellNumber);
//    //Col<int> neighborPosition = zeros<Col<int> >(3);
//    int neiPosX=0, neiPosY=0, neiPosZ=0;
//    vector<int> neighborNumbers;// = zeros<Col<int> >(26);
//    if(myCells.size()<=1){
//        return neighborNumbers;
//    }
//    int index=0;
//    for(int z=0;z<2;z++){
//        for(int y=0;y<2;y++){
//            for(int x=0;x<2;x++){
//                if(x==0 && y==0 && z==0){

//                }else{
//                    neiPosX = (myPosition[0] + x + Nx)%Nx;
//                    neiPosY = (myPosition[1] + y + Ny)%Ny;
//                    neiPosZ = (myPosition[2] + z + Nz)%Nz;
//                    neighborNumbers.push_back(findCellNumberForCell(neiPosX, neiPosY, neiPosZ));
//                }
//            }
//        }
//    }
//    return neighborNumbers;

//}



void CellContainer::writeVMDfile(string filename, string comment, string element){
    /*
      This method writes a VMD-formatted file called filename.
      */
    ofstream writeToFile;

    const char* filename_ch = filename.c_str();
    writeToFile.open(filename_ch);
    writeToFile<<numberOfAtoms<<endl;
    writeToFile<<comment<<endl;


    /////////////////////////////////////////////////////////////////////////////////////////////////////
    //runs through the cells, and write the position and velocity of all the atoms in each cell to file//
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    Cell currentCell;
    Atom * a;

    for(int j=0;j<myCells.size();j++){
        currentCell = myCells[j];

        if(!currentCell.isEmpty()){
            for(int i=0;i<currentCell.numberOfAtoms;i++){
                a = currentCell.myAtoms[i];


                writeToFile<<element<<" ";
                for(int p=0;p<3;p++){
                    //AAngstrom
                    writeToFile<<a->position[p]<<" ";
                }
                for(int k=0;k<3;k++){
                     //aangstrom per femtosekund
                    writeToFile<<a->velocity[k]<<" ";
                }
                writeToFile<<"\n";

            }

        }

    }
}
