#include "cell.h"


Cell::Cell()
{

    numberOfAtoms = 0;


}


void Cell::addAtom(Atom *a){
//    cout<<numberOfAtoms<<endl;
//    if(numberOfAtoms>0){
//        first->addAtom(a, numberOfAtoms);
//    }else{
//        first = new AtomNode(a); //add the first element to the cell

//    }
//    numberOfAtoms++;

    myAtoms.push_back(a);
    numberOfAtoms++;
}

bool Cell::isEmpty(){
    return numberOfAtoms==0;
}

//bool Cell::removeAtom(int index){
//    if(index==0){
//        if(first==NULL){
//            return false;
//        }
//        first = first->nextAtom;
//        numberOfAtoms--;
//        return true;
//    }else{
//        AtomNode* tmp = first;
//        for(int i=0; i<index;i++){
//            if(tmp == NULL){
//                return false;
//            }
//            tmp = tmp->nextAtom;
//         }
//        tmp->removeNodeAfterMe();
//        numberOfAtoms--;
//        return true;
//    }


//}
