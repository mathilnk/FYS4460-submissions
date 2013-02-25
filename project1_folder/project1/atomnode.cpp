#include "atomnode.h"

AtomNode::AtomNode(Atom *atom)
{
    this->myAtom = atom;
}


bool AtomNode::isLast(){
    bool answer = false;
    if(nextAtom==NULL){
        answer = true;
    }
    return answer;
}

AtomNode* AtomNode::findLast(){
    AtomNode* tmp = this;
    while(tmp->nextAtom!=NULL){//!tmp->isLast()){
        tmp = tmp->nextAtom;
    }
    return tmp;
}

void AtomNode::addAtom(Atom *a, int numberOfAtomsInCell){

    AtomNode* tmp = this;
    AtomNode* last = this;
    int count = 0;
    while(tmp!=NULL){
        count++;
        last = tmp;
        tmp = tmp->nextAtom;
        cout<<count<<" "<<numberOfAtomsInCell<<endl;
    }
//    for(int i=0;i<numberOfAtomsInCell-1;i++){
//        tmp = tmp->nextAtom;
//        if(tmp == NULL){
//            cout<<"nullpointer"<<endl;
//        }
//    }

    last->nextAtom = new AtomNode(a);

}

void AtomNode::removeNodeAfterMe(){
    AtomNode* toBeRemoved = nextAtom;
    if(toBeRemoved!= NULL){
        AtomNode* newNext = toBeRemoved->nextAtom;
        nextAtom = newNext;
    }
}


Atom* AtomNode::getMyAtom(){
    return myAtom;
}
