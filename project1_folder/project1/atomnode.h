#ifndef ATOMNODE_H
#define ATOMNODE_H
#include<atom.h>

class AtomNode
{
public:
    AtomNode();
    AtomNode(Atom * atom);
    AtomNode* nextAtom;
    Atom * myAtom;
    Atom* getMyAtom();
    bool isLast();
    AtomNode* findLast();
    void addAtom(Atom * a, int numberOfAtomsInCell);
    void removeNodeAfterMe();
};

#endif // ATOMNODE_H
