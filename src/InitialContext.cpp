#include "InitialContext.hpp"
#include "InitialContextRep.hpp"

// InitialContext constructor

InitialContext::InitialContext(){
    rep = new InitialContextRep;
}

// InitialContext destructor

InitialContext::~InitialContext(){
    delete rep;
}

// Insert an atom

void InitialContext::addAtom(int number){
    rep->addAtom(number);
}

// Simple print function

void InitialContext::Print(void){
    rep->Print();
}


