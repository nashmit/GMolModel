#include "Context.hpp"
#include "ContextRep.hpp"

// Context constructor

Context::Context(){
    rep = new ContextRep;
}

// Context destructor

Context::~Context(){
    delete rep;
}

// Simple print function

void Context::Print(void){
    rep->Print();
}


