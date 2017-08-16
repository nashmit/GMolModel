#include "Context.hpp"

// Constructor

Context::Context(World *inp_p_world, Sampler * inp_p_sampler){
    p_world = inp_p_world;
    p_sampler = inp_p_sampler;
}


// Destructor

Context::~Context(){
    //
}

// Get world

World * Context::getWorld(void) const{
    return p_world;
}

// Get mutable world

World * Context::updWorld(void){
    return p_world;
}

// Get sampler

Sampler * Context::getSampler(void){
    return p_sampler;
}
