#include "Context.hpp"

// Includes to get the structure of additional classes

#include "World.hpp"
#include "Sampler.hpp"

// Default constructor
Context::Context(void){
}

// Constructor
Context::Context(World *inp_p_world){
    worlds.push_back(inp_p_world);
}

// Add another world and a sampler to the context
World * Context::AddWorld(World *inp_p_world){
    worlds.push_back(inp_p_world);
    return worlds.back();
}

// Destructor
Context::~Context(){
    worlds.clear();
}

// Get world
World * Context::getWorld(void) const{
    return worlds.back();
}

// Get a specific world
World * Context::getWorld(int which) const{
    return worlds[which];
}

// Get the last mutable world
World * Context::updWorld(void){
    return worlds.back();
}

// Get a mutable specific world
World * Context::updWorld(int which){
    return worlds[which];
}

