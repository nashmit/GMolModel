/**@file
Implementation of HamiltonianMonteCarloSampler class. **/

#include "Robo.hpp"
#include "ParaMolecularDecorator.hpp"

using namespace SimTK;
ParaMolecularDecorator::ParaMolecularDecorator(void)
{
    ;
}

void ParaMolecularDecorator::generateDecorations(const State& state,
        Array_<DecorativeGeometry>& geometry) 
{
    geometry.push_back( DecorativeLine( Vec3(0,0,0), Vec3(1,0,0) ) );
    geometry.push_back( DecorativeSphere() );
}

ParaMolecularDecorator::~ParaMolecularDecorator()
{
    ;
}
