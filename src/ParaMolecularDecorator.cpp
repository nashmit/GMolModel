/**@file
Implementation of HamiltonianMonteCarloSampler class. **/

#include "Robo.hpp"
#include "ParaMolecularDecorator.hpp"

using namespace SimTK;
ParaMolecularDecorator::ParaMolecularDecorator(void)
{
    ;
}

void ParaMolecularDecorator::loadPoint(const Vec3 point)
{
    points.push_back(point);
}

void ParaMolecularDecorator::loadLine(const Vec3 p1, const Vec3 p2)
{
    lines.push_back(std::pair<Vec3, Vec3>(p1, p2));
}

void ParaMolecularDecorator::clearPoints(void)
{
    points.clear();
}

void ParaMolecularDecorator::clearLines(void)
{
    lines.clear();
}


void ParaMolecularDecorator::generateDecorations(const State& state,
        Array_<DecorativeGeometry>& geometry) 
{
    for (auto p = points.begin(); p != points.end(); p++) {
        geometry.push_back(DecorativeSphere(0.1));
    }
    geometry.push_back( DecorativeLine( Vec3(0,0,0), Vec3(1,0,0) ) );
    geometry.push_back( DecorativeSphere(0.5) );
}

ParaMolecularDecorator::~ParaMolecularDecorator()
{
    points.clear();
    lines.clear();
}
