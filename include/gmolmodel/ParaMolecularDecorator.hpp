#ifndef __PARAMOLECULARDECORATOR_HPP__
#define __PARAMOLECULARDECORATOR_HPP__

/* -------------------------------------------------------------------------- *
 *                               Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

#include "Robo.hpp"

using namespace SimTK;

class ParaMolecularDecorator : public DecorationGenerator {
public:
    ParaMolecularDecorator(void);

    void generateDecorations(const State& state,
        Array_<DecorativeGeometry>& geometry);

    ~ParaMolecularDecorator(void);

private:
    Array_< Vec3 >  points;
    Array_< std::pair< Vec3, Vec3 > > lines;
};

#endif // __PARAMOLECULARDECORATOR_HPP__
