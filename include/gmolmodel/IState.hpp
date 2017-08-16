#ifndef __ISTATE_HPP__
#define __ISTATE_HPP__

#include "Robo.hpp"

class IState : public SimTK::State
{
public:
    virtual ~IState()=0;
};

#endif // __ISTATE_HPP__
