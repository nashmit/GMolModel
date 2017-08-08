#ifndef __IINITIALCONTEXT_HPP__
#define __IINITIALCONTEXT_HPP__

#include "Robo.hpp"
#include "ICurrentState.hpp"
//#include "MidVVIntegrator.hpp"

class IInitialContext
{
    public:
        //IInitialContext()=0;

        virtual ~IInitialContext() {}

        virtual void PrintAtomList()=0;

        //virtual ICurrentState* GetState()=0;
        //virtual MidVVIntegrator* GetIntegrator()=0;
};

#endif // __IINITIALCONTEXT_HPP__
