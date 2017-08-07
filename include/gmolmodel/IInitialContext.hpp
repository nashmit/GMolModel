#ifndef __IINITIALCONTEXT_HPP__
#define __IINITIALCONTEXT_HPP__

class IInitialContext
{
    public:
        virtual IInitialContext()=0;
        virtual ICurrentState* GetState()=0;
        virtual IIntegrator* GetIntegrator()=0;
        virtual ~IInitialContext()=0;
};

#endif // __IINITIALCONTEXT_HPP__
