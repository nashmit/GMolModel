#ifndef __ICURRENTSTATE_HPP__
#define __ICURRENTSTATE_HPP__

class ICurrentState : public SimTK::State
{
    public:
        virtual ~ICurrentState()=0;
};

#endif // __ICURRENTSTATE_HPP__
