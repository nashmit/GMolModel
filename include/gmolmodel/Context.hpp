#include "Robo.hpp"
#include "World.hpp"
#include "Sampler.hpp"

class Context{
private:
    World *p_world;
    Sampler *p_sampler;

public:
    Context(World *, Sampler *);
    ~Context();
    World * getWorld(void) const;
    World * updWorld(void);
    Sampler * getSampler(void);
};

