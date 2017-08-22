#include "Robo.hpp"

class Sampler;
class World;

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

