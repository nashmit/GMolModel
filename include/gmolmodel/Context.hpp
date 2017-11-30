#include "Robo.hpp"

class Sampler;
class World;

class Context{
private:
    World *p_world;
    Sampler *p_sampler;

public:
    Context(World *, Sampler *);
    Context();
    ~Context();

    World * AddWorld(World *, Sampler *);

    World * getWorld(void) const;
    World * getWorld(int which) const;

    World * updWorld(void);
    World * updWorld(int which);

    Sampler * getSampler(void) const;
    Sampler * getSampler(int which) const;

private:
    std::vector<World *> worlds;
    std::vector<Sampler *> samplers;

};

