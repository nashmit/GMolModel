#include "Robo.hpp"

class Sampler;
class World;

class Context{
//private:
//    World *p_world;

public:
    Context(World *);
    Context();
    ~Context();

    World * AddWorld(World *);

    World * getWorld(void) const;
    World * getWorld(int which) const;

    World * updWorld(void);
    World * updWorld(int which);

private:
    std::vector<World *> worlds;

};

