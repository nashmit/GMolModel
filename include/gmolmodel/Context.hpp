#ifndef __CONTEXT_HPP__
#define __CONTEXT_HPP__

#include "Robo.hpp"
//#include "SetupReader.hpp"

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

    // Use a SetupReader Object to read worlds information from a file
    void LoadWorldsFromSetup(SetupReader&);

private:
    std::vector<World *> worlds;
    std::vector<int> worldIndexes;

};

#endif //__CONTEXT_HPP__
