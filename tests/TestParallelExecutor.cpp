#include "Simbody.h"

#ifndef DEBUG
#define DEBUG 1
#endif

#ifdef DEBUG
#define TRACE(STR) printf("%s", STR); 
#else
#define TRACE(STR) 
#endif

class NonbondedForceTask : public SimTK::Parallel2DExecutor::Task {
public:
    NonbondedForceTask() {
    }

    void initialize() {
        TRACE("NonbondedForceTask::initialize BEGIN\n");
        TRACE("NonbondedForceTask::initialize END\n");
    }

    void finish() {
        TRACE("NonbondedForceTask::finish BEGIN\n");
        TRACE("NonbondedForceTask::finish END\n");
    }

    void execute(int body1, int body2) {
        TRACE("NonbondedForceTask::execute BEGIN\n");
        TRACE("NonbondedForceTask::execute END\n");
    }

private:
    int temp;
};

int main (int argc, char **argv)
{
    int numThreadsInUse = 1;
    int nofIncludedBodies = 10;

    NonbondedForceTask task;

    SimTK::ParallelExecutor *executor = new SimTK::ParallelExecutor(numThreadsInUse);

    SimTK::Parallel2DExecutor *nonbondedExecutor =  new SimTK::Parallel2DExecutor(nofIncludedBodies, *executor);

    nonbondedExecutor->execute(task, SimTK::Parallel2DExecutor::HalfMatrix);

    return 0;
}
