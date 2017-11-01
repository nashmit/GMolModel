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
    NonbondedForceTask(int& number) {
	this->number = number;
    }

    void initialize() {
        //TRACE("NonbondedForceTask::initialize BEGIN\n");
        //TRACE("NonbondedForceTask::initialize END\n");
	TRACE("NonbondedForceTask::initialize\n");
    }

    void finish() {
        TRACE("NonbondedForceTask::finish\n");
        //TRACE("NonbondedForceTask::finish END\n");
    }

    void execute(int body1, int body2) {
        TRACE("NonbondedForceTask::execute\n");
        //TRACE("NonbondedForceTask::execute END\n");
    }

private:
	int number;
    
};

int main (int argc, char **argv)
{
    int numThreadsInUse = 3;
    int nofIncludedBodies = 4;

    int number = 0;

    NonbondedForceTask task(number);

    SimTK::ParallelExecutor *executor = new SimTK::ParallelExecutor(numThreadsInUse);

    SimTK::Parallel2DExecutor *nonbondedExecutor =  new SimTK::Parallel2DExecutor(nofIncludedBodies, *executor);

    nonbondedExecutor->execute(task, SimTK::Parallel2DExecutor::HalfMatrix);

    return 0;
}
