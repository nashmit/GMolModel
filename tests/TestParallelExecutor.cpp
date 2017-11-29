#include "Simbody.h"

#ifndef DEBUG
#define DEBUG 1
#endif

#ifdef DEBUG
#define TRACE(STR) printf("%s", STR); 
#else
#define TRACE(STR) 
#endif

#include <stdio.h>

class NonbondedForceTask : public SimTK::Parallel2DExecutor::Task {
public:
    NonbondedForceTask(int mat4x4[4][4], int& _total) 
	: total(_total) 
    {
	for(int i = 0 ; i < 4 ; i++ )
		for( int j = 0 ; j < 4 ; j++ )
			this->mat4x4[i][j] = mat4x4[i][j];
    }

    void initialize() {
        //TRACE("NonbondedForceTask::initialize BEGIN\n");
        //TRACE("NonbondedForceTask::initialize END\n");
	TRACE("NonbondedForceTask::initialize\n");

	sum.upd() = 0;

    }

    void finish() {
        TRACE("NonbondedForceTask::finish\n");
        //TRACE("NonbondedForceTask::finish END\n");
	printf("sum = %d \n", sum.get() );
	total += sum.get();
    }

    void execute(int body1, int body2) {
        TRACE("NonbondedForceTask::execute\n");
        //TRACE("NonbondedForceTask::execute END\n");
	
	printf("body1, body2, sum before : %d %d %d \n", body1, body2, sum.get() );

	sum.upd() += mat4x4[ body1 ][ body2 ];
    }

private:
	int mat4x4[4][4];
	SimTK::ThreadLocal<int> sum;
	int& total;
    
};

int main (int argc, char **argv)
{
	
    int mat4x4[][4] = { { 1, 1, 1, 1 } ,
		   	{ 1, 1, 1, 1 } ,
			{ 1, 1, 1, 1 } ,
			{ 1, 1, 1, 1 } };

    int numThreadsInUse = 3;
    int nofIncludedBodies = 4;

    int total = 0;

    NonbondedForceTask task( mat4x4, total);

    SimTK::ParallelExecutor *executor = new SimTK::ParallelExecutor(numThreadsInUse);

    SimTK::Parallel2DExecutor *nonbondedExecutor =  new SimTK::Parallel2DExecutor(nofIncludedBodies, *executor);

    nonbondedExecutor->execute(task, SimTK::Parallel2DExecutor::HalfMatrix);

    printf("Before finish: %d\n",total);

    return 0;
}
