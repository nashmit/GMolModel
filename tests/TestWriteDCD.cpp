#include "format/PDBObject.hpp"
#include "format/TrajectoryObject.hpp"


using namespace std;

int main(int argc, char *argv[])
{


// no unitcell

int unitcell = 0;

PDBObject MOL;
TrajectoryObject MOL_traj;

MOL.readPDBfile(argv[1]);

printf("NumberAtoms %i \n", MOL.NumberAtoms);
printf("AtomsXcoord[0] %f \n", MOL.AtomsXcoord[0]);


MOL_traj.createTrajectory(argv[2], "dcd", MOL.NumberAtoms, unitcell);

printf("natoms %i \n", MOL_traj.natoms);

for(int i=0; i<10; i++)
{
  // change coordinates:
  for(int j =0; j<MOL_traj.natoms; j++)
  {
     MOL.AtomsXcoord[j] = MOL.AtomsXcoord[j] + i;
     MOL.AtomsYcoord[j] = MOL.AtomsYcoord[j] + i;
     MOL.AtomsZcoord[j] = MOL.AtomsZcoord[j] + i;
  }

  MOL_traj.appendTimestep("dcd", MOL.AtomsXcoord, MOL.AtomsYcoord, MOL.AtomsZcoord);
}

printf("nsets %i \n", MOL_traj.nsets);



return 0;

}
