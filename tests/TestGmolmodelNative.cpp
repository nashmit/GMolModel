#include "simmain.hpp"

int main(int argc, char **argv)
{
  SymSystem *sys;
  TARGET_TYPE *shm;
  int nosteps;
  int ntrials;
  int steps_per_trial;
  double temperature;
  double delta_t;
  int trouble;
  int sweep;
  int SHMSZ;
  TARGET_TYPE **coords;
  std::string ligdir;
  std::string gaff_fn;
  temperature = 300.0;
  delta_t = 0.0015;
  trouble = 0;
 
  vector3 *x;
  vector3 *f;

  int parser_argc = 6;
  ligdir += '/';

  /* Options are:
    IC: Internal Coordinates: fully flexible
    TD: Torsional Dynamics: only torsions are flexible 
    RR: Rigid Rings: torsional dynamics with rigid rings
   */
  const char *parser_argv[6] = {
  "-ligdir", ligdir.c_str(),
  "-gaff",  gaff_fn.c_str(),
  "-ictd", "TD"
  };

  //+++++++ Simbody PART ++++++++++++
  bArgParser parser(parser_argc, parser_argv);
  parser.Print();

  // * Read number of atoms  * //
  std::string line;
  ifstream in_stream;
  in_stream.open(parser.mol2F);
  while(in_stream >> line){
      std::cout<<"mol2 lines: "<<line<<std::endl;
  }
  in_stream.close();
 
//  //int natoms = configuration->dimensions[0]; // WATCHOUT
//  int natoms;
//  int order[natoms+2];
//  int acceptance;
//  order[natoms] = 1;
//  order[natoms+1] = 1945;
//  acceptance = order[natoms];
//
//  TARGET_TYPE **indexMap = NULL;
//  TARGET_TYPE *PrmToAx_po = NULL;
//  TARGET_TYPE *MMTkToPrm_po = NULL;
//
//  indexMap = new TARGET_TYPE*[(natoms)];
//  int _indexMap[natoms][3];
//  PrmToAx_po = new TARGET_TYPE[natoms];
//  MMTkToPrm_po = new TARGET_TYPE[natoms];
//
//  int natoms3 = 3*(natoms);
//  int arrays_cut = 2 + 4*natoms3;
//
//  SHMSZ = (
//    2*sizeof(TARGET_TYPE) +       // Counter and flag
//    natoms3*sizeof(TARGET_TYPE) + // Positions
//    natoms3*sizeof(TARGET_TYPE) + // Velocities
//    natoms3*sizeof(TARGET_TYPE) + // indexMap
//    natoms3*sizeof(TARGET_TYPE) + // Gradients
//    5*sizeof(TARGET_TYPE) +       // ac + 0 -> 4  // step, nosteps, temperature, timestep, trouble
//    1*sizeof(TARGET_TYPE) +       // ac + 5       // potential energy
//    2*sizeof(TARGET_TYPE) +       // ac + 6->7    // DAE step done; Move accepted
//    1*sizeof(TARGET_TYPE) +       // ac + 8       // KE
//    1*sizeof(TARGET_TYPE) +       // ac + 9       // steps_per_trial
//    1*sizeof(TARGET_TYPE) //+     // ac + 10      // trial
//  );
//  shm = new TARGET_TYPE[SHMSZ];
//
//  std::cout<<"parser.mol2F "<<parser.mol2F<<std::endl<<std::flush;
//  std::cout<<"parser.rbF "<<parser.rbF<<std::endl<<std::flush;
//  std::cout<<"parser.gaffF "<<parser.gaffF<<std::endl<<std::flush;
//  std::cout<<"parser.frcmodF "<<parser.frcmodF<<std::endl<<std::flush;
//  std::cout<<"parser.ictdF "<<parser.ictdF<<std::endl<<std::flush;
//
//  sys = new SymSystem(
//    parser.mol2F, parser.rbF, parser.gaffF, parser.frcmodF,
//    parser.ictdF, 
//    PrmToAx_po, MMTkToPrm_po,
//    pyFFEvaluatorObject,
//    p_energy_po,
//    configuration,
//    universe_spec,
//    shm
//  );
//  for(int i=0; i<natoms; i++){
//    _indexMap[i][2] = order[i];
//  }
//
//  // * Memory alloc for convinient arrays * //
//  coords = new TARGET_TYPE*[sys->mr->natms];
//  TARGET_TYPE **vels = new TARGET_TYPE*[sys->mr->natms];
//  TARGET_TYPE **inivels = new TARGET_TYPE*[sys->mr->natms];
//  TARGET_TYPE **grads = new TARGET_TYPE*[sys->mr->natms];
//
//  // * Seed the random number generator * //
//  srand (time(NULL));
//
//  //+++++++ SERVER PART +++++++++
//  float mysweep = 0;
//  shm[0] = CNT_START;
//
//  int cli;
//  sweep = 0;
//  bool system_initialized = false;
//  int big_loop_i;
//
//  shm[1] = CLIENT_FLAG;
//        int a;
//        int start, start_1, stop;
//        int i = natoms*2*3 + 2;
//
//        // * Assign convenient pointers for order * /
//        start = 2 + natoms3 + natoms3; 
//        for(a=0; a<natoms; a++){
//          indexMap[a] = &(shm[a*3 + start]);
//        }
//        // * Assign order in shm[][2] * //
//        start = 2 + natoms3 + natoms3;
//        cli = start;
//        for (a=0; a<natoms; a++){ // Order
//          cli += 2;
//          shm[cli++] = order[a];
//        }
//
//        cli = 2;
//        for (a=0; a<natoms; a++){ // Positions
//          shm[cli++] = x[a][0]; shm[cli++] = x[a][1]; shm[cli++] = x[a][2];
//        }
//        for (a=0; a<natoms; a++){ // Velocities
//          shm[cli++] = .0; //vels[a][0];//shm[cli++] = v[a][0];
//          shm[cli++] = .0; //vels[a][1];//shm[cli++] = v[a][1];
//          shm[cli++] = .0; //vels[a][2];//shm[cli++] = v[a][2];
//        }
//        for (a=0; a<natoms; a++){ // Order
//          cli += 2;
//          shm[cli++] = order[a];
//        }
//        for (a=0; a<natoms; a++){ // Forces
//          shm[cli++] = .0; shm[cli++] = .0; shm[cli++] = .0;
//        }
//        shm[cli++] = big_loop_i;          // Step
//        shm[cli++] = 50.0; //(TARGET_TYPE)nosteps;    // Number of steps
//        shm[cli++] = temperature; // Temeperature
//        shm[cli++] = delta_t;
//        shm[cli++] = trouble;
//        shm[cli++] = p_energy_po->energy;
//        cli++; // DAE set only by the server
//        shm[cli] = acceptance + 0.0;
//
//        //Set the client flag
//        shm[1] = CLIENT_FLAG;
//
//      shm[arrays_cut + 6] = 13.0; // set DAE to done
//
//      // * Assign convenient pointers for coordinates * //
//      int shm_coords_sup = (natoms3)+2;
//      int j=-1, k=0;
//      start = 2; 
//      for(j=0; j<natoms; j++){
//        coords[j] = &(shm[j*3 + start]);
//      }
//
//      mysweep = shm[0] - CNT_START;
//
//      // * Assign convenient pointers for velocities * //
//      start = 2 + natoms3;
//      start_1 = start - 1;
//      stop = 2 + natoms3 + natoms3;
//      for(j=0; j<natoms; j++){
//        vels[j]    = &(shm[j*3 + start]);
//        inivels[j] = &(shm[j*3 + start]);
//      }
//
//
//      // * Assign convenient pointers for gradients * //
//      start = 2 + natoms3 + natoms3 + natoms3;
//      start_1 = start - 1;
//      stop = 2 + natoms3 + natoms3 + natoms3 + natoms3;
//      for(j=0; j<natoms; j++){
//        grads[j] = &(shm[j*3 + start]);
//      }
//
//      // * Rewrite indexMap to memory * /
//      start = 2 + natoms3 + natoms3;
//      start_1 = start - 1;
//      stop = 2 + natoms3 + natoms3 + natoms3;
//
//      // * Initialize Simulation * /
//      TARGET_TYPE timestep, mytimestep;
//      mytimestep = shm[arrays_cut + 3];
//
//      // * Build bMainResidue and fill indexMap * /
//      sys->InitSimulation(coords, vels, inivels, indexMap, grads, mytimestep, true);
//      system_initialized = true;
}

//boost::python::object GCHMCIntegrator::Call(
//  int nosteps_arg,
//  int steps_per_trial_arg,
//  TARGET_TYPE temp_arg,
//  TARGET_TYPE ts,
//  int pyseed,
//  int _massMatNumOpt, // EU
//  int _metroFixmanOpt, // EU
//  double _lj14sf //--
//){
//  #ifdef DEBUG_TIME
//  //boost::timer boost_timer;
//  #endif
//
//  int ntrials_arg = 0;
//  if(nosteps_arg%steps_per_trial_arg){
//    fprintf(stderr, 
//      "GCHMCIntegrator::Call(): Incorrect nosteps/steps_per_trial combination: %d/%d\n",
//       nosteps_arg, steps_per_trial_arg);
//    exit(1);
//  }
//  ntrials_arg = nosteps_arg/steps_per_trial_arg;
//
//  int tx, tshm;
//  boost::python::list booConfList;
//  int i, j;
//  double **retConfsPois = new double* [ntrials_arg];
//  double *retPotEsPoi = new double[ntrials_arg];
//  double *accs = new double;
//  *accs = 0.0;
//  vector3 *xCall;
//
//  PyObject *SrcConfObj = PyObject_CallMethod(this->universe, "copyConfiguration", NULL);
//  // START GET CURR CONF
//  PyObject *CallConfArr = PyObject_GetAttrString(SrcConfObj, "array");
//  PyArrayObject *CallConfiguration = (PyArrayObject *)CallConfArr;
//  xCall = (vector3 *)CallConfiguration->data;
//  // STOP  GET CURR CONF
//  PyObject *ConfObjs[ntrials_arg];
//  PyObject *ConfArrs[ntrials_arg];
//  PyArrayObject *Confs[ntrials_arg];
//  boost::python::object booConfObjs[ntrials_arg];
//  boost::python::list booPotEs;
//  
//  for(i=0; i<ntrials_arg; i++){
//    ConfObjs[i] = PyObject_CallMethod(this->universe, "copyConfiguration", NULL);
//    ConfArrs[i] = PyObject_GetAttrString(ConfObjs[i], "array");
//    Confs[i] = (PyArrayObject *)(ConfArrs[i]);
//    retConfsPois[i] = (double *)(Confs[i])->data;
//  }
//
//  // Put curr conf from universe into shm
//  for(int a=0; a<sys->mr->natms; a++){
//    tx = (int)(sys->indexMap[ a ][2]);
//    tshm = ((int)(sys->indexMap[ a ][1]) * 3) + 2;
//    shm[ tshm +0] = xCall[tx][0];
//    shm[ tshm +1] = xCall[tx][1];
//    shm[ tshm +2] = xCall[tx][2];
//  }
//
//  *sys->pyseed = pyseed;
//  sys->massMatNumOpt = _massMatNumOpt; // EU
//  sys->metroFixmanOpt = _metroFixmanOpt; // EU
//  sys->lj14sf = _lj14sf; //--
//  sys->sysRetConfsPois = retConfsPois;
//  sys->sysRetPotEsPoi = retPotEsPoi;
//  sys->sysAccs = accs;
//
//  if(nosteps_arg % ntrials_arg){
//    fprintf(stderr, "nosteps_arg % ntrials_arg not 0\n");
//    exit(1);
//  }
//  this->nosteps = nosteps_arg;
//  this->ntrials = ntrials_arg;
//  this->temperature = temp_arg;
//  this->delta_t = ts;
//  ++(this->sweep);
//
//  int arrays_cut = 2 + 4*3*(sys->mr->natms);
//  shm[arrays_cut + 0] = 0.0; // step (will be incremented in MidVV
//  shm[arrays_cut + 1] = (TARGET_TYPE)(this->nosteps);
//  shm[arrays_cut + 2] = this->temperature;
//  shm[arrays_cut + 3] = delta_t; // picosec
//  shm[arrays_cut + 7] = 1.0;     // acceptance always 1 !!
//  shm[arrays_cut + 9] = nosteps/ntrials; // steps_per_trial
//  shm[arrays_cut + 10] = 0.0;    // initialize trial
//
//  // * Get <Universe>'s pyo_evaluator //
//  PyObject *pyo_evaluator;
//  PyObject *cpyo_evaluator;
//  pyo_evaluator = PyObject_CallMethod(this->universe, "energyEvaluator", NULL);
//  cpyo_evaluator = PyObject_CallMethod(pyo_evaluator, "CEvaluator", NULL);
//
//  this->pyFFEvaluatorObject = (PyFFEvaluatorObject *)cpyo_evaluator;
//  sys->pyFFEvaluatorObject = this->pyFFEvaluatorObject;
//
//  #ifdef DEBUG_TIME
//  //std::cout<<"Time simmain nosteps"<<this->nosteps<<" time "<<boost_timer.elapsed()<<std::endl;
//  #endif
//  sys->Advance(this->nosteps); 
//  #ifdef DEBUG_TIME
//  //std::cout<<"Time simmain nosteps"<<this->nosteps<<" time "<<boost_timer.elapsed()<<std::endl;
//  #endif
//
//  xCall = (vector3 *)configuration->data;
//
//  for(int a=0; a<sys->mr->natms; a++){
//    tx = (int)(sys->indexMap[ a ][2]);
//    tshm = ((int)(sys->indexMap[ a ][1]) * 3) + 2;
//    xCall[tx][0] = shm[ tshm +0];
//    xCall[tx][1] = shm[ tshm +1];
//    xCall[tx][2] = shm[ tshm +2];
//  }
//
//  for(i=0; i<ntrials_arg; i++){
//    booConfObjs[i] = boost::python::object(boost::python::handle<>( boost::python::borrowed((PyObject *)(ConfObjs[i])) ));
//    booConfList.append(booConfObjs[i]);
//
//    booPotEs.append(retPotEsPoi[i]);  // potential energies
//  }
//
//  double ntrials = (shm[arrays_cut + 1]) / (shm[arrays_cut + 9]);
//  double booAccPerTrials = *accs / ntrials;  //float(acc)/float(ntrials)
//  boost::python::list booRetList;
//  booRetList.append(booConfList);
//  booRetList.append(booPotEs);
//  //booRetList.append(booAccPerTrials); // RESTORE
//  booRetList.append(int(*accs));
//  booRetList.append(int(ntrials));
//  booRetList.append(ts); // delta_t
//
//  delete accs;
//
//  return booRetList; // (xs, energies, acc, ntrials, delta_t)
//}
//
//
//BOOST_PYTHON_MODULE(GCHMC)
//{
//  class_<GCHMCIntegrator>("GCHMCIntegrator", init<PyObject *, std::string, std::string>())
//    .def_readonly("nosteps",     &GCHMCIntegrator::nosteps)
//    .def_readonly("ntrials",     &GCHMCIntegrator::ntrials)
//    .def_readonly("steps_per_trial",     &GCHMCIntegrator::steps_per_trial)
//    .def_readonly("temperature", &GCHMCIntegrator::temperature)
//    .def_readonly("delta_t",     &GCHMCIntegrator::delta_t)
//    .def_readonly("trouble",     &GCHMCIntegrator::trouble)
//    .def_readonly("sweep",       &GCHMCIntegrator::sweep)
//    .def("Call", &GCHMCIntegrator::Call)
//    .def("Clear", &GCHMCIntegrator::Clear)
//  ;
//}



