#ifndef __FIXMANTORQUE_HPP__
#define __FIXMANTORQUE_HPP__

/* -------------------------------------------------------------------------- *
 *                               Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

//==============================================================================
//                           CLASS FixmanTorque
//==============================================================================
/**
 **/

class FixmanTorque : public SimTK::Force::Custom::Implementation {
 public:
  SimTK::CompoundSystem *compoundSystem;
  int *fassno;
  int *flag;
  //World *Caller;

  FixmanTorque(SimTK::CompoundSystem *compoundSystem, SimTK::SimbodyMatterSubsystem& matter
            , int *fassno
            //, World *Caller
            );

  void calcForce(const SimTK::State& state, SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
    SimTK::Vector_<SimTK::Vec3>& particleForces, SimTK::Vector& mobilityForces) const;

  SimTK::Real calcPotentialEnergy(const SimTK::State& state) const;

  bool dependsOnlyOnPositions() const;

 private:
  SimTK::SimbodyMatterSubsystem& matter;
};

#endif //__FIXMANTORQUE_HPP__
