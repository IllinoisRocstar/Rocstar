//
// Created by agondolo on 9/14/18.
//

#ifndef _ROCSTAR_SURFDIVER_HPP_
#define _ROCSTAR_SURFDIVER_HPP_

#include "RocstarAction.h"

class FluidAgent;
class SolidAgent;

// stop and run surfdiver
class SurfDiver : public RocstarAction {
 public:
  SurfDiver(FluidAgent *fag, SolidAgent *sag);
  void init(double t);
  virtual void run(double t, double dt, double alpha);
 protected:
  //void read_file( const char *fname, const string &wname, double alpha);
  FluidAgent *fagent;
  SolidAgent *sagent;
  std::string outdir;
  std::string fluid_mesh_str, solid_mesh_str;
  int fluid_mesh, solid_mesh;
  int RFC_transfer, RFC_interpolate, RFC_readcntr, RFC_overlay;
  int RFC_write, RFC_read;
};

// run surfdiver if overlay mesh is missing
class SurfDiverAfterRemeshing : public SurfDiver {
 public:
  SurfDiverAfterRemeshing(FluidAgent *fag, SolidAgent *sag) :
      SurfDiver(fag, sag) {}
  void run(double t, double dt, double alpha);
};

#endif //_ROCSTAR_SURFDIVER_HPP_
