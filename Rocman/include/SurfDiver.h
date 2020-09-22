#ifndef _ROCSTAR_SURFDIVER_HPP_
#define _ROCSTAR_SURFDIVER_HPP_

#include "Action.h"

class FluidAgent;
class SolidAgent;

// stop and run surfdiver
class SurfDiver : public Action {
 public:
  SurfDiver(FluidAgent *fag, SolidAgent *sag);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

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
  void run(double t, double dt, double alpha) override;
};

#endif //_ROCSTAR_SURFDIVER_HPP_
