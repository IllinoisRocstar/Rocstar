/* *******************************************************************
 * Rocstar Simulation Suite                                          *
 * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
 *                                                                   *
 * Illinois Rocstar LLC                                              *
 * Champaign, IL                                                     *
 * www.illinoisrocstar.com                                           *
 * sales@illinoisrocstar.com                                         *
 *                                                                   *
 * License: See LICENSE file in top level of distribution package or *
 * http://opensource.org/licenses/NCSA                               *
 *********************************************************************/
/* *******************************************************************
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
 * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
 * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
 *********************************************************************/
// $Id: basic_actions.h,v 1.53 2009/08/27 14:04:54 mtcampbe Exp $

#ifndef _BASIC_ACTIONS_H_
#define _BASIC_ACTIONS_H_

#include "rocman.h"

class FluidAgent;
class SolidAgent;
class BurnAgent;

class DummyAction : public Action {
public:
  DummyAction() : Action(0, (const char **)NULL, NULL, NULL, (char *)"DummyAction") {}
  void init(double t) {}
  void run(double t, double dt, double alpha) {}
};

class SetValueDouble : public Action {
 public:
  SetValueDouble(const std::string at, const double val);
  void init(double t);
  void run(double t, double dt, double alpha);
private:
  int    attr_hdl;
  double v;
};

class SetZero : public SetValueDouble {
 public:
  SetZero(const std::string at);
};

class CopyValue: public Action {
 public:
  CopyValue(const std::string from, const std::string to, int *cond=NULL);
  void init(double t);
  void run(double t, double dt, double alpha);
private:
  int    from_hdl, to_hdl;
  int * condition;
};

class DummyPrint: public Action {
 public:
  DummyPrint(BurnAgent *bag, SolidAgent *sag, FluidAgent *fag, const std::string l);
  void init(double t);
  void run(double t, double dt, double alpha);
private:
  BurnAgent *bagent;
  SolidAgent *sagent;
  FluidAgent *fagent;
  std::string label;
};

class BCInvoker: public Action {
 public:
  BCInvoker(Agent *ag, int l=1);
  void init(double t);
  void run(double t, double dt, double alpha);
private:
  Agent *agent;
  int level;
};

class GMInvoker: public Action {
 public:
  GMInvoker(Agent *ag);
  void init(double t);
  void run(double t, double dt, double alpha);
private:
  Agent *agent;
};

class BCInitInvoker: public Action {
 public:
  BCInitInvoker(Agent *ag);
  void init(double t);
  void run(double t, double dt, double alpha);
private:
  Agent *agent;
};

// stop and run surfdiver
class SurfDiver: public Action {
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
class SurfDiverAfterRemeshing: public SurfDiver {
 public:
  SurfDiverAfterRemeshing(FluidAgent *fag, SolidAgent *sag):
         SurfDiver(fag, sag) {}
  virtual void run(double t, double dt, double alpha);
};

// POST_UPDATE_FLUID
class ComputeFluidLoad_ALE: public Action {
  public:								
    explicit ComputeFluidLoad_ALE( FluidAgent *fag, SolidAgent *sag, const std::string f_pf, const std::string fb_mdot, const std::string b_rb, const std::string f_ts);
    void init(double t);
    void run(double t, double dt, double alpha);
  private:
    int traction_mode;
    FluidAgent *fagent;
    SolidAgent *sagent;
    int f_pf_hdl, f_tf_hdl, f_ts_hdl;
    int fb_ts_hdl, fb_pf_hdl;
    int fb_mdot_hdl, fb_rhof_alp_hdl, fb_mdot_tmp_hdl;
    int b_rb_hdl, fb_nf_alp_hdl, fb_tf_hdl;
};

// UPDATE_INBUF_GM_FLUID
class ComputeMeshMotion : public Action {
 public:
  ComputeMeshMotion(FluidAgent *ag, const std::string a_vm, const std::string f_du_alp, double z);
  void init(double t);
  void run(double t, double dt, double alpha);
 private:
  FluidAgent *fagent;
  double zoom;
  int a_vm_hdl, f_du_alp_hdl, f_zoom_hdl;
};

class ComputeFaceCenters : public Action {
 public:
  ComputeFaceCenters(BurnAgent *ag, const std::string b_nc, const std::string b_cnts,
		     const std::string b_nrml="");
  void init(double t);
  void run(double t, double dt, double alpha);
 private:
  BurnAgent *agent;
  int b_nc_hdl, b_cnts_hdl, b_nrml_hdl;
  int SURF_n2f, SURF_fn;
};

class FluidPropagateSurface : public Action {
 public:
  FluidPropagateSurface(FluidAgent *ag, BurnAgent *bag, const std::string b_rb, const std::string a_vm, double z);
  void init(double t);
  void run(double t, double dt, double alpha);
 private:
  FluidAgent *fagent;
  BurnAgent *bagent;
  double zoom;
  int p_rb_hdl, b_rb_hdl, fb_rb_hdl, p_cnstr_type, p_pmesh_hdl, p_vm_hdl, a_vm_hdl;
  int p_cflag_hdl, p_pos_hdl, p_bflag_hdl;
  int MAP_reduce_maxabs;
  int MAP_reduce_minabs;
  int PROP_set_cnstr, PROP_propagate;
  int PROPCON_find_intersections, PROPCON_constrain_displacements;
  int PROPCON_burnout, PROPCON_burnout_filter;
};

class MassTransfer: public Action {
 public:
  MassTransfer(FluidAgent *ag, BurnAgent *bag, const std::string b_rhos, const std::string b_rb, const std::string f_mdot);
  void init(double t);
  void run(double t, double dt, double alpha);
 private:
  FluidAgent *fagent;
  BurnAgent *bagent;
  int f_mdot_hdl, b_rhos_hdl, b_rb_hdl, fb_mdot_hdl;
};

class ZoomInterface : public Action {
 public:
  ZoomInterface(FluidAgent *ag, BurnAgent *bag, const std::string fb_mdot_alp, double z);
  void init(double t);
  void run(double t, double dt, double alpha);
 private:
  FluidAgent *fagent;
  BurnAgent *bagent;
  int b_rb_alp_hdl, fb_rhof_alp_hdl, fb_mdot_alp_hdl, rhos_hdl;
  double zoom;
};

class ComputeRhofvf : public Action {
 public:
  ComputeRhofvf(FluidAgent *ag, std::string f_vs_alp, std::string f_rhof_alp, std::string f_rhofvf_alp);
  void init(double t);
  void run(double t, double dt, double alpha);
 private:
  FluidAgent *fagent;
  std::string f_vs_alp_str, f_rhof_alp_str, f_rhofvf_alp_str;
  int f_vs_alp_hdl, f_rhof_alp_hdl, f_rhofvf_alp_hdl;
};

class ComputeBurnPane : public Action {
 public:
  ComputeBurnPane(FluidAgent *ag, BurnAgent *bag, SolidAgent *sag, const std::string fb_mdot_alp, const std::string rhofvf_alp, double z);
  void init(double t);
  void run(double t, double dt, double alpha);
 private:
  FluidAgent *fagent;
  BurnAgent *bagent;
  SolidAgent *sagent;
  double zoom;
  int b_rb_alp_hdl, fb_rhof_alp_hdl, fb_mdot_alp_hdl;
  int fb_nf_alp_hdl, fb_rhofvf_alp_hdl;
};

// copy mesh from parent window to fluid face
class CopyBurnFromParentMesh : public Action {
 public:
  CopyBurnFromParentMesh(BurnAgent *bag, FluidAgent *fag);
  void init(double t);
  void run(double t, double dt, double alpha);
 private:
  BurnAgent *bagent;
  FluidAgent *fagent;
  std::string burn_mesh, parent_mesh;
};

// copy bflag from rocburn if not at time 0
class CopyBflagFromBurn : public Action {
 public:
  CopyBflagFromBurn(BurnAgent *bag);
  void init(double t);
  void run(double t, double dt, double alpha);
 private:
  BurnAgent *bagent;
  std::string  parent_bflag, burn_bflag;
};

// compute pconn for the whole surface, which will be used by surface propagation
class ComputePconn : public Action {
 public:
  ComputePconn(Agent *ag, std::string a_mesh, std::string a_pconn, std::string p_pmesh, int *cond=NULL);
  void init(double t);
  void run(double t, double dt, double alpha);
 private:
  Agent *agent;
  std::string a_mesh_str, a_pconn_str, p_pmesh_str;
  int a_mesh_hdl, a_pconn_hdl, p_pmesh_hdl;
  int MAP_compute_pconn;
  int PROP_initialize;
  int *cond_addr;
};

/*
class InitBurnBuffer : public Action {
 public:
  InitBurnBuffer(FluidAgent *ag, BurnAgent *bag, double z);
  void init(double t);
  void run(double t, double dt, double alpha);
 private:
  FluidAgent *fagent;
  BurnAgent *bagent;
  double zoom;
  int b_rb_alp_hdl, fb_rhof_alp_hdl, fb_mdot_alp_hdl;
};
*/

class SolidPropagateSurface_ALE : public Action {
 public:
  SolidPropagateSurface_ALE(SolidAgent *ag, const std::string p_rb, const std::string a_vbar, double z);
  void init(double t);
  void run(double t, double dt, double alpha);
 private:
  SolidAgent *sagent;
  double zoom;
  int p_rb_hdl, b_rb_hdl, p_pmesh_hdl, p_vm_hdl, p_vbar_hdl, a_vbar_hdl;
  int p_cnstr_hdl, p_pos_hdl;
  int MAP_reduce_maxabs;
  int PROP_propagate;
  int PROPCON_find_intersections, PROPCON_constrain_displacements;
};

class Reset_du_alp : public Action {
 public:
  Reset_du_alp(FluidAgent *fag);
  void init(double t);
  void run(double t, double dt, double alpha);
private:
  FluidAgent *fagent;
  int    du_alp_hdl;
};

#endif








