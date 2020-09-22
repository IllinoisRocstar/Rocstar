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

#include "Action.h"

#include "SurfDiver.h"

class RocstarAgent;
class FluidAgent;
class SolidAgent;
class BurnAgent;

class DummyAction : public Action {
public:
  DummyAction() : Action("DummyAction") {}
  void init(double t) override {}
  void run(double t, double dt, double alpha) override {}
  void finalize() override {}
};

class SetValueDouble : public Action {
public:
  SetValueDouble(const std::string &at, double val);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  double v;
  int attr_hdl{-1};
};

class SetZero : public SetValueDouble {
public:
  explicit SetZero(const std::string &at);
};

class CopyValue : public Action {
public:
  CopyValue(const std::string &from, const std::string &to, bool cond = false);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  bool condition;
  int from_hdl{-1}, to_hdl{-1};
};

class DummyPrint : public Action {
public:
  DummyPrint(BurnAgent *bag, SolidAgent *sag, FluidAgent *fag, std::string l);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  BurnAgent *bagent;
  SolidAgent *sagent;
  FluidAgent *fagent;
  std::string label;
};

/* AEG: Invokers removed.
class BCInvoker : public Action {
 public:
  BCInvoker(RocstarAgent *ag, int l = 1);
  void init(double t);
  void run(double t, double dt, double alpha);
  void finalize() override {}

 private:
  RocstarAgent *agent;
  int level;
};

class GMInvoker : public Action {
 public:
  GMInvoker(RocstarAgent *ag);
  void init(double t);
  void run(double t, double dt, double alpha);
  void finalize() override {}

 private:
  RocstarAgent *agent;
};

class BCInitInvoker : public Action {
 public:
  BCInitInvoker(RocstarAgent *ag);
  void init(double t);
  void run(double t, double dt, double alpha);
  void finalize() override {}

 private:
  RocstarAgent *agent;
};
*/

// POST_UPDATE_FLUID
class ComputeFluidLoad_ALE : public Action {
public:
  explicit ComputeFluidLoad_ALE(FluidAgent *fag, SolidAgent *sag,
                                const std::string &f_pf,
                                const std::string &fb_mdot,
                                const std::string &b_rb,
                                const std::string &f_ts);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  int traction_mode;
  FluidAgent *fagent;
  SolidAgent *sagent;
  int f_pf_hdl{-1}, f_tf_hdl{-1}, f_ts_hdl{-1};
  int fb_ts_hdl{-1}, fb_pf_hdl{-1};
  int fb_mdot_hdl{-1}, fb_rhof_alp_hdl{-1}, fb_mdot_tmp_hdl{-1};
  int b_rb_hdl{-1}, fb_nf_alp_hdl{-1}, fb_tf_hdl{-1};
};

// UPDATE_INBUF_GM_FLUID
class ComputeMeshMotion : public Action {
public:
  ComputeMeshMotion(FluidAgent *ag, const std::string &a_vm,
                    const std::string &f_du_alp, double z);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  FluidAgent *fagent;
  double zoom;
  int a_vm_hdl{-1}, f_du_alp_hdl{-1}, f_zoom_hdl{-1};
};

class ComputeFaceCenters : public Action {
public:
  ComputeFaceCenters(BurnAgent *ag, const std::string &b_nc,
                     const std::string &b_cnts, const std::string &b_nrml = "");
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  BurnAgent *agent;
  int b_nc_hdl{-1}, b_cnts_hdl{-1}, b_nrml_hdl{-1};
  int SURF_n2f{-1};
  int SURF_fn{-1};
};

class FluidPropagateSurface : public Action {
public:
  FluidPropagateSurface(FluidAgent *ag, BurnAgent *bag, const std::string &b_rb,
                        const std::string &a_vm, double z);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  FluidAgent *fagent;
  BurnAgent *bagent;
  double zoom;
  int p_rb_hdl{-1}, b_rb_hdl{-1}, fb_rb_hdl{-1}, p_cnstr_type{-1};
  int p_pmesh_hdl{-1}, p_vm_hdl{-1}, a_vm_hdl{-1};
  int p_cflag_hdl{-1}, p_pos_hdl{-1}, p_bflag_hdl{-1};
  int MAP_reduce_maxabs{-1};
  int MAP_reduce_minabs{-1};
  int PROP_set_cnstr{-1};
  int PROP_propagate{-1};
  int PROPCON_find_intersections{-1};
  int PROPCON_constrain_displacements{-1};
  int PROPCON_burnout{-1};
  int PROPCON_burnout_filter{-1};
};

class MassTransfer : public Action {
public:
  MassTransfer(FluidAgent *ag, BurnAgent *bag, const std::string &b_rhos,
               const std::string &b_rb, const std::string &f_mdot);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  FluidAgent *fagent;
  BurnAgent *bagent;
  int f_mdot_hdl{-1}, b_rhos_hdl{-1}, b_rb_hdl{-1}, fb_mdot_hdl{-1};
};

class ZoomInterface : public Action {
public:
  ZoomInterface(FluidAgent *ag, BurnAgent *bag, const std::string &fb_mdot_alp,
                double z);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  FluidAgent *fagent;
  BurnAgent *bagent;
  int b_rb_alp_hdl{-1}, fb_rhof_alp_hdl{-1}, fb_mdot_alp_hdl{-1}, rhos_hdl{-1};
  double zoom;
};

class ComputeRhofvf : public Action {
public:
  ComputeRhofvf(FluidAgent *ag, const std::string &f_vs_alp,
                const std::string &f_rhof_alp, const std::string &f_rhofvf_alp);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  FluidAgent *fagent;
  int f_vs_alp_hdl{-1}, f_rhof_alp_hdl{-1}, f_rhofvf_alp_hdl{-1};
};

class ComputeBurnPane : public Action {
public:
  ComputeBurnPane(FluidAgent *ag, BurnAgent *bag, SolidAgent *sag,
                  const std::string &fb_mdot_alp, const std::string &rhofvf_alp,
                  double z);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  FluidAgent *fagent;
  BurnAgent *bagent;
  SolidAgent *sagent;
  double zoom;
  int b_rb_alp_hdl{-1}, fb_rhof_alp_hdl{-1}, fb_mdot_alp_hdl{-1};
  int fb_nf_alp_hdl{-1}, fb_rhofvf_alp_hdl{-1};
};

// copy mesh from parent window to fluid face
class CopyBurnFromParentMesh : public Action {
public:
  CopyBurnFromParentMesh(BurnAgent *bag, FluidAgent *fag);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  BurnAgent *bagent;
  FluidAgent *fagent;
  std::string burn_mesh, parent_mesh;
};

// copy bflag from rocburn if not at time 0
class CopyBflagFromBurn : public Action {
public:
  explicit CopyBflagFromBurn(BurnAgent *bag);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  BurnAgent *bagent;
  std::string parent_bflag, burn_bflag;
};

// compute pconn for the whole surface, which will be used by surface
// propagation
class ComputePconn : public Action {
public:
  ComputePconn(RocstarAgent *ag, const std::string &a_mesh,
               const std::string &a_pconn, const std::string &p_pmesh,
               bool cond = false);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  RocstarAgent *agent;
  bool condition;
  int a_mesh_hdl{-1}, a_pconn_hdl{-1}, p_pmesh_hdl{-1};
  int MAP_compute_pconn{-1};
  int PROP_initialize{-1};
};

/*
class InitBurnBuffer : public Action {
 public:
  InitBurnBuffer(FluidAgent *ag, BurnAgent *bag, double z);
  void init(double t);
  void run(double t, double dt, double alpha);
  void finalize() override {}

 private:
  FluidAgent *fagent;
  BurnAgent *bagent;
  double zoom;
  int b_rb_alp_hdl{-1}, fb_rhof_alp_hdl{-1}, fb_mdot_alp_hdl{-1};
};
*/

class SolidPropagateSurface_ALE : public Action {
public:
  SolidPropagateSurface_ALE(SolidAgent *ag, const std::string &p_rb,
                            const std::string &a_vbar, double z);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  SolidAgent *sagent;
  double zoom;
  int p_rb_hdl{-1}, p_pmesh_hdl{-1}, p_vbar_hdl{-1}, a_vbar_hdl{-1};
  int MAP_reduce_maxabs{-1};
  int PROP_propagate{-1};
};

class Reset_du_alp : public Action {
public:
  explicit Reset_du_alp(FluidAgent *fag);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;
  void finalize() override {}

private:
  FluidAgent *fagent;
  int du_alp_hdl{-1};
};

#endif //_BASIC_ACTIONS_H_
