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
// $Id: transfer_actions.h,v 1.23 2008/12/06 08:45:22 mtcampbe Exp $

#ifndef _TRANSFER_ACTIONS_H_
#define _TRANSFER_ACTIONS_H_

#include "Action.h"

class FluidAgent;
class SolidAgent;
class BurnAgent;

// base class for mesh transfer
class InterMeshTransfer : public Action {
public:
  InterMeshTransfer(ActionDataList adl, FluidAgent *fag, SolidAgent *sag,
                    const std::string &name = "");
  void finalize() override{};

protected:
  void load_rocface(int rfc_verb);
  FluidAgent *fagent;
  SolidAgent *sagent;
};

// Fluid => Solid (pressure)
class LoadTransfer_FS : public InterMeshTransfer {
public:
  explicit LoadTransfer_FS(FluidAgent *fag, SolidAgent *sag,
                           const std::string &f_ts, const std::string &s_ts,
                           const std::string &s_pf);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;

private:
  int traction_mode;
  int size_ts{0};
  int f_tf_hdl{-1}, f_ts_hdl{-1}, s_ts_hdl{-1}, s_pf_hdl{-1};
  int f_pf_hdl{-1};
  int RFC_transfer{-1};
  int SURF_compute_face_normals{-1};
};

// Fluid => Solid (pressure)  with Burn
class LoadTransfer_FSc_ALE : public InterMeshTransfer {
public:
  explicit LoadTransfer_FSc_ALE(FluidAgent *fag, SolidAgent *sag,
                                BurnAgent *bag, const std::string &f_pf,
                                const std::string &fb_mdot,
                                const std::string &b_rb,
                                const std::string &s_ts,
                                const std::string &s_pf);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;

private:
  int traction_mode;
  int size_ts{0};
  BurnAgent *bagent;
  int f_tf_hdl{-1}, f_ts_hdl{-1}, s_ts_hdl{-1}, s_pf_hdl{-1};
  int f_pf_hdl{-1}, fb_ts_hdl{-1}, fb_pf_hdl{-1};
  int fb_mdot_hdl{-1}, fb_rhof_alp_hdl{-1}, fb_mdot_tmp_hdl{-1};
  int b_rb_hdl{-1}, fb_nf_alp_hdl{-1}, fb_tf_hdl{-1};
  int RFC_transfer{-1};
  int SURF_compute_face_normals{-1};
};

// Fluid => Solid (pressure)  with Burn : part 2 only
class LoadTransferOnly_FSc_ALE : public InterMeshTransfer {
public:
  explicit LoadTransferOnly_FSc_ALE(FluidAgent *fag, SolidAgent *sag,
                                    BurnAgent *bag, const std::string &f_ts,
                                    const std::string &s_ts,
                                    const std::string &s_pf);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;

private:
  int traction_mode;
  int size_ts{0};
  BurnAgent *bagent;
  int f_ts_hdl{-1}, s_ts_hdl{-1}, s_pf_hdl{-1};
  int RFC_transfer{-1};
  int SURF_compute_face_normals{-1};
};

// Solid => Fluid
// POST_UPDATE_SOLID solid_agent.f90
class GetDeformedMesh : public InterMeshTransfer {
public:
  explicit GetDeformedMesh(FluidAgent *fag, SolidAgent *sag,
                           const std::string &s_x, const std::string &uhat,
                           const std::string &s_y);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;

private:
  int s_x_hdl{-1}, s_uhat_hdl{-1}, s_y_hdl{-1};
};

// Solid => Fluid
// POST_UPDATE_SOLID solid_agent.f90
class GetDeformedMesh_ALE : public InterMeshTransfer {
public:
  explicit GetDeformedMesh_ALE(FluidAgent *fag, SolidAgent *sag,
                               const std::string &s_x, const std::string &uhat,
                               const std::string &s_y, double z);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;

private:
  double zoom;
  bool with_ALE{false};
  int s_x_hdl{-1}, s_uhat_hdl{-1}, s_y_hdl{-1};
  int s_rb_hdl{-1}, s_mdot_hdl{-1}, s_areas_hdl{-1}, s_rhos_hdl{-1};
  int SURF_compute_bounded_volumes{-1}, SURF_compute_face_areas{-1};
};

// Solid => Fluid
// INIT_INBUFF_FLUID fluid_agent.f90
class MeshMotionTransfer_SF : public InterMeshTransfer {
public:
  MeshMotionTransfer_SF(FluidAgent *fag, SolidAgent *sag,
                        const std::string &s_u, const std::string &f_total_disp,
                        const std::string &f_vm);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;

private:
  int s_u_hdl{-1}, f_total_disp_hdl{-1}, f_vm_hdl{-1};
  int f_nc_hdl{-1}, f_nc_t0_hdl{-1};
  int RFC_interpolate{-1};
};

// Solid => Fluid
// INIT_INBUFF_FLUID fluid_agent.f90
class DeformationVelTransfer_SF : public InterMeshTransfer {
public:
  DeformationVelTransfer_SF(FluidAgent *fag, SolidAgent *sag,
                            const std::string &s_vs, const std::string &f_vs);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;

private:
  int s_vs_hdl{-1}, f_vs_hdl{-1};
  int RFC_transfer{-1};
};

// Solid => Fluid
class MeshMotionTransferISS : public InterMeshTransfer {
public:
  explicit MeshMotionTransferISS(FluidAgent *fag, SolidAgent *sag,
                                 const std::string &s_u,
                                 const std::string &s_vs,
                                 const std::string &f_vm); // fluid agent
  void init(double t) override;
  void run(double t, double dt, double alpha) override;

private:
  int s_u_hdl{-1}, s_vs_hdl{-1}, f_vm_hdl{-1};
  int f_u_hdl{-1}, f_vs_hdl{-1};
  int f_uold_hdl{-1}, f_vsold_hdl{-1}, f_utmp_hdl{-1}, f_vtmp_hdl{-1};
  int RFC_interpolate{-1};
};

// Burn
// INIT_INBUFF_FLUID fluid_agent.f90
class TransferSolidDensity : public InterMeshTransfer {
public:
  explicit TransferSolidDensity(FluidAgent *fag, SolidAgent *sag,
                                const std::string &s_rhos,
                                const std::string &f_rhos);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;

private:
  int s_rhos_hdl{-1}, f_rhos_hdl{-1};
  int RFC_transfer{-1};
};

// Burn
// INIT_INBUFF_FLUID fluid_agent.f90
class TransferBurnRate_FS_ALE : public InterMeshTransfer {
public:
  explicit TransferBurnRate_FS_ALE(FluidAgent *fag, SolidAgent *sag,
                                   const std::string &b_rb,
                                   const std::string &s_rb);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;

private:
  int b_rb_hdl{-1}, s_rb_hdl{-1};
  int p_rb_hdl{-1}, f_mdot_tmp_hdl{-1}, fb_mdot_tmp_hdl{-1};
  int RFC_transfer{-1};
};

class MassTransfer_SF_ALE : public InterMeshTransfer {
public:
  MassTransfer_SF_ALE(FluidAgent *ag, SolidAgent *sag, BurnAgent *bag,
                      const std::string &f_mdot);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;

private:
  BurnAgent *bagent;
  int f_mdot_hdl{-1};
  int b_rhos_hdl{-1}, b_rb_hdl{-1}, fb_mdot_hdl{-1};
  int f_total_disp_hdl{-1}, f_nc_hdl{-1}, f_nc_t0_hdl{-1}, s_mdot_hdl{-1};
  int RFC_transfer{-1};
};

class TemperatureTransfer_SF : public InterMeshTransfer {
public:
  TemperatureTransfer_SF(SolidAgent *sag, FluidAgent *ag,
                         const std::string &s_Tf,
                         const std::string &fb_Tflm_alp,
                         const std::string &fn_Tb);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;

private:
  int s_Ts_hdl{-1}, fb_Tflm_hdl{-1}, fn_Tb_hdl{-1};
  int f_Ts_hdl{-1}, bcflag_hdl{-1}; // temporary
  int RFC_transfer{-1};
};

class HeatTransfer_FS : public InterMeshTransfer {
public:
  HeatTransfer_FS(FluidAgent *ag, SolidAgent *sag, BurnAgent *bag,
                  const std::string &f_qc, const std::string &f_qr,
                  const std::string &fb_qev, const std::string &s_qs);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;

private:
  BurnAgent *bagent;
  int f_qc_hdl{-1}, f_qr_hdl{-1}, b_qev_hdl{-1}, s_qs_hdl{-1};
  int s_qc_hdl{-1}, s_qr_hdl{-1}, f_qev_hdl{-1}, s_qev_hdl{-1};
  int RFC_transfer{-1};
};

class RemeshInit : public InterMeshTransfer {
public:
  explicit RemeshInit(FluidAgent *fag, SolidAgent *sag, const std::string &s_u,
                      const std::string &f_total_disp, const std::string &f_nc,
                      const std::string &f_nc_t0);
  void init(double t) override;
  void run(double t, double dt, double alpha) override;

private:
  int s_u_hdl{-1}, f_total_disp_hdl{-1};
  int f_nc_hdl{-1}, f_nc_t0_hdl{-1};
  int RFC_interpolate{-1};
};

#endif //_TRANSFER_ACTIONS_H_
