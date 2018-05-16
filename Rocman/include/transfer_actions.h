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

#include "rocman.h"

class FluidAgent;
class SolidAgent;
class BurnAgent;

// base class for mesh transfer
class InterMeshTransfer : public Action {
  public:
    explicit InterMeshTransfer( FluidAgent *fag, SolidAgent *sag, char *name=NULL);
  protected:
    void load_rocface(const RocmanControl_parameters *param);
    FluidAgent *fagent;
    SolidAgent *sagent;
};

// Fluid => Solid (pressure)
class LoadTransfer_FS: public InterMeshTransfer {
  public:								
    explicit LoadTransfer_FS( FluidAgent *fag, SolidAgent *sag, const std::string f_ts, const std::string s_ts, const std::string s_pf);
    void init(double t);
    void run(double t, double dt, double alpha);
  private:
    int traction_mode;
    int size_ts;
    std::string f_ts_str, s_ts_str, s_pf_str;
    int f_tf_hdl, f_ts_hdl, s_ts_hdl, s_pf_hdl;
    int f_pf_hdl; 
    int RFC_transfer;
    int SURF_compute_face_normals;
};

// Fluid => Solid (pressure)  with Burn
class LoadTransfer_FSc_ALE: public InterMeshTransfer {
  public:								
    explicit LoadTransfer_FSc_ALE( FluidAgent *fag, SolidAgent *sag, BurnAgent *bag, const std::string f_pf, const std::string fb_mdot, const std::string b_rb, const std::string s_ts, const std::string s_pf);
    void init(double t);
    void run(double t, double dt, double alpha);
  private:
    int traction_mode;
    int size_ts;
    BurnAgent *bagent;
    int f_tf_hdl, f_ts_hdl, s_ts_hdl, s_pf_hdl;
    int f_pf_hdl, fb_ts_hdl, fb_pf_hdl;
    int fb_mdot_hdl, fb_rhof_alp_hdl, fb_mdot_tmp_hdl;
    int b_rb_hdl, fb_nf_alp_hdl, fb_tf_hdl;
    int RFC_transfer;
    int SURF_compute_face_normals;
};

// Fluid => Solid (pressure)  with Burn : part 2 only
class LoadTransferOnly_FSc_ALE: public InterMeshTransfer {
  public:								
    explicit LoadTransferOnly_FSc_ALE( FluidAgent *fag, SolidAgent *sag, BurnAgent *bag, const std::string f_ts, const std::string s_ts, const std::string s_pf);
    void init(double t);
    void run(double t, double dt, double alpha);
  private:
    int traction_mode;
    int size_ts;
    BurnAgent *bagent;
    int f_ts_hdl, s_ts_hdl, s_pf_hdl;
    int RFC_transfer;
    int SURF_compute_face_normals;
};

// Solid => Fluid
// POST_UPDATE_SOLID solid_agent.f90
class GetDeformedMesh: public InterMeshTransfer {
public:
  explicit GetDeformedMesh( FluidAgent *fag, SolidAgent *sag, const std::string s_x, const std::string uhat, const std::string s_y);
  void init(double t);
  void run(double t, double dt, double alpha);
private:
  std::string s_x_str, s_uhat_str, s_y_str;
  int s_x_hdl, s_uhat_hdl, s_y_hdl;
};

// Solid => Fluid
// POST_UPDATE_SOLID solid_agent.f90
class GetDeformedMesh_ALE: public InterMeshTransfer {
public:
  explicit GetDeformedMesh_ALE( FluidAgent *fag, SolidAgent *sag, const std::string s_x, const std::string uhat, const std::string s_y, double z);
  void init(double t);
  void run(double t, double dt, double alpha);
private:
  std::string s_x_str, s_uhat_str, s_y_str;
  int s_x_hdl, s_uhat_hdl, s_y_hdl;
  int s_rb_hdl, s_mdot_hdl, s_areas_hdl, s_rhos_hdl;
  int withALE;
  double zoom;
  int SURF_compute_bounded_volumes, SURF_compute_face_areas;
};

// Solid => Fluid
// INIT_INBUFF_FLUID fluid_agent.f90
class MeshMotionTransfer_SF: public InterMeshTransfer {
public:
  explicit MeshMotionTransfer_SF( FluidAgent *fag, SolidAgent *sag, const std::string s_u, const std::string f_total_disp, const std::string f_vm);
  void init(double t);
  void run(double t, double dt, double alpha);
private:
  int s_u_hdl, f_total_disp_hdl, f_vm_hdl;
  int f_nc_hdl, f_nc_t0_hdl;
  int RFC_interpolate;
};

// Solid => Fluid
// INIT_INBUFF_FLUID fluid_agent.f90
class DeformationVelTransfer_SF: public InterMeshTransfer {
public:
  explicit DeformationVelTransfer_SF( FluidAgent *fag, SolidAgent *sag, const std::string s_vs, const std::string f_vs);	            // fluid agent
  void init(double t);
  void run(double t, double dt, double alpha);
private:
  std::string s_vs_str, f_vs_str;
  int s_vs_hdl, f_vs_hdl;
  int RFC_transfer;
};

// Solid => Fluid
class MeshMotionTransferISS: public InterMeshTransfer {
public:
  explicit MeshMotionTransferISS( FluidAgent *fag, SolidAgent *sag, const std::string s_u, const std::string s_vs, const std::string f_vm);	// fluid agent
  void init(double t);
  void run(double t, double dt, double alpha);
private:
  std::string s_u_str, s_vs_str, f_vm_str;
  int s_u_hdl, s_vs_hdl, f_vm_hdl;
  int f_u_hdl, f_vs_hdl;
  int f_uold_hdl, f_vsold_hdl, f_utmp_hdl, f_vtmp_hdl;
  int RFC_interpolate;
};

// Burn
// INIT_INBUFF_FLUID fluid_agent.f90
class TransferSolidDensity: public InterMeshTransfer {
public:
  explicit TransferSolidDensity( FluidAgent *fag, SolidAgent *sag, const std::string s_rhos, const std::string f_rhos);
  void init(double t);
  void run(double t, double dt, double alpha);
private:
  std::string s_rhos_str, f_rhos_str;
  int s_rhos_hdl, f_rhos_hdl;
  int RFC_transfer;
};

// Burn
// INIT_INBUFF_FLUID fluid_agent.f90
class TransferBurnRate_FS_ALE: public InterMeshTransfer {
public:
  explicit TransferBurnRate_FS_ALE( FluidAgent *fag, SolidAgent *sag, const std::string b_rb, const std::string s_rb);
  void init(double t);
  void run(double t, double dt, double alpha);
private:
  std::string b_rb_str, s_rb_str;
  int b_rb_hdl, p_rb_hdl, s_rb_hdl, f_mdot_tmp_hdl, fb_mdot_tmp_hdl;
  int RFC_transfer;
};

class MassTransfer_SF_ALE : public InterMeshTransfer {
 public:
  MassTransfer_SF_ALE(FluidAgent *ag, SolidAgent *sag, BurnAgent *bag, const std::string f_mdot);
  void init(double t);
  void run(double t, double dt, double alpha);
 private:
  BurnAgent *bagent;
  int f_mdot_hdl, b_rhos_hdl, b_rb_hdl, fb_mdot_hdl;
  int f_total_disp_hdl, f_nc_hdl, f_nc_t0_hdl, s_mdot_hdl;
  int RFC_transfer;
};

class TemperatureTransfer_SF : public InterMeshTransfer {
 public:
  TemperatureTransfer_SF(SolidAgent *sag, FluidAgent *ag, const std::string s_Tf, const std::string fb_Tflm_alp, const std::string fn_Tb);
  void init(double t);
  void run(double t, double dt, double alpha);
 private:
  int s_Ts_hdl, fb_Tflm_hdl, fn_Tb_hdl;
  int f_Ts_hdl, bcflag_hdl;              // temporary
  int RFC_interpolate;
  int RFC_transfer;
};

class HeatTransfer_FS : public InterMeshTransfer {
 public:
  HeatTransfer_FS(FluidAgent *ag, SolidAgent *sag, BurnAgent *bag, const std::string f_qc, const std::string f_qr, const std::string fb_qev, const std::string s_qs);
  void init(double t);
  void run(double t, double dt, double alpha);
 private:
  BurnAgent *bagent;
  int f_qc_hdl, f_qr_hdl, b_qev_hdl, s_qs_hdl;
  int s_qc_hdl, s_qr_hdl, f_qev_hdl, s_qev_hdl;
  int RFC_interpolate;
  int RFC_transfer;
};

class RemeshInit : public InterMeshTransfer {
public:
  explicit RemeshInit( FluidAgent *fag, SolidAgent *sag, const std::string s_u, const std::string f_total_disp, const std::string f_nc, const std::string f_nc_t0);
  void init(double t);
  void run(double t, double dt, double alpha);
private:
  int s_u_hdl, f_total_disp_hdl ;
  int f_nc_hdl, f_nc_t0_hdl;
  int RFC_interpolate;
};

#define DECLARE_NEW_ACTION( ActionName)					\
  class ActionName : public InterMeshTransfer {					\
  public:								\
    explicit ActionName( FluidAgent *fag, SolidAgent *sag, const char *at[], int *is=NULL, void *p=NULL);	\
    void init(double t);						\
    void run(double t, double dt, double alpha);			\
  };

// This action set the dataitem to 0s.
DECLARE_NEW_ACTION( LoadTransfer_SF);

#endif






