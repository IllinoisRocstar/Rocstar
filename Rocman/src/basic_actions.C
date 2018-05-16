/* *******************************************************************
 * Rocstar Simulation Suite                                          *
 *                                                                   *
 * License: See LICENSE file in top level of distribution package or *
 * http://opensource.org/licenses/NCSA                               *
 *********************************************************************/
/* *******************************************************************
 * Rocstar Simulation Suite                                          *
 * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
 *                                                                   *
 * Illinois Rocstar LLC                                              *
 * Champaign, IL                                                     *
 * www.illinoisrocstar.com                                           *
 * sales@illinoisrocstar.com                                         *
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
// $Id: basic_actions.C,v 1.82 2009/09/14 14:18:27 mtcampbe Exp $

#include "basic_actions.h"
#include "FluidAgent.h"
#include "SolidAgent.h"
#include "BurnAgent.h"

COM_EXTERN_MODULE(SurfUtil);
COM_EXTERN_MODULE(Rocprop);
COM_EXTERN_MODULE(SurfMap);
COM_EXTERN_MODULE(Rocon);

SetValueDouble::SetValueDouble( const std::string at, const double val) :
  Action( 1, &at, NULL, NULL, (char *)"SetValueDouble"), v(val)
{
  // Define I/O of the arguments.
  int io[] = {OUT};
  set_io( 1, io);
}

// Obtain the dataitem handles
void SetValueDouble::init( double t) {
  attr_hdl = get_dataitem_handle( 0);
}

void SetValueDouble::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: SetValueDouble::run() %s with %e.\n", attr[0], v));
  COM_call_function( RocBlas::copy_scalar, &v, &attr_hdl);
}

SetZero::SetZero(const std::string at) : SetValueDouble( at, 0)
{ action_name = (char *)"SetZero"; }

CopyValue::CopyValue( const std::string from, const std::string to, int *cond) :
  Action( 0, (const char**)NULL, NULL, NULL, (char *)"CopyValue"), condition(cond)
{
  // Define I/O of the arguments.
  int io[] = {IN, OUT};
  set_io( 2, io);

  std::string atts[2];
  atts[0] = from;
  atts[1] = to;
  set_attr(2, atts);
}

// Obtain the dataitem handles
void CopyValue::init( double t) {
  from_hdl = get_dataitem_handle( 0);
  to_hdl = get_dataitem_handle( 1);
}

void CopyValue::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: CopyValue::run() %s => %s.\n", attr[0], attr[1]));
  if (!condition || *condition)
  COM_call_function( RocBlas::copy, &from_hdl, &to_hdl);
}

DummyPrint::DummyPrint( BurnAgent *bag, SolidAgent *sag, FluidAgent *fag, const std::string l): 
  Action( 0, (const char**)NULL, NULL, NULL, (char *)"DummyPrint"),
      bagent(bag), sagent(sag), fagent(fag), label(l)
{
  set_io( 0, NULL);

  set_attr(0, (const char**)NULL);
}

// Obtain the dataitem handles
void DummyPrint::init( double t) {
  MAN_DEBUG(3, ("Rocstar: DummyPrint::init() with %s.\n", label.c_str()));
}

void DummyPrint::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: DummyPrint::run() with %s.\n", label.c_str()));

  //debug_print(sagent->solidBuf+".x", 202, 1, "DummyPrint");
  //debug_print(bagent->iburn_ng+".Tflm_alp", 102, 0);
  //debug_print(bagent->iburn_ng+".rb", 102, 0);
  //debug_print(bagent->iburn_ng+".rb_alp", 102, 0);
  //debug_print(bagent->iburn_ng+".pf_alp", 102, 0);
  //debug_print(fagent->fluidBufB+".mdot", 2, 0);
  //debug_print(fagent->fluidBufB+".mdot_old", 2, 0);
  //debug_print(fagent->fluidBufB+".mdot_alp", 2, 0);
  //debug_print(fagent->fluidBufNG+".mdot_alp", 2, 0);
}

BCInvoker::BCInvoker( Agent *ag, int l): 
  Action( 0, (const char**)NULL, NULL, NULL, (char *)"BCInvoker"),
      agent(ag), level(l)
{
  set_io( 0, NULL);
  set_attr(0, (const char**)NULL);
}

void BCInvoker::init( double t) {
    //  init bcactions
  agent->init_bcactions(t);
}

void BCInvoker::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: BCInvoker::run() agent: %s level: %d.\n", agent->get_agent_name().c_str(), level));

//RAF pass alpha, not time!!
//  agent->obtain_bc(&t, &level);
  agent->obtain_bc(&alpha, &level);

  MAN_DEBUG(3, ("Rocstar: BCInvoker::run() agent: %s level: %d DONE.\n", agent->get_agent_name().c_str(), level));
}

GMInvoker::GMInvoker( Agent *ag): 
  Action( 0, (const char**)NULL, NULL, NULL, (char *)"GMInvoker"),
      agent(ag)
{
  set_io( 0, NULL);
  set_attr(0, (const char**)NULL);
}

void GMInvoker::init( double t) {
    //  init bcactions
  agent->init_gmactions(t);
}

void GMInvoker::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: GMInvoker::run() agent: %s with alpha=%e.\n", agent->get_agent_name().c_str(), alpha));

  agent->obtain_gm(&alpha);

  MAN_DEBUG(3, ("Rocstar: GMInvoker::run() agent: %s with alpha=%e DONE.\n", agent->get_agent_name().c_str(), alpha));
}

BCInitInvoker::BCInitInvoker( Agent *ag): 
  Action( 0, (const char**)NULL, NULL, NULL, (char *)"BCInitInvoker"),
      agent(ag)
{
  set_io( 0, NULL);
  set_attr(0, (const char**)NULL);
}

void BCInitInvoker::init( double t) {
  MAN_DEBUG(3, ("Rocstar: BCInitInvoker::init() agent: %s .\n", agent->get_agent_name().c_str()));
    //  init bcactions
  agent->init_bcinitaction(t);
  MAN_DEBUG(3, ("Rocstar: BCInitInvoker::init() agent: %s DONE.\n", agent->get_agent_name().c_str()));
}

void BCInitInvoker::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: BCInitInvoker::run() agent: %s t:%e dt:%e alpha:%e.\n", agent->get_agent_name().c_str(), t, dt, alpha));

  agent->run_bcinitaction(t, dt);

  MAN_DEBUG(3, ("Rocstar: BCInitInvoker::run() agent: %s DONE.\n", agent->get_agent_name().c_str()));
}

// POST_UPDATE_FLUID
// first part of LoadTransfer_FSc_ALE: output f_ts
ComputeFluidLoad_ALE::ComputeFluidLoad_ALE( FluidAgent *fag, SolidAgent *sag, 
                         const std::string f_pf, const std::string fb_mdot, 
                         const std::string b_rb, 
                         const std::string f_ts) : 
  Action(0, (const char**)NULL, NULL, NULL, (char *)"ComputeFluidLoad_ALE"),
    fagent(fag), sagent(sag) 
{
  int io[] = {IN, IN, IN, OUT};
  set_io( 4, io); 
  
  std::string atts[4];
  atts[0] = f_pf;
  atts[1] = fb_mdot;
  atts[2] = b_rb;
  atts[3] = f_ts;
  set_attr(4, atts);

  traction_mode = fagent->get_coupling()->get_rocmancontrol_param()->traction_mode;

  // create_dataitem
  // for SolidAgent ??????????????????
//  if (traction_mode == NO_SHEER && size_ts == 3) {
    sagent->register_new_dataitem( sagent->solidBufBase, ".pf", 'e', COM_DOUBLE, 1, "Pa");
//  }

  // for FluidAgent
  if (traction_mode == NO_SHEER)  {
    fagent->register_clone_dataitem( 0, fagent->ifluid_i, ".ts", fagent->get_surface_window(), ".pf");
  }
  else
    fagent->register_clone_dataitem( 0, fagent->ifluid_i, ".ts", fagent->get_surface_window(), ".tf");
}

// Create buffer data
void ComputeFluidLoad_ALE::init( double t) 
{
  if (traction_mode == NO_SHEER) {
    f_pf_hdl = get_dataitem_handle_const( 0);
    fb_pf_hdl = COM_get_dataitem_handle( fagent->fluidBufB+".pf");
  }
  else {
    f_pf_hdl = -1;
    fb_pf_hdl = -1;
  }
  fb_mdot_hdl = get_dataitem_handle( 1);
  b_rb_hdl = get_dataitem_handle( 2);
  f_ts_hdl = get_dataitem_handle( 3);

  if (traction_mode != NO_SHEER) {
    f_tf_hdl = COM_get_dataitem_handle_const( fagent->fluidBufNG+".tf");
    fb_tf_hdl = COM_get_dataitem_handle_const( fagent->fluidBufB+".tf");
    fb_nf_alp_hdl = COM_get_dataitem_handle( fagent->fluidBufB+".nf_alp");
  }
  else {
    f_tf_hdl = -1;
    fb_tf_hdl = -1;
    fb_nf_alp_hdl = -1;
  }
  fb_ts_hdl = COM_get_dataitem_handle( fagent->fluidBufB+".ts");
  fb_rhof_alp_hdl = COM_get_dataitem_handle ( fagent->fluidBufB+".rhof_alp");
  fb_mdot_tmp_hdl = COM_get_dataitem_handle ( fagent->fluidBufB+".mdot_tmp");

  MAN_DEBUG(3, ("ComputeFluidLoad_ALE::init() called - traction_mode: %d .\n", traction_mode));
}

void ComputeFluidLoad_ALE::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: ComputeFluidLoad_ALE::run() with t:%e dt:%e.\n", t, dt));
  // ts = tf + ( mdot*Vs

  // part 1  (POST_UPDATE_FLUID in fluid_agent.f90)
  if (traction_mode != NO_SHEER) {
    COM_call_function( RocBlas::copy, &f_tf_hdl, &f_ts_hdl);
                                                                                
      // Compute tfmts = t_s-t_f (see developers guide).
      // Here fb_mdot_tmp is used as buffer space.
    COM_call_function( RocBlas::div, &fb_mdot_hdl, &fb_rhof_alp_hdl, &fb_mdot_tmp_hdl);
                                                                                
    if ( sagent->withALE) {
      COM_call_function( RocBlas::sub, &fb_mdot_tmp_hdl, &b_rb_hdl, &fb_mdot_tmp_hdl);
    }      
                                                                                
    COM_call_function( RocBlas::mul, &fb_mdot_tmp_hdl, &fb_mdot_hdl, &fb_mdot_tmp_hdl);
    COM_call_function( RocBlas::mul, &fb_nf_alp_hdl, &fb_mdot_tmp_hdl, &fb_ts_hdl);
                                                                                
       // Compute ts = tf - (tf-ts)
    COM_call_function( RocBlas::sub, &fb_tf_hdl, &fb_ts_hdl, &fb_ts_hdl);
  }
  else {   // ts is a scalar

    COM_call_function( RocBlas::copy, &f_pf_hdl, &f_ts_hdl);

    //    Compute tsmpf = p_f-t_s (see developers guide).
    //    Here fb_mdot_tmp is used as buffer space.
    COM_call_function( RocBlas::div, &fb_mdot_hdl, &fb_rhof_alp_hdl, &fb_mdot_tmp_hdl);

    if ( sagent->withALE) {
        COM_call_function( RocBlas::sub, &fb_mdot_tmp_hdl, &b_rb_hdl, &fb_mdot_tmp_hdl);
    }

    COM_call_function( RocBlas::mul, &fb_mdot_tmp_hdl, &fb_mdot_hdl, &fb_ts_hdl);

     //    Compute ts = pf - (pf-ts)
    COM_call_function( RocBlas::sub, &fb_pf_hdl, &fb_ts_hdl, &fb_ts_hdl);

     // Subtract from P_ambient if not zeros 
    double P_ambient = fagent->get_coupling()->get_rocmancontrol_param()->P_ambient;
    if ( P_ambient != 0.0) {
       COM_call_function( RocBlas::sub_scalar, &f_ts_hdl, &P_ambient, &f_ts_hdl);
    }
  }

  // debug_print(fagent->fluidBufNG+".ts", 102, 0, "LOAD");
}

// UPDATE_INBUFF_GM_FLUID in fluid_agent.f90
ComputeMeshMotion::ComputeMeshMotion(FluidAgent *ag, 
                const std::string a_vm, const std::string f_du_alp, double z): 
  Action(0, (const char**)NULL, NULL, NULL, (char *)"ComputeMeshMotion"),
                fagent(ag), zoom(z)
{
  int io[] = {IN, OUT};
  set_io( 2, io); 
  
  std::string atts[2];
  atts[0] = a_vm;
  atts[1] = f_du_alp;
  set_attr(2, atts);

  fagent->register_new_dataitem( fagent->propBufAll, ".vm", 'n', COM_DOUBLE, 3, "m/s");
}

void ComputeMeshMotion::init( double t) {
  a_vm_hdl = get_dataitem_handle(0);
  f_du_alp_hdl = get_dataitem_handle(fagent->get_surface_window()+attr[1]);

    // can not be set earlier because surface window was unkown
  //f_du_alp_hdl = COM_get_dataitem_handle(fagent->get_surface_window()+".du_alp");
  f_zoom_hdl = COM_get_dataitem_handle(fagent->fluidBufNG+".zoomFact");  // used internally
}

// gm callback
void ComputeMeshMotion::run( double t_dummy, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: ComputeMeshMotion::run() with dt:%e alpha:%e.\n", dt, alpha));
  
  COM_assertion_msg(alpha>=0.0, "ERROR: ComputeMeshMotion called with invalid alpha!");

//  debug_print(fagent->get_surface_window()+".du_alp", 102, 0);

  double time = alpha * dt;
  COM_call_function( RocBlas::mul_scalar, &a_vm_hdl, &time, &f_du_alp_hdl);
  
//  debug_print(fagent->propBufAll+".vm", 202, 1);
//  debug_print(fagent->get_surface_window()+".du_alp", 102, 0);

  if (f_zoom_hdl>0) {
    COM_call_function( RocBlas::copy_scalar, &zoom, &f_zoom_hdl);
  }
}

static inline void load_rocsurf()
{
  if (COM_get_window_handle("SURF") <= 0)
    COM_LOAD_MODULE_STATIC_DYNAMIC( SurfUtil, "SURF");
}

static inline void load_rocon(const RocmanControl_parameters *param, MPI_Comm comm)
{
  
  if(param->PROPCON_enabled && (COM_get_window_handle("PROPCON") <= 0)){
    //    std::cout << "Rocstar: LOADING ROCON" << std::endl;
    COM_LOAD_MODULE_STATIC_DYNAMIC(Rocon,"PROPCON");
    int PROPCON_init_hndl = COM_get_function_handle("PROPCON.init_from_file");
    std::string rocon_config_file("Rocon/RoconControl.txt");
    COM_call_function(PROPCON_init_hndl,rocon_config_file.c_str(),&(param->PROPCON_ndiv));
    //    std::cout << "Rocstar: LOADED ROCON" << std::endl;
  }
}

static inline void load_rocprop(const RocmanControl_parameters *param, MPI_Comm comm)
{
  if (COM_get_window_handle("PROP") <= 0) {
    MAN_DEBUG(3, ("Rocstar: Loading module RocProp...\n"));
    COM_LOAD_MODULE_STATIC_DYNAMIC( Rocprop, "PROP");

    int PROP_set_option = COM_get_function_handle( "PROP.set_option");
    if (param->PROP_fom)
      COM_call_function( PROP_set_option, "method", "fo");
    else
      COM_call_function( PROP_set_option, "method", "mp");
    char val[256];
    sprintf(val, "%d", param->PROP_rediter);
    COM_call_function( PROP_set_option, "rediter", val);
    sprintf(val, "%f", param->PROP_fangle);
    COM_call_function( PROP_set_option, "fangle", val);

    int PROP_initialize = COM_get_function_handle("PROP.initialize");
    int PROP_propagate = COM_get_function_handle("PROP.propagate");
    COM_set_profiling_barrier( PROP_initialize, comm);
    COM_set_profiling_barrier( PROP_propagate, comm);
  }
}

static inline void load_rocmap()
{
  if (COM_get_window_handle("MAP") <= 0) {
    MAN_DEBUG(3, ("Rocstar: Loading module RocMap...\n"));
    COM_LOAD_MODULE_STATIC_DYNAMIC( SurfMap, "MAP");
  }
}

ComputeFaceCenters::ComputeFaceCenters(BurnAgent *ag, const std::string b_nc, const std::string b_cnts,
				       const std::string b_nrml): 
  Action(0, (const char**)NULL, NULL, NULL, (char *)"ComputeFaceCenters"), agent(ag)
{
  int io[] = {IN, OUT};
  set_io( 2, io); 
  
  std::string atts[3];
  atts[0] = b_nc;
  atts[1] = b_cnts;
  atts[2] = b_nrml;
  set_attr(3, atts);
}

void ComputeFaceCenters::init( double t) {
/*
  b_nc_hdl = get_dataitem_handle(0);
  if (COM_get_dataitem_handle(agent->get_surface_window()+".centers") > 0)
    b_cnts_hdl = get_dataitem_handle(1);
  else
    b_cnts_hdl = -1;

  load_rocsurf();
  SURF_n2f = COM_get_function_handle( "SURF.interpolate_to_centers");
*/
}

void ComputeFaceCenters::run( double t, double dt, double alpha_dummy) {
  
  // ???????????
/*
  double a = 0.0;
  int level = 1;
  agent->obtain_bc(&a, &level);
*/

  b_nc_hdl = get_dataitem_handle(attr[0]);
  b_cnts_hdl = -1;
  if (COM_get_dataitem_handle(agent->get_surface_window()+".centers") > 0){
    b_cnts_hdl = get_dataitem_handle(1);
  }
  b_nrml_hdl = -1;
  if(COM_get_dataitem_handle(agent->get_surface_window()+".normals") > 0){
    b_nrml_hdl = get_dataitem_handle(2);
  }
  
  // ?????????????
//  double zero = 0.0;
//  agent->obtain_bc(&zero);

  if (b_cnts_hdl >= 0 || b_nrml_hdl >= 0) {
    load_rocsurf();
    MAN_DEBUG(3, ("Rocstar: ComputeFaceCenters::run() with t:%e dt:%e b_cnts_hdl:%d b_nrml_hdl:%d.\n", t, dt, b_cnts_hdl, b_nrml_hdl));
    if(b_cnts_hdl >= 0){
      SURF_n2f = COM_get_function_handle( "SURF.interpolate_to_centers");
      COM_call_function( SURF_n2f, &b_nc_hdl, &b_cnts_hdl);
    }
    if(b_nrml_hdl >= 0){
      SURF_fn  = COM_get_function_handle( "SURF.compute_element_normals");
      COM_call_function( SURF_fn, &b_nc_hdl,  &b_nrml_hdl);
    }
  }
}

// zoom > 0 effective
FluidPropagateSurface::FluidPropagateSurface(FluidAgent *fag, BurnAgent *bag, 
             const std::string b_rb,
             const std::string a_vm,
             double z): 
  Action(0, (const char**)NULL, NULL, NULL, (char *)"FluidPropagateSurface"), fagent(fag), bagent(bag), zoom(z)
{
  int io[] = {IN, OUT};
  set_io( 2, io); 
  
  std::string atts[2];
  atts[0] = b_rb;
  atts[1] = a_vm;
  set_attr(2, atts);

  fagent->register_new_dataitem( fagent->propBufAll, ".vm", 'n', COM_DOUBLE, 3, "m/s");
}

void FluidPropagateSurface::init( double t) {
  b_rb_hdl = get_dataitem_handle( 0);
  a_vm_hdl = get_dataitem_handle( 1);

  std::string propBuf = fagent->propBufAll;
  int PROP_fom = fagent->get_coupling()->get_rocmancontrol_param()->PROP_fom;
  if (!PROP_fom) propBuf = fagent->propBuf;
  p_rb_hdl     = COM_get_dataitem_handle( propBuf+".rb");
  p_pmesh_hdl  = COM_get_dataitem_handle( propBuf+".pmesh");
  p_vm_hdl     = COM_get_dataitem_handle( propBuf+".vm");
  p_bflag_hdl  = COM_get_dataitem_handle( propBuf+".bflag");
  p_cnstr_type = COM_get_dataitem_handle( propBuf+".cnstr_type");
  p_cflag_hdl  = COM_get_dataitem_handle( propBuf+".cflag");
  p_pos_hdl    = COM_get_dataitem_handle( propBuf+".positions");
  fb_rb_hdl    = COM_get_dataitem_handle( fagent->fluidBufB+".rb");

  load_rocprop(fagent->get_coupling()->get_rocmancontrol_param(), fagent->get_communicator());
  load_rocon(fagent->get_coupling()->get_rocmancontrol_param(),fagent->get_communicator());
  PROP_set_cnstr = COM_get_function_handle("PROP.set_constraints");
  PROP_propagate = COM_get_function_handle("PROP.propagate");
  if(COM_get_window_handle("PROPCON") >= 0) {
    PROPCON_find_intersections = COM_get_function_handle("PROPCON.find_intersections");
    PROPCON_constrain_displacements = COM_get_function_handle("PROPCON.constrain_displacements");
    PROPCON_burnout = COM_get_function_handle("PROPCON.burnout");
    PROPCON_burnout_filter = COM_get_function_handle("PROPCON.burnout_filter");
  }
  else {
    PROPCON_find_intersections = -1;
    PROPCON_constrain_displacements = -1;
    PROPCON_burnout = -1;
  }
  load_rocmap();
  MAP_reduce_maxabs = COM_get_function_handle("MAP.reduce_maxabs_on_shared_nodes");
  MAP_reduce_minabs = COM_get_function_handle("MAP.reduce_minabs_on_shared_nodes");
}

void FluidPropagateSurface::run( double t, double dt, double alpha_dummy) 
{
  if (zoom <= 0.0 ||  fagent->get_coupling()->initial_start() ) return;

  int rank = 0;
  MPI_Comm mycomm = fagent->get_communicator();
  MPI_Comm_rank(mycomm,&rank);
  if(rank == 0)
    MAN_DEBUG(2, ("Rocstar: FluidPropagateSurface::run() with t:%e dt:%e zoom:%f p_cnstr_type:%d.\n", t, dt, zoom, p_cnstr_type));

  double zero = 0.0;
  COM_call_function( RocBlas::copy_scalar, &zero, &p_rb_hdl);

  COM_call_function( RocBlas::copy, &b_rb_hdl, &fb_rb_hdl);

  // Burnout strategy:
  //
  // integer array (bflag) lives on interacting boundary 
  // faces.  bflag indicates whether the boundary face
  // is burning/reacting/injecting, etc.  In the current
  // implementation, bflag can be 1 or 0 for burning or
  // not-burning.
  //
  // integer array (cflag) lives on boundary nodes. 
  // cflag indicates the status of one or more geometrical
  // or other imposed constraints on the deformation or
  // solution space of the domain.
  // In the current implementation cflag = 1 for nodes
  // which have been constrained from further motion by
  // the physical boundary.  This happens when the domain
  // boundary and the physical boundary come in contact.
  //
  // Based on the value of bflag, certain models may
  // need to be used to determine physical quantities
  // on the boundary.   In the current implementation,
  // the burnrate solver (Rocburn) can determine whether
  // the propellant material should ignite, and if so
  // it sets bflag to 1 and sets a burning rate
  // for each burning face based on (aP^n), where
  // [P = fluid pressure], and the temperature of the
  // reacting surface.
  //
  //
  // This function sets bflag to 0 for every burned out
  // face.  Burned out faces are determined by faces
  // in which every node has cflag=1
  // BurnOut(cflag,bflag);
  if(PROPCON_burnout >= 0){
    //    std::cout << "Rocstar> Burning out surface elements with Rocon." 
    //	      << std::endl;
    double data2 = 20000*dt;
    COM_call_function( PROP_propagate, &p_pmesh_hdl, &p_rb_hdl, &data2, &p_vm_hdl);
    MPI_Barrier(mycomm);
    COM_call_function(PROPCON_find_intersections,&p_pmesh_hdl,&p_vm_hdl,&p_pos_hdl,&p_cflag_hdl);
    COM_call_function( MAP_reduce_maxabs, &p_cflag_hdl);
    COM_call_function(PROPCON_burnout,&p_pmesh_hdl,&p_cflag_hdl,&p_bflag_hdl);
    COM_call_function( MAP_reduce_maxabs, &p_cflag_hdl); 
    COM_call_function(PROPCON_burnout,&p_pmesh_hdl,&p_cflag_hdl,&p_bflag_hdl);
    //    COM_call_function( MAP_reduce_maxabs,
  }

  //
  // This function sets the burning rate (rb) to 0.
  // The propagation code (Rocprop) sees rb as the 
  // (speed function).  If the speed function is 0
  // then Rocprop will not propagate the face.  This will
  // help Rocprop and Rocon to work together in keeping
  // the propagation constrained to the constraint surface.
  // p_rb_hdl = burning_rate
  if(PROPCON_burnout_filter >= 0)
    COM_call_function(PROPCON_burnout_filter,&p_bflag_hdl,&p_rb_hdl);
  //
  
  // Now original code to propagate the mesh:
  if ( p_cnstr_type > 0) 
    COM_call_function( PROP_set_cnstr, &p_cnstr_type);

  double data = dt * zoom;

  MPI_Barrier(mycomm);
  if(!rank && man_verbose > 2)
    std::cout << "Rocstar: Calling Rocprop"  << std::endl;
  // p_rb_hdl is the "speed function"
  COM_call_function( PROP_propagate, &p_pmesh_hdl, &p_rb_hdl, &data, &p_vm_hdl);
  MPI_Barrier(mycomm);
  if(!rank && man_verbose > 2)
    std::cout << "Rocstar: Rocprop done." << std::endl;
  // Rocon determines if the domain boundary and the physical boundary have
  // contacted and sets cflag (p_cflag_hdl) if so, and also constrains the 
  // nodes to not breach the physical boundary.
  if(PROPCON_constrain_displacements >= 0){
    MPI_Barrier(mycomm);
    if(!rank && man_verbose > 2)
      std::cout << "Rocstar: Finding constraint intersections." << std::endl;
    //    std::cout << "Constraining displacements with Rocon" << std::endl;
    //    COM_call_function(PROPCON_find_intersections,&p_pmesh_hdl,&p_vm_hdl,&p_pos_hdl,&p_cflag_hdl);
    MPI_Barrier(mycomm);
    if(!rank && man_verbose > 2)
      std::cout << "Rocstar: Reducing intersections." << std::endl;
    //    COM_call_function( MAP_reduce_maxabs, &p_cflag_hdl);
    //    COM_call_function( MAP_reduce_minabs, &p_cflag_hdl);
    MPI_Barrier(mycomm);
    if(!rank && man_verbose > 2)
      std::cout << "Rocstar: Constraining displacements." << std::endl;
    //    COM_call_function(PROPCON_constrain_displacements,&p_pmesh_hdl,&p_vm_hdl,&p_pos_hdl,&p_cflag_hdl);
    MPI_Barrier(mycomm);
    if(!rank && man_verbose > 2)
      std::cout << "Rocstar: Displacements constrained." << std::endl;
  }

  COM_call_function( RocBlas::div_scalar, &p_vm_hdl, &dt, &p_vm_hdl);

  const RocmanControl_parameters *param = fagent->get_coupling()->get_rocmancontrol_param();

  // Let's try this to help fix Rocflo symmetry plane motion
  COM_call_function( MAP_reduce_maxabs, &a_vm_hdl);
  //  COM_call_function(MAP_reduce_minabs,&a_vm_hdl);

  if ( ! param->PROP_fom) {
     // Move shared nodes between burning and nonburning panes.
     COM_call_function( MAP_reduce_maxabs, &a_vm_hdl);
  }

  //
  // Other functions that go elsewhere that are associated 
  // with burnout that I can think of are:
  // 
  // Sets mass flux to 0.0 on burned out faces
  // BurnOutFilter(bflag,0.0,mass_flux);
  // or BurnOutFilter(bflag,mass_flux_insulation,mass_flux);
  //
  // Sets Tflame to Tgas for a perfect insulator, or
  // something else if maybe the insulation chars or
  // something.
  // BurnOutFilter(bflag,gas_temp,tflame);
  // or maybe BurOutFilter(bflag,ins_flame_temp,tflame);
  //
  // 

//debug_print(fagent->propBufAll+".vm", 102, 0);
}


// If solid has no ALE, mass conservation is not possible
// and the above approach will generate zero mdot, so we have to compute
// mdot as rhos*rb in this case.
MassTransfer::MassTransfer(FluidAgent *fag, BurnAgent *bag,
                  const std::string b_rhos, const std::string b_rb,
                  const std::string f_mdot): 
  Action(0, (const char**)NULL, NULL, NULL, (char *)"MassTransfer"), 
                   fagent(fag), bagent(bag)
{
  int io[] = {IN, IN, OUT};
  set_io( 3, io); 
  
  std::string atts[3];
  atts[0] = b_rhos;
  atts[1] = b_rb;
  atts[2] = f_mdot;
  set_attr(3, atts);

  bagent->register_use_dataitem( bagent->iburn_all, ".rhos", bagent->parentWin, ".rhos");
}

void MassTransfer::init( double t) {
  b_rhos_hdl = get_dataitem_handle( 0);
  b_rb_hdl = get_dataitem_handle( 1);
  f_mdot_hdl = get_dataitem_handle( 2);

  fb_mdot_hdl = COM_get_dataitem_handle( fagent->fluidBufB+".mdot");
}

void MassTransfer::run( double t, double dt, double alpha_dummy) {
  MAN_DEBUG(3, ("Rocstar: MassTransfer::run() with t:%e dt:%e.\n", t, dt));

  double zero = 0.0;
  COM_call_function( RocBlas::copy_scalar, &zero, &f_mdot_hdl);

  COM_call_function( RocBlas::mul, &b_rhos_hdl, &b_rb_hdl, &fb_mdot_hdl);

  // Sets mass flux to 0.0 on burned out faces
  // BurnOutFilter(bflag,0.0,mass_flux);
  // or BurnOutFilter(bflag,mass_flux_insulation,mass_flux);

//  debug_print(fagent->fluidBufB+".mdot", 102, 0);
//  debug_print(bagent->iburn_ng+".rb", 102, 0);
}

ZoomInterface::ZoomInterface(FluidAgent *fag, BurnAgent *bag, 
               const std::string fb_mdot_alp,
               double z): 
  Action(0, (const char**)NULL, NULL, NULL, (char *)"ZoomInterface"), 
               fagent(fag), bagent(bag), zoom(z)
{
  int io[] = {INOUT};
  set_io( 1, io); 
  
  std::string atts[1];
  atts[0] = fb_mdot_alp;
  set_attr(1, atts);
}

void ZoomInterface::init( double t) {
  fb_mdot_alp_hdl = get_dataitem_handle( 0);
  rhos_hdl        = COM_get_dataitem_handle(bagent->iburn_ng+".rhos");
  //  b_rb_alp_hdl = COM_get_dataitem_handle( bagent->iburn_ng+".rb_alp");
  b_rb_alp_hdl = COM_get_dataitem_handle( bagent->iburn_ng+".rb");
  fb_rhof_alp_hdl = COM_get_dataitem_handle( fagent->fluidBufB+".rhof_alp");
}

void ZoomInterface::run( double t, double dt, double alpha_dummy) {

  if (zoom <= 1.0) return;

  MAN_DEBUG(3, ("Rocstar: ZoomInterface::run() with t:%e dt:%e zoom:%e.\n", t, dt, zoom));

  //  double z = zoom - 1;

  //  COM_call_function( RocBlas::mul_scalar, &b_rb_alp_hdl, &z, &b_rb_alp_hdl);
 
  //  Rocman should just pass in the normal mdot  
  //  COM_call_function( RocBlas::axpy, &b_rb_alp_hdl, &fb_rhof_alp_hdl, &fb_mdot_alp_hdl, &fb_mdot_alp_hdl); (zoom)
  COM_call_function( RocBlas::mul, &b_rb_alp_hdl, &rhos_hdl, &fb_mdot_alp_hdl);
  // Sets mass flux to 0.0 on burned out faces
  // BurnOutFilter(bflag,0.0,mass_flux);
  // or BurnOutFilter(bflag,mass_flux_insulation,mass_flux);

  //  z = zoom / z;

  //  COM_call_function( RocBlas::mul_scalar, &b_rb_alp_hdl, &z, &b_rb_alp_hdl);
}

// UPDATE_INBUFF_BC_FLUID
// compute rhofvf
ComputeRhofvf::ComputeRhofvf(FluidAgent *fag, std::string f_vs_alp, std::string f_rhof_alp, std::string f_rhofvf_alp): 
  Action(0, (const char**)NULL, NULL, NULL, (char *)"ComputeRhofvf"), 
               fagent(fag)
{
  int io[] = {IN, IN, OUT};
  set_io( 3, io); 
  
  std::string atts[3];
  atts[0] = f_vs_alp;
  atts[1] = f_rhof_alp;
  atts[2] = f_rhofvf_alp;
  set_attr(3, atts);
}

void ComputeRhofvf::init( double t) {
  f_vs_alp_hdl = COM_get_dataitem_handle(attr[0]);
  f_rhof_alp_hdl = COM_get_dataitem_handle(attr[1]);
  f_rhofvf_alp_hdl = COM_get_dataitem_handle(attr[2]);
}

void ComputeRhofvf::run( double t, double dt, double alpha_dummy) {

  MAN_DEBUG(3, ("Rocstar: ComputeRhofvf::run() with t:%e dt:%e.\n", t, dt));

  COM_call_function( RocBlas::mul, &f_vs_alp_hdl, &f_rhof_alp_hdl, &f_rhofvf_alp_hdl);
}

// UPDATE_INBUFF_BC_FLUID fluid_agent.f90
// used in zomm>=1 and standalone
ComputeBurnPane::ComputeBurnPane(FluidAgent *fag, BurnAgent *bag, 
             SolidAgent *sag,
             const std::string fb_mdot_alp, const std::string rhofvf_alp,
             double z): 
  Action(0, (const char**)NULL, NULL, NULL, (char *)"ComputeBurnPane"),
             fagent(fag), bagent(bag), sagent(sag), zoom(z)
{
  int io[] = {IN, OUT};
  set_io( 2, io); 
  
  std::string atts[2];
  atts[0] = fb_mdot_alp;
  atts[1] = rhofvf_alp;
  set_attr(2, atts);
}

void ComputeBurnPane::init( double t) {
  fb_mdot_alp_hdl = get_dataitem_handle( 0);
  fb_rhofvf_alp_hdl = get_dataitem_handle( 1);

  b_rb_alp_hdl = COM_get_dataitem_handle( bagent->iburn_ng+".rb_alp");
  fb_rhof_alp_hdl = COM_get_dataitem_handle( fagent->fluidBufB+".rhof_alp");
  fb_nf_alp_hdl = COM_get_dataitem_handle( fagent->fluidBufB+".nf_alp");
  fb_rhofvf_alp_hdl = COM_get_dataitem_handle( fagent->fluidBufB+".rhofvf_alp");
}

void ComputeBurnPane::run( double t, double dt, double alpha_dummy) {
  MAN_DEBUG(3, ("Rocstar: ComputeBurnPane::run() with t:%e dt:%e zoom:%f.\n", t, dt, zoom));

  // The rest of computation are on burning panes only
  // Note that b_rb_alp is used here as a temporary variable.
  //  b_rb_alp is not used elsewhere.
  if (zoom>=1.0 && sagent==NULL || sagent && sagent->withALE) {                      
    COM_call_function( RocBlas::mul, &fb_rhof_alp_hdl, &b_rb_alp_hdl, &b_rb_alp_hdl);
    //  debug_print(fagent->fluidBufB+".mdot_alp", 102, 0);
    //  debug_print(fagent->fluidBufNG+".mdot_alp", 102, 0);
    COM_call_function( RocBlas::sub, &b_rb_alp_hdl, &fb_mdot_alp_hdl, &b_rb_alp_hdl);  
  }                                                                                   
  else {                                                              
    COM_call_function( RocBlas::neg, &fb_mdot_alp_hdl, &b_rb_alp_hdl);              
  }                                                                               
  
  // debug_print(bagent->iburn_ng+".rb_alp", 102, 0);
  //   debug_print(fagent->fluidBufB+".rhofvf_alp", 102, 0);
  
  COM_call_function( RocBlas::axpy, &b_rb_alp_hdl, &fb_nf_alp_hdl, &fb_rhofvf_alp_hdl, &fb_rhofvf_alp_hdl); 
  
  //debug_print(fagent->fluidBufB+".rhofvf_alp", 102, 0);
  //debug_print(fagent->fluidBufB+".mdot_alp", 102, 0);
  //debug_print(fagent->fluidBufB+".rhof_alp", 102, 0);
  //debug_print(fagent->fluidBufB+".nf_alp", 102, 0);
  //debug_print(bagent->iburn_ng+".rb_alp", 102, 0);
}

// INIT_INBUFF_BURN in burn_agent.f90
CopyBurnFromParentMesh::CopyBurnFromParentMesh(BurnAgent *bag, FluidAgent *fag): 
  Action(0, (const char**)NULL, NULL, NULL, (char *)"CopyBurnFromParentMesh"),
            bagent(bag), fagent(fag)
{
  set_io( 0, (int *)NULL); 
  
  set_attr(0, (const char**)NULL);
}

void CopyBurnFromParentMesh::init( double t) {
  burn_mesh = bagent->get_surface_window()+".mesh";
  parent_mesh = bagent->parentWin+".mesh";
  MAN_DEBUG(3, ("Rocstar: CopyBurnFromParentMesh::init() with t:%e %s  %s.\n", t, burn_mesh.c_str(), parent_mesh.c_str()));
}

void CopyBurnFromParentMesh::run( double t, double dt, double alpha_dummy) {
  MAN_DEBUG(3, ("Rocstar: CopyBurnFromParentMesh::run() with t:%e dt:%e %s %s.\n", t, dt, burn_mesh.c_str(), parent_mesh.c_str()));
  COM_copy_dataitem( burn_mesh.c_str(), parent_mesh.c_str());

/*
    // ????????
  double v = 0.0;
  int f_du_alp_hdl = COM_get_dataitem_handle(fagent->get_surface_window()+".du_alp");
  COM_call_function( RocBlas::copy_scalar, &v, &f_du_alp_hdl);
*/
}

CopyBflagFromBurn::CopyBflagFromBurn(BurnAgent *bag): 
  Action(0, (const char**)NULL, NULL, NULL, (char *)"CopyBflagFromBurn"),
                   bagent(bag)
{
  set_io( 0, (int *)NULL); 
  
  set_attr(0, (const char**)NULL);
}

void CopyBflagFromBurn::init( double t) {
  burn_bflag = bagent->get_surface_window()+".bflag";
  parent_bflag = bagent->parentWin+".bflag";
}

void CopyBflagFromBurn::run( double t, double dt, double alpha_dummy) {
  MAN_DEBUG(3, ("Rocstar: CopyBflagFromBurn::run() with t:%e dt:%e.\n", t, dt));

  
  if (!bagent->get_coupling()->initial_start() && 
      (bagent->ignmodel || t == bagent->get_coupling()->get_control_param()->init_time)) {
    COM_copy_dataitem( parent_bflag.c_str(), burn_bflag.c_str());
  }
  
}

// POST_INIT_FLUID
ComputePconn::ComputePconn(Agent *ag, std::string a_mesh, std::string a_pconn, std::string p_pmesh, int *cond):
  Action(0, (const char**)NULL, NULL, NULL, (char *)"ComputePconn"),
             agent(ag), a_mesh_str(a_mesh), a_pconn_str(a_pconn), p_pmesh_str(p_pmesh), cond_addr(cond)
{
  int io[3] = {OUT, OUT, OUT};
  set_io( 3, io); 
  
  std::string atts[3];
  atts[0] = a_mesh_str;
  atts[1] = a_pconn_str;
  atts[2] = p_pmesh_str;
  set_attr(3, atts);
}

void ComputePconn::init(double t)
{
  if (cond_addr && *cond_addr == 0) return;

  a_mesh_hdl = get_dataitem_handle( 0);
  a_pconn_hdl = get_dataitem_handle( 1);
  p_pmesh_hdl = get_dataitem_handle( 2);

  load_rocmap();
  MAP_compute_pconn = COM_get_function_handle("MAP.compute_pconn");

  load_rocprop(agent->get_coupling()->get_rocmancontrol_param(), agent->get_communicator());
  PROP_initialize = COM_get_function_handle("PROP.initialize");
}

void ComputePconn::run( double t, double dt, double alpha_dummy) {

  if (cond_addr && *cond_addr == 0) return;

  MAN_DEBUG(3, ("Rocstar: (%s) ComputePconn::run() with t:%e dt:%e.\n", agent->get_agent_name().c_str(), t, dt));

  // Compute pconn for the whole surface, which will also be used
  // by surface propagation.
  COM_call_function( MAP_compute_pconn, &a_mesh_hdl, &a_pconn_hdl);

  COM_call_function( PROP_initialize, &p_pmesh_hdl);

  // ??????????
/*
  double v = 0.0;
  int f_du_alp_hdl = COM_get_dataitem_handle(fagent->get_surface_window()+".du_alp");
  COM_call_function( RocBlas::copy_scalar, &v, &f_du_alp_hdl);

  double a = 0.0;
  int level = 1;
  fagent->obtain_bc(&a, &level);
  level = 2;
  fagent->obtain_bc(&a, &level);
*/
  // ??????????????
}

// zoom > 0 effective
SolidPropagateSurface_ALE::SolidPropagateSurface_ALE(SolidAgent *fag, 
            const std::string p_rb, const std::string a_vbar, double z): 
  Action(0, (const char**)NULL, NULL, NULL, (char *)"SolidPropagateSurface"),
            sagent(fag), zoom(z)
{
  int io[] = {IN, OUT};
  set_io( 2, io); 
  
  std::string atts[2];
  atts[0] = p_rb;
  atts[1] = a_vbar;
  set_attr(2, atts);

/*
  fagent->register_new_dataitem( fagent->propBufAll, ".vm", 'n', COM_DOUBLE, 3, "m/s");
  fagent->register_new_dataitem( fagent->propBufAll, ".rb", 'e', COM_DOUBLE, 1, "m/s");
*/
}

void SolidPropagateSurface_ALE::init( double t) {

  if (sagent->withALE == 0) return;

  p_rb_hdl = get_dataitem_handle( 0);
  a_vbar_hdl = get_dataitem_handle( 1);

  std::string propBuf = sagent->propBufAll;
  int PROP_fom = sagent->get_coupling()->get_rocmancontrol_param()->PROP_fom;
  if (!PROP_fom) propBuf = sagent->propBuf;
  p_vbar_hdl = COM_get_dataitem_handle( propBuf+".vbar");
  p_pmesh_hdl = COM_get_dataitem_handle( propBuf+".pmesh");
  p_vm_hdl = COM_get_dataitem_handle( propBuf+".vm");
  p_cnstr_hdl  = COM_get_dataitem_handle( propBuf+".cflag");
  p_pos_hdl    = COM_get_dataitem_handle( propBuf+".positions");
  //p_cnstr_type = COM_get_dataitem_handle( propBuf+".cnstr_type");

  load_rocprop(sagent->get_coupling()->get_rocmancontrol_param(), sagent->get_communicator());
  load_rocon(sagent->get_coupling()->get_rocmancontrol_param(), sagent->get_communicator());
  PROP_propagate = COM_get_function_handle("PROP.propagate");
  if(COM_get_window_handle("PROPCON") >= 0){
    PROPCON_find_intersections = COM_get_function_handle("PROPCON.find_intersections");
    PROPCON_constrain_displacements = COM_get_function_handle("PROPCON.constrain_displacements");
  }
  else {
    PROPCON_find_intersections = -1;
    PROPCON_constrain_displacements = -1;
  }
  load_rocmap();
  MAP_reduce_maxabs = COM_get_function_handle("MAP.reduce_maxabs_on_shared_nodes");
}

void SolidPropagateSurface_ALE::run( double t, double dt, double alpha_dummy) 
{
  if (sagent->withALE == 0) return;

  int rank = 0;
  MPI_Comm mycomm = sagent->get_communicator();
  MPI_Comm_rank(mycomm,&rank);
  if(!rank && man_verbose > 1)
    std::cout << "Rocstar: SolidPropagateSurface_ALE::run() with t:"<< t 
	      <<  "dt:" << dt << " zoom:" << zoom << std::endl;

  if ( !sagent->get_coupling()->initial_start() && zoom>0) 
  {
    COM_call_function( RocBlas::neg, &p_rb_hdl, &p_rb_hdl);   // Flip the sign of burn rate
    double dtz = dt*zoom;

    COM_call_function( PROP_propagate, &p_pmesh_hdl, &p_rb_hdl, &dtz, &p_vbar_hdl);

  if(!rank && man_verbose > 2)
    std::cout << "Rocstar: Rocprop done" << std::endl;

    COM_call_function( RocBlas::div_scalar, &p_vbar_hdl, &dt, &p_vbar_hdl);

      // Update shared nodes for the whole surface
    COM_call_function( MAP_reduce_maxabs, 1, &a_vbar_hdl);
  } 
  else {
    double dtz = 0;
    COM_call_function( RocBlas::copy_scalar, &dtz, &a_vbar_hdl);
  }
}

Reset_du_alp::Reset_du_alp( FluidAgent *fag) :
  Action( 0, (const char**)NULL, NULL, NULL, (char *)"Reset_du_alp"), fagent(fag)
{
  // Define I/O of the arguments.
  set_io( 0, NULL);
  set_attr(0, (const char**)NULL);
}

// Obtain the dataitem handles
void Reset_du_alp::init( double t) {
  du_alp_hdl = COM_get_dataitem_handle( fagent->get_surface_window()+".du_alp");
}

void Reset_du_alp::run( double t, double dt, double alpha) {
  MAN_DEBUG(3, ("Rocstar: Reset_du_alp::run().\n"));
  double zero = 0.0;
  COM_call_function( RocBlas::copy_scalar, &zero, &du_alp_hdl);
}






