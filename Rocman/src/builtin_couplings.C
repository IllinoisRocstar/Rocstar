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
// $Id: builtin_couplings.C,v 1.41 2010/02/18 21:47:40 juzhang Exp $

/** \file builtin_couplings.C
 *  * Contains the implementation of builtin coupling schemes. 
 *       */

/* Author: Gengbin Zheng */

#include "rocman.h"
#include "builtin_couplings.h"
#include "FluidAgent.h"
#include "SolidAgent.h"
#include "BurnAgent.h"
#include "basic_actions.h"
#include "transfer_actions.h"
#include "Interpolate.h"

/**************************************************************************
            Fluid Alone no Burn
**************************************************************************/

// fluid alone
static void declare_fluid_actions( FluidAgent *fluid_agent) {
  const std::string ifluid_i = fluid_agent->ifluid_i;
  const std::string fluidBufB = fluid_agent->fluidBufB;
  const std::string fluidBufNB = fluid_agent->fluidBufNB;
  const std::string fluidBufNG = fluid_agent->fluidBufNG;

  // No action for initial call back routine
  fluid_agent->add_icaction( new DummyAction());

  // UPDATE_INBUFF_BC_FLUID
  // Level-1 boundary condition: Set bflag and Tflm_alp to 0
  fluid_agent->add_bcaction( new SetZero( fluidBufB+".Tflm_alp"), 1);
  //fluid_agent->add_bcaction( new SetZero( fluidBufNB+".Tb_alp"), 1);
  //  fluid_agent->add_bcaction( new SetZero( fluidBufB+".bflag"), 1);

  // Level-2 boundary condition: Set mdot and rhofvf to 0
  fluid_agent->add_bcaction( new SetZero( fluidBufB+".mdot_alp"), 2);
  fluid_agent->add_bcaction( new SetZero( fluidBufNG+".rhofvf_alp"), 2);

  // Set displacement to zero.
  //fluid_agent->add_gmaction( new SetZero( fluidBufNG+".du_alp"));
  fluid_agent->add_gmaction( new Reset_du_alp( fluid_agent));
}

FluidAlone::FluidAlone( const char *module, MPI_Comm com, Control_parameters *p, const RocmanControl_parameters *mp): Coupling("FluidAloneNoBurn", module, p, mp)
{
//  p->zoomFactor = 1;    // ignore zoom factor
  maxPredCorr = 1;

  // Create agents
  FluidAgent *fluid_agent = new FluidAgent( this, module, module, com);
  add_agent( fluid_agent);

  // Create and register actions
  scheduler.add_action( fluid_agent->get_main_action());

  // Declare the actions for fluid-agent
  declare_fluid_actions( fluid_agent);
}

/**************************************************************************
            Solid Alone without Burn
**************************************************************************/

// solid alone
static void declare_solid_actions( SolidAgent *solid_agent) {
  std::string solidBuf = solid_agent->solidBuf;

    // No action for initial call back routine
  solid_agent->add_icaction( new DummyAction());

  double pressure = solid_agent->get_coupling()->get_rocmancontrol_param()->pressure;

    // INIT_INBUFF_SOLID() in "solid_agent.f90"
  solid_agent->add_bcinitaction( new SetValueDouble(solidBuf+".ts", pressure));

    // INIT_INTERP_HANDLES() in "solid_agent.f90"
  solid_agent->add_bcaction( new Extrapolate_Linear(solid_agent, solid_agent,
                                                    solidBuf+".ts"));
}

SolidAlone::SolidAlone( const char *module, MPI_Comm com, Control_parameters *p, const RocmanControl_parameters *mp): Coupling("SolidAlone", module, p, mp) 
{
  maxPredCorr = 1;

  // Create agents
  SolidAgent *solid_agent = new SolidAgent( this, module, module, com);
  add_agent( solid_agent);

  // Create and register actions
  scheduler.add_action( solid_agent->get_main_action());

  // Declare the actions for solid-agent
  declare_solid_actions( solid_agent);
}

/**************************************************************************
            Solid Fluid Coupling without Burn
**************************************************************************/

SolidFluidSPC::SolidFluidSPC( 
                      const char *fluidmodule, const char *solidmodule, 
                      MPI_Comm com, 
                      Control_parameters *p, 
                      const RocmanControl_parameters *mp): 
                 Coupling("FluidSolidSPC", fluidmodule, solidmodule, p, mp)
{
  maxPredCorr = param->maxNumPredCorrCycles;

  // Create agents
  SolidAgent *solid_agent = new SolidAgent(this, normalize_modname(solidmodule), solidmodule, com, 1);
  add_agent( solid_agent);
  const std::string solidBufBase = solid_agent->solidBufBase;
  const std::string solidBuf = solid_agent->solidBuf;

  FluidAgent *fluid_agent = new FluidAgent(this, normalize_modname(fluidmodule), fluidmodule, com, 1);
  add_agent( fluid_agent);
  const std::string propBufAll = fluid_agent->propBufAll;
  const std::string fluidBufNG = fluid_agent->fluidBufNG;

  double zoom = 1;    // ignore zoom

    // ###########  SOLID  ############

    // INIT_INBUFF_SOLID() 
  solid_agent->add_bcinitaction( new LoadTransfer_FS(fluid_agent, solid_agent, 
                            fluidBufNG+".ts", solidBuf+".ts", solidBuf+".pf"));

    // INIT_INTERP_HANDLES() in solid_agent.f90
  solid_agent->add_bcaction( new Extrapolate_Linear(solid_agent, solid_agent,
                                           solidBuf+".ts"), 1);

    // solid main physics routines
  scheduler.add_action( solid_agent->get_main_action());

    // ###########  FLUID  ############

    // POST_UPDATE_SOLID 
  fluid_agent->add_bcinitaction( new GetDeformedMesh(fluid_agent, solid_agent, 
                                     solidBuf+".x", solidBuf+".uhat", 
                                     solidBuf+".nc"));

    // INIT_INBUFF_FLUID 
  fluid_agent->add_bcinitaction(new MeshMotionTransfer_SF(fluid_agent, 
                                    solid_agent,
                                    solidBuf+".u", fluidBufNG+".total_disp",
                                    fluidBufNG+".vm"));

  fluid_agent->add_bcinitaction( new DeformationVelTransfer_SF(
                                     fluid_agent, solid_agent, 
                                     solidBuf+".vs", fluidBufNG+".vs"));

    // fluid main Physics routine
  scheduler.add_action( fluid_agent->get_main_action());

    // fluid boundary conditions
  fluid_agent->add_bcaction( new DummyAction(), 1);

    // INIT_INTERP_HANDLES() in "fluid_agent.f90"
  fluid_agent->add_bcaction( new Interpolate_Linear(fluid_agent, fluid_agent,
                                           fluidBufNG+".vs"), 2);

    // UPDATE_INBUFF_GM_FLUID() 
  fluid_agent->add_gmaction( new ComputeMeshMotion( fluid_agent, 
                                    propBufAll+".vm", 
                                    fluid_agent->get_surface_window()+".du_alp",
                                    zoom));
}



/**************************************************************************
            Fluid alone with Burn
**************************************************************************/

FluidBurnAlone::FluidBurnAlone( const char *fluidmodule, const char *burnmodule, 
            MPI_Comm com, Control_parameters *p, 
            const RocmanControl_parameters *mp): 
            Coupling("FluidBurnAlone", fluidmodule, burnmodule, p, mp)
{
  maxPredCorr = 1;

  // Create agents
  FluidAgent *fluid_agent = new FluidAgent( this, fluidmodule, fluidmodule, com);
  add_agent( fluid_agent);
  const std::string propBufAll = fluid_agent->propBufAll;
  const std::string fluidBufB = fluid_agent->fluidBufB;
  const std::string fluidBufNG = fluid_agent->fluidBufNG;

  BurnAgent *burn_agent = new BurnAgent( this, "Rocburn", burnmodule, com, fluidBufB);
  add_agent( burn_agent);
  const std::string iburn_ng = burn_agent->iburn_ng;

  const double rhoc = mp->rhoc;
  double zoom = burn_agent->get_coupling()->get_control_param()->zoomFactor;

    // ###########  INITIALIZATION  ############

   //  POST_INIT_FLUID
  if (initial_start()) {
    if (zoom > 0) {
      std::string propBuf = fluid_agent->propBufAll;
      int PROP_fom = mp->PROP_fom;
      if (!PROP_fom) propBuf = fluid_agent->propBuf;
      init_scheduler.add_action( new ComputePconn(fluid_agent, 
                                               propBufAll+".mesh", 
                                               propBufAll+".pconn",
                                               propBuf+".pmesh"));
    }
      // INIT_INBUFF_FLUID
    init_scheduler.add_action( new BCInitInvoker( fluid_agent));
    init_scheduler.add_action( new Reset_du_alp( fluid_agent));
    init_scheduler.add_action( new BCInvoker(fluid_agent, 1));
    init_scheduler.add_action( new BCInvoker(fluid_agent, 2));
  }

    // ###########  BURN  ############

//RAF  if (initial_start()) {
  if (initial_start()) {
    burn_agent->add_icaction( new BCInitInvoker(burn_agent));
      // UPDATE_INBUFF_BC_BURN called in POST_INIT_BURN
    burn_agent->add_icaction( new BCInvoker(burn_agent));
    burn_agent->add_icaction( new ComputeFaceCenters(burn_agent, 
						     iburn_ng+".nc", 
						     iburn_ng+".centers", 
						     iburn_ng+".normals"));
  }

    // POST_INIT_BURN and INIT_BUFF_BURN
    // INIT_INBUFF_BURN
  burn_agent->add_bcinitaction( new CopyBurnFromParentMesh(burn_agent, fluid_agent));
  burn_agent->add_bcinitaction( new SetValueDouble(iburn_ng+".rhos", rhoc));

  burn_agent->add_bcaction( new Extrapolate_Linear(burn_agent, fluid_agent,
                                   iburn_ng+".pf"));
  burn_agent->add_bcaction( new Extrapolate_Central(burn_agent, fluid_agent,
                                   iburn_ng+".qc", "_alp", 1));
  burn_agent->add_bcaction( new Extrapolate_Central(burn_agent, fluid_agent,
                                   iburn_ng+".qr", "_alp", 1));
  burn_agent->add_bcaction( new Extrapolate_Linear(burn_agent, fluid_agent,
                                   iburn_ng+".Tf", "_alp", 1));
  burn_agent->add_bcaction( new Extrapolate_Linear(burn_agent, fluid_agent,
                                   iburn_ng+".Tv", "_alp", 1));
  burn_agent->add_bcaction( new Extrapolate_Linear(burn_agent, fluid_agent,
                                   iburn_ng+".dn", "_alp", 1));
  burn_agent->add_bcaction( new Extrapolate_Linear(burn_agent, burn_agent,
                                   iburn_ng+".rhos", "_alp", 1));

    // burn main Physics routine
  scheduler.add_action( burn_agent->get_main_action());

    // ###########  FLUID  ############

    // INIT_INBUFF_FLUID
//  fluid_agent->add_bcinitaction( new DummyPrint(burn_agent, NULL, fluid_agent, "fluid"));
  fluid_agent->add_bcinitaction( new SetZero( propBufAll+".vm"));
  fluid_agent->add_bcinitaction( new FluidPropagateSurface( fluid_agent, 
                                   burn_agent, 
                                   iburn_ng+".rb", propBufAll+".vm", zoom));
  fluid_agent->add_bcinitaction( new CopyBflagFromBurn( burn_agent));
  fluid_agent->add_bcinitaction( new MassTransfer( fluid_agent, burn_agent,
                                   iburn_ng+".rhos", iburn_ng+".rb",
                                   fluidBufNG+".mdot"));

  // UPDATE_INBUFF_BC_FLUID
  // Level-1 boundary condition: 
  fluid_agent->add_bcaction( new Interpolate_Linear( burn_agent, burn_agent,
                                   iburn_ng+".Tflm"), 1);

  // Level-2 boundary condition
  fluid_agent->add_bcaction( new Interpolate_Central( fluid_agent, fluid_agent,
                                   fluidBufB+".mdot"), 2);
  //fluid_agent->add_bcaction( new DummyPrint(burn_agent, NULL, fluid_agent, "MDOT"), 2);
  fluid_agent->add_bcaction( new Interpolate_Central( burn_agent, burn_agent,
                                   iburn_ng+".rb"), 2);
  fluid_agent->add_bcaction( new ZoomInterface( fluid_agent, burn_agent, 
                                     fluidBufB+".mdot_alp", zoom), 2);
  fluid_agent->add_bcaction( new SetZero( fluidBufNG+".rhofvf_alp"), 2);
  fluid_agent->add_bcaction( new ComputeBurnPane( fluid_agent, burn_agent, NULL,
                                     fluidBufB+".mdot_alp",
                                     fluidBufB+".rhofvf_alp",
                                     zoom), 2);

    // UPDATE_INBUFF_GM_FLUID() 
  fluid_agent->add_gmaction( new ComputeMeshMotion( fluid_agent, 
                                    propBufAll+".vm", 
                                    fluid_agent->get_surface_window()+".du_alp",
                                    zoom));

  // Create and register actions
  scheduler.add_action( fluid_agent->get_main_action());

}

/**************************************************************************
            Solid with Burn 
**************************************************************************/

SolidBurn::SolidBurn( const char *solidmodule, const char *burnmodule, 
           MPI_Comm com, 
           Control_parameters *p, const RocmanControl_parameters *mp): 
           Coupling("SolidBurn", solidmodule, burnmodule, p, mp)
{
#if 0
  maxPredCorr = param->maxNumPredCorrCycles;

    // Create agents
  SolidAgent *solid_agent = new SolidAgent( this, normalize_modname(solidmodule), solidmodule, com, 1);
  add_agent( solid_agent);
  const std::string solidBuf = solid_agent->solidBuf;
  const std::string solid_propBufAll = solid_agent->propBufAll;
  const std::string solid_propBuf = solid_agent->propBuf;

  BurnAgent *burn_agent = new BurnAgent( this, "Rocburn", burnmodule, com, fluidBufB);
  add_agent( burn_agent);
  const std::string iburn_ng = burn_agent->iburn_ng;

  double zoom = 1.0;       // ignore zoom for coupled simulation
  int PROP_fom = mp->PROP_fom;

    // ###########  INITIALIZATION  ############

    // POST_INIT_SOLID
  init_scheduler.add_action( new GetDeformedMesh(fluid_agent, solid_agent, 
                                   solidBuf+".x", solidBuf+".uhat", 
                                   solidBuf+".nc"));
    // withALE POST_INIT_SOLID
  std::string propBuf = solid_agent->propBufAll;
  if (!PROP_fom) propBuf = solid_agent->propBuf;
  init_scheduler.add_action( new ComputePconn(solid_agent, 
                                   solid_agent->propBufAll+".mesh", 
                                   solid_agent->propBufAll+".pconn",
                                   propBuf+".pmesh", &solid_agent->withALE));
  if (initial_start())
    init_scheduler.add_action( new BCInvoker(solid_agent, 1));

  //if (param->remeshed) 
  //  init_scheduler.add_action( new SurfDiverAfterRemeshing( fluid_agent, solid_agent));

    // POST_INIT_FLUID
  if (initial_start()) {
    init_scheduler.add_action( new CopyValue(fluidBufNG+".nc", fluidBufNG+".nc_t0"));
      // INIT_INBUFF_FLUID called in POST_INIT_FLUID
    init_scheduler.add_action( new BCInitInvoker( fluid_agent));
      // UPDATE_INBUFF_GM_FLUID
    init_scheduler.add_action( new GMInvoker( fluid_agent));
     // UPDATE_INBUFF_BC_FLUID called in POST_INIT_FLUID
    init_scheduler.add_action( new BCInvoker(fluid_agent, 1));
    init_scheduler.add_action( new BCInvoker(fluid_agent, 2));
  }

    // POST_UPDATE_FLUID
  init_scheduler.add_action( new BCInitInvoker( solid_agent));
/*
  init_scheduler.add_action( new ComputeFluidLoad_ALE( fluid_agent, solid_agent,
                                    fluidBufNG+".pf", fluidBufB+".mdot", 
                                    iburn_ng+".rb", 
                                    fluidBufNG+".ts"));
*/
  init_scheduler.add_action( new LoadTransfer_FSc_ALE(fluid_agent, 
                                   solid_agent, burn_agent, 
                                   fluidBufNG+".pf", fluidBufB+".mdot", 
                                   iburn_ng+".rb", solidBuf+".ts", 
                                   solidBuf+".pf"));

    // INIT_INBUFF_SOLID
  init_scheduler.add_action( new BCInitInvoker( solid_agent));

    // POST_UPDATE_SOLID 
  init_scheduler.add_action( new GetDeformedMesh_ALE(
                                   fluid_agent, solid_agent, 
                                   solidBuf+".x", solidBuf+".uhat", 
                                   solidBuf+".nc", zoom));

    // ###########  BURN  ############

    // POST_INIT_BURN (IC) 
//RAF  if (initial_start()) {
  if (1) {
    burn_agent->add_icaction( new BCInitInvoker(burn_agent));
    burn_agent->add_icaction( new BCInvoker(burn_agent));
    burn_agent->add_icaction( new ComputeFaceCenters(burn_agent, 
                                   iburn_ng+".nc", iburn_ng+".centers"));
  }

  burn_agent->add_bcinitaction( new CopyBurnFromParentMesh(burn_agent, fluid_agent));
  burn_agent->add_bcinitaction( new TransferSolidDensity(fluid_agent, 
                                   solid_agent, 
                                   solidBuf+".rhos", fluidBufNG+".rhos"));

  burn_agent->add_bcaction( new Extrapolate_Linear(burn_agent, fluid_agent,
                                   iburn_ng+".pf"));
  burn_agent->add_bcaction( new Extrapolate_Central(burn_agent, fluid_agent,
                                   iburn_ng+".qc", "_alp", 1));
  burn_agent->add_bcaction( new Extrapolate_Central(burn_agent, fluid_agent,
                                   iburn_ng+".qr", "_alp", 1));
  burn_agent->add_bcaction( new Extrapolate_Linear(burn_agent, fluid_agent,
                                   iburn_ng+".Tf", "_alp", 1));
  burn_agent->add_bcaction( new Extrapolate_Linear(burn_agent, burn_agent,
                                   iburn_ng+".rhos", "_alp", 1));

    // burn main Physics routine
  scheduler.add_action( burn_agent->get_main_action());

    // ###########  SOLID  ############

    // INIT_INBUFF_SOLID() 
/*
  solid_agent->add_bcinitaction( new LoadTransfer_FSc_ALE(fluid_agent, 
                                   solid_agent, burn_agent, 
                                   fluidBufNG+".pf", fluidBufB+".mdot", 
                                   iburn_ng+".rb", solidBuf+".ts", 
                                   solidBuf+".pf"));
*/
  solid_agent->add_bcinitaction( new TransferBurnRate_FS_ALE(
                                   fluid_agent, solid_agent,
                                   iburn_ng+".rb", solidBuf+".rb"));
  solid_agent->add_bcinitaction( new SolidPropagateSurface_ALE( solid_agent, 
                                   solid_propBuf+".rb",
                                   solid_propBufAll+".vbar", zoom));
  solid_agent->add_bcinitaction( new CopyValue(solidBuf+".x", solidBuf+".nc", 
                                   &solid_agent->withALE));

    // INIT_INTERP_HANDLES() in "solid_agent.f90" // TODO: Change _alp as optional
  solid_agent->add_bcaction( new Extrapolate_Linear(solid_agent, solid_agent,
                                    solidBuf+".ts"), 1);
    // for ALE
  solid_agent->add_bcaction( new Extrapolate_Central(solid_agent, solid_agent,
                                    solid_propBufAll+".vbar", "_alp", 1), 1);

    // solid main physics routines
  scheduler.add_action( solid_agent->get_main_action());

#endif
}

/**************************************************************************
            Fully coupled (Fluid/solid with Burn)
**************************************************************************/

SolidFluidBurnSPC::SolidFluidBurnSPC( 
           const char *fluidmodule, const char *solidmodule, 
           const char *burnmodule, MPI_Comm com, 
           Control_parameters *p, const RocmanControl_parameters *mp): 
           FullyCoupling("FluidSolidBurnSPC", fluidmodule, solidmodule, burnmodule, p, mp)
{
  maxPredCorr = param->maxNumPredCorrCycles;

    // Create agents
  fluid_agent = new FluidAgent( this, normalize_modname(fluidmodule), fluidmodule, com, 1);
  add_agent( fluid_agent);
  const std::string ifluid_all = fluid_agent->ifluid_all;
  const std::string fluid_propBufAll = fluid_agent->propBufAll;
  const std::string fluidBufB = fluid_agent->fluidBufB;
  const std::string fluidBufNG = fluid_agent->fluidBufNG;

  solid_agent = new SolidAgent( this, normalize_modname(solidmodule), solidmodule, com, 1);
  add_agent( solid_agent);
  const std::string solidBuf = solid_agent->solidBuf;
  const std::string solid_propBufAll = solid_agent->propBufAll;
  const std::string solid_propBuf = solid_agent->propBuf;

  burn_agent = new BurnAgent( this, "Rocburn", burnmodule, com, fluidBufB);
  add_agent( burn_agent);
  const std::string iburn_ng = burn_agent->iburn_ng;

  double zoom = burn_agent->get_coupling()->get_control_param()->zoomFactor;
  //  double zoom = 1.0;       // ignore zoom for coupled simulation
  int PROP_fom = mp->PROP_fom;

    // ###########  INITIALIZATION  ############

    // POST_INIT_SOLID
  init_scheduler.add_action( new GetDeformedMesh(fluid_agent, solid_agent, 
                                   solidBuf+".x", solidBuf+".uhat", 
                                   solidBuf+".nc"));
    // withALE POST_INIT_SOLID
  std::string propBuf = solid_agent->propBufAll;
  if (!PROP_fom) propBuf = solid_agent->propBuf;
  init_scheduler.add_action( new ComputePconn(solid_agent, 
                                   solid_agent->propBufAll+".mesh", 
                                   solid_agent->propBufAll+".pconn",
                                   propBuf+".pmesh", &solid_agent->withALE));
  if (initial_start() || mp->remeshed)
    init_scheduler.add_action( new BCInvoker(solid_agent, 1));

    // POST_INIT_FLUID
  if (initial_start() || mp->remeshed) {
    if (!mp->remeshed)
      init_scheduler.add_action( new CopyValue(fluidBufNG+".nc", fluidBufNG+".nc_t0"));
    else    // remeshed, restore nc_t0 from estimate
      init_scheduler.add_action(new RemeshInit(fluid_agent, 
                                     solid_agent,
                                     solidBuf+".u", fluidBufNG+".total_disp",
                                     fluidBufNG+".nc", fluidBufNG+".nc_t0"));
      // INIT_INBUFF_FLUID called in POST_INIT_FLUID
    init_scheduler.add_action( new BCInitInvoker( fluid_agent));
      // UPDATE_INBUFF_GM_FLUID
    init_scheduler.add_action( new GMInvoker( fluid_agent));
      // UPDATE_INBUFF_BC_FLUID called in POST_INIT_FLUID
    init_scheduler.add_action( new BCInvoker(fluid_agent, 1));
    init_scheduler.add_action( new BCInvoker(fluid_agent, 2));
  }

  if (initial_start() || mp->remeshed) {
      // POST_UPDATE_FLUID
    init_scheduler.add_action( new ComputeFluidLoad_ALE( 
                                    fluid_agent, solid_agent,
                                    fluidBufNG+".pf", fluidBufB+".mdot", 
                                    iburn_ng+".rb", 
                                    fluidBufNG+".ts"));
      // INIT_INBUFF_SOLID
/*
    init_scheduler.add_action( new LoadTransfer_FSc_ALE(fluid_agent, 
                                   solid_agent, burn_agent, 
                                   fluidBufNG+".pf", fluidBufB+".mdot", 
                                   iburn_ng+".rb", solidBuf+".ts", 
                                   solidBuf+".pf"));
*/
    init_scheduler.add_action( new BCInitInvoker( solid_agent));

      // POST_UPDATE_SOLID 
    init_scheduler.add_action( new GetDeformedMesh_ALE(
                                   fluid_agent, solid_agent, 
                                   solidBuf+".x", solidBuf+".uhat", 
                                   solidBuf+".nc", zoom));
  }

    // ###########  BURN  ############

    // POST_INIT_BURN (IC) 
//RAF  if (initial_start() || mp->remeshed) {
  if (1) {
    burn_agent->add_icaction( new BCInitInvoker(burn_agent));
    burn_agent->add_icaction( new BCInvoker(burn_agent));
    burn_agent->add_icaction( new ComputeFaceCenters(burn_agent, 
                                   iburn_ng+".nc", iburn_ng+".centers"));
  }

  burn_agent->add_bcinitaction( new CopyBurnFromParentMesh(burn_agent, fluid_agent));
  burn_agent->add_bcinitaction( new TransferSolidDensity(fluid_agent, 
                                   solid_agent, 
                                   solidBuf+".rhos", fluidBufNG+".rhos"));

  burn_agent->add_bcaction( new Extrapolate_Linear(burn_agent, fluid_agent,
                                   iburn_ng+".pf"));
  burn_agent->add_bcaction( new Extrapolate_Central(burn_agent, fluid_agent,
                                   iburn_ng+".qc", "_alp", 1));
  burn_agent->add_bcaction( new Extrapolate_Central(burn_agent, fluid_agent,
                                   iburn_ng+".qr", "_alp", 1));
  burn_agent->add_bcaction( new Extrapolate_Linear(burn_agent, fluid_agent,
                                   iburn_ng+".Tf", "_alp", 1));
  burn_agent->add_bcaction( new Extrapolate_Linear(burn_agent, burn_agent,
                                   iburn_ng+".rhos", "_alp", 1));

    // burn main Physics routine
  scheduler.add_action( burn_agent->get_main_action());

    // ###########  SOLID  ############

    // INIT_INBUFF_SOLID() 
  solid_agent->add_bcinitaction( new LoadTransferOnly_FSc_ALE(fluid_agent, 
                                   solid_agent, burn_agent, 
                                   fluidBufNG+".ts", 
                                   solidBuf+".ts", solidBuf+".pf"));
  solid_agent->add_bcinitaction( new TransferBurnRate_FS_ALE(
                                   fluid_agent, solid_agent,
                                   iburn_ng+".rb", solidBuf+".rb"));
  solid_agent->add_bcinitaction( new SolidPropagateSurface_ALE( solid_agent, 
                                   solid_propBuf+".rb",
                                   solid_propBufAll+".vbar", zoom));
  solid_agent->add_bcinitaction( new CopyValue(solidBuf+".x", solidBuf+".nc", 
                                   &solid_agent->withALE));

    // INIT_INTERP_HANDLES() in "solid_agent.f90" // TODO: Change _alp as optional
  solid_agent->add_bcaction( new Extrapolate_Linear(solid_agent, solid_agent,
                                    solidBuf+".ts"), 1);
    // for ALE
  solid_agent->add_bcaction( new Extrapolate_Central(solid_agent, solid_agent,
                                    solid_propBufAll+".vbar", "_alp", 1), 1);

    // solid main physics routines
  scheduler.add_action( solid_agent->get_main_action());

    // ###########  FLUID  ############

    // POST_UPDATE_SOLID 
  fluid_agent->add_bcinitaction( new GetDeformedMesh_ALE(
                                     fluid_agent, solid_agent, 
                                     solidBuf+".x", solidBuf+".uhat", 
                                     solidBuf+".nc", zoom));

    // INIT_INBUFF_FLUID 
  fluid_agent->add_bcinitaction( new SetZero( fluid_propBufAll+".vm"));
  fluid_agent->add_bcinitaction(new MeshMotionTransfer_SF(fluid_agent, 
                                     solid_agent,
                                     solidBuf+".u", fluidBufNG+".total_disp",
                                     fluidBufNG+".vm"));
  fluid_agent->add_bcinitaction( new DeformationVelTransfer_SF(
                                     fluid_agent, solid_agent, 
                                     solidBuf+".vs", fluidBufNG+".vs"));
  fluid_agent->add_bcinitaction( new MassTransfer_SF_ALE(fluid_agent, 
                                     solid_agent, burn_agent, 
                                     fluidBufNG+".mdot"));
  fluid_agent->add_bcinitaction( new CopyBflagFromBurn(burn_agent));

  // UPDATE_INBUFF_BC_FLUID
  // Level-1 boundary condition: 
  fluid_agent->add_bcaction( new Interpolate_Linear(burn_agent, burn_agent,
                                     iburn_ng+".Tflm"), 1);

  // Level-2 boundary condition
  fluid_agent->add_bcaction( new Interpolate_Central(fluid_agent, fluid_agent,
                                     fluidBufB+".mdot"), 2);
  fluid_agent->add_bcaction( new Interpolate_Central(burn_agent, burn_agent,
                                     iburn_ng+".rb"), 2);
  fluid_agent->add_bcaction( new ZoomInterface( fluid_agent, burn_agent, 
                                     fluidBufB+".mdot_alp", zoom), 2);
  fluid_agent->add_bcaction( new Interpolate_Linear(fluid_agent, fluid_agent,
                                     fluidBufNG+".vs"), 2);
  fluid_agent->add_bcaction( new ComputeRhofvf( fluid_agent, 
                                     fluidBufNG+".vs_alp", 
                                     fluidBufNG+".rhof_alp",
                                     fluidBufNG+".rhofvf_alp"), 2);
  fluid_agent->add_bcaction( new ComputeBurnPane( fluid_agent, burn_agent, 
                                     solid_agent, 
                                     fluidBufB+".mdot_alp",
                                     fluidBufB+".rhofvf_alp",
                                     zoom), 2);

    // UPDATE_INBUFF_GM_FLUID() 
  fluid_agent->add_gmaction( new ComputeMeshMotion( fluid_agent, 
                                    fluid_propBufAll+".vm", 
                                    fluid_agent->get_surface_window()+".du_alp",
                                    zoom));

    // fluid main Physics routine
  scheduler.add_action( fluid_agent->get_main_action());

    // POST_UPDATE_FLUID
  scheduler.add_action( new ComputeFluidLoad_ALE( fluid_agent, solid_agent, 
                                    fluidBufNG+".pf", fluidBufB+".mdot", 
                                    iburn_ng+".rb", 
                                    fluidBufNG+".ts"));

/*
  scheduler.add_action( new LoadTransfer_FSc_ALE(fluid_agent, 
                                   solid_agent, burn_agent, 
                                   fluidBufNG+".pf", fluidBufB+".mdot", 
                                   iburn_ng+".rb", solidBuf+".ts", 
                                   solidBuf+".pf"));
*/

//  scheduler.add_action( new SurfDiver(fluid_agent,  solid_agent));
}







