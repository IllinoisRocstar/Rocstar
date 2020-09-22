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
 * Contains the implementation of builtin coupling schemes.
 */

/* Author: Gengbin Zheng */

#include "builtin_couplings.h"

#include "Control_parameters.h"
#include "Interpolate.h"
#include "basic_actions.h"
#include "rocman.h"
#include "transfer_actions.h"

/**************************************************************************
            Fluid Alone no Burn
**************************************************************************/

FluidAlone::FluidAlone(MPI_Comm com, const Control_parameters *p,
                       const RocmanControl_parameters *mp,
                       const std::string &module)
    : RocstarCoupling("FluidAlone", com, p, mp, module) {
  maxPredCorr = param->maxNumPredCorrCycles;

  // Create agents
  auto *fluid_agent = new FluidAgent(this, module, module, com);
  add_agent(fluid_agent);
  const std::string fluidBufB = fluid_agent->fluidBufB;
  const std::string fluidBufNB = fluid_agent->fluidBufNB;
  const std::string fluidBufNG = fluid_agent->fluidBufNG;

  // Declare the actions for fluid-agent

  // No action for initial call back routine
  fluid_agent->add_icaction(new DummyAction());

  // UPDATE_INBUFF_BC_FLUID
  // Level-1 boundary condition: Set bflag and Tflm_alp to 0
  fluid_agent->add_bcaction(new SetZero(fluidBufB + ".Tflm_alp"), 1);
  // fluid_agent->add_bcaction(new SetZero(fluidBufNB + ".Tb_alp"), 1);
  // fluid_agent->add_bcaction(new SetZero(fluidBufB + ".bflag"), 1);

  // Level-2 boundary condition: Set mdot and rhofvf to 0
  fluid_agent->add_bcaction(new SetZero(fluidBufB + ".mdot_alp"), 2);
  fluid_agent->add_bcaction(new SetZero(fluidBufNG + ".rhofvf_alp"), 2);

  // Set displacement to zero.
  // fluid_agent->add_gmaction(new SetZero(fluidBufNG + ".du_alp"));
  fluid_agent->add_gmaction(new Reset_du_alp(fluid_agent));

  // Register main physics action
  scheduler->add_action(fluid_agent->get_main_action());
}

/**************************************************************************
            Fluid Alone with Burn
**************************************************************************/

FluidBurnAlone::FluidBurnAlone(MPI_Comm com, const Control_parameters *p,
                               const RocmanControl_parameters *mp,
                               const std::string &fluidmodule,
                               const std::string &burnmodule)
    : RocstarCoupling("FluidBurnAlone", com, p, mp, fluidmodule, burnmodule) {
  maxPredCorr = param->maxNumPredCorrCycles;

  const double rhoc = rocmanparam->rhoc;
  const double zoom = param->zoomFactor;
  const int order = rocmanparam->order;
  const int PROP_fom = rocmanparam->PROP_fom;

  // Create agents
  auto *fluid_agent = new FluidAgent(this, fluidmodule, fluidmodule, com);
  add_agent(fluid_agent);
  const std::string propBufAll = fluid_agent->propBufAll;
  const std::string fluidBufB = fluid_agent->fluidBufB;
  const std::string fluidBufNG = fluid_agent->fluidBufNG;
  const std::string propBuf = (PROP_fom ? propBufAll : fluid_agent->propBuf);

  auto *burn_agent = new BurnAgent(this, "Rocburn", burnmodule, com, fluidBufB);
  add_agent(burn_agent);
  const std::string iburn_ng = burn_agent->iburn_ng;

  // ###########  INITIALIZATION  ############

  // POST_INIT_FLUID
  if (is_initial_start()) {
    if (zoom > 0) {
      init_scheduler->add_action(
          new ComputePconn(fluid_agent, propBufAll + ".mesh",
                           propBufAll + ".pconn", propBuf + ".pmesh"));
    }
    // INIT_INBUFF_FLUID
    /* AEG: Invokers removed.
    init_scheduler->add_action(new BCInitInvoker(fluid_agent));
    */
    init_scheduler->add_action(new Reset_du_alp(fluid_agent));
    /* AEG: Invokers removed.
    init_scheduler->add_action(new BCInvoker(fluid_agent, 1));
    init_scheduler->add_action(new BCInvoker(fluid_agent, 2));
    */
  }

  // ###########  BURN  ############

  if (is_initial_start()) {
    /* AEG: Invokers removed.
    burn_agent->add_icaction(new BCInitInvoker(burn_agent));
    // UPDATE_INBUFF_BC_BURN called in POST_INIT_BURN
    burn_agent->add_icaction(new BCInvoker(burn_agent));
    */
    burn_agent->add_icaction(
        new ComputeFaceCenters(burn_agent, iburn_ng + ".nc",
                               iburn_ng + ".centers", iburn_ng + ".normals"));
  } else {
    // No action for initial call back routine
    burn_agent->add_icaction(new DummyAction());
  }

  // POST_INIT_BURN and INIT_BUFF_BURN
  // INIT_INBUFF_BURN
  burn_agent->add_bcinitaction(
      new CopyBurnFromParentMesh(burn_agent, fluid_agent));
  burn_agent->add_bcinitaction(new SetValueDouble(iburn_ng + ".rhos", rhoc));

  burn_agent->add_bcaction(new Extrapolate_Linear(
      burn_agent, fluid_agent, iburn_ng + ".pf", false, order));
  burn_agent->add_bcaction(new Extrapolate_Central(
      burn_agent, fluid_agent, iburn_ng + ".qc", true, order));
  burn_agent->add_bcaction(new Extrapolate_Central(
      burn_agent, fluid_agent, iburn_ng + ".qr", true, order));
  burn_agent->add_bcaction(new Extrapolate_Linear(
      burn_agent, fluid_agent, iburn_ng + ".Tf", true, order));
  burn_agent->add_bcaction(new Extrapolate_Linear(
      burn_agent, fluid_agent, iburn_ng + ".Tv", true, order));
  burn_agent->add_bcaction(new Extrapolate_Linear(
      burn_agent, fluid_agent, iburn_ng + ".dn", true, order));
  burn_agent->add_bcaction(new Extrapolate_Linear(
      burn_agent, burn_agent, iburn_ng + ".rhos", true, order));

  // Register main physics action
  scheduler->add_action(burn_agent->get_main_action());

  // ###########  FLUID  ############

  // No action for initial call back routine
  fluid_agent->add_icaction(new DummyAction());

  // INIT_INBUFF_FLUID
  // fluid_agent->add_bcinitaction(
  //     new DummyPrint(burn_agent, nullptr, fluid_agent, "fluid"));
  // MS : starts debugging
  fluid_agent->add_bcinitaction(new SetZero(propBufAll + ".vm"));
  fluid_agent->add_bcinitaction(new FluidPropagateSurface(
      fluid_agent, burn_agent, iburn_ng + ".rb", propBufAll + ".vm", zoom));
  fluid_agent->add_bcinitaction(new CopyBflagFromBurn(burn_agent));
  fluid_agent->add_bcinitaction(
      new MassTransfer(fluid_agent, burn_agent, iburn_ng + ".rhos",
                       iburn_ng + ".rb", fluidBufNG + ".mdot"));

  // UPDATE_INBUFF_BC_FLUID
  // Level-1 boundary condition:
  fluid_agent->add_bcaction(
      new Interpolate_Linear(burn_agent, burn_agent, iburn_ng + ".Tflm", order),
      1);

  // Level-2 boundary condition
  fluid_agent->add_bcaction(new Interpolate_Central(fluid_agent, fluid_agent,
                                                    fluidBufB + ".mdot", order),
                            2);
  // fluid_agent->add_bcaction(
  //     new DummyPrint(burn_agent, nullptr, fluid_agent, "MDOT"), 2);
  fluid_agent->add_bcaction(
      new Interpolate_Central(burn_agent, burn_agent, iburn_ng + ".rb", order),
      2);
  fluid_agent->add_bcaction(
      new ZoomInterface(fluid_agent, burn_agent, fluidBufB + ".mdot_alp", zoom),
      2);
  fluid_agent->add_bcaction(new SetZero(fluidBufNG + ".rhofvf_alp"), 2);
  fluid_agent->add_bcaction(
      new ComputeBurnPane(fluid_agent, burn_agent, nullptr,
                          fluidBufB + ".mdot_alp", fluidBufB + ".rhofvf_alp",
                          zoom),
      2);

  // UPDATE_INBUFF_GM_FLUID()
  fluid_agent->add_gmaction(new ComputeMeshMotion(
      fluid_agent, propBufAll + ".vm",
      fluid_agent->get_surf_win() + ".du_alp", zoom));

  // Register main physics action
  scheduler->add_action(fluid_agent->get_main_action());
}

/**************************************************************************
            Solid Alone without Burn
**************************************************************************/

SolidAlone::SolidAlone(MPI_Comm com, const Control_parameters *p,
                       const RocmanControl_parameters *mp,
                       const std::string &module)
    : RocstarCoupling("SolidAlone", com, p, mp, module) {
  maxPredCorr = param->maxNumPredCorrCycles;

  const double pressure = rocmanparam->pressure;
  const int order = rocmanparam->order;

  // Create agents
  auto *solid_agent = new SolidAgent(this, module, module, com);
  add_agent(solid_agent);
  const std::string solidBuf = solid_agent->solidBuf;

  // No action for initial call back routine
  solid_agent->add_icaction(new DummyAction());

  // INIT_INBUFF_SOLID() in "solid_agent.f90"
  solid_agent->add_bcinitaction(new SetValueDouble(solidBuf + ".ts", pressure));

  // INIT_INTERP_HANDLES() in "solid_agent.f90"
  solid_agent->add_bcaction(new Extrapolate_Linear(
      solid_agent, solid_agent, solidBuf + ".ts", false, order));

  // Register main physics action
  scheduler->add_action(solid_agent->get_main_action());
}

/**************************************************************************
            Solid Alone with Burn
**************************************************************************/

SolidBurnAlone::SolidBurnAlone(MPI_Comm com, const Control_parameters *p,
                               const RocmanControl_parameters *mp,
                               const std::string &solidmodule,
                               const std::string &burnmodule)
    : RocstarCoupling("SolidBurnAlone", com, p, mp, solidmodule, burnmodule) {
#if 0
  maxPredCorr = param->maxNumPredCorrCycles;

  const double zoom = 1.0; // ignore zoom for coupled simulation
  const int PROP_fom = rocmanparam->PROP_fom;
  const int order = rocmanparam->order;
  const int remeshed = rocmanparam->remeshed;

  // Create agents
  auto *solid_agent = new SolidAgent(this, solidmodule, solidmodule, com, true);
  add_agent(solid_agent);
  const std::string solidBuf = solid_agent->solidBuf;
  const std::string solid_propBufAll = solid_agent->propBufAll;
  const std::string solid_propBuf = solid_agent->propBuf;
  const std::string propBuf = PROP_fom ? solid_propBufAll : solid_propBuf;

  auto *burn_agent = new BurnAgent(this, "Rocburn", burnmodule, com, fluidBufB);
  add_agent(burn_agent);
  const std::string iburn_ng = burn_agent->iburn_ng;

  // ###########  INITIALIZATION  ############

  // POST_INIT_SOLID
  init_scheduler->add_action(
      new GetDeformedMesh(fluid_agent, solid_agent, solidBuf + ".x",
                          solidBuf + ".uhat", solidBuf + ".nc"));
  // with_ALE POST_INIT_SOLID
  init_scheduler->add_action(
      new ComputePconn(solid_agent, solid_agent->propBufAll + ".mesh",
                       solid_agent->propBufAll + ".pconn", propBuf + ".pmesh",
                       solid_agent->with_ALE));
  /* AEG: Invokers removed.
  if (is_initial_start())
    init_scheduler->add_action(new BCInvoker(solid_agent, 1));
  */

  /*
  if (remeshed)
    init_scheduler->add_action(
        new SurfDiverAfterRemeshing(fluid_agent, solid_agent));
  */

  // POST_INIT_FLUID
  if (is_initial_start()) {
    init_scheduler->add_action(
        new CopyValue(fluidBufNG + ".nc", fluidBufNG + ".nc_t0"));
    /* AEG: Invokers removed.
    // INIT_INBUFF_FLUID called in POST_INIT_FLUID
    init_scheduler->add_action(new BCInitInvoker(fluid_agent));
    // UPDATE_INBUFF_GM_FLUID
    init_scheduler->add_action(new GMInvoker(fluid_agent));
    // UPDATE_INBUFF_BC_FLUID called in POST_INIT_FLUID
    init_scheduler->add_action(new BCInvoker(fluid_agent, 1));
    init_scheduler->add_action(new BCInvoker(fluid_agent, 2));
    */
  }

  /* AEG: Invokers removed.
  // POST_UPDATE_FLUID
  init_scheduler->add_action(new BCInitInvoker(solid_agent));
  */

  /*
  init_scheduler->add_action(new ComputeFluidLoad_ALE(
      fluid_agent, solid_agent, fluidBufNG + ".pf", fluidBufB + ".mdot",
      iburn_ng + ".rb", fluidBufNG + ".ts"));
  */
  init_scheduler->add_action(new LoadTransfer_FSc_ALE(
      fluid_agent, solid_agent, burn_agent, fluidBufNG + ".pf",
      fluidBufB + ".mdot", iburn_ng + ".rb", solidBuf + ".ts",
      solidBuf + ".pf"));

  /* AEG: Invokers removed.
  // INIT_INBUFF_SOLID
  init_scheduler->add_action(new BCInitInvoker(solid_agent));
  */

  // POST_UPDATE_SOLID
  init_scheduler->add_action(
      new GetDeformedMesh_ALE(fluid_agent, solid_agent, solidBuf + ".x",
                              solidBuf + ".uhat", solidBuf + ".nc", zoom));

  // ###########  BURN  ############

  // POST_INIT_BURN (IC)
  // RAF  if (is_initial_start()) {
  if (true) {
    /* AEG: Invokers removed.
    burn_agent->add_icaction(new BCInitInvoker(burn_agent));
    burn_agent->add_icaction(new BCInvoker(burn_agent));
    */
    burn_agent->add_icaction(new ComputeFaceCenters(
        burn_agent, iburn_ng + ".nc", iburn_ng + ".centers"));
  } else {
    // No action for initial call back routine
    burn_agent->add_icaction(new DummyAction());
  }

  burn_agent->add_bcinitaction(
      new CopyBurnFromParentMesh(burn_agent, fluid_agent));
  burn_agent->add_bcinitaction(new TransferSolidDensity(
      fluid_agent, solid_agent, solidBuf + ".rhos", fluidBufNG + ".rhos"));

  burn_agent->add_bcaction(new Extrapolate_Linear(
      burn_agent, fluid_agent, iburn_ng + ".pf", false, order));
  burn_agent->add_bcaction(new Extrapolate_Central(
      burn_agent, fluid_agent, iburn_ng + ".qc", true, order));
  burn_agent->add_bcaction(new Extrapolate_Central(
      burn_agent, fluid_agent, iburn_ng + ".qr", true, order));
  burn_agent->add_bcaction(new Extrapolate_Linear(
      burn_agent, fluid_agent, iburn_ng + ".Tf", true, order));
  burn_agent->add_bcaction(new Extrapolate_Linear(
      burn_agent, burn_agent, iburn_ng + ".rhos", true, order));

  // burn main Physics routine
  scheduler->add_action(burn_agent->get_main_action());

  // ###########  SOLID  ############

  // No action for initial call back routine
  solid_agent->add_icaction(new DummyAction());

  // INIT_INBUFF_SOLID()
  /*
  solid_agent->add_bcinitaction(new LoadTransfer_FSc_ALE(
      fluid_agent, solid_agent, burn_agent, fluidBufNG + ".pf",
      fluidBufB + ".mdot", iburn_ng + ".rb", solidBuf + ".ts",
      solidBuf + ".pf"));
  */
  solid_agent->add_bcinitaction(new TransferBurnRate_FS_ALE(
      fluid_agent, solid_agent, iburn_ng + ".rb", solidBuf + ".rb"));
  solid_agent->add_bcinitaction(new SolidPropagateSurface_ALE(
      solid_agent, solid_propBuf + ".rb", solid_propBufAll + ".vbar", zoom));
  solid_agent->add_bcinitaction(
      new CopyValue(solidBuf + ".x", solidBuf + ".nc", solid_agent->with_ALE));

  // INIT_INTERP_HANDLES() in "solid_agent.f90"
  solid_agent->add_bcaction(new Extrapolate_Linear(solid_agent, solid_agent,
                                                   solidBuf + ".ts", false,
                                                   order),
                            1);
  // for ALE
  solid_agent->add_bcaction(new Extrapolate_Central(solid_agent, solid_agent,
                                                    solid_propBufAll + ".vbar",
                                                    true, order),
                            1);

  // solid main physics routines
  scheduler->add_action(solid_agent->get_main_action());
#endif
}

/**************************************************************************
            Solid Fluid Coupling without Burn
**************************************************************************/

SolidFluidSPC::SolidFluidSPC(MPI_Comm com, const Control_parameters *p,
                             const RocmanControl_parameters *mp,
                             const std::string &fluidmodule,
                             const std::string &solidmodule)
    : FullyCoupling("SolidFluidSPC", com, p, mp, fluidmodule, solidmodule) {
  maxPredCorr = param->maxNumPredCorrCycles;

  const double zoom = 1.0; // ignore zoom for coupled simulation
  const int order = rocmanparam->order;

  // Create agents
  solid_agent = new SolidAgent(this, solidmodule, solidmodule, com, true);
  add_agent(solid_agent);
  const std::string solidBufBase = solid_agent->solidBufBase;
  const std::string solidBuf = solid_agent->solidBuf;

  fluid_agent = new FluidAgent(this, fluidmodule, fluidmodule, com, true);
  add_agent(fluid_agent);
  const std::string propBufAll = fluid_agent->propBufAll;
  const std::string fluidBufNG = fluid_agent->fluidBufNG;

  // ###########  SOLID  ############

  // No action for initial call back routine
  solid_agent->add_icaction(new DummyAction());

  // INIT_INBUFF_SOLID()
  solid_agent->add_bcinitaction(
      new LoadTransfer_FS(fluid_agent, solid_agent, fluid_agent->get_surf_win_i() + ".ts",
                          solidBuf + ".ts", solidBuf + ".pf"));

  // INIT_INTERP_HANDLES() in solid_agent.f90
  solid_agent->add_bcaction(new Extrapolate_Linear(solid_agent, solid_agent,
                                                   solidBuf + ".ts", false,
                                                   order),
                            1);

  // solid main physics routines
  scheduler->add_action(solid_agent->get_main_action());

  // ###########  FLUID  ############

  // No action for initial call back routine
  fluid_agent->add_icaction(new DummyAction());

  // POST_UPDATE_SOLID
  fluid_agent->add_bcinitaction(
      new GetDeformedMesh(fluid_agent, solid_agent, solidBuf + ".x",
                          solidBuf + ".uhat", solidBuf + ".nc"));

  // INIT_INBUFF_FLUID
  fluid_agent->add_bcinitaction(new MeshMotionTransfer_SF(
      fluid_agent, solid_agent, solidBuf + ".u", fluid_agent->get_surf_win_i() + ".total_disp",
      fluid_agent->get_surf_win_i() + ".vm"));

  fluid_agent->add_bcinitaction(new DeformationVelTransfer_SF(
      fluid_agent, solid_agent, solidBuf + ".vs", fluid_agent->get_surf_win_i() + ".vs"));
/*
  fluid_agent->add_bcinitaction(new MeshMotionTransfer_SF(
      fluid_agent, solid_agent, solidBuf + ".u", fluidBufNG + ".total_disp",
      fluidBufNG + ".vm"));

  fluid_agent->add_bcinitaction(new DeformationVelTransfer_SF(
      fluid_agent, solid_agent, solidBuf + ".vs", fluidBufNG + ".vs"));
*/

  // fluid boundary conditions
  fluid_agent->add_bcaction(new DummyAction(), 1);

  // INIT_INTERP_HANDLES() in "fluid_agent.f90"
  fluid_agent->add_bcaction(new Interpolate_Linear(fluid_agent, fluid_agent,
                                                   fluid_agent->get_surf_win_i() + ".vs", order),
                            2);

  // UPDATE_INBUFF_GM_FLUID()
  fluid_agent->add_gmaction(new ComputeMeshMotion(
      fluid_agent, propBufAll + ".vm",
      fluid_agent->get_surf_win() + ".du_alp", zoom));

  // fluid main Physics routine
  scheduler->add_action(fluid_agent->get_main_action());
}

/**************************************************************************
            Fully coupled (Fluid/solid with Burn)
**************************************************************************/

SolidFluidBurnSPC::SolidFluidBurnSPC(MPI_Comm com, const Control_parameters *p,
                                     const RocmanControl_parameters *mp,
                                     const std::string &fluidmodule,
                                     const std::string &solidmodule,
                                     const std::string &burnmodule)
    : FullyCoupling("SolidFluidBurnSPC", com, p, mp, fluidmodule, solidmodule,
                    burnmodule) {
  maxPredCorr = param->maxNumPredCorrCycles;

  const double zoom = 1.0; // ignore zoom for coupled simulation
  const int PROP_fom = rocmanparam->PROP_fom;
  const int order = rocmanparam->order;
  const int remeshed = rocmanparam->remeshed;

  // Create agents
  fluid_agent = new FluidAgent(this, fluidmodule, fluidmodule, com, true);
  add_agent(fluid_agent);
  const std::string fluid_propBufAll = fluid_agent->propBufAll;
  const std::string fluidBufB = fluid_agent->fluidBufB;
  const std::string fluidBufNG = fluid_agent->fluidBufNG;

  solid_agent = new SolidAgent(this, solidmodule, solidmodule, com, true);
  add_agent(solid_agent);
  const std::string solidBuf = solid_agent->solidBuf;
  const std::string solid_propBufAll = solid_agent->propBufAll;
  const std::string solid_propBuf = solid_agent->propBuf;
  const std::string propBuf = PROP_fom ? solid_propBufAll : solid_propBuf;
  const bool with_ALE = solid_agent->with_ALE;

  burn_agent = new BurnAgent(this, "Rocburn", burnmodule, com, fluidBufB);
  add_agent(burn_agent);
  const std::string iburn_ng = burn_agent->iburn_ng;

  // ###########  INITIALIZATION  ############

  // POST_INIT_SOLID
  init_scheduler->add_action(
      new GetDeformedMesh(fluid_agent, solid_agent, solidBuf + ".x",
                          solidBuf + ".uhat", solidBuf + ".nc"));
  // with_ALE POST_INIT_SOLID
  init_scheduler->add_action(
      new ComputePconn(solid_agent, solid_agent->propBufAll + ".mesh",
                       solid_agent->propBufAll + ".pconn", propBuf + ".pmesh",
                       solid_agent->with_ALE));
  /* AEG: Invokers removed.
  if (is_initial_start() || remeshed)
    init_scheduler->add_action(new BCInvoker(solid_agent, 1));
  */

  if (is_initial_start() || remeshed) {
    // POST_INIT_FLUID
    if (!remeshed) {
      init_scheduler->add_action(
          new CopyValue(fluidBufNG + ".nc", fluidBufNG + ".nc_t0"));
    } else { // remeshed, restore nc_t0 from estimate
      init_scheduler->add_action(new RemeshInit(
          fluid_agent, solid_agent, solidBuf + ".u", fluidBufNG + ".total_disp",
          fluidBufNG + ".nc", fluidBufNG + ".nc_t0"));
    }
    /* AEG: Invokers removed.
    // INIT_INBUFF_FLUID called in POST_INIT_FLUID
    init_scheduler->add_action(new BCInitInvoker(fluid_agent));
    // UPDATE_INBUFF_GM_FLUID
    init_scheduler->add_action(new GMInvoker(fluid_agent));
    // UPDATE_INBUFF_BC_FLUID called in POST_INIT_FLUID
    init_scheduler->add_action(new BCInvoker(fluid_agent, 1));
    init_scheduler->add_action(new BCInvoker(fluid_agent, 2));
    */

    // POST_UPDATE_FLUID
    init_scheduler->add_action(new ComputeFluidLoad_ALE(
        fluid_agent, solid_agent, fluidBufNG + ".pf", fluidBufB + ".mdot",
        iburn_ng + ".rb", fluidBufNG + ".ts"));
    /*
    // INIT_INBUFF_SOLID
    init_scheduler->add_action(new LoadTransfer_FSc_ALE(
        fluid_agent, solid_agent, burn_agent, fluidBufNG + ".pf",
        fluidBufB + ".mdot", iburn_ng + ".rb", solidBuf + ".ts",
        solidBuf + ".pf"));
    */
    /* AEG: Invokers removed.
    init_scheduler->add_action(new BCInitInvoker(solid_agent));
    */

    // POST_UPDATE_SOLID
    init_scheduler->add_action(
        new GetDeformedMesh_ALE(fluid_agent, solid_agent, solidBuf + ".x",
                                solidBuf + ".uhat", solidBuf + ".nc", zoom));
  }

  // ###########  BURN  ############

  // POST_INIT_BURN (IC)
  // RAF if (is_initial_start() || remeshed) {
  if (true) {
    /* AEG: Invokers removed.
    burn_agent->add_icaction(new BCInitInvoker(burn_agent));
    burn_agent->add_icaction(new BCInvoker(burn_agent));
    */
    burn_agent->add_icaction(new ComputeFaceCenters(
        burn_agent, iburn_ng + ".nc", iburn_ng + ".centers"));
  } else {
    // No action for initial call back routine
    burn_agent->add_icaction(new DummyAction());
  }

  burn_agent->add_bcinitaction(
      new CopyBurnFromParentMesh(burn_agent, fluid_agent));
  burn_agent->add_bcinitaction(new TransferSolidDensity(
      fluid_agent, solid_agent, solidBuf + ".rhos", fluidBufNG + ".rhos"));

  burn_agent->add_bcaction(new Extrapolate_Linear(
      burn_agent, fluid_agent, iburn_ng + ".pf", false, order));
  burn_agent->add_bcaction(new Extrapolate_Central(
      burn_agent, fluid_agent, iburn_ng + ".qc", true, order));
  burn_agent->add_bcaction(new Extrapolate_Central(
      burn_agent, fluid_agent, iburn_ng + ".qr", true, order));
  burn_agent->add_bcaction(new Extrapolate_Linear(
      burn_agent, fluid_agent, iburn_ng + ".Tf", true, order));
  burn_agent->add_bcaction(new Extrapolate_Linear(
      burn_agent, burn_agent, iburn_ng + ".rhos", true, order));

  // burn main Physics routine
  scheduler->add_action(burn_agent->get_main_action());

  // ###########  SOLID  ############

  // No action for initial call back routine
  solid_agent->add_icaction(new DummyAction());

  // INIT_INBUFF_SOLID()
  solid_agent->add_bcinitaction(new LoadTransferOnly_FSc_ALE(
      fluid_agent, solid_agent, burn_agent, fluidBufNG + ".ts",
      solidBuf + ".ts", solidBuf + ".pf"));
  solid_agent->add_bcinitaction(new TransferBurnRate_FS_ALE(
      fluid_agent, solid_agent, iburn_ng + ".rb", solidBuf + ".rb"));
  solid_agent->add_bcinitaction(new SolidPropagateSurface_ALE(
      solid_agent, solid_propBuf + ".rb", solid_propBufAll + ".vbar", zoom));
  solid_agent->add_bcinitaction(
      new CopyValue(solidBuf + ".x", solidBuf + ".nc", with_ALE));

  // INIT_INTERP_HANDLES() in "solid_agent.f90"
  solid_agent->add_bcaction(new Extrapolate_Linear(solid_agent, solid_agent,
                                                   solidBuf + ".ts", false,
                                                   order),
                            1);
  // for ALE
  solid_agent->add_bcaction(new Extrapolate_Central(solid_agent, solid_agent,
                                                    solid_propBufAll + ".vbar",
                                                    true, order),
                            1);

  // solid main physics routines
  scheduler->add_action(solid_agent->get_main_action());

  // ###########  FLUID  ############

  // No action for initial call back routine
  fluid_agent->add_icaction(new DummyAction());

  // POST_UPDATE_SOLID
  fluid_agent->add_bcinitaction(
      new GetDeformedMesh_ALE(fluid_agent, solid_agent, solidBuf + ".x",
                              solidBuf + ".uhat", solidBuf + ".nc", zoom));

  // INIT_INBUFF_FLUID
  fluid_agent->add_bcinitaction(new SetZero(fluid_propBufAll + ".vm"));
  fluid_agent->add_bcinitaction(new MeshMotionTransfer_SF(
      fluid_agent, solid_agent, solidBuf + ".u", fluidBufNG + ".total_disp",
      fluidBufNG + ".vm"));
  fluid_agent->add_bcinitaction(new DeformationVelTransfer_SF(
      fluid_agent, solid_agent, solidBuf + ".vs", fluidBufNG + ".vs"));
  fluid_agent->add_bcinitaction(new MassTransfer_SF_ALE(
      fluid_agent, solid_agent, burn_agent, fluidBufNG + ".mdot"));
  fluid_agent->add_bcinitaction(new CopyBflagFromBurn(burn_agent));

  // UPDATE_INBUFF_BC_FLUID
  // Level-1 boundary condition:
  fluid_agent->add_bcaction(
      new Interpolate_Linear(burn_agent, burn_agent, iburn_ng + ".Tflm", order),
      1);

  // Level-2 boundary condition
  fluid_agent->add_bcaction(new Interpolate_Central(fluid_agent, fluid_agent,
                                                    fluidBufB + ".mdot", order),
                            2);
  fluid_agent->add_bcaction(
      new Interpolate_Central(burn_agent, burn_agent, iburn_ng + ".rb", order),
      2);
  fluid_agent->add_bcaction(
      new ZoomInterface(fluid_agent, burn_agent, fluidBufB + ".mdot_alp", zoom),
      2);
  fluid_agent->add_bcaction(new Interpolate_Linear(fluid_agent, fluid_agent,
                                                   fluidBufNG + ".vs", order),
                            2);
  fluid_agent->add_bcaction(
      new ComputeRhofvf(fluid_agent, fluidBufNG + ".vs_alp",
                        fluidBufNG + ".rhof_alp", fluidBufNG + ".rhofvf_alp"),
      2);
  fluid_agent->add_bcaction(
      new ComputeBurnPane(fluid_agent, burn_agent, solid_agent,
                          fluidBufB + ".mdot_alp", fluidBufB + ".rhofvf_alp",
                          zoom),
      2);

  // UPDATE_INBUFF_GM_FLUID()
  fluid_agent->add_gmaction(new ComputeMeshMotion(
      fluid_agent, fluid_propBufAll + ".vm",
      fluid_agent->get_surf_win() + ".du_alp", zoom));

  // fluid main Physics routine
  scheduler->add_action(fluid_agent->get_main_action());

  // ###########  INTERACTION  ############

  // POST_UPDATE_FLUID
  scheduler->add_action(new ComputeFluidLoad_ALE(
      fluid_agent, solid_agent, fluidBufNG + ".pf", fluidBufB + ".mdot",
      iburn_ng + ".rb", fluidBufNG + ".ts"));

  /*
  scheduler->add_action(new LoadTransfer_FSc_ALE(
      fluid_agent, solid_agent, burn_agent, fluidBufNG + ".pf",
      fluidBufB + ".mdot", iburn_ng + ".rb", solidBuf + ".ts",
      solidBuf + ".pf"));
  */

  /*
  scheduler->add_action(new SurfDiver(fluid_agent, solid_agent));
  */
}
