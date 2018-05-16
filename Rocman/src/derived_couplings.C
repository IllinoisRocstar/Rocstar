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
// $Id: derived_couplings.C,v 1.74 2008/12/06 08:45:22 mtcampbe Exp $

#include "rocman.h"
#include "derived_couplings.h"
#include "FluidAgent.h"
#include "SolidAgent.h"
#include "BurnAgent.h"
#include "basic_actions.h"
#include "transfer_actions.h"
#include "Interpolate.h"


/**************************************************************************
            Fluid Solid Coupling with Farhat's Improved Staggered Scheme
                          (ISS method) without Burn
**************************************************************************/

FluidSolidISS::FluidSolidISS( 
                         const char *fluidmodule, const char *solidmodule, 
                         MPI_Comm com, 
                         Control_parameters *p, 
                         const RocmanControl_parameters *mp): 
         Coupling("FluidSolidISS", fluidmodule, solidmodule, p, mp)
{
  maxPredCorr = param->maxNumPredCorrCycles;

  // Create agents
  FluidAgent *fluid_agent = new FluidAgent(this, normalize_modname(fluidmodule), fluidmodule, com, 1);
  add_agent( fluid_agent);
  const std::string fluidBufNG = fluid_agent->fluidBufNG;
  const std::string fluid_propBufAll = fluid_agent->propBufAll;

  SolidAgent *solid_agent = new SolidAgent(this, normalize_modname(solidmodule), solidmodule, com, 1);
  add_agent( solid_agent);
  const std::string solidBuf = solid_agent->solidBuf;

  double zoom = 1.0;

    // ###########  FLUID  ############

  fluid_agent->add_bcinitaction( new GetDeformedMesh(fluid_agent, solid_agent, 
                                     solidBuf+".x", solidBuf+".uhat", 
                                     solidBuf+".nc"));

  fluid_agent->add_bcinitaction( new MeshMotionTransferISS(
                                     fluid_agent, solid_agent, 
                                     solidBuf+".u", solidBuf+".vs", 
				     fluidBufNG+".vm"));

    // fluid boundary conditions
  fluid_agent->add_bcaction( new DummyAction(), 1);

//  fluid_agent->add_bcaction( new Extrapolate_Linear(fluid_agent, 
//                                           ifluid_i+".vs", "_alp"), 2);

  fluid_agent->add_gmaction( new ComputeMeshMotion( fluid_agent, 
                                    fluid_propBufAll+".vm", 
                                    fluid_agent->get_surface_window()+".du_alp",
                                    zoom));

    // fluid main Physics routine
  scheduler.add_action( fluid_agent->get_main_action());

    // ###########  SOLID  ############

  solid_agent->add_bcinitaction( new LoadTransfer_FS(fluid_agent, solid_agent, 
                             fluidBufNG+".ts", solidBuf+".ts", solidBuf+".pf"));

  solid_agent->add_bcaction( new Interpolate_Constant(solid_agent, solid_agent,
                                                    solidBuf+".ts"), 1);
  
    // solid main physics routines
  scheduler.add_action( solid_agent->get_main_action());
}


/**************************************************************************
            Fully coupled with Burn Energy
**************************************************************************/

SolidFluidBurnEnergySPC::SolidFluidBurnEnergySPC( 
           const char *fluidmodule, const char *solidmodule, 
           const char *burnmodule, MPI_Comm com, 
           Control_parameters *p, const RocmanControl_parameters *mp): 
           FullyCoupling("SolidFluidBurnEnergySPC", fluidmodule, solidmodule, burnmodule, p, mp)
{
  maxPredCorr = param->maxNumPredCorrCycles;

    // Create agents

  fluid_agent = new FluidAgent( this, normalize_modname(fluidmodule), fluidmodule, com, 1);
  add_agent( fluid_agent);
  const std::string ifluid_all = fluid_agent->ifluid_all;
  const std::string fluid_propBufAll = fluid_agent->propBufAll;
  const std::string fluidBufB = fluid_agent->fluidBufB;
  const std::string fluidBufNB = fluid_agent->fluidBufNB;
  const std::string fluidBufNG = fluid_agent->fluidBufNG;

  solid_agent = new SolidAgent( this, normalize_modname(solidmodule), solidmodule, com, 1);
  add_agent( solid_agent);
  const std::string solidBuf = solid_agent->solidBuf;
  const std::string solid_propBufAll = solid_agent->propBufAll;
  const std::string solid_propBuf = solid_agent->propBuf;

  burn_agent = new BurnAgent( this, "Rocburn", burnmodule, com, fluidBufB);
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

  if (initial_start()) {
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
  if (initial_start()) {
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
  solid_agent->add_bcinitaction( new HeatTransfer_FS(fluid_agent, solid_agent,
				   burn_agent,
                                   fluidBufNG+".qc", fluidBufNG+".qr", 
                                   iburn_ng+".qev",
                                   solidBuf+".qs"));
/*
  solid_agent->add_bcinitaction( new LoadTransfer_FSc_ALE(fluid_agent, 
                                   solid_agent, burn_agent, 
                                   fluidBufNG+".pf", fluidBufB+".mdot", 
                                   iburn_ng+".rb", solidBuf+".ts", 
                                   solidBuf+".pf"));
*/
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
  fluid_agent->add_bcinitaction( new TemperatureTransfer_SF(solid_agent,
                                     fluid_agent, 
                                     solidBuf+".Ts", 
                                     iburn_ng+".Tflm", fluidBufNB+".Tb"));
  fluid_agent->add_bcinitaction( new MeshMotionTransfer_SF(fluid_agent, 
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
  fluid_agent->add_bcaction( new Interpolate_Linear(fluid_agent, fluid_agent,
                                     fluidBufNB+".Tb"), 1);

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
}







