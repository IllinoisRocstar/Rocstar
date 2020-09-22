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

#ifndef _FLUIDAGENT_H_
#define _FLUIDAGENT_H_

#include "RocstarAgent.h"

class FluidAgent : public RocstarAgent {
public:
  FluidAgent(RocstarCoupling *coup, std::string mod, std::string obj,
             MPI_Comm com, bool withSolid = false);

  void read_restart_data() override;
  void output_visualization_files(double t) override;

  void init_convergence(int iPredCorr) override;
  bool check_convergence() const override;

  int compute_integrals() override;

private:
  bool with_plag;
  std::string plag_window;

  const bool with_solid;

private:
  std::string fluid_plag;

  std::string fluidPlagIn;
  std::string fluidVPIn;

public:
  std::string propBufAll;
  std::string fluidBufNG;
  std::string propBuf;
  std::string fluidBufB;
  std::string fluidBufNB;

private:
  std::string fluidBufPC;
  std::string fluidBufBak;
  std::string fluidBufPRE;

  std::string fluidVolBak;
  std::string fluidPlagBak;

  /**
   * @name Predictor-corrector iteration data
   */
  ///@{
  double tolerTrac;   ///< Traction tolerance for convergence check
  double tolerVelo;   ///< Mesh motion velocity tolerance for convergence check
  double tolerMass;   ///< Mass flux tolerance for convergence check
  int f_vm_hdl;       ///< COM handle to mesh motion velocity
  int f_vm_pre_hdl;   ///< COM handle to previous mesh motion velocity
  int f_mdot_hdl;     ///< COM handle to mass flux
  int f_mdot_pre_hdl; ///< COM handle to previous mass flux
  int f_ts_hdl;       ///< COM handle to traction
  int f_ts_pre_hdl;   ///< COM handle to previous traction
  ///@}

  // for update distances
  int nc_hdl, nc_tmp_hdl, sq_dist_hdl;

private:
  void finalize_windows() override;

  void create_buffer_all() override;
};

#endif
