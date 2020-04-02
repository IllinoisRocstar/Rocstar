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

#ifndef _SOLIDAGENT_H_
#define _SOLIDAGENT_H_

#include "RocstarAgent.h"

const int NO_SHEER = 1;
const int WITH_SHEER = 2;

class SolidAgent : public RocstarAgent {
public:
  SolidAgent(RocstarCoupling *coup, std::string mod, std::string obj,
             MPI_Comm com, bool withFluid = false);

  void read_restart_data() override;
  void output_visualization_files(double t) override; // not implemented

  int compute_integrals() override;

public:
  std::string solidBufBase; // a window for intermediate buffers for fluid-solid
                            // interface
  std::string solidBuf;     // a window for solid/fluid interaction
  std::string propBufAll;   // for surface propagation (containing all panes)
  std::string propBuf;      // for surface propagation

private:
  std::string solidBufBak;
  std::string solidVolBak;

  bool with_fluid; // flag: coupled
public:
  bool with_ALE; // TODO init it
  int rhos_mode;
  int size_ts;
  int traction_mode;

private:
  int y_hdl; // for compute distances

private:
  void finalize_windows() override;

  void create_buffer_all() override;
};

#endif
