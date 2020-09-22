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

#ifndef _BURNAGENT_H_
#define _BURNAGENT_H_

#include "RocstarAgent.h"

class BurnAgent : public RocstarAgent {
public:
  BurnAgent(RocstarCoupling *coup, std::string mod, std::string obj,
            MPI_Comm com, std::string parent);

  void read_restart_data() override;
  void output_visualization_files(double t) override;

private:
  int tbl_flag;

public:
  std::string parentWin;

public:
  std::string iburn_ng;

private:
  std::string burnIntBak;

public:
  bool ignmodel;

private:
  void finalize_windows() override;

  void parse_ic_options(void *option_data) override;
  void input_fallback_surf(const std::string &surface_window_in) override;
  void input_fallback_vol(const std::string &volume_window_in) override;
  void create_buffer_all() override;
};

#endif
