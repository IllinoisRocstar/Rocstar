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

#include "Agent.h"

class BurnAgent: public Agent {
public:

  BurnAgent(Coupling *coup, std::string mod, std::string obj, MPI_Comm com, const std::string parent);

  virtual void input( double t);

  virtual void load_module();
  virtual void unload_module();
  virtual void init_module(double t, double dt);
  virtual void finalize();

  virtual void read_restart_data();
  virtual void output_restart_files( double t);
  virtual void output_visualization_files( double t);

  virtual void create_buffer_all();
protected:
  // Window name.
  static const char *window_name;  
  int  tbl_flag;
public:
  std::string parentWin;  
//  std::string iburn_i;

  std::string iburn_all;		// for registering dataitems
  std::string iburn_ng;		

    // for input
  std::string iburn;
  std::string burn;

   // Rocin windows
  std::string burnSurfIN;
  std::string burnVolIN;

  std::string burnIntBak;

  std::string burnBufOUT;
  bool ignmodel;
};

#endif






