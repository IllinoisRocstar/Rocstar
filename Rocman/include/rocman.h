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

#ifndef _ROCMAN_H_
#define _ROCMAN_H_

#include <string>

#include "com.h"

enum {
  MAN_INTEG_VOL = 1,     ///< Volume
  MAN_INTEG_MASS = 2,    ///< Mass
  MAN_INTEG_XMOM = 3,    ///< Momentum x
  MAN_INTEG_YMOM = 4,    ///< Momentum y
  MAN_INTEG_ZMOM = 5,    ///< Momentum z
  MAN_INTEG_ENER = 6,    ///< Energy
  MAN_INTEG_IBAREA = 7,  ///< Area of burning f-s interface
  MAN_INTEG_INBAREA = 8, ///< Area of non-burning f-s interface
  MAN_INTEG_VOL_UND = 9  ///< Undeformed volume
};

class RocmanControl_parameters {
public:
  int verbose;      // verbosity of rocman
  int separate_out; // write to individual directory of each proc

  int order; // order of interpolation scheme
  int traction_mode;
  double P_ambient;
  // for solid
  double rhoc; // Solid density
  double pressure;
  double burn_rate;

  // Data transfer parameters
  int rfc_verb;
  int rfc_order;  // order of quadrature rules
  int rfc_iter;   // max iterations
  double rfc_tol; // tolerance for iterative solver

  char PROP_fom; // Using face-offsetting method for surface propagation - true
                 // or false 'T' or 'F'
  int PROP_rediter;     // number of smoothing iterations
  double PROP_fangle;   // Threshold for feature angles
  char PROPCON_enabled; // 'T' or the default, 'F' (Rocon)
  int PROPCON_ndiv; // number of divisions for propagation constraints (Rocon)
  char async_in;    //
  char async_out;

  int remeshed;

public:
  RocmanControl_parameters();
  void read(MPI_Comm comm, int commrank);
  void print();
};

extern int man_verbose;

#define MAN_DEBUG(l, x)                                                        \
  if (man_verbose >= l)                                                        \
  printf x

void debug_print(const std::string &str, int pane, int pe, MPI_Comm comm,
                 const char *memo = nullptr);
void debug_int_print(const std::string &str, int pane, int pe, MPI_Comm comm,
                     const char *memo = nullptr);

#endif //_ROCMAN_H_
