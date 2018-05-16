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

using namespace std;
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <cmath>

#include "commpi.h"

#define MAN_INTEG_SIZE             9

#include "Action.h"

#include "Scheduler.h"

#include "Coupling.h"

#include "Agent.h"

#include "RocBlas.h"

#define CI_DATE "$Date: 2009/11/04 15:15:15 $ GMT"
#define VERSION "3.0.0"
#define BUILD_NUM  2

enum { MAN_INTEG_VOL=1, MAN_INTEG_MASS=2,                     // Volume and mass
       MAN_INTEG_XMOM=3, MAN_INTEG_YMOM=4, MAN_INTEG_ZMOM=5,  // Momemtum
       MAN_INTEG_ENER=6,             // Energy
       MAN_INTEG_IBAREA=7,           // Area of burning f-s interface
       MAN_INTEG_INBAREA=8,          // Area of non-burning f-s interface
       MAN_INTEG_VOL_UND=9           // Undeformed volume
};

class RocmanControl_parameters {
public:
  int    verbose;	     // verbosity of rocman
  int    separate_out;	     // write to individual directory of each proc

  int    order;              // order of interpolation scheme
  int    traction_mode;
  double P_ambient;
    // for solid
  double rhoc;			// Solid density
  double pressure;
  double burn_rate;

    // Data transfer parameters
  int    rfc_verb;
  int    rfc_order;		// order of quadrature rules
  int    rfc_iter;		// max iterations
  double rfc_tol;		// tolerance for iterative solver

  char    PROP_fom;		// Using face-offsetting method for surface propagation - true or false 'T' or 'F'
  int     PROP_rediter;		// number of smoothing iterations
  double  PROP_fangle;          // Threshold for feature angles
  char    PROPCON_enabled;      // 'T' or the default, 'F' (Rocon)
  int     PROPCON_ndiv;         // number of divisions for propagation constraints (Rocon)
  char    async_in;		// 
  char    async_out;

  int     remeshed;
public:
  RocmanControl_parameters();
  void read(MPI_Comm comm, int commrank);
  void print();
};

extern int man_verbose;

#define MAN_DEBUG(l, x)   if (man_verbose>=l) printf x;

void debug_print(const std::string str, int pane, int pe, MPI_Comm comm, const char *memo=NULL);
void debug_int_print(const std::string str, int pane, int pe, MPI_Comm comm, const char *memo=NULL);

#endif






