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
// $Id: Remesher_Simmetrix.h,v 1.4 2008/12/06 08:45:27 mtcampbe Exp $

#ifndef _REMESHING_SIMMETRIX_H_
#define _REMESHING_SIMMETRIX_H_

#include "Rocprop.h"

PROP_BEGIN_NAMESPACE

// A serial remesher for triangular surface meshes based on Simmetrix tools.
class Remesher_Simmetrix : Remesher_base {

public:
  // Initialize the remesher (in particular, MeshSim).
  explicit Remesher_Simmetrix(const char *logfile=NULL);

  // Finalize the remesher (in particular, MeshSim).
  virtual ~Remesher_Simmetrix();

  // Main entry for remeshing. 
  virtual void remesh_serial( Manifold *wm, COM::DataItem *mesh_out,
			      double lave, double fangle);

  // Set global mesh size. See MS_setGlobalMeshSize() of MeshSim for 
  // the semantics of the parameters.
  void set_global_mesh_size( int type, double val) 
  { _size_type = type; _size_val = val; }

  // Set feature angle (in degrees). If 0 (the default), 
  // then no feature detection is performed.
  void set_fangle( double fangle) { _fangle = fangle; }

  // Set short ratio. See DM_modifySurfaceMesh() of GeomSim for
  // the semantics of the parameter.
  void set_short_ratio( double r) { _shrtRatio = r; }

protected:  
  // Convert a window object into a MeshSim object
  void *window_to_simmesh( const COM::Window *outwin);

  // Convert a MeshSim object into a window
  void simmesh_to_window( void *mesh, COM::Window *outwin);   

protected:
  static int instances; // Number of instances of the remesher
  int    _size_type;    // Type of mesh size
  double _size_val;     // Value of mesh size 
  double _fangle;       // Feature angle
  double _shrtRatio;    // Collapse an edge if it is incident on an edge 
                        // that is at least _shrtRatio times as long
};

PROP_END_NAMESPACE

#endif






