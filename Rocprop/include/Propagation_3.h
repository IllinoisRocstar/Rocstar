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
// $Id: Propagation_3.h,v 1.17 2008/12/06 08:45:27 mtcampbe Exp $

#ifndef _PROPAGATION_3_H_
#define _PROPAGATION_3_H_

#include <iostream>
#include <vector>
#include <map>
#include "propbasic.h"

PROP_BEGIN_NAMESPACE

const int MAX_LINES = 6;
const int BOUND_LEN = 15;

/** This class provides some common primitives used by various
 *  propagation algorithms. It is a pure virtual class. */
class Propagation_3 {
public: 
  Propagation_3()
    : _surf(NULL), _mode( SURF::ACROSS_PANE), _buf(NULL), _rank(0), 
      _verb(0), _cnstr_set(0), _cnstr_nodes(0), _cnstr_faces(0), 
      _cnstr_bndry_edges(0), _cnstr_bndry_nodes(), _cnstr_bound(0) {}

  virtual ~Propagation_3() {}

  /// Construct an object from a window manifold
  Propagation_3( Manifold *wm, COM::Window *buf);

  /** Main entry for invoking a surface propagation algorithm.
   *  Note that disps serves as both input and output. At input, it
   *  contains the (normalized) constrained nodal direction for the
   *  vertices with a customized directional or planar constraint 
   *  (\sa set_constraints()). At output, it contains the displacement
   *   vector for all vertices.
   *   
   *  \param spd nodal or facial speed function
   *  \param dt  time step
   *  \param disps nodal constraints at input and nodal displacements at output
   */
  virtual double time_stepping( const COM::DataItem *spd, double dt,
				COM::DataItem *disps, int *smoothed=NULL)=0;

  /** Set the types and directions of nodal constraints. 
   *  \param cnstr_types is an integer-type attribute. If it is a nodal
   *         attribute, then the value for each node is one of the following:
   *         0: No constraint.
   *         1: Move along the direction to be given in the
   *            displacement vector when calling \sa time_stepping().
   *        -1: Move orthogonal to the direction to be given in the 
   *            displacement vector when calling \sa time_stepping().
   *         2: Fixed point.
   *       'x': Move along x direction.
   *       'y': Move along y direction.
   *       'z': Move along z direction.
   *      -'x': Move orthogonal to x direction (in the yz plane).
   *      -'y': Move orthogonal to y direction (in the xz plane).
   *      -'z': Move orthogonal to z direction (in the xy plane).
   *       't': Move orthogonal to normal direction.
   *
   *  If the attribute is facial, panel, or windowed, then its value 
   *  cannot be 1 or -1 but can be any other value as above. These constraints
   *  will be converted into nodal constraints by \sa convert_constraints().
   */
  virtual void set_constraints( const COM::DataItem *cnstr_types);

  /// Set the bounds
  void set_bounds( const COM::DataItem *bnd);

  // Bound a face based on nodal constraints. If all the nodes are beyond 
  // the bounded, then set the value to zero.
  void bound_facial_speed( COM::DataItem *fa);

  /// Enforces the nodal constraints by projecting motion onto given direction.
  virtual void enforce_nodal_constraints( COM::DataItem *du);

  virtual void bound_nodal_motion( COM::DataItem *disps);

  /// Set the verbose level. 
  void set_verbose( bool b) { _verb = b; }

protected:
  /// Convert facial constraints to nodal constraints
  void convert_constraints( const COM::DataItem *ctypes_faces,
			    COM::DataItem *ctypes_nodes);
  /// Convert facial or panel constraints to nodal constraints
  void determine_constraint_boundary( const COM::DataItem *ctypes_faces,
				      COM::DataItem *ctypes_bndry_edges,
				      COM::DataItem *ctypes_bndry_nodes);

  /// Get orthonormals of the constraints. If nnrms is 0, then not constrained.
  /// If nnrms is 3, then the point is fixed. Otherwise, the point is 
  /// constrained either in a plane or a line.
  static void get_constraint_directions( int type, int &ndirs, 
					 Vector_3 dirs[2]);

  /// Enforce constraint for a specific vector.
  static void enforce_nodal_constraint( int type, Vector_3 &du);

  static void bound_nodal_motion( const Point_3 &pnt, const double *bnd, 
				  Vector_3 &du, double eps=0);

  static bool check_spherical_bound( const Point_3 &pnt, const Point_3 &org, 
				     const double rad_min, const double rad_max,
				     double eps=0);
  static void bound_spherical_disp( const Point_3 &pnt, const Point_3 &org, 
				    const double rad_min, const double rad_max, 
				    Vector_3 &disp, double eps=0);

  static bool check_radial_bound( const double x, const double y,
				  const double bnd_min, const double bnd_max,
				  double eps=0);
  static void bound_radial_disp( const double x, const double y, 
				 const double bnd_min, const double bnd_max,
				 double &dx, double &dy, double eps=0);

  static bool check_axial_bound( const double x, const double bnd_min, 
				 const double bnd_max, double eps=0);
  static void bound_axial_disp( const double x, const double bnd_min, 
				const double bnd_max, double &dx, double eps=0);

  static bool reached_nodal_bound( const Point_3 &pnt, const double *bnd, 
				   double eps=0);

  static bool in_bounding_box( const Point_3 &pnt, 
			       const Point_3 &lb, const Point_3 &ub) {
    return pnt[0] >= lb[0] && pnt[0] <= ub[0] &&
      pnt[1] >= lb[1] && pnt[1] <= ub[1] &&
      pnt[2] >= lb[2] && pnt[2] <= ub[2];
  }

protected:
  static double    square(double x) { return x*x; }
  Manifold         *_surf;
  SURF::Access_Mode _mode;	// Mode for accessing manifold
  
  COM::Window      *_buf;	// The buffer window created by Rocprop.

  int             _rank;        // Process rank.
  bool            _verb;        // Whether to constraints were set
  bool            _cnstr_set;   // Whether to constraints were set
  bool            _bnd_set;     // Whether to bound was set
  COM::DataItem *_cnstr_nodes; // Stores types of nodal constaints.
  COM::DataItem *_cnstr_faces;	// Stores facial constaints
  COM::DataItem *_cnstr_bndry_edges;	// Stores edges along constraint boundaries
  COM::DataItem *_cnstr_bndry_nodes;	// Stores nodes along constraint boundaries
  COM::DataItem *_cnstr_bound;	// Stores cylindrical constraint

  std::vector< COM::Pane*>      _panes; // panes of buffer window
};

PROP_END_NAMESPACE

#endif






