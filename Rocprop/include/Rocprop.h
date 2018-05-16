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
// $Id: Rocprop.h,v 1.21 2008/12/06 08:45:27 mtcampbe Exp $

#ifndef _ROCPROP_H_
#define _ROCPROP_H_

#include "propbasic.h"
#include "com.h"
#include "com_devel.hpp"
#include "Propagation_3.h"
#include "Rocsurf.h"

PROP_BEGIN_NAMESPACE

class Remesher_base {
public:
  virtual void remesh_serial( Manifold *wm, COM::DataItem *mesh_out,
			      double lave, double fangle)=0;

  virtual ~Remesher_base() {}
};

class Rocprop : public SURF::Rocsurf {
public:
  enum { PROP_FO, PROP_MP};

  /// Defautl constructor.
  Rocprop() : _parent(NULL), _prop(NULL), _win(NULL), _buf(NULL), _rem(NULL),
	      _prop_method(PROP_FO), _wf_expn(-1), _wght_scheme(-1),
	      _nrm_dfsn(-1), _feature_layer(-1), _eig_thres(-1), 
	      _courant( -1), _fangle_strong(-1), 
	      _fangle_weak(-1), _fangle_turn(-1), _verb(0), _cnstr_types(NULL), 
	      _cnstr_bound(NULL), _time_lb(1./128), _smoother( -1), _rediter(4), 
	      _conserv(-1), _cookie(PROP_COOKIE) {}
  
  virtual ~Rocprop();

  /// Initialize Rocprop with given mesh.
  void initialize( const COM::DataItem *pmesh,
		   SURF::Rocsurf *rsurf=NULL);

  /// Perturb the given mesh by alpha times shortest edge length.
  void perturb_mesh( COM::DataItem *pmesh, const double &alpha);

  /** Set the types and directions of constraints
   *  \seelaso Propagation_3::set_constraints
   */
  void set_constraints( const COM::DataItem *cnstr_types);

  /** Set the bounds
   *  \seelaso Propagation_3::set_bounds
   */
  void set_bounds( const COM::DataItem *bnd);

  /**  Propagates the interface. This function subcycles until it reachs
   *   given time step or sub-time-step becomes too small. If the time
   *   step became too small, then abort if code is NULL, or return
   *   -1 if code is not NULL.
   */
  void propagate( const COM::DataItem *pmesh,
		  COM::DataItem *vel,
		  const double *dt, 
		  COM::DataItem *du,
		  double *dt_elapsed=NULL,
		  int *code=NULL);

  /// Register remesher. Rocprop does not own the remesher.
  void set_remesher( void *rem, int *owner=0) { 
    if ( _rem && _rem_owner) delete _rem;
    _rem = (Remesher_base*)rem; 
    _rem_owner = owner && *owner;
  }

  /// Invoke serial remeshing
  void remesh_serial( COM::DataItem *mesh_out, 
		      double *lave=NULL, double *fangle=NULL);

  /** Set options for propagation. The options supported are:
   *    "method": "mp" (marker particles, default), 
   *              "fo" (face offsetting)
   *    "wavefrontal": 0 or 1 (for true or false)
   *    "weight":    weighting scheme
   *    "eigthres":  encoding relative threshold for null space.
   *    "dirthres":  encoding relative threshold for tangent space.
   *    "courant":   encoding courant constant (between 0 and 1)
   *    "fangle":  feature face angle, between 0 and 180.
   *    "fmode":   "auto" (propagate calls f-d automatically), and "manual"
   *    "smoother": mesh smoother.
   *    "rediter": redistribution iterations.
   *    "verbose": verbose level.
   */
  void set_option( const char *opt, const char *val);

  static void load( const std::string &mname);

  static void unload( const std::string &mname);

  // Obtain a reference to the manifold.
  virtual Manifold *manifold() 
  { return _wm?_wm:(_parent?_parent->manifold():NULL); }

  int validate_object() const {
    if ( _cookie != PROP_COOKIE) return -1;
    else return Rocsurf::validate_object();
  }
  
protected:
  enum { PROP_COOKIE=7627873};

  SURF::Rocsurf          *_parent;   // Parent Rocsurf object
  Propagation_3          *_prop;     // Propagation object
  COM::Window            *_win;      // User window.
  COM::Window            *_buf;      // A buffer window that clones the 
                                     // coordinates but use nc and pconn.
  Remesher_base          *_rem;      // Base for remeshing

  int                     _rem_owner;// Whether Rocprop owns the remesher
  int                     _prop_method;

  char                    _wf_expn;        //< Whether expansion is wavefrontal
  char                    _wght_scheme;    //< Weighting scheme
  char                    _nrm_dfsn;       //< Whether diffuse normals
  char                    _feature_layer;  //< Whether to generate feature layer
  double                  _eig_thres;  //< Threshold for primary-null psace
  double                  _courant;    //< Constant in front of CFL condition (0<c<1)
  double                  _fangle_strong;  //< Feature face angle.
  double                  _fangle_weak;    //< Feature face angle.
  double                  _fangle_turn;    //< Feature turning angle.

  int                     _verb;       // Verbose level
  const COM::DataItem   *_cnstr_types;// Types of nodal constraints
  const COM::DataItem   *_cnstr_bound;// cylindrical bound

  double                  _time_lb;    // Lowerbound of relative time step
  int                     _smoother;   // Choice of smoother
  int                     _rediter;    // Redistribution iterations
  int                     _conserv;    // Whether or not to conserve mass
  int                     _rank;       // process rank.
  
  int                     _cookie;
};

PROP_END_NAMESPACE

#endif // _ROCPROP_H_






