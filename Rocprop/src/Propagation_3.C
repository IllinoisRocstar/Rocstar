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
// $Id: Propagation_3.C,v 1.17 2008/12/06 08:45:27 mtcampbe Exp $

#include "Propagation_3.h"
#include "Rocblas.h"
#include "Generic_element_2.h"

PROP_BEGIN_NAMESPACE

Propagation_3::Propagation_3( Manifold *wm, COM::Window *buf) 
  : _surf(wm), _buf(buf), _cnstr_set(false), _bnd_set(0), 
    _cnstr_nodes(NULL), _cnstr_bound(NULL) {

  // Depending on whether we have pconn for ghost nodes and elements,
  // run in either WHOLE_PANE or ACROSS_PANE mode.
  if ( wm->pconn_nblocks() >=3)
    _mode = SURF::WHOLE_PANE;
  else
    _mode = SURF::ACROSS_PANE;

  if ( _buf) {
    _cnstr_nodes = _buf->new_dataitem( "PROP_cnstr_nodes", 'n', COM_INT, 1, "");
    _cnstr_faces = _buf->new_dataitem( "PROP_cnstr_faces", 'e', COM_INT, 1, "");
    _cnstr_bndry_edges = _buf->new_dataitem( "PROP_cnstr_bndry_edges", 'e', COM_CHAR, 1, "");
    _cnstr_bndry_nodes = _buf->new_dataitem( "PROP_cnstr_bndry_nodes", 'n', COM_CHAR, 1, "");
    _cnstr_bound = _buf->new_dataitem( "PROP_cnstr_bound", 'p', COM_DOUBLE, BOUND_LEN, "");

    _buf->panes( _panes); 
  }

  MPI_Comm comm = _buf->get_communicator();
  if ( COMMPI_Initialized())  _rank = COMMPI_Comm_rank(comm);
  else _rank = 0;
}

void Propagation_3::
set_bounds( const COM::DataItem *bnd) {
  if ( _buf == NULL) return;

  if ( bnd == NULL) {
    // Deallocate cnstr_bound if constraints are unset
    _buf->dealloc_array( "PROP_cnstr_bound", 0);
    _bnd_set = false;
  }
  else {
    COM_assertion_msg( COM_compatible_types( bnd->data_type(), COM_DOUBLE) &&
		       bnd->size_of_components()==BOUND_LEN && 
		       bnd->is_panel() || bnd->is_windowed(),
		       "Propagation_3::set_bounds expects double-precision 3-vector");

    _buf->inherit( const_cast<COM::DataItem*>(bnd), "PROP_cnstr_bound",
		   COM::Pane::INHERIT_CLONE, true, NULL, 0);
    _bnd_set = true;
  }
  
  _buf->init_done( false);
}

void Propagation_3::
set_constraints( const COM::DataItem *ct_in) {
  if ( _buf == NULL) return;

  if ( ct_in == NULL) {
    // Deallocate cnstr_types if constraints are unset
    _buf->dealloc_array( "PROP_cnstr_nodes", 0);
    _buf->dealloc_array( "PROP_cnstr_faces", 0);
    _buf->dealloc_array( "PROP_cnstr_bndry_edges", 0);
    _buf->dealloc_array( "PROP_cnstr_bndry_nodes", 0);
    _cnstr_set = false;
  }
  else {
    COM_assertion_msg( ct_in->is_panel() || ct_in->is_elemental(), 
		       "Constraints must be associated with panes or faces.");

    // Copy the constraints onto the window
    _buf->resize_array( _cnstr_faces, 0);
    Rocblas::copy( ct_in, _cnstr_faces);   // Set facal flags

    // Determine the boundary of patches with constraints
    _buf->resize_array( _cnstr_bndry_edges, 0);
    _buf->resize_array( _cnstr_bndry_nodes, 0);
    determine_constraint_boundary( _cnstr_faces, _cnstr_bndry_edges, _cnstr_bndry_nodes);

    // Determine nodal constraints.
    _buf->resize_array( _cnstr_nodes, 0);

    // Convert constraints from facial/panel to nodal
    convert_constraints( _cnstr_faces, _cnstr_nodes);
    _cnstr_set = true;
  }
}

inline bool matchall( int a, int b) { return (a & b) == b; }
inline bool matchany( int a, int b) { return (a & b); }

void Propagation_3::
convert_constraints( const COM::DataItem *ctypes_faces,
		     COM::DataItem *ctypes_nodes) {
  int zero = 0;
  Rocblas::copy_scalar( &zero, ctypes_nodes);

  // Loop through all the nodes to compute their displacements
  enum { FIXED=1, XDIR = 2, YDIR = 4, ZDIR = 8, 
	 NOX = 16, NOY = 32, NOZ = 64};

  std::vector< COM::Pane*>::iterator it = _panes.begin(), iend=_panes.end();
  // Loop through the panes and its real nodes
  for ( it=_panes.begin(); it!=iend; ++it) {
    COM::Pane *pane = *it;

    COM_assertion( pane->dataitem(COM_NC)->stride()==3);
    const Point_3 *pnts = 
      reinterpret_cast<const Point_3*>(pane->coordinates()); 
    
    const int *val_faces=(const int*)pane->dataitem(ctypes_faces->id())->pointer();
    int *val_nodes=(int*)pane->dataitem(ctypes_nodes->id())->pointer();

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1); 
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      int type = val_faces[j];

      int ne = ene.size_of_edges();
      // Loop through all vertices of current face.
      for ( int k=0; k<ne; ++k) {
	int vindex = ene[k]-1;

	switch ( type) {
	case 0: break;
	case 2:    val_nodes[vindex] |= FIXED; break;
	case 'x':  val_nodes[vindex] |= XDIR; break;
	case -'x': val_nodes[vindex] |= NOX;   break;
	case 'y':  val_nodes[vindex] |= YDIR; break;
	case -'y': val_nodes[vindex] |= NOY;   break;
	case 'z':  val_nodes[vindex] |= ZDIR; break;
	case -'z': val_nodes[vindex] |= NOZ;   break;
	case 't':
	case -'t': {
	  SURF::Generic_element_2 e(ene.size_of_edges(), ene.size_of_nodes());
	  COM::Element_node_vectors_k_const<Point_3 > ps;
	  ps.set( pnts, ene, 1);

	  Vector_2 nc(0.5,0.5);
	  Vector_3 J[2], nrm;
	  e.Jacobian( ps, nc, J);
	  nrm = Vector_3::cross_product( J[0], J[1]).normalize();
	  
	  double eps = 1.e-1;
	  if ( std::abs(std::abs(nrm[0])-1) < eps) {
	    if ( type > 0) val_nodes[vindex] |= NOX;
	    else val_nodes[vindex] |= XDIR;
	  }
	  else if ( std::abs(std::abs(nrm[1])-1) < eps) {
	    if ( type > 0) val_nodes[vindex] |= NOY;
	    else val_nodes[vindex] |= YDIR;
	  }
	  else if ( std::abs(std::abs(nrm[2])-1) < eps) {
	    if ( type > 0) val_nodes[vindex] |= NOZ;
	    else val_nodes[vindex] |= ZDIR;
	  }
	  else
	    COM_assertion_msg(false, "Plane must be along axial direction");
	  break;
	}
	default: 
	  COM_assertion_msg( false, "Unknown/unsupported type of constraint");
	} // switch
      } // for k
    } // for j
  } // for it

  _surf->reduce_on_shared_nodes( ctypes_nodes, Manifold::OP_BOR);

  // Loop through the panes and its real nodes
  for ( it=_panes.begin(); it!=iend; ++it) {
    COM::Pane *pane = *it;
    int *val_nodes = (int*)pane->dataitem( ctypes_nodes->id())->pointer();

    for ( int j=0, nj=pane->size_of_real_nodes(); j<nj; ++j) {
      int &type = val_nodes[j];

      if ( type == 0) continue; // no constrained

      if ( (type & FIXED) || matchall(type, NOX|NOY|NOZ) ||
	   matchall(type, XDIR|YDIR) || matchall(type, XDIR|ZDIR) ||
	   matchall(type, YDIR|ZDIR) )
	type = 2;
      else if ( matchall(type, NOX|NOY)) {
	if ( matchany(type, XDIR|YDIR)) { type = 2; }
	else { type = 'z'; }
      }
      else if ( matchall(type, NOX|NOZ)) {
	if ( matchany(type, XDIR|ZDIR)) { type = 2; }
	else { type = 'y'; }
      }
      else if ( matchall(type, NOY|NOZ)) {
	if ( matchany(type, YDIR|ZDIR)) { type = 2; }
	else { type = 'x'; }
      }
      else if ( type & NOX) {
	if ( type & XDIR) { type = 2; }
	else if ( type & YDIR) { type = 'y'; }
	else if ( type & ZDIR) { type = 'z'; }
	else type = -'x';
      }
      else if ( type & NOY) {
	if ( type & YDIR) { type = 2; }
	else if ( type & XDIR) { type = 'x'; }
	else if ( type & ZDIR) { type = 'z'; }
	else type = -'y';
      }
      else if ( type & NOZ) {
	if ( type & ZDIR) { type = 2; }
	else if ( type & XDIR) { type = 'x'; }
	else if ( type & YDIR) { type = 'y'; }
	else type = -'z';
      }
      else if ( type & XDIR)
	type = 'x';
      else if ( type & YDIR)
	type = 'y';
      else {
	COM_assertion( type & ZDIR);
	type = 'z'; 
      }
    } // for j
  } // for it
}

/// Convert facial or panel constraints to nodal constraints
void Propagation_3::
determine_constraint_boundary( const COM::DataItem *ctypes_faces,
			       COM::DataItem *ctypes_bndry_edges,
			       COM::DataItem *ctypes_bndry_nodes) {
  _surf->update_bd_flags( ctypes_faces); // Update facal flags along boundaries

  char zero = 0;
  Rocblas::copy_scalar( &zero, ctypes_bndry_edges);
  Rocblas::copy_scalar( &zero, ctypes_bndry_nodes);

  Manifold::PM_iterator pm_it=_surf->pm_begin();
  std::vector< COM::Pane*>::iterator it = _panes.begin(), iend=_panes.end();

  // Loop through the panes and its real elements
  for ( it=_panes.begin(); it!=iend; ++it, ++pm_it) {
    COM::Pane *pane = *it;

    const int *val_faces=(const int*)pane->dataitem(ctypes_faces->id())->pointer();
    char *val_bndry_edges=(char*)pane->dataitem(ctypes_bndry_edges->id())->pointer();
    char *val_bndry_nodes=(char*)pane->dataitem(ctypes_bndry_nodes->id())->pointer();

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1);
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      for ( int k=0, ne=ene.size_of_edges(); k<ne; ++k) {
	Edge_ID eid(j+1,k), eid_opp = (*pm_it)->get_opposite_real_edge( eid);
	
	if (!(*pm_it)->is_physical_border_edge( eid_opp)) {
	  int flag_opp = (eid_opp.is_border() ? 
			  (*pm_it)->get_bd_flag( eid_opp) : val_faces[eid_opp.eid()-1]);
	  if ( val_faces[j] != flag_opp) {
	    val_bndry_edges[j] |= (1<<k);
	    val_bndry_nodes[ene[k]-1] = 1;
	  }
	}
      } // for k
    } // for j
  } // for it

  _surf->reduce_on_shared_nodes( ctypes_bndry_nodes, Manifold::OP_MAX);
}

bool Propagation_3::
check_spherical_bound( const Point_3 &pnt, const Point_3 &org, 
		       const double rad_min, const double rad_max, 
		       double eps) {
  double sqnrm = (pnt - org).squared_norm();

  return sqnrm >= square( rad_max-eps) || 
    rad_min >=0 && sqnrm <= square( rad_min+eps);
}

void Propagation_3::
bound_spherical_disp( const Point_3 &pnt, const Point_3 &org, 
		      const double rad_min, const double rad_max, 
		      Vector_3 &disp, double eps) {
  Vector_3 dir = pnt - org;
  double sqnrm = dir.squared_norm(), sqnrm1 = (dir+disp).squared_norm(), t;

  if ( std::max( sqnrm, sqnrm1) >= square((t=rad_max)-eps) && 
       std::min( sqnrm, sqnrm1) <= square(rad_max+eps) ||
       (t=rad_min) >= 0 && std::min( sqnrm, sqnrm1) <= square(rad_min+eps) && 
       std::max( sqnrm, sqnrm1) >= square(rad_min-eps)) {

    if ( dir.is_null()) 
      disp[0] = disp[1] = disp[2] = 0;
    else {
      Vector_3 dir_unit; (dir_unit = dir).normalize();
      disp += ( t - disp*dir_unit) * dir_unit - dir;
    }
  }
  else if ( sqnrm >= square(rad_max) && sqnrm1 > sqnrm || 
	    rad_min >= 0 && sqnrm <= square(rad_min) && sqnrm1 < sqnrm) {
    if ( dir.is_null()) 
      disp[0] = disp[1] = disp[2] = 0;
    else
      // Only move tangentially
      disp -= disp * dir / dir.squared_norm();
  }
}

bool Propagation_3::check_radial_bound( const double x, const double y,
					const double bnd_min, 
					const double bnd_max,
					double eps) {

  double sq_xy = square( x)+square( y);
  return  sq_xy >= square(bnd_max-eps) || 
    bnd_min >=0 && sq_xy <= square(bnd_min+eps);
}

void Propagation_3::bound_radial_disp( const double x, const double y,
                                       const double bnd_min, const double bnd_max,
                                       double &dx, double &dy, double eps) {

  double sq_xy = square( x)+square( y);
  double sq_rad_actual=square( x+dx)+square( y+dy);
  double t;

  if ( std::max( sq_xy, sq_rad_actual) >= square((t=bnd_max)-eps) && 
       std::min( sq_xy, sq_rad_actual) <= square(bnd_max+eps) ||
       std::min( sq_xy, sq_rad_actual) <= square((t=bnd_min)+eps) && 
       std::max( sq_xy, sq_rad_actual) >= square(bnd_min-eps)) {

    if ( sq_xy == 0)
      dx = dy = 0;
    else {
      // Try to propagate toward the cylinder
      double ratio = t / std::sqrt( sq_xy);
      dx = ratio*x - x;
      dy = ratio*y - y;
    }
  }
  else if ( sq_xy >= square(bnd_max-eps) && sq_rad_actual >= sq_xy || 
	    bnd_min>=0 && sq_xy <= square(bnd_min+eps) && sq_rad_actual<=sq_xy) {
    dx = dy = 0; // Do not alter dx if moves inwards
  }

}

bool Propagation_3::check_axial_bound( const double x, const double bnd_min, 
				       const double bnd_max, double eps) {
  return x >= bnd_max-eps || x <= bnd_min+eps;
}

void Propagation_3::bound_axial_disp( const double x, const double bnd_min,
                                      const double bnd_max, double &dx, 
				      double eps) {

  double xnew = x+dx;
  if ( std::max( x, xnew) >= bnd_max-eps && std::min(x, xnew) <= bnd_max+eps) {
    dx = bnd_max - x;
  }
  else if ( std::min( x, xnew) <= bnd_min+eps && std::max(x, xnew) >= bnd_min-eps) {
    dx = bnd_min - x;
  }
  else if ( x >= bnd_max+eps && xnew > x || x<=bnd_min-eps && xnew < x)
    dx = 0; // Do not alter dx if moves inwards
}

void Propagation_3::
bound_nodal_motion( const Point_3 &pnt, const double *bnd, 
		    Vector_3 &du, double eps) {
  if ( !in_bounding_box( pnt, (const Point_3&)bnd[6], (const Point_3&)bnd[9]))
    return;

  switch ( (int)bnd[0]) {
  case 's':  bound_spherical_disp( pnt, (const Point_3&)bnd[3],
				   bnd[1], bnd[2], du, eps); break;

  case 'x':  bound_axial_disp( pnt[0]-bnd[3], bnd[1], bnd[2], du[0], eps); break;
  case 'y':  bound_axial_disp( pnt[1]-bnd[4], bnd[1], bnd[2], du[1], eps); break;
  case 'z':  bound_axial_disp( pnt[2]-bnd[5], bnd[1], bnd[2], du[2], eps); break;

  case -'x':  bound_radial_disp( pnt[1]-bnd[4], pnt[2]-bnd[5], bnd[1], bnd[2], 
				 du[1], du[2], eps); break;
  case -'y':  bound_radial_disp( pnt[0]-bnd[3], pnt[2]-bnd[5], bnd[1], bnd[2], 
				 du[0], du[2], eps); break;
  case -'z':  bound_radial_disp( pnt[0]-bnd[3], pnt[1]-bnd[4], bnd[1], bnd[2], 
				 du[0], du[1], eps); break;
    
  default: COM_assertion_msg( false, "Unknown type of bnd");
  }
}

bool Propagation_3::
reached_nodal_bound( const Point_3 &pnt, const double *bnd, double eps) {
  if ( !in_bounding_box( pnt, (const Point_3&)bnd[6], (const Point_3&)bnd[9]))
    return false;

  switch ( (int)bnd[0]) {
  case 's':  return check_spherical_bound( pnt, (const Point_3&)bnd[3],
					   bnd[1], bnd[2], eps);

  case 'x':  return check_axial_bound( pnt[0]-bnd[3], bnd[1], bnd[2], eps); 
  case 'y':  return check_axial_bound( pnt[1]-bnd[4], bnd[1], bnd[2], eps); 
  case 'z':  return check_axial_bound( pnt[2]-bnd[5], bnd[1], bnd[2], eps); 

  case -'x': return check_radial_bound( pnt[1]-bnd[4], pnt[2]-bnd[5], bnd[1], bnd[2], eps); 
  case -'y': return check_radial_bound( pnt[0]-bnd[3], pnt[2]-bnd[5], bnd[1], bnd[2], eps); 
  case -'z': return check_radial_bound( pnt[0]-bnd[3], pnt[1]-bnd[4], bnd[1], bnd[2], eps); 
    
  default: COM_assertion_msg( false, "Unknown type of bnd");
  }

  return false;
}

void Propagation_3::
enforce_nodal_constraint( int type, Vector_3 &du) {
  switch ( type) {
  case 0: break;
  case 2:    du = Vector_3(0.); break;
  case 'x':  du[1] = du[2] = 0.; break;
  case -'x': du[0] = 0.; break;
  case 'y':  du[0] = du[2] = 0.; break;
  case -'y': du[1] = 0.; break;
  case 'z':  du[0] = du[1] = 0.; break;
  case -'z': du[2] = 0.; break;
  default: COM_assertion_msg( false, "Unknown type of constraint");
  }
}

void Propagation_3::
get_constraint_directions( int type, int &ndirs, Vector_3 dirs[2]) {
  switch ( type) {
  case 0: 
    ndirs = 0;
    break;
  case 2: 
    ndirs = 3;
    break;
  case 'x':  
    ndirs = 1; 
    dirs[0] = Vector_3(1,0,0); 
    break;
  case -'x': 
    ndirs = 2; 
    dirs[0] = Vector_3(0,1,0); 
    dirs[1] = Vector_3(0,0,1); 
    break;
  case 'y': 
    ndirs = 1; 
    dirs[0] = Vector_3(0,1,0); 
    break;
  case -'y':
    ndirs = 2;
    dirs[0] = Vector_3(1,0,0); 
    dirs[1] = Vector_3(0,0,1); 
    break;
  case 'z': 
    ndirs = 1; 
    dirs[0] = Vector_3(0,0,1); 
    break;
  case -'z':
    ndirs = 2;
    dirs[0] = Vector_3(1,0,0); 
    dirs[1] = Vector_3(0,1,0); 
    break;
  default: COM_assertion_msg( false, "Unknown type of constraint");
  }
}

// Enforce nodal constraints.
void Propagation_3::
enforce_nodal_constraints( COM::DataItem *du) {
  if ( !_cnstr_set) return;

  std::vector< COM::Pane*>::iterator it = _panes.begin(), iend=_panes.end();
  // Loop through the panes and its real nodes
  for ( it=_panes.begin(); it!=iend; ++it) {
    COM::Pane *pane = *it;
    
    const int *types = reinterpret_cast<const int*>
      (pane->dataitem(_cnstr_nodes->id())->pointer());
    Vector_3 *ds = reinterpret_cast<Vector_3*>
      (pane->dataitem( du->id())->pointer());

    for (int i=0, n=pane->size_of_real_nodes(); i<n; ++i) {
      enforce_nodal_constraint( types[i], ds[i]);
    }
  }
}

// Enforce nodal constraints.
void Propagation_3::
bound_nodal_motion( COM::DataItem *du) {
  if ( !_bnd_set) return;

  std::vector< COM::Pane*>::iterator it = _panes.begin(), iend=_panes.end();
  // Loop through the panes and its real nodes
  for ( it=_panes.begin(); it!=iend; ++it) {
    COM::Pane *pane = *it;
    
    Vector_3 *ds = reinterpret_cast<Vector_3*>
      (pane->dataitem( du->id())->pointer());

    const Point_3 *pnts = reinterpret_cast<const Point_3*>
      (pane->dataitem(COM_NC)->pointer());
    const double *bnd = _bnd_set ? reinterpret_cast<const double*>
      ( pane->dataitem( _cnstr_bound->id())->pointer()) : NULL;
    int nbnd = _bnd_set ? pane->dataitem( _cnstr_bound->id())->size_of_items() : 0;

    // Bound all nodes
    for ( int i=0; i<nbnd; ++i, bnd+=BOUND_LEN) {
      double eps = 1.e-5 * std::min( std::min( bnd[9]-bnd[6], bnd[10]-bnd[7]),
				     bnd[11]-bnd[8]);
      if ( eps>1 || eps <= 0) eps = 1.e-10;

      for (int i=0, n=pane->size_of_real_nodes(); i<n; ++i)
	bound_nodal_motion( pnts[i], bnd, ds[i]);
    }
  }
}

// Bound a face based on nodal constraints. If all the nodes are beyond 
// the bounded, then set the value to zero.
void Propagation_3::
bound_facial_speed( COM::DataItem *fa) {
  COM_assertion( fa->window() == _buf); 

  COM_assertion_msg( _bnd_set, "Bound must be set first before calling bound_facial_speed");

  COM_assertion_msg( fa->is_elemental() && fa->size_of_components()==1 &&
		     COM_compatible_types( fa->data_type(), COM_DOUBLE),
		     "Propagation_3::bound_facial_speed expects a scalar double-precision elemental dataitem.");

  std::vector< COM::Pane*>::iterator it = _panes.begin(), iend=_panes.end();
  // Loop through the panes and its real nodes
  for ( it=_panes.begin(); it!=iend; ++it) {
    COM::Pane *pane = *it;
    
    double *val_dbl = (double*)pane->dataitem(fa->id())->pointer();
    // Skip the panes without the dataitem.
    if ( val_dbl==NULL) continue;

    COM_assertion( pane->dataitem(COM_NC)->stride()==3);
    const Point_3 *pnts = reinterpret_cast<const Point_3*>(pane->coordinates()); 
    int nbnd = _bnd_set ? pane->dataitem( _cnstr_bound->id())->size_of_items() : 0;
    if ( nbnd==0) continue;

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1);
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      // Skip the faces that are already initialized to 0.
      if ( val_dbl[j]==0) continue;
      int ne = ene.size_of_edges();

      const double *bnd = reinterpret_cast<const double*>
	( pane->dataitem( _cnstr_bound->id())->pointer());
      bool out[] = {false, false, false, false};
      // Bound all nodes
      for ( int i=0; i<nbnd; ++i, bnd+=BOUND_LEN) {
	double eps = 1.e-5 * std::min( std::min( bnd[9]-bnd[6], bnd[10]-bnd[7]),
				       bnd[11]-bnd[8]);
	if ( eps>1 || eps <=0 ) eps = 1.e-10;

	for ( int k=0; k<ne; ++k) {
	  out[k] |= reached_nodal_bound( pnts[ene[k]-1], bnd, eps);
	}
      }
      
      if ( out[0] && out[1] && out[2] && (ne==3 || out[3])) val_dbl[j] = 0;
    }
  }
}

PROP_END_NAMESPACE






