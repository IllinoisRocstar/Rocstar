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
// $Id: FaceOffset_3.C,v 1.34 2008/12/06 08:45:27 mtcampbe Exp $

#include "FaceOffset_3.h"
#include "Rocblas.h"
#include "Generic_element_2.h"
#include "Rocsurf.h"
#include <algorithm>

PROP_BEGIN_NAMESPACE

//============== Debug features =======
// Export debugging information to the user window, if user has
// the registered "back-door" dataitems to get debug info such 
// as normals. It should be set to true for development purpose.
#define EXPORT_DEBUG_INFO     1

//=================================================

const double FaceOffset_3::pi= 3.14159265358979;

// Set parameters to default values
void FaceOffset_3::reset_default_parameters( bool feature_preserve) {

  // Two default feature modes:
  if ( feature_preserve) {
    // Feature-preserving mode: (Recommended for structured meshes)
    _wght_scheme = SURF::E2N_ANGLE;
    _tol_angle_weak = 18.*pi/180;
  }
  else {
    // Feature-smearing mode: (Recommended for unstructured meshes)
    _wght_scheme = SURF::E2N_AREA;
    _tol_angle_weak = 30.*pi/180.;
  }

  _dir_thres_weak = square(tan(0.5*_tol_angle_weak)); // 0.003;
  _tol_angle_strong = 65.*pi/180;
  _dir_thres_strong = 0.7;

  _tol_kstrong    = 5;
  _tol_kangle     = (49./180.)*pi;
  _tol_eangle     = (25/180.)*pi;

  _tol_turn_angle = pi/4.5;
  _tol_ad = pi/3;
  _tol_cos_osta   = std::cos(pi/9.);

  _nrm_dfsn = 1; _eig_thres = 0.003;
  _eig_weak_lb = 1.e-4; _eig_weak_ub = 0.025;

  _courant = 0.5; _inv_sqrt_c = 1/std::sqrt(_courant);
  _smoother = SMOOTHER_LAPLACIAN; _conserv = 0; 
  _wf_expn = true; _feature_layer=0;
}

// Construct an object from a window manifold and a buffer window.
// Extend the buffer window to store eigenvalues, eigenvectors, 
// tangent ranks, rescaling factors, and weights.
FaceOffset_3::FaceOffset_3( Manifold *wm, COM::Window *buf)  
  : Propagation_3( wm, buf)
{
  // Determine whether the mesh is structured
  _is_strtd = true;
  std::vector< COM::Pane*>::iterator it = _panes.begin();
  for (int i=0, local_npanes = _panes.size(); i<local_npanes; ++i, ++it) { 
    COM::Pane *pane = *it;
    if ( !pane->is_structured()) _is_strtd = false;
  }
  if ( COMMPI_Initialized()) {
    int s = _is_strtd;
    MPI_Allreduce( &s, &_is_strtd, 1, MPI_INT, MPI_MIN, 
		   _buf->get_communicator());
  }

  // Set parameters to default values. 
  reset_default_parameters( _is_strtd);

  // Extend buffer window
  _facenormals = _buf->new_dataitem( "facenormals", 'e', COM_DOUBLE, 3, "");
  _buf->resize_array( _facenormals, 0);

  _facecenters = _buf->new_dataitem( "facecenters", 'e', COM_DOUBLE, 3, "");
  _buf->resize_array( _facecenters, 0);

  _faceareas = _buf->new_dataitem( "faceareas", 'e', COM_DOUBLE, 1, "");
  _buf->resize_array( _faceareas, 0);

  // This is used to store v-center for smoothing offset.
  _vcenters = _buf->new_dataitem( "vcenters", 'n', COM_DOUBLE, 3, "");
  _buf->resize_array( _vcenters, 0);

  _eigvalues = _buf->new_dataitem( "lambda", 'n', COM_DOUBLE, 3, "");
  _buf->resize_array( _eigvalues, 0);

  _vnormals = _buf->new_dataitem( "vnormals", 'n', COM_DOUBLE, 3, "");
  _buf->resize_array( _vnormals, 0);

  _ridges = _buf->new_dataitem( "ridges", 'p', COM_INT, 2, "");
  
  _ridgeneighbors = _buf->new_dataitem( "ridgeneighbors", 'n', COM_INT, 4, "");
  _buf->resize_array( _ridgeneighbors, 0);

  _weak = _buf->new_dataitem( "weakfeature", 'n', COM_CHAR, 1, "");
  _buf->resize_array( _weak, 0);

  _strong = _buf->new_dataitem( "strongfeature", 'n', COM_CHAR, 1, "");
  _buf->resize_array( _strong, 0);

  _As = _buf->new_dataitem( "As", 'n', COM_DOUBLE, 9, "");
  _buf->resize_array( _As, 0);

  _boffset = _buf->new_dataitem( "boffset", 'n', COM_DOUBLE, 3, "");
  _buf->resize_array( _boffset, 0);

  _bmedial = _buf->new_dataitem( "bmedial", 'n', COM_DOUBLE, 3, "");
  _buf->resize_array( _bmedial, 0);

  _eigvecs = _buf->new_dataitem( "eigvecs", 'n', COM_DOUBLE, 9, "");
  _buf->resize_array( _eigvecs, 0);

  _tangranks = _buf->new_dataitem( "tangranks", 'n', COM_CHAR, 1, "");
  _buf->resize_array( _tangranks, 0);

  _ctangranks = _buf->new_dataitem( "ctangranks", 'n', COM_CHAR, 1, "");
  _buf->resize_array( _ctangranks, 0);

  _scales = _buf->new_dataitem( "scales", 'n', COM_DOUBLE, 1, "");
  _buf->resize_array( _scales, 0);

  _weights = _buf->new_dataitem( "weights", 'n', COM_DOUBLE, 1, "");
  _buf->resize_array( _weights, 0);

  _buf->init_done(false);
}

// Main entry of the algorithm
double FaceOffset_3::
time_stepping( const COM::DataItem *spds, double dt,
	       COM::DataItem *disps, int *smoothed) {

  COM_assertion_msg( spds || dt==0, 
		     "If no speed function, then time step must be zero.");

  // Inherit speeds and displacement onto the buffer window
  const COM::DataItem *spds_buf = spds ?
    ( (spds->window() == _buf) ? spds : 
      _buf->inherit( const_cast<COM::DataItem*>(spds), "speeds_inherited", 
		     COM::Pane::INHERIT_USE, true, NULL, 0)) : NULL;
  COM::DataItem *disps_buf = 
    _buf->inherit( disps, "disps_inherited", COM::Pane::INHERIT_USE, 
		   true, NULL, 0);
  COM::DataItem *disps_tmp =
    _buf->inherit( disps, "disps_inherited2", COM::Pane::INHERIT_CLONE, 
		   true, NULL, 0);
  _buf->init_done( false);

  // Compute the face areas and store them.
  SURF::Rocsurf::compute_element_areas( _faceareas);

  // First, compute the quadrics to determine the offset intersection 
  // with given time step. 
  bool with_nodal_velo = dt>0 && spds->is_nodal() && spds->size_of_components()==3;

  // If wavefrontal, then perform two setps.
  for ( int i=0; i<1+(_wf_expn && dt>0); ++i) {
    compute_quadrics( dt, spds_buf, _boffset, i ? disps_buf : NULL);

    if ( with_nodal_velo)
      Rocblas::mul_scalar( spds, &dt, disps_buf);
    else
      Rocblas::copy( _boffset, disps_buf);

    // Compute mean normals and offset directions during first iteration
    if ( i==0 && _wf_expn && dt>0) {
      compute_directions( !with_nodal_velo ? disps_buf : NULL, false);

      // Enforce constraints.
      obtain_constrained_directions( disps_buf, NULL);
    }
  }
  
  // Filter out isolated ridge vertices and identify ridge edges.
  filter_and_identify_ridge_edges(!_is_strtd && _wght_scheme==SURF::E2N_ANGLE && dt==0);
  mark_weak_vertices();

  // Recompute mean normals and offset directions with only primary component
  compute_directions( (dt!=0 && !with_nodal_velo) ? disps_buf : NULL, true);

  if ( dt==0)
    denoise_surface( NULL, disps_buf);
  else {
    double zero=0.;
    Rocblas::copy_scalar( &zero, disps_tmp);
    denoise_surface( disps_buf, disps_tmp);
    Rocblas::add( disps_buf, disps_tmp, disps_buf);
  }

  bool needs_smoothing=(!smoothed || *smoothed);
  if ( smoothed) *smoothed = true;

  if ( !needs_smoothing) {
    Rocblas::copy( disps_buf, _vcenters);
    if ( smoothed && *smoothed) *smoothed=false;
  }
  else if ( _is_strtd) {
    update_vertex_centers();  // Update vertex centers for ridge vertices
    nulaplacian_smooth( _vnormals, _tangranks, disps_buf,
			_vcenters, disps_tmp, _weights);
  }
  else if ( _smoother == SMOOTHER_ANISOTROPIC) {
    // Compute anisotropic vertex centers.
    compute_anisotropic_vertex_centers( disps_buf);
  }
  update_vertex_centers();  // Update vertex centers for ridge vertices

  // Recompute face normals if propagation under nodal velocity
  if ( with_nodal_velo) {
    _surf->compute_normals( _facenormals);
    _surf->update_bd_normals( _facenormals, false);
  }

  // Constrain the displacements and vertex-centers.
  obtain_constrained_directions( disps_buf, _vcenters);
  if (_bnd_set) {
    bound_nodal_motion( disps_buf);
    _surf->reduce_on_shared_nodes( disps_buf, Manifold::OP_MAXABS);
  }

#if EXPORT_DEBUG_INFO
  COM::DataItem *ctypes = disps->window()->dataitem( "cnstr_types_nodal");
  if ( _cnstr_set && ctypes) Rocblas::copy( _cnstr_nodes, ctypes);

  COM::DataItem *tranks = disps->window()->dataitem( "tangranks");
  if ( tranks) Rocblas::copy( _tangranks, tranks);

  COM::DataItem *facenormals = disps->window()->dataitem( "facenormals");
  if ( facenormals) Rocblas::copy( _facenormals, facenormals);
  
  COM::DataItem *facecenters = disps->window()->dataitem( "facecenters");
  if ( facecenters) Rocblas::copy( _facecenters, facecenters);
    
  COM::DataItem *eigvecs = disps->window()->dataitem( "eigvecs");
  if ( eigvecs) Rocblas::copy( _eigvecs, eigvecs);
  
  COM::DataItem *lambdas = disps->window()->dataitem( "lambdas");
  if ( lambdas) Rocblas::copy( _eigvalues, lambdas);

  COM::DataItem *ridges = disps->window()->dataitem( "ridges");
  if ( ridges) 
    disps->window()->inherit( _ridges, "ridges", COM::Pane::INHERIT_CLONE, 
			      true, NULL, 0);
 
  if ( dt>0) {
    COM::DataItem *disps_novis = disps->window()->dataitem( "disps_novis");
    if ( disps_novis) Rocblas::copy( disps, disps_novis);
  }
#endif

  // Second, reduce time step based on stability limit, and rescale
  // displacements by the reduced time step and wavefrontal expansion.
  double dt_sub = dt; 
  if ( _courant>0) {
    for ( ;;) {
      double scale = rescale_displacements( disps_buf, _vcenters);
      
      if ( scale==1) break;

      if ( dt==0 && _verb && _rank==0)
	std::cout << "Rocprop: Scaled back displacement with a factor of " 
		  << scale << std::endl;

      dt_sub *= scale;
      
      if ( !with_nodal_velo || _smoother == SMOOTHER_LAPLACIAN) break;
      // If nodal velocity
      Rocblas::mul_scalar( spds, &dt_sub, disps_buf);
      if (smoothed) *smoothed=false;
    }
  }
  else if ( dt==0) {
    Rocblas::copy( _vcenters, disps_buf);
  }

#if EXPORT_DEBUG_INFO
  if ( dt>0) {
    COM::DataItem *scales = disps->window()->dataitem( "scales");
    if ( scales) Rocblas::copy( _scales, scales);
  }
#endif

  _surf->reduce_on_shared_nodes( disps_buf, Manifold::OP_MAXABS);

  // Deallocate spds_buf and disps_buf
  _buf->delete_dataitem( disps_tmp->name());
  _buf->delete_dataitem( disps_buf->name());
  if ( spds_buf && spds_buf!=spds)
    _buf->delete_dataitem( spds_buf->name());
  _buf->init_done( false);

  // Finally, return the current time step.
  return dt_sub;
}

// This function computes the normal and center of the offset face. The 
// speed can be a nodal dataitem (scalar or 3-vector), or be an elemental
// dataitem. (scalar associated with face center, 3-vector associated
// with face center, or 3-vectors associated with three quadrature points.
// Only triangular meshes supported for latter two cases.)
void FaceOffset_3::
obtain_face_offset( const Point_3 *pnts, const double *spd_ptr,
		    const DataItem *spd, double dt, const Vector_3 *dirs, 
		    Element_node_enumerator &ene, 
		    Vector_3 &ns_nz, Point_3 &cnt, 
		    Vector_3 *disps, Vector_3 *ns) {

  const Vector_3 null_vec(0,0,0);
  int ne=ene.size_of_edges();

  // Evaluate current face normal 
  Element_node_vectors_k_const<Point_3> ps; ps.set( pnts, ene, 1);
  SURF::Generic_element_2 e(ne);
  Point_3 pnts_buf[9];

  // Compute new position for all vertices 
  int ncomp = spd?spd->size_of_components():0;
  COM_assertion_msg( !ncomp || ncomp==1 && spd->is_elemental() ||
		     ncomp==3 && spd->is_nodal(), 
		     "Speed must be scalar-elemental or 3-vector-nodal");

  // Copy coordinates into buffer.
  for ( int i=0; i<ne; ++i) pnts_buf[i] = ps[i];

  // If dirs is present, then compute the displacement
  if ( dt!=0 && ncomp==3) { // Add dt*speed to each point
    for ( int i=0; i<ne; ++i)
      pnts_buf[i] += disps[i] = dt*((const Vector_3*)spd_ptr)[ene[i]-1];
  }
  else if ( dt==0) {
    for ( int i=0; i<ne; ++i) disps[i] = null_vec;
  }

  // Offset face normal equal to current face normal 
  e.interpolate_to_center((Vector_3*)pnts_buf, (Vector_3*)&cnt);

  // Compute normal directions.
  Vector_2 nc(0.5,0.5);
  Vector_3 J[2];
  e.Jacobian( pnts_buf, nc, J);
  ns_nz = Vector_3::cross_product( J[0], J[1]).normalize();

  if ( ne==3) {
    ns[0] = ns[1] = ns[2] = ns_nz;
  }
  else {
    for ( int k=0; k<ne; ++k) {
      Vector_3 ts[] = { pnts_buf[(k+1)%ne]-pnts_buf[k], 
			pnts_buf[(k+ne-1)%ne]-pnts_buf[k]};
      ns[k] = Vector_3::cross_product( ts[0], ts[1]);

      // If the angle is far from 180 or 0 degrees, then use triangle's normal
      double l = ns[k].norm();
      if ( l > std::abs(0.05*ts[0]*ts[1])) ns[k] /= l;
      else ns[k] = ns_nz;
    }
  }

  // Compute points for normal motion.
  if (dt!=0 && ncomp==1) {
    double l = dt*spd_ptr[ene.id()-1];
    Vector_3 d = l*ns_nz;
    (Vector_3&)cnt += l*ns_nz; 

    if ( ne==3) {
      for ( int k=0; k<ne; ++k) pnts_buf[k] += (disps[k] = d);
    }
    else {
      for ( int k=0; k<ne; ++k) pnts_buf[k] += (disps[k] = l * ns[k]);
    }
  }

  // If wave-frontal motion, recalculate the directions and displacements.
  if ( dirs) {
    for ( int i=0; i<ne; ++i) {
      int vindex= ene[i]-1;
      
      if ( dirs[vindex].is_null()) {
	disps[i] = null_vec;
	pnts_buf[i] = ps[i] + disps[i];
      }
      else {
	Vector_3 dir = dirs[vindex];
	dir.normalize();

	Vector_3 tng = ( pnts_buf[(i+ne-1)%ne] - pnts_buf[i]) +
	  (pnts_buf[(i+1)%ne] - pnts_buf[i]);	

	if ( tng * dir < 0) { // Expansion
	  double l = (disps[i] * ns[i]);
	  if ( dir * ns[i]<0) l=-l;

	  (disps[i] = dir) *= l;

	  pnts_buf[i] = ps[i] + disps[i];
	}
      }
    }

    // Reevaluate normal directions.
    e.Jacobian( pnts_buf, nc, J);
    ns_nz = Vector_3::cross_product( J[0], J[1]).normalize();

    if ( ne==3) {
      ns[0] = ns[1] = ns[2] = ns_nz;
    }
    else {
      for ( int k=0; k<ne; ++k) {
	ns[k] = Vector_3::cross_product
	  ( pnts_buf[(k+1)%ne]-pnts_buf[k],
	    pnts_buf[(k+ne-1)%ne]-pnts_buf[k]).normalize();
      }
    }
  }
}

void FaceOffset_3::compute_directions( COM::DataItem *bo_attr, bool princ) {
  // Loop through the panes and its real vertices
  std::vector< COM::Pane*>::iterator it = _panes.begin();
  for (int i=0, local_npanes = _panes.size(); i<local_npanes; ++i, ++it) { 
    COM::Pane *pane = *it;

    // Obtain pointers
    Vector_3 *evs = reinterpret_cast<Vector_3*>
      ( pane->dataitem(_eigvecs->id())->pointer());
    const Vector_3 *lambdas = reinterpret_cast<const Vector_3*>
      ( pane->dataitem(_eigvalues->id())->pointer());
    Vector_3 *vnrms = reinterpret_cast<Vector_3*>
      ( pane->dataitem(_vnormals->id())->pointer());
    Vector_3 *voffset = bo_attr ? reinterpret_cast<Vector_3*>
      ( pane->dataitem(bo_attr->id())->pointer()) : NULL;
    const char *trank = princ ? reinterpret_cast<char*>
      ( pane->dataitem(_tangranks->id())->pointer()) : NULL;

    // Loop through all real nodes of the pane
    for ( int j=0, jn=pane->size_of_real_nodes(); j<jn; ++j) {
      obtain_directions( &evs[3*j], lambdas[j], princ ? 3-trank[j] : 3, 
			 vnrms[j], voffset?voffset+j:NULL);
    }
  }
}

// Obtain mean normal and offset direction at a vertex.
void FaceOffset_3::
obtain_directions( Vector_3 es[3], const Vector_3 &lambdas, int nrank, 
		   Vector_3 &b_medial, Vector_3 *b_offset) const
{
  // Use smaller threshold for acute angles.
  double eps = lambdas[0]*_eig_thres;
  
  if ( b_offset) {
    // Solve for x within primary space
    Vector_3 d_offset(0,0,0);
    for ( int k=0; k<nrank; ++k) {
      if ( lambdas[k]>eps) {
	d_offset += (es[k]* *b_offset/-lambdas[k])*es[k];
      }
    }

    *b_offset = d_offset;
  }

  // Solve for x within primary space
  Vector_3 d_medial(0,0,0);
  for ( int k=0; k<nrank; ++k) {
    if ( lambdas[k]>eps) {
      d_medial += (es[k] * b_medial/-lambdas[k])*es[k];
    }
  }

  b_medial = d_medial.normalize();
}

// Decompose propagation directions based on constraints
bool FaceOffset_3::
obtain_constrained_directions( COM::DataItem *disps_buf,
			       COM::DataItem *vcenters) {

  // Loop through the panes and its real faces
  std::vector< COM::Pane*>::iterator it = _panes.begin();
  Manifold::PM_iterator pm_it=_surf->pm_begin();
  
  const Vector_3 null_vec(0,0,0);

  // Loop through the panes and its real vertices
  for (int i=0, local_npanes = _panes.size(); i<local_npanes; ++i, ++it) { 
    COM::Pane *pane = *it;

    // Obtain pointers 
    Vector_3 *evs = reinterpret_cast<Vector_3*>
      ( pane->dataitem(_eigvecs->id())->pointer());
    Vector_3 *lambdas = reinterpret_cast<Vector_3*>
      ( pane->dataitem(_eigvalues->id())->pointer());
    Vector_3 *ds = disps_buf ? reinterpret_cast<Vector_3*>
      ( pane->dataitem(disps_buf->id())->pointer()) : NULL;
    const Vector_3 *bs = reinterpret_cast<Vector_3*>
      ( pane->dataitem(_boffset->id())->pointer());
    const Vector_3 *bms = reinterpret_cast<Vector_3*>
      ( pane->dataitem(_bmedial->id())->pointer());

    // Feature information
    const char *tranks = reinterpret_cast<char*>
      ( pane->dataitem(_tangranks->id())->pointer());
    Vector_3 *vcnts = vcenters ? reinterpret_cast<Vector_3*>
      ( pane->dataitem(vcenters->id())->pointer()) : NULL;
    const char *val_bndry_nodes = reinterpret_cast<const char*>
      ( pane->dataitem(_cnstr_bndry_nodes->id())->pointer());

    Vector_3 *As = reinterpret_cast<Vector_3*>
      ( pane->dataitem(_As->id())->pointer());
    const int *cnstrs = _cnstr_nodes ? reinterpret_cast<int*>
      ( pane->dataitem(_cnstr_nodes->id())->pointer()) : NULL;

    // Loop through all real nodes of the pane
    for ( int j=0, jn=pane->size_of_real_nodes(); j<jn; ++j) {
      // Skip vertices with zero right-hand side, such as those at
      // edge- and face-centers for quadratic elements
      if ( lambdas[j][0]==0) continue;

      Vector_3 *es_v = &evs[3*j];   // eigen-vectors of current vertex.
      int nrank = 3-tranks[j];

      Vector_3 V[2];
      int     ndirs = 0;
      if ( cnstrs) get_constraint_directions( cnstrs[j], ndirs, V);
      
      if ( ndirs==3 ) { // Fixed point
	lambdas[j] = es_v[0] = es_v[1] = es_v[2] = null_vec;

	if (ds) ds[j] = null_vec;
	if (vcnts) vcnts[j] = null_vec;
	continue;
      }

      const Vector_3 *A_v = &As[3*j];
      double eps = _eig_thres*lambdas[j][0];
      // Propagate along the constrained direction.
      double cos_fangle = 0.98481, sin_fangle = 0.17365; // (10 deg)
      
      if ( ds) {
	Vector_3 b(0, 0, 0); 
	for ( int k=0; k<nrank; ++k) {
	  if ( lambdas[j][k]>eps)
	    b += ds[j]*es_v[k]*es_v[k];
	}
	ds[j] = b;
      }

      if ( vcnts) { // Project tangent motion within null space
	Vector_3 t(0, 0, 0);
	for ( int k=1; k<3; ++k) {
	  if ( k>=nrank || lambdas[j][k]<eps)
	    // Note that vcnts may have component in domain space.
	    t += vcnts[j]*es_v[k]*es_v[k];
	}
	vcnts[j] = t;
      }
      if ( ndirs == 0) continue;

      // Evaluate surface normal
      Vector_3 nrm_v;
      if ( ds && ds[j]!=null_vec) {
	nrm_v = ds[j]; nrm_v.normalize(); 
      }
      else {
	if ( ds) ds[j] = null_vec;

	nrm_v = bms[j]; 
	obtain_directions( es_v, lambdas[j], tranks[j], nrm_v); 
      }

      // Propagate along the constrained direction.
      if ( ndirs == 2) {
	// Split planar constraints
	Vector_3 nrm = Vector_3::cross_product( V[0],V[1]);
	double 	det = std::abs(nrm*nrm_v);

      	if ( det < cos_fangle || val_bndry_nodes[j]) {
	  // Split into two directions.
	  V[0] = nrm_v; V[0] -= V[0]*nrm*nrm; V[0].normalize();
	  V[1] = Vector_3::cross_product( V[0], nrm);
	
	  if (vcnts) vcnts[j] = vcnts[j]*V[1]*V[1];
	}
	else if (vcnts) { // Plane is tangent space
	  vcnts[j] -= vcnts[j]*nrm*nrm;

	  if ( ds) ds[j] = null_vec; 
	  continue; 
	}
      }
      else if ( vcnts && !val_bndry_nodes[j] &&
		( tranks[j]==1 && std::abs(V[0]*es_v[2]) > cos_fangle ||
		  tranks[j]==2 && std::abs(V[0]*nrm_v) < sin_fangle)) {
	// ndir == 1
	vcnts[j] = vcnts[j]*V[0]*V[0];
	if ( ds) ds[j] = null_vec;
	continue;
      }

      if ( ds && !ds[j].is_null()) {
	double VAV=Vector_3(V[0]*A_v[0],V[0]*A_v[1],V[0]*A_v[2])*V[0];
	double Vb = V[0]*bs[j];

	if ( VAV<eps) 
	  ds[j] = null_vec;
	else {
	  ds[j] = (-Vb/VAV)*V[0];
	  if ( ndirs==1 && vcnts) vcnts[j] = null_vec;
	}
      }
    }
  }

  if ( _conserv && !disps_buf) balance_mass();

  return true;
}

// Reduce time step based on stability limit, and rescale displacements 
// by the reduced time step. Introduce dissipation.
// This function returns the reduced time step.
double FaceOffset_3::
rescale_displacements( COM::DataItem *disps, 
		       COM::DataItem *vcenters, int depth) {
  
  double alpha = HUGE_VAL;

  // Loop through the panes and its real faces
  std::vector< COM::Pane*>::iterator it = _panes.begin();
  Manifold::PM_iterator pm_it=_surf->pm_begin();
  for ( int i=0, local_npanes = _panes.size(); i<local_npanes; 
	++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;
    
    const Point_3 *pnts = reinterpret_cast<Point_3*>
      (pane->dataitem(COM_NC)->pointer());
    const Vector_3 *ds = reinterpret_cast<const Vector_3*>
      ( pane->dataitem(disps->id())->pointer());    
    const Vector_3 *vcnts = reinterpret_cast<const Vector_3*>
      ( pane->dataitem(vcenters->id())->pointer());
    const Vector_3 *fnrms = reinterpret_cast<Vector_3*>
      ( pane->dataitem(_facenormals->id())->pointer());
    const char      *tranks = reinterpret_cast<char*>
      ( pane->dataitem(_tangranks->id())->pointer());

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1); 

    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      int ne = ene.size_of_edges();
      int uindex=ene[ne-1]-1, vindex=ene[0]-1, windex=0;

      // Loop through all vertices of current face.
      for ( int k=0; k<ne; ++k, uindex=vindex, vindex=windex ) {
	windex = ene[(k+1==ne)?0:k+1]-1;

	if ( k!=0 && ne!=4 && (tranks[uindex] == 2 || tranks[vindex] == 2)) 
	  continue; 

	Vector_3 tngs[2] = { pnts[uindex]-pnts[vindex],
			     pnts[windex]-pnts[vindex]};

	Vector_3 disp_v  = ds[vindex]+vcnts[vindex];
	Vector_3 dsps[2] = { ds[uindex]+vcnts[uindex]-disp_v, 
			     ds[windex]+vcnts[windex]-disp_v };

	Vector_3 nrm = Vector_3::cross_product( tngs[0], tngs[1]);
	Vector_3 dt  = Vector_3::cross_product( dsps[0], tngs[1])+
	  Vector_3::cross_product( tngs[0], dsps[1]);
	Vector_3 dd = Vector_3::cross_product( dsps[0], dsps[1]);
  
	double a = dd*nrm;
	double b = dt*nrm;
	double c = nrm*nrm;

	std::pair<double, double>  roots=solve_quadratic_func(a,b,c);
	if ( roots.first<alpha && roots.first>0) alpha = roots.first;
	if ( roots.second<alpha && roots.second>0) alpha = roots.second;

	// Solve for time-step constraints at edges.
	Edge_ID eid(j+1,k), eid_opp = (*pm_it)->get_opposite_real_edge( eid);
	bool is_border_edge = (*pm_it)->is_physical_border_edge( eid_opp);

	if ( tranks[uindex] != 2 && tranks[vindex] != 2 && !is_border_edge) {
	  const Vector_3 &n1 = fnrms[j];
	  const Vector_3 &n2 = eid_opp.is_border() ? 
	    (*pm_it)->get_bd_normal( eid_opp) : fnrms[eid_opp.eid()-1];

	  nrm = n1+n2;
	  double a = dd*nrm;
	  double b = dt*nrm;
	  double c = nrm*nrm;

	  roots=solve_quadratic_func(a,b,c);
	  if ( roots.first<alpha && roots.first>0) alpha = roots.first;
	  if ( roots.second<alpha && roots.second>0) alpha = roots.second;
	}
      }
    }
  }

  // Reduce alpha on all processors.
  MPI_Comm comm = _buf->get_communicator();
  int  needs_rescale = alpha*_courant<=1;

  if ( COMMPI_Initialized() && COMMPI_Comm_size( comm)>0) {
    double local = alpha;
    MPI_Allreduce(&local, &alpha, 1, MPI_DOUBLE, MPI_MIN, comm);

    int local_nr = needs_rescale;
    MPI_Allreduce( &local_nr, &needs_rescale, 1, MPI_INT, MPI_LOR, comm);    
  }

  // We need to scale back displacement
  if ( _courant>0 && alpha*_courant<=1) {
    double scale = 0.9*_courant*alpha;

    Rocblas::mul_scalar( disps, &scale, disps);
    Rocblas::mul_scalar( vcenters, &scale, vcenters);

    const int maxdepth=6;
    if ( depth>maxdepth) {
      std::cerr << "Rocprop: Mesh too distorted. Stopping" << std::endl;
      COM_assertion_msg( depth<=maxdepth, "Too many sub-iterations.");
      abort();
    }

    // Recursively call the function to recompute the scale.
    return scale*rescale_displacements( disps, vcenters, depth+1);
  }
  else {
    // Add tangential components to normal components.
    Rocblas::add( disps, vcenters, disps);
    
    return 1.;
  }
}

// Solve for equation a*t^2+b*t+c=0 and find the roots between 0 and 1.
std::pair<double,double>
FaceOffset_3::solve_quadratic_func( double a, double b, double c) {
  double max_abc = std::max( std::max(std::abs(a), std::abs(b)), std::abs(c));

  if ( max_abc > 0) 
  { a /= max_abc; b /= max_abc; c /= max_abc; }

  if ( a == 0) {
    if ( b == 0) {
      return std::make_pair(HUGE_VAL, HUGE_VAL);
    }

    return std::make_pair(-c/b, HUGE_VAL);   
  }
  else {
    double det = b*b - 4*a*c;
    if ( det < 0) { 
      return std::make_pair(HUGE_VAL, HUGE_VAL);
    }
    
    double sqrt_d = std::sqrt( det);

    if ( b>0) {
      double temp = -b - sqrt_d;
      return std::make_pair(2*c/temp, 0.5*temp/a);
    }
    else {
      double temp = -b + sqrt_d;
      return std::make_pair( 0.5*temp/a, 2*c/temp);
    }
  }
}

PROP_END_NAMESPACE






