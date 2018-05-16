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
// $Id: quadric_analysis.C,v 1.15 2008/12/06 08:45:28 mtcampbe Exp $

#include "FaceOffset_3.h"
#include "Rocblas.h"
#include "Rocsurf.h"
#include <algorithm>

PROP_BEGIN_NAMESPACE

// If predicted is present, then the motion is wavefrontal.
void FaceOffset_3::
compute_quadrics( double dt, const COM::DataItem *spds, 
		  COM::DataItem *bo_attr, COM::DataItem *predicted) {

  COM_assertion_msg( spds || dt==0, 
		     "If no speed function, then time step must be zero.");

  const double zero = 0., eps = 1.e-100;
  if ( bo_attr) Rocblas::copy_scalar( &zero, bo_attr);
  if ( dt==0)   bo_attr = NULL;

  // Use the following buffer spaces: _As for matrix A; _eigvalues to 
  // store bm; _sumangles to store sum of angles
  COM::DataItem *A_attr = _As;
  COM::DataItem *bm_attr = _eigvalues;
  COM::DataItem *w_attr = _weights;

  Rocblas::copy_scalar( &zero, A_attr);
  Rocblas::copy_scalar( &zero, bm_attr);
  Rocblas::copy_scalar( &eps, w_attr);
  Rocblas::copy_scalar( &zero, _vcenters);

  COM::DataItem *sa_attr = 
    _buf->new_dataitem( "sumangles", 'n', COM_DOUBLE, 1, "");
  _buf->resize_array( sa_attr, 0);
  
  bool shear= _smoother==SMOOTHER_LAPLACIAN || _is_strtd ||
    _smoother!=SMOOTHER_ANISOTROPIC && (dt==0 || bo_attr==0);

  // Loop through the panes and its real faces
  std::vector< COM::Pane*>::iterator it = _panes.begin();
  Manifold::PM_iterator pm_it=_surf->pm_begin();
  for ( int i=0, local_npanes = _panes.size(); 
	i<local_npanes; ++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;

    // Obtain nodal coordinates of current pane, assuming contiguous layout
    const Point_3 *pnts = reinterpret_cast<Point_3*>
      (pane->dataitem(COM_NC)->pointer());
    Vector_3 *fnrms = reinterpret_cast<Vector_3*>
      ( pane->dataitem(_facenormals->id())->pointer());
    Point_3 *fcnts = reinterpret_cast<Point_3*>
      ( pane->dataitem(_facecenters->id())->pointer());

    double   *ws = reinterpret_cast<double*>
      ( pane->dataitem(w_attr->id())->pointer());
    Vector_3 *vcnts = reinterpret_cast<Vector_3*>
      ( pane->dataitem(_vcenters->id())->pointer());

    // Obtain the pointer to speed function
    const double *spds_ptr = spds?reinterpret_cast<const double*>
      (pane->dataitem( spds->id())->pointer()):NULL;

    // Obtain pointers
    Vector_3 *As = reinterpret_cast<Vector_3*>
      ( pane->dataitem(A_attr->id())->pointer());
    Vector_3 *bs_m = reinterpret_cast<Vector_3*>
      ( pane->dataitem(bm_attr->id())->pointer());
    Vector_3 *bs_o = bo_attr ? reinterpret_cast<Vector_3*>
      ( pane->dataitem(bo_attr->id())->pointer()) : NULL;
    Vector_3 *pred_dirs = predicted ? reinterpret_cast<Vector_3*>
      ( pane->dataitem(predicted->id())->pointer()) : NULL;
    double   *sa = reinterpret_cast<double*>
      ( pane->dataitem(sa_attr->id())->pointer());
    double   *areas = reinterpret_cast<double*>
      ( pane->dataitem(_faceareas->id())->pointer());

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1); 
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      // If speed is uniform within a face, then unit normal is the same
      // after propagation.
      Point_3  &cnt = fcnts[j];

      int ne = ene.size_of_edges();
      int uindex=ene[ne-1]-1, vindex=ene[0]-1, windex=0;

      Vector_3 disps[9], nrms[9];

      // Compute signed distance from current face.
      obtain_face_offset( pnts, spds_ptr, spds, dt, pred_dirs,
			  ene, fnrms[j], cnt, disps, nrms);

      double   delta=0;
      // Loop through all vertices of current face.
      for ( int k=0; k<ne; ++k, uindex=vindex, vindex=windex) {
	windex = ene[(k+1==ne)?0:k+1]-1;

	// Evaluate angle at each vertex. 
	Vector_3 ts[2] = { pnts[uindex]-pnts[vindex], 
			   pnts[windex]-pnts[vindex]};
	Vector_3 ns_nz = nrms[k];

	if ( bo_attr) { // Update delta
	  delta = -disps[k]*ns_nz;
	}

	double angle = eval_angle( ts[0], ts[1]);
	Vector_3 nrm = ns_nz;
	double w=1;
	if ( _wght_scheme!=SURF::E2N_ONE) {
	  switch ( _wght_scheme) {
	  case SURF::E2N_ANGLE: {
	    double s = std::sqrt((ts[0]*ts[0])*(ts[1]*ts[1]));
	    if ( s==0) w = 0; 
	    else w = angle; 

	    break;
	  }
	  case SURF::E2N_AREA: {
	    w = areas[j]; break;
	  }
	  default: ; // w=1;
	  }
	  nrm *= w;
	}

	// Obtain the addresses of the linear system and weight
	Vector_3 *A_v = &As[3*vindex];

	// Update A, b_o, b_m, and as_v for the current vertex
	A_v[0] += ns_nz[0]*nrm;
	A_v[1] += ns_nz[1]*nrm;
	A_v[2] += ns_nz[2]*nrm;

	if ( delta!=0) bs_o[vindex] += delta*nrm;
	bs_m[vindex] -= nrm;
	sa[vindex] += angle;

	if ( shear) {
	  Vector_3 diff = cnt-pnts[vindex];
	  ws[vindex] += 1.;
	  vcnts[vindex] += diff;
	}
      } // for k (edges)
    } // for j (elements)
  } // for i (panes)
  
  // Reduce on shared nodes for A, b_m, b_o, and r
  _surf->reduce_on_shared_nodes( A_attr, Manifold::OP_SUM);
  _surf->reduce_on_shared_nodes( bm_attr, Manifold::OP_SUM);
  if ( bo_attr) _surf->reduce_on_shared_nodes( bo_attr, Manifold::OP_SUM);
  _surf->reduce_on_shared_nodes( sa_attr, Manifold::OP_SUM);

  // Copy bm_attr into _vnormals and _bmedial
  Rocblas::copy( bm_attr, _vnormals);
  Rocblas::copy( bm_attr, _bmedial);

  // Reduce to maxabs so that shared nodes have exactly the same value
  _surf->reduce_on_shared_nodes( A_attr, Manifold::OP_MAXABS);
  _surf->reduce_on_shared_nodes( sa_attr, Manifold::OP_MAXABS);

  if ( shear) {
    _surf->reduce_on_shared_nodes( w_attr, Manifold::OP_SUM);
    _surf->reduce_on_shared_nodes( _vcenters, Manifold::OP_SUM);
    // Obtain angle-weighted average of mass center at each vertex.
    Rocblas::div( _vcenters, w_attr, _vcenters);
  }

  // Loop through the panes and its real vertices to determine trank
  it = _panes.begin(); 
  for (int i=0, local_npanes = _panes.size(); i<local_npanes; ++i, ++it) { 
    COM::Pane *pane = *it;

    // Obtain pointers
    Vector_3 *As = reinterpret_cast<Vector_3*>
      ( pane->dataitem(A_attr->id())->pointer());
    Vector_3 *evs = reinterpret_cast<Vector_3*>
      ( pane->dataitem(_eigvecs->id())->pointer());
    Vector_3 *bs_m = reinterpret_cast<Vector_3*>
      ( pane->dataitem(bm_attr->id())->pointer());
    double   *sa = reinterpret_cast<double*>
      ( pane->dataitem(sa_attr->id())->pointer());

    // Feature information
    char      *tranks = reinterpret_cast<char*>
      ( pane->dataitem(_tangranks->id())->pointer());

    // Loop through all real nodes of the pane
    for ( int j=0, jn=pane->size_of_real_nodes(); j<jn; ++j) {
      // Skip vertices with zero right-hand side, such as those at
      // edge- and face-centers for quadratic elements
      if ( bs_m[j].is_null()) continue;

      Vector_3 *A_v = &As[3*j];
      Vector_3 *es = &evs[3*j];    // eigen-vectors of current vertex.
      std::copy( A_v, A_v+3, es);
      
      // Perform eigen-analysis for the vertex using medial quadric
      // nrank is the rank of the normal (primary) space of A.
      char nrank = eigen_analyze_vertex( es, bs_m[j], 2*pi-sa[j]);

      // Disable features if f-angle is greater than 1.
      if ( _dir_thres_weak>0.99) nrank=1;
      tranks[j] = 3-nrank;
    }
  }

  _buf->delete_dataitem( sa_attr->name());
  _buf->init_done( false);

  // Reduce tangranks to ensure shared nodes have the same value across panes
  _surf->reduce_on_shared_nodes( _tangranks, Manifold::OP_MIN);
}

// Solve the minimization problem (x^T)Ax+2b^Tx with eigen-decomposition.
int FaceOffset_3::
eigen_analyze_vertex( Vector_3 A_io[3], Vector_3 &b_io, 
		      double angle_defect) {

  // A_io will be replaced by its eigenvectors.
  Vector_3 *es = A_io;

  // Create a backup of b_io, as b_io will be replaced by eigenvalues.
  const Vector_3 b = b_io;
  Vector_3 &lambdas = b_io;

  // Compute the eigenvalues and eigenvectors of the matrix. 
  // Note that lambdas will be in decreasing order!
  compute_eigenvectors( es, lambdas);
  
  // Classify offset intersection based on acos(-b.norm()/rv) and lambdas
  if ( std::abs(angle_defect) >= _tol_ad)
    return 3; // Sharp corner
  else if (lambdas[1] < _dir_thres_weak*lambdas[0])
    return 1; // Ridge
  else if ( lambdas[2] < _dir_thres_strong*lambdas[1])
    return 2; // Flat
  else        // Vertices that are neither sharp nor flat, but have
    return 4; // ambiguous ridge directions
}

// Update tangential motion of ridge vertices
void FaceOffset_3::
update_vertex_centers() {
  COM::DataItem *vcenters_buf =
    _buf->new_dataitem( "FO_vcenters_buf", 'n', COM_DOUBLE, 3, "");
  COM::DataItem *vecdiff_buf =
    _buf->new_dataitem( "FO_vecdiff_buf", 'n', COM_DOUBLE, 3, "");
  _buf->resize_array( vcenters_buf, 0);
  _buf->resize_array( vecdiff_buf, 0);

  // Loop through the panes and ridge edges
  std::vector< COM::Pane*>::iterator it = _panes.begin();
  Manifold::PM_iterator pm_it=_surf->pm_begin();

  for ( int i=0, local_npanes = _panes.size(); i<local_npanes; 
	++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;
    const std::set< Edge_ID> &eset_pn = _edges[i];  // List of ridge edges

    const Point_3 *pnts = reinterpret_cast<Point_3*>
      (pane->dataitem(COM_NC)->pointer());
    Vector_3 *vcs = reinterpret_cast<Vector_3*>
      (pane->dataitem(vcenters_buf->id())->pointer());
    Vector_3 *vdiff = reinterpret_cast<Vector_3*>
      (pane->dataitem(vecdiff_buf->id())->pointer());

    for ( std::set< Edge_ID>::const_iterator eit=eset_pn.begin(),
	    eend=eset_pn.end(); eit!=eend; ++eit) {
      // Make sure that eid_opp is not border
      Edge_ID eid = *eit, eid_opp = (*pm_it)->get_opposite_real_edge( eid);
      bool is_border_edge = (*pm_it)->is_physical_border_edge( eid_opp);

      Element_node_enumerator ene( pane, eid.eid());
      int lid = eid.lid(), neighbor = (lid+1)%ene.size_of_edges();
      int vindex = ene[lid]-1, windex = ene[neighbor]-1;
      
      // For ridge and flat vertices, compute based on edge centers
      Vector_3 diff = pnts[ windex] - pnts[vindex];
      vcs[vindex] += diff;
      if ( vdiff[vindex].is_null()) vdiff[vindex] = diff; 
      else vdiff[vindex] -= diff;

      if ( is_border_edge) {
	vcs[windex] -= diff;
	if ( vdiff[windex].is_null()) vdiff[windex] = -diff; 
	else vdiff[windex] += diff;
      }
    }
  }
  
  _surf->reduce_on_shared_nodes( vcenters_buf, Manifold::OP_SUM);
  _surf->reduce_on_shared_nodes( vecdiff_buf, Manifold::OP_DIFF);

  const Vector_3 null_vec(0,0,0);
  it = _panes.begin();
  // Finally, copy from vcenters_buf to vcenters
  for ( int i=0, local_npanes = _panes.size(); 
	i<local_npanes; ++i, ++it) {
    COM::Pane *pane = *it;

    Vector_3 *vcs_buf = reinterpret_cast<Vector_3*>
      (pane->dataitem(vcenters_buf->id())->pointer());
    Vector_3 *vcs = reinterpret_cast<Vector_3*>
      (pane->dataitem(_vcenters->id())->pointer());
    Vector_3 *vdiff = reinterpret_cast<Vector_3*>
      (pane->dataitem(vecdiff_buf->id())->pointer());
    const char *tranks = reinterpret_cast<const char*>
      ( pane->dataitem(_tangranks->id())->pointer());
    const char *cranks = reinterpret_cast<const char*>
      ( pane->dataitem(_ctangranks->id())->pointer());
    const char *val_bndry_nodes = reinterpret_cast<const char*>
      ( pane->dataitem(_cnstr_bndry_nodes->id())->pointer());
    
    for (int j=0, jn=pane->size_of_real_nodes(); j<jn; ++j) if (cranks[j]!=2) {
      if ( cranks[j]==0)
	vcs[ j] = null_vec;
      else if ( val_bndry_nodes && val_bndry_nodes[j] && !vdiff[j].is_null()) {
	vcs_buf[j] = (vcs_buf[j]*vdiff[j]*vdiff[j])/vdiff[j].squared_norm();

	vcs[ j] = (1./4.)*vcs_buf[j];
      }
    }
  }

  _buf->delete_dataitem( vecdiff_buf->name());
  _buf->delete_dataitem( vcenters_buf->name());
  _buf->init_done( false);
}

// Mark weak vertices, i.e. the vertices that are not on or incident 
// on features but has non-trivial second-largest eigenvalues.
void FaceOffset_3::mark_weak_vertices() {
  char zero_c=0;
  Rocblas::copy_scalar( &zero_c, _weak);

  std::vector< COM::Pane*>::iterator it = _panes.begin();
  Manifold::PM_iterator pm_it=_surf->pm_begin();
  for ( int i=0, local_npanes = _panes.size(); 
	i<local_npanes; ++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;
    char *tranks = reinterpret_cast<char*>
      ( pane->dataitem(_tangranks->id())->pointer());
    char *weak = reinterpret_cast<char*>
      ( pane->dataitem(_weak->id())->pointer());
    const Vector_3 *lambdas = reinterpret_cast<const Vector_3*>
      ( pane->dataitem(_eigvalues->id())->pointer());

    for ( int j=0, nj=pane->size_of_real_nodes(); j<nj; ++j) {
      weak[j] = tranks[j]==2 && 
	lambdas[j][1] >= _eig_weak_lb*lambdas[j][0] &&
	lambdas[j][1] <= _eig_weak_ub*lambdas[j][0] ||
	tranks[j]<=1 && _eig_thres*lambdas[j][0] > lambdas[j][1];
    }
  }

  _surf->reduce_on_shared_nodes( _weak, Manifold::OP_MIN);
}

inline double SQR( double x) { return (x*x); }

// ----------------------------------------------------------------------------
void dsytrd3(double A[3][3], double Q[3][3], double d[3], double e[2])
// ----------------------------------------------------------------------------
// Reduces a symmetric 3x3 matrix to tridiagonal form by applying
// (unitary) Householder transformations:
//            [ d[0]  e[0]       ]
//    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
//            [       e[1]  d[2] ]
// The function accesses only the diagonal and upper triangular parts of
// A. The access is read-only.
// ----------------------------------------------------------------------------
// Note: Based on implementat of Joachim Kopp (www.mpi-hd.mpg.de/~jkopp/3x3/)
// ---------------------------------------------------------------------------
{
  double u[3], q[3];
  double omega, f;
  double K, h, g;
  
  // Initialize Q to the identitity matrix
  for (int i=0; i < 3; i++)
  {
    Q[i][i] = 1.0;
    for (int j=0; j < i; j++)
      Q[i][j] = Q[j][i] = 0.0;
  }

  // Bring first row and column to the desired form 
  h = SQR(A[0][1]) + SQR(A[0][2]);
  if (A[0][1] > 0)
    g = -sqrt(h);
  else
    g = sqrt(h);
  e[0] = g;
  f    = g * A[0][1];
  u[1] = A[0][1] - g;
  u[2] = A[0][2];
  
  omega = h - f;
  if (omega > 0.0)
  {
    omega = 1.0 / omega;
    K     = 0.0;
    for (int i=1; i < 3; i++)
    {
      f    = A[1][i] * u[1] + A[i][2] * u[2];
      q[i] = omega * f;                  // p
      K   += u[i] * f;                   // u* A u
    }
    K *= 0.5 * SQR(omega);

    for (int i=1; i < 3; i++)
      q[i] = q[i] - K * u[i];
    
    d[0] = A[0][0];
    d[1] = A[1][1] - 2.0*q[1]*u[1];
    d[2] = A[2][2] - 2.0*q[2]*u[2];
    
    // Store inverse Householder transformation in Q
    for (int j=1; j < 3; j++)
    {
      f = omega * u[j];
      for (int i=1; i < 3; i++)
        Q[i][j] = Q[i][j] - f*u[i];
    }

    // Calculate updated A[1][2] and store it in e[1]
    e[1] = A[1][2] - q[1]*u[2] - u[1]*q[2];
  }
  else
  {
    for (int i=0; i < 3; i++)
      d[i] = A[i][i];
    e[1] = A[1][2];
  }
}

// ----------------------------------------------------------------------------
int dsyevq3(double A[3][3], double Q[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using the QL algorithm with implicit shifts, preceded by a
// Householder reduction to tridiagonal form.
// The function accesses only the diagonal and upper triangular parts of A.
// The access is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
// Dependencies:
//   dsytrd3()
// ----------------------------------------------------------------------------
// Note: Based on implementat of Joachim Kopp (www.mpi-hd.mpg.de/~jkopp/3x3/)
// ----------------------------------------------------------------------------
{
  double e[3];                   // The third element is used only as temporary workspace
  double g, r, p, f, b, s, c, t; // Intermediate storage
  int nIter;
  int m;

  // Transform A to real tridiagonal form by the Householder method
  dsytrd3(A, Q, w, e);
  
  // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
  // with the QL method
  //
  // Loop over all off-diagonal elements
  for (int l=0; l < 2; l++)
  {
    nIter = 0;
    for (;;)
    {
      // Check for convergence and exit iteration loop if off-diagonal
      // element e(l) is zero
      for (m=l; m < 2; m++)
      {
        g = fabs(w[m])+fabs(w[m+1]);
        if (fabs(e[m]) + g == g)
          break;
      }
      if (m == l)
        break;
      
      if (nIter++ >= 30)
        return -1;

      // Calculate g = d_m - k
      g = (w[l+1] - w[l]) / (e[l] + e[l]);
      r = sqrt(SQR(g) + 1.0);
      if (g > 0)
        g = w[m] - w[l] + e[l]/(g + r);
      else
        g = w[m] - w[l] + e[l]/(g - r);

      s = c = 1.0;
      p = 0.0;
      for (int i=m-1; i >= l; i--)
      {
        f = s * e[i];
        b = c * e[i];
        if (fabs(f) > fabs(g))
        {
          c      = g / f;
          r      = sqrt(SQR(c) + 1.0);
          e[i+1] = f * r;
          c     *= (s = 1.0/r);
        }
        else
        {
          s      = f / g;
          r      = sqrt(SQR(s) + 1.0);
          e[i+1] = g * r;
          s     *= (c = 1.0/r);
        }
        
        g = w[i+1] - p;
        r = (w[i] - g)*s + 2.0*c*b;
        p = s * r;
        w[i+1] = g + p;
        g = c*r - b;

        // Form eigenvectors
        for (int k=0; k < 3; k++)
        {
          t = Q[k][i+1];
          Q[k][i+1] = s*Q[k][i] + c*t;
          Q[k][i]   = c*Q[k][i] - s*t;
        }
      }
      w[l] -= p;
      e[l]  = g;
      e[m]  = 0.0;
    }
  }

  return 0;
}

// Compute eigenvalues and eigenvectors of a 3x3 matrix A.
// At output, the eigenvalues are saved in lambdas in decending order.
// The columns of A are replaced by the orthonormal eigenvectors. 
void FaceOffset_3::
compute_eigenvectors( Vector_3 A[3], Vector_3 &lambdas) {

  double abuf[3][3] = { {A[0][0], A[0][1], A[0][2]},
			{A[1][0], A[1][1], A[1][2]},
			{A[2][0], A[2][1], A[2][2]}};
  double ebuf[3][3];
  
  int info = dsyevq3( abuf, ebuf, &lambdas[0]);
  COM_assertion_msg( info==0, "Computation of eigenvectos failed");

  std::swap( ebuf[0][1], ebuf[1][0]);
  std::swap( ebuf[0][2], ebuf[2][0]);
  std::swap( ebuf[1][2], ebuf[2][1]);

  const Vector_3 *buf = (const Vector_3*)&ebuf[0][0];

  if ( lambdas[0]<lambdas[1]) {
    if ( lambdas[1]<lambdas[2]) { // l2>l1>l0
      std::swap( lambdas[0], lambdas[2]);
      A[0] = buf[2]; A[1] = buf[1]; A[2] = buf[0];
    }
    else { 
      double t = lambdas[0];
      lambdas[0] = lambdas[1]; 
      A[0] = buf[1];

      if ( t<lambdas[2]) { // l1>=l2>l0
	lambdas[1] = lambdas[2]; lambdas[2] = t;
	A[1] = buf[2]; A[2] = buf[0];
      }
      else { // l1>l0>=l2
	lambdas[1] = t;
	A[1] = buf[0]; A[2] = buf[2];
      }
    }
  }
  else if ( lambdas[0]<lambdas[2]) { // l2>l0>=l1
    double t = lambdas[0];
    lambdas[0] = lambdas[2]; lambdas[2] = lambdas[1]; lambdas[1] = t;
    A[0] = buf[2]; A[1] = buf[0]; A[2] = buf[1];
  }
  else {
    A[0] = buf[0];
    if ( lambdas[1]<lambdas[2]) { // l0>=l2>l1
      std::swap( lambdas[1], lambdas[2]);
      A[1] = buf[2]; A[2] = buf[1];
    }
    else {// l0>=l1>=l2
      A[1] = buf[1]; A[2] = buf[2]; 
    }
  }
}

PROP_END_NAMESPACE






