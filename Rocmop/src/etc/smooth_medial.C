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
#include "Rocmop.h"
#include "roccom.h"
#include "Pane.h"
#include "Rocblas.h"
#include "Rocmap.h"
#include "geometry.h"
#include "Rocsurf.h"
#include "Generic_element_2.h"

MOP_BEGIN_NAMESPACE

using MAP::Rocmap;

//============== Debug features =======
// Export debugging information to user window, if user has registered
// "back-door" attributes to get debug info such as normals. 
// Should be set to true for development purpose.
#define EXPORT_DEBUG_INFO     1

static const double half_pi = 1.5707963267949;
static const double pi      = 3.14159265358979;

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

void Rocmop::smooth_surf_medial(){
  if(_verb > 1) 
    std::cout << "Entering Rocmop::smooth_medial" << std::endl;

  COM_assertion_msg(_wrk_window, "Unexpected NULL pointer encountered.");

  // Calculate the face normals
  Rocmop::evaluate_face_normals();

  // Compute the medial quadric
  Rocmop::compute_medial_quadric();

#if EXPORT_DEBUG_INFO
 COM::Attribute *vnormals = const_cast<COM::Attribute*>(_usr_window->attribute("vnormals"));
 COM::Attribute *vnrms = _wrk_window->attribute("vnormals");
 COM::Attribute *awnormals = const_cast<COM::Attribute*>(_usr_window->attribute("angleweighted"));
 COM::Attribute *awnrms = _wrk_window->attribute("awnormals");
 COM::Attribute *uwnormals = const_cast<COM::Attribute*>(_usr_window->attribute("unitweighted"));
 COM::Attribute *uwnrms = _wrk_window->attribute("uwnormals");
 COM::Attribute *tangranks = const_cast<COM::Attribute*>(_usr_window->attribute("tangranks"));
 COM::Attribute *tranks = _wrk_window->attribute("tangranks");
 COM::Attribute *evals = _wrk_window->attribute("lambda");
 COM::Attribute *eigenvalues = const_cast<COM::Attribute*>(_usr_window->attribute("eigenvalues"));
 COM::Attribute *ws = _wrk_window->attribute("weights");
 COM::Attribute *weights = const_cast<COM::Attribute*>(_usr_window->attribute("adefects"));

 if ( vnormals)  Rocblas::copy( vnrms, vnormals);
 if ( awnormals)  Rocblas::copy( awnrms, awnormals);
 if ( uwnormals)  Rocblas::copy( uwnrms, uwnormals);
 if ( tangranks) Rocblas::copy( tranks, tangranks);
 if ( eigenvalues) Rocblas::copy(evals, eigenvalues);
 if ( weights) Rocblas::copy(ws,weights);
#endif
 
 // Identify ridge edges
 //identify_ridge_edges();
 
 // Redistribute vertices within their tangent spaces
 // FIXME.... MODIFY THESE TO USE PN TRIANGLES
 for ( int i=0; i<_rediter; ++i) {
   redistribute_vertices_ridge();
   redistribute_vertices_smooth();
 }
 
 // Add the displacements onto the nodal coords
 // FIXME
 const std::string att1("disps");
 COM::Attribute *disps = _wrk_window->attribute(att1);
 COM::Attribute *nc = _wrk_window->attribute( COM::COM_NC);
 // Rocblas::add(disps,nc,nc);
 
 if(_verb > 2) 
   std::cout << "Exiting Rocmop::smooth_medial" << std::endl;
}

void Rocmop::
compute_medial_quadric() {
 if(_verb > 1) 
   std::cout << "Entering Rocmop::compute_medial_quadric" << std::endl;
  // Use the following buffer spaces: _eigvecs for matrix A; 
  // _eigvalues to store bm; disps to store bo, _weights to store r
 const std::string att1("eigvecs"), att2("lambda"), att3("weights"),
    att4("facenormals"), att5("facecenters"), att6("tangranks"),
   att7("cntnranks"), att8("cntnvecs"), att9("vnormals"), att10("disps"),
   att11("awnormals"), att12("uwnormals");
  COM::Attribute *A_attr =  _wrk_window->attribute(att1);
  COM::Attribute *bm_attr = _wrk_window->attribute(att2);
  COM::Attribute *aw_attr = _wrk_window->attribute(att11);
  COM::Attribute *uw_attr = _wrk_window->attribute(att12);
  COM::Attribute *r_attr = _wrk_window->attribute(att3);

  int facenormals_id = _wrk_window->attribute(att4)->id();
  int facecenters_id = _wrk_window->attribute(att5)->id();
  int tangranks_id = _wrk_window->attribute(att6)->id();
  int cntnranks_id = _wrk_window->attribute(att7)->id();
  int cntnvecs_id = _wrk_window->attribute(att8)->id();
  int vnormals_id = _wrk_window->attribute(att9)->id();
  int awnormals_id = _wrk_window->attribute(att11)->id();
  int uwnormals_id = _wrk_window->attribute(att12)->id();
  int disps_id = _wrk_window->attribute(att10)->id();

  double zero = 0.;
  Rocblas::copy_scalar( &zero, A_attr);
  Rocblas::copy_scalar( &zero, bm_attr);
  Rocblas::copy_scalar( &zero, r_attr);

  // Initialize disps to 0
  std::vector< COM::Pane*> allpanes;
  _wrk_window->panes(allpanes);
  for(int i =0, local_npanes = allpanes.size();
      i<local_npanes; ++i){
    Vector_3<double> *ds = reinterpret_cast<Vector_3<double>*>
      ( allpanes[i]->attribute(disps_id)->pointer());
    for(int i =0, ni = allpanes[i]->size_of_real_nodes(); i<ni; ++i){
      ds[i] = Vector_3<double>(0.0,0.0,0.0);
    }
  }
  
  // Loop through the panes and its real faces, calculate A, b, and c = sum_i(wi)
  std::vector< COM::Pane*>::iterator it = allpanes.begin();
  SURF::Window_manifold_2::PM_iterator pm_it=_wm->pm_begin();
  for ( int i=0, local_npanes = allpanes.size(); 
	i<local_npanes; ++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;

    // Obtain nodal coordinates of current pane, assuming contiguous layout
    const Vector_3<double> *pnts = reinterpret_cast<Vector_3<double>*>
      (pane->attribute(COM_NC)->pointer());
    Vector_3<double> *nrms = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(facenormals_id)->pointer());

    // Obtain pointers to A_attr, bm_attr, and r_attr 
    Vector_3<double> *As = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(A_attr->id())->pointer());
    Vector_3<double> *bs_m = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(bm_attr->id())->pointer());
    Vector_3<double> *aws = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(aw_attr->id())->pointer());
    Vector_3<double> *uws = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(uw_attr->id())->pointer());
    double   *rs = reinterpret_cast<double*>
      ( pane->attribute(r_attr->id())->pointer());
    Vector_3<double> *cnts = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(facecenters_id)->pointer());

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1); 
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      // If speed is uniform within a face, then unit normal is the same
      // after propagation.
      const Vector_3<double> &ns_nz = nrms[j];
      Vector_3<double> nrm = ns_nz;

      int ne = ene.size_of_edges();
      int uindex=ene[ne-1]-1, vindex=ene[0]-1;

      // Loop through all vertices of current face.
      for ( int k=0; k<ne; ++k) {
	int windex = ene[(k+1==ne)?0:k+1]-1;
	
	Vector_3<double> cnt = cnts[j];

	// Calculate the weight of the current face for this node.
	double w=1;

	
	Vector_3<double> ts[2] = { pnts[uindex]-pnts[vindex], 
				   pnts[windex]-pnts[vindex]};
	
	double s = std::sqrt((ts[0]*ts[0])*(ts[1]*ts[1]));
	double angle=0;
	if ( s>=0) {
	  double cosw = ts[0]*ts[1]/s; 
	  if (cosw>1) cosw=1; else if ( cosw<-1) cosw=-1;
	  angle = std::acos( cosw);
	}
	
	if ( _wght_scheme!=SURF::E2N_ONE) {	    
	  
	  switch ( _wght_scheme) {
	  case SURF::E2N_ANGLE: 
	  case SURF::E2N_SPHERE: {
	    
	    if ( _wght_scheme==SURF::E2N_SPHERE)
	      w = std::sin(angle)/s; 
	    else if ( w>half_pi) 
	      w = std::max(0.,pi-angle); 
	    else
	      w = angle;
	    break;
	  }
	  case SURF::E2N_AREA: {
	    w = Vector_3<double>::cross_product(ts[0],ts[1]).norm()/2.; break;
	  }
	  }
	  nrm = w*ns_nz;
	}

	// Obtain the addresses of the linear system and weight
	Vector_3<double> *A_v = &As[3*vindex];
	Vector_3<double> &b_mv = bs_m[vindex];
	Vector_3<double> &aw = aws[vindex];
	Vector_3<double> &uw = uws[vindex];
	double   &r_v  = rs[vindex];

	// Update A, b_m, and r for the current vertex
	A_v[0] += ns_nz[0]*nrm;
	A_v[1] += ns_nz[1]*nrm;
	A_v[2] += ns_nz[2]*nrm;
	b_mv   -= nrm;
	aw     += (std::sin(angle)/s)*ns_nz;
	uw     += ns_nz;
	r_v    += angle;

	// Update indices for next node
	uindex=vindex; vindex=windex;
      }

#if defined(BOUNDARY_TREATMENT) && BOUNDARY_TREATMENT
      // Loop through the halfedges of the face to handle border edges
      Halfedge h( &*pm_it, Edge_ID(j+1,0), SURF::REAL_PANE), h0=h;
      do {
	if ( h.is_border()) {
	  // Incorporate normal plane passing through the border edge into 
	  // its incident vertices with weight pi/2.
	  Vector_3<double> nrm_nz = Vector_3<double>::cross_product( ns_nz, h.tangent()).normalize();
	  Vector_3<double> nrm = half_pi*nrm_nz;

	  // Update the origin and destination of h with the orthogonal pane
	  for ( int k=0; k<2; ++k) {
	    int vindex = ene[ (k?h:h.next()).id().lid()]-1;
	    // Obtain the addresses of the linear system and weight
	    Vector_3<double> *A_v = &As[3*vindex];
	    Vector_3<double> &b_mv = bs_m[vindex]; 
	    double   &r_v  = rs[vindex];

	    // Update A, b_m, and r for the current vertex. 
	    // No need to update as sigma is 0 for it.
	    A_v[0] += nrm_nz[0]*nrm;
	    A_v[1] += nrm_nz[1]*nrm;
	    A_v[2] += nrm_nz[2]*nrm;
	    b_mv   -= nrm;
	    r_v    += half_pi;
	  }
	}
      } while ( (h=h.next()) != h0);
#endif
    }
  }

  // Reduce on shared nodes for A, b_m, and r
  _wm->reduce_on_shared_nodes( A_attr, SURF::Window_manifold_2::OP_SUM);
  _wm->reduce_on_shared_nodes( bm_attr, SURF::Window_manifold_2::OP_SUM);
  _wm->reduce_on_shared_nodes( aw_attr, SURF::Window_manifold_2::OP_SUM);
  _wm->reduce_on_shared_nodes( uw_attr, SURF::Window_manifold_2::OP_SUM);
  _wm->reduce_on_shared_nodes( r_attr, SURF::Window_manifold_2::OP_SUM);

  // Loop through the panes and its real vertices
  it = allpanes.begin(); 
  for (int i=0, local_npanes = allpanes.size(); i<local_npanes; ++i, ++it) { 

    COM::Pane *pane = *it;

    // Obtain pointers to A_attr, bo_attr, bm_attr, and r_attr 
    Vector_3<double> *As = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(A_attr->id())->pointer());
    Vector_3<double> *bs_m = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(bm_attr->id())->pointer());
    double   *rs = reinterpret_cast<double*>
      ( pane->attribute(r_attr->id())->pointer());
    int      *tranks = reinterpret_cast<int*>
      ( pane->attribute(tangranks_id)->pointer());
    int      *cranks = reinterpret_cast<int*>
      ( pane->attribute(cntnranks_id)->pointer());
    Vector_3<double> *cvecs = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(cntnvecs_id)->pointer());
    int      *cnstrs = _cnstr_types?reinterpret_cast<int*>
      ( pane->attribute(_cnstr_types->id())->pointer()):NULL;
    Vector_3<double> *cnstr_dirs = _cnstr_dirs?reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(_cnstr_dirs->id())->pointer()):NULL;
    Vector_3<double>  *vnrms = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(vnormals_id)->pointer());
    Vector_3<double>  *awnrms = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(awnormals_id)->pointer());
    Vector_3<double>  *uwnrms = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(uwnormals_id)->pointer());


    // Loop through all real nodes of the pane
    //std::cout << "Number of nodes = " << pane->size_of_real_nodes() << std::endl;
    for ( int j=0, jn=pane->size_of_real_nodes(); j<jn; ++j) {

      // Skip vertices with zero weight, such as those at
      // edge- and face-centers for quadratic elements
      if ( rs[j]==0) continue;
      Vector_3<double> *es = &As[3*j];    // eigen-vectors of current vertex.
      Vector_3<double> *cs = &cvecs[2*j]; // cnstr-tangent of current vertex.
      // Make copy of matrix A as it will be overwritten by eigenvectors
      Vector_3<double> A_v[3] = {es[0], es[1], es[2]};
      // Likewise, b will be overwritten by eigenvalues
      Vector_3<double> b_v = bs_m[j];
      awnrms[j].normalize();
      uwnrms[j].normalize();
      

      // Perform eigen-analysis for the vertex using medial quadric
      // orank is the rank of the orthogonal (primar) space of A.
      // The non-constrained solution is stored in vnrm.
      int orank = eigen_analyze_vertex( es, bs_m[j], &vnrms[j], 2*pi-rs[j]);
      // EDIT
      //orank = (orank == 3) ? 2 : orank;

      //COM_assertion_msg((orank <=3) && (orank >= 1), "Orank is invalid\n");

      // Store the rank of the tangent space.
      tranks[j] = 3-orank;

      // Copy eigenvectors into cvecs.
      for ( int k=tranks[j]-1; k>=0; --k)
	cs[k] = es[orank+k];

      // Compute normal estimate using medial quadric, 
      // and save it into vnrms (overwrite the non-constrained solution)
      Vector_3<double> d(0,0,0);
      // Set the tolerance to avoid perturbation due to small eigenvalues
      double tol = bs_m[j][0]*5.e-3;
      for ( int k=0; k<orank; ++k)
	d += (es[k]*b_v/-std::max(bs_m[j][k],tol))*es[k];
      vnrms[j] = d;
    }
  }
  if(_verb > 1) 
    std::cout << "Exiting Rocmop::compute_medial_quadric" << std::endl;
}

// Regularize all nodes in their tangent spaces
void Rocmop::
redistribute_vertices_ridge() {
  // Use the following buffer spaces: _eigvalues for vector c (only first one)

  const std::string att1("disps"), att2("lambda"), att3("weights"),
    att4("cntnvecs"), att5("cntnranks");
  int disps_id = _wrk_window->attribute(att1)->id();
  // c_attr are the edge centers
  // w_attr is the nodal safe factor 
  COM::Attribute *c_attr = _wrk_window->attribute(att2);
  COM::Attribute *w_attr = _wrk_window->attribute(att3);
  int cntnvecs_id = _wrk_window->attribute(att4)->id();
  int cntnranks_id = _wrk_window->attribute(att5)->id();

  double zero = 0.;
  Rocblas::copy_scalar( &zero, c_attr);

  // Loop through the panes and its real faces
  std::vector< COM::Pane*> allpanes;
  _wrk_window->panes(allpanes);
  std::vector< COM::Pane*>::iterator it = allpanes.begin();
  Window_manifold_2::PM_iterator pm_it=_wm->pm_begin();
  for ( int i=0, local_npanes = allpanes.size(); 
	i<local_npanes; ++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;
    
    Vector_3<double> *es = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(cntnvecs_id)->pointer());
    const Point_3<double> *pnts = reinterpret_cast<Point_3<double>*>
      (pane->attribute(COM_NC)->pointer());
    Vector_3<double> *cs = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(c_attr->id())->pointer());
    int *cranks = reinterpret_cast<int*>
      ( pane->attribute(cntnranks_id)->pointer());

    std::set< Edge_ID> &eset_pn = _edges[i];
    for ( std::set< Edge_ID>::const_iterator eit=eset_pn.begin(),
	    eend=eset_pn.end(); eit!=eend; ++eit) {
      Edge_ID eid = *eit;
      bool    is_border = eid.is_border();
      if ( is_border) eid = pm_it->get_opposite_real_edge( eid);

      // Loop through real elements of the current pane
      Element_node_enumerator ene( pane, eid.eid()); 
      int ne = ene.size_of_edges();
      
      int vindex=ene[eid.lid()]-1, windex=ene[(eid.lid()+1)%ne]-1;
      if ( is_border) std::swap( vindex, windex);

      if ( cranks[vindex]>0) {
	Vector_3<double> diff_wv = pnts[windex]-pnts[vindex];
	// OLD:: Vector_3<double> diff_wv = pnts[windex]-pnts[vindex] + ds[windex] - ds[vindex];

	// get average of edge centers
	cs[vindex][0] += 0.5*(es[2*vindex]*diff_wv);
      }
    }
  }

  _wm->reduce_on_shared_nodes( c_attr, Window_manifold_2::OP_SUM);

  // Loop through all vertices to computer center
  it = allpanes.begin();  pm_it=_wm->pm_begin();
  for ( int i=0, local_npanes = allpanes.size(); 
	i<local_npanes; ++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;

    Vector_3<double> *es = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(cntnvecs_id)->pointer());
    Vector_3<double> *cs = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute( c_attr->id())->pointer());
    int *cranks = reinterpret_cast<int*>
      ( pane->attribute(cntnranks_id)->pointer());

    // Loop through all real nodes of the pane
    std::set< Edge_ID> &eset_pn = _edges[i];
    for ( std::set< Edge_ID>::const_iterator eit=eset_pn.begin(),
	    eend=eset_pn.end(); eit!=eend; ++eit) {
      Edge_ID eid = *eit;
      bool    is_border = eid.is_border();
      if ( is_border) eid = pm_it->get_opposite_real_edge( eid);

      // Loop through real elements of the current pane
      Element_node_enumerator ene( pane, eid.eid()); 
      int ne = ene.size_of_edges();
      
      int vindex=ene[eid.lid()]-1, windex=ene[(eid.lid()+1)%ne]-1;
      if ( is_border) std::swap( vindex, windex);

      if ( cranks[vindex]>0)
	cs[vindex] = 0.5*cs[vindex][0]*es[2*vindex];
    }
  }

  // Loop through the panes and its real faces to compute safe factors.
  get_redist_safe_factor( c_attr, w_attr, 1);

  // Loop through all vertices to update displacement vector 
  it = allpanes.begin();  pm_it=_wm->pm_begin();
  for ( int i=0, local_npanes = allpanes.size(); 
	i<local_npanes; ++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;

    Vector_3<double> *cs = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute( c_attr->id())->pointer());
    Vector_3<double> *ds = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute( disps_id)->pointer());
    double   *ws = reinterpret_cast<double*>
      ( pane->attribute(w_attr->id())->pointer());

    // Loop through all real nodes of the pane
    std::set< Edge_ID> &eset_pn = _edges[i];
    for ( std::set< Edge_ID>::const_iterator eit=eset_pn.begin(),
	    eend=eset_pn.end(); eit!=eend; ++eit) {
      Edge_ID eid = *eit;
      bool    is_border = eid.is_border();
      if ( is_border) eid = pm_it->get_opposite_real_edge( eid);

      // Loop through real elements of the current pane
      Element_node_enumerator ene( pane, eid.eid()); 
      int ne = ene.size_of_edges();
      
      int vindex=ene[eid.lid()]-1, windex=ene[(eid.lid()+1)%ne]-1;
      if ( is_border) std::swap( vindex, windex);

      double w=0.5;
      if ( ws[vindex]>0) w = std::min( w, 1./ws[vindex]);
      // displacement = w * average of two edge segments
      ds[vindex] += w*cs[vindex];
      // do PN triangle stuff here
    }
  }
}

// Regularize all nodes in their tangent spaces
void Rocmop::
redistribute_vertices_smooth() {

  // Use the following buffer spaces: _eigvalues for vector c 
  // (only lsast two components); _weights to store weight
  const std::string att1("disps"), att2("lambda"), att3("weights"),
    att4("cntnvecs"), att5("cntnranks"), att6("tangranks");
  int disps_id = _wrk_window->attribute(att1)->id();
  COM::Attribute *c_attr = _wrk_window->attribute(att2);
  COM::Attribute *w_attr = _wrk_window->attribute(att3);
  int cntnvecs_id = _wrk_window->attribute(att4)->id();
  int cntnranks_id = _wrk_window->attribute(att5)->id();
  int tangranks_id = _wrk_window->attribute("tangranks")->id();

  double zero = 0.;
  Rocblas::copy_scalar( &zero, c_attr);
  Rocblas::copy_scalar( &zero, w_attr);

  std::vector< COM::Pane*> allpanes;
  _wrk_window->panes(allpanes);
  // Loop through the panes and its real faces
  std::vector< COM::Pane*>::iterator it = allpanes.begin();
  for (int i=0, local_npanes = allpanes.size(); i<local_npanes; ++i, ++it) { 
    COM::Pane *pane = *it;
    
    Vector_3<double> *es = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(cntnvecs_id)->pointer());
    const Point_3<double> *pnts = reinterpret_cast<Point_3<double>*>
      (pane->attribute(COM_NC)->pointer());
    Vector_3<double> *cs = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(c_attr->id())->pointer());
    double   *ws = reinterpret_cast<double*>
      ( pane->attribute(w_attr->id())->pointer());
    Vector_3<double> *ds = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(disps_id)->pointer());
    int *tranks = reinterpret_cast<int*>
      ( pane->attribute(tangranks_id)->pointer());
    int *cranks = reinterpret_cast<int*>
      ( pane->attribute(cntnranks_id)->pointer());

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1); 
    Element_node_vectors_k_const<Point_3<double> > ps; ps.set( pnts, ene, 1);

    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      int ne = ene.size_of_edges();
      int uindex=ene[ne-1]-1, vindex=ene[0]-1;

      // Loop through all vertices of current face.
      for ( int k=0; k<ne; ++k) {
	int windex = ene[(k+1==ne)?0:k+1]-1;

	// Only redistribute smooth vertices 
	if ( tranks[vindex]==2) {
	  std::cout << "cranks[vindex]-1 = " << cranks[vindex]-1 << std::endl;
	  Vector_3<double> tngs[2] = {pnts[uindex]-pnts[vindex], 
			      pnts[windex]-pnts[vindex]};
	  Vector_3<double> disp_v = ds[vindex];
	  Vector_3<double> diff_uv = ds[uindex]-disp_v+tngs[0];
	  Vector_3<double> diff_wv = ds[windex]-disp_v+tngs[1];

	  // Compute angle at vertex
	  double theta = std::acos
	    ( diff_uv*diff_wv/std::sqrt((diff_wv*diff_wv)*(diff_uv*diff_uv)));
	  std::cout << "theta = " << theta << std::endl;

	  if ( theta>half_pi) theta = std::max(0.,pi-theta);
	  std::cout << "theta 2 = " << theta << std::endl;
  
	  for ( int k=cranks[vindex]-1; k>=0; --k) {
	    cs[vindex][k] += theta*(es[2*vindex+k]*(diff_uv+diff_wv));
	    //	    std::cout << "cs[" << vindex << "][" << k << "] += " 
	    //      << theta << "*(" << es[2*vindex+k] << " * (" 
	    //      << diff_uv << " + " << diff_wv << "))\n";
	  }
	  ws[vindex] += theta;
	  //	  std::cout << "ws[" << vindex << "] = " << ws[vindex]  << std::endl;
	}

	uindex=vindex; vindex=windex;
      }
    }
  }
  _wm->reduce_on_shared_nodes( c_attr, Window_manifold_2::OP_SUM);
  _wm->reduce_on_shared_nodes( w_attr, Window_manifold_2::OP_SUM);

  // Loop through all vertices to compute center of stars.
  it = allpanes.begin();
  for (int i=0, local_npanes = allpanes.size(); i<local_npanes; ++i, ++it) { 
    COM::Pane *pane = *it;

    Vector_3<double> *cs = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute( c_attr->id())->pointer());
    double   *ws = reinterpret_cast<double*>
      ( pane->attribute( w_attr->id())->pointer());
    int *tranks = reinterpret_cast<int*>
      ( pane->attribute(tangranks_id)->pointer());
    Vector_3<double> *es = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(cntnvecs_id)->pointer());
    int *cranks = reinterpret_cast<int*>
      ( pane->attribute(cntnranks_id)->pointer());

    // Loop through all real nodes of the pane
    for ( int j=0, nj=pane->size_of_real_nodes(); j<nj; ++j) {
      if ( tranks[j]==2) {
	Vector_3<double> s(0,0,0);
	for ( int k=cranks[j]-1; k>=0; --k) {
	  s += (cs[j][k]*es[2*j+k])/ws[j];
	}
	cs[j] = s;
      }
    }
  }

  // Loop through the panes and its real faces to compute safe factors.
  get_redist_safe_factor( c_attr, w_attr, 2);

  // Loop through all vertices to update displacement vector 
  it = allpanes.begin();
  for (int i=0, local_npanes = allpanes.size(); i<local_npanes; ++i, ++it) { 
    COM::Pane *pane = *it;

    Vector_3<double> *cs = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute( c_attr->id())->pointer());
    double   *ws = reinterpret_cast<double*>
      ( pane->attribute( w_attr->id())->pointer());
    Vector_3<double> *ds = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute( disps_id)->pointer());
    int *tranks = reinterpret_cast<int*>
      ( pane->attribute( tangranks_id)->pointer());

    // Loop through all real nodes of the pane
    for ( int j=0, nj=pane->size_of_real_nodes(); j<nj; ++j) {
      if ( tranks[j]==2) {
	double w=0.5;
	if ( ws[j]>0) w = std::min( w, 1./ws[j]);
	ds[j] += w*cs[j];
	std::cout << "ds[" << j << "] += " << w << " * " << cs[j] << std::endl;
      }
    }
  }
}

void Rocmop::
get_redist_safe_factor( COM::Attribute *c_attr, COM::Attribute *w_attr, 
			int rank) {
  double zero=0;
  Rocblas::copy_scalar( &zero, w_attr);

  const std::string att1("tangranks"), att2("weights2"),
    att3("barycrds"), att4("PNelemids");
  int tangranks_id = _wrk_window->attribute(att1)->id(),
    weights2_id = _wrk_window->attribute(att2)->id(),
    barycrds_id = _wrk_window->attribute(att3)->id(),
    PNelemid_id = _wrk_window->attribute(att4)->id();

  std::vector< COM::Pane*> allpanes;
  _wrk_window->panes(allpanes);
  std::vector< COM::Pane*>::iterator it = allpanes.begin();
  for (int i=0, local_npanes = allpanes.size(); i<local_npanes; ++i, ++it) { 
    COM::Pane *pane = *it;
    
    const Point_3<double> *pnts = reinterpret_cast<Point_3<double>*>
      (pane->attribute(COM_NC)->pointer());
    Vector_3<double> *cs = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(c_attr->id())->pointer());
    double   *ws = reinterpret_cast<double*>
      ( pane->attribute(w_attr->id())->pointer());
    int *tranks = reinterpret_cast<int*>
      ( pane->attribute(tangranks_id)->pointer());
    double *ws2 = reinterpret_cast<double*>
      ( pane->attribute(weights2_id)->pointer());
    double *bcs = reinterpret_cast<double*>
      ( pane->attribute(barycrds_id)->pointer());
    int *PNids = reinterpret_cast<int*>
      ( pane->attribute(PNelemid_id)->pointer());

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1); 

    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      int ne = ene.size_of_edges();
      int uindex=ene[ne-1]-1, vindex=ene[0]-1;
      Vector_3<double> nrm_nz;

      // Loop through all vertices of current face.
      for ( int k=0; k<ne; ++k) {
	int windex = ene[(k+1==ne)?0:k+1]-1;

	// Only redistribute correct classification of vertices 
	if ( tranks[vindex]==rank) {
	  Vector_3<double> tngs[2] = {pnts[uindex]-pnts[vindex], 
			      pnts[windex]-pnts[vindex]};

	  // Compute normal of face offset
	  if ( k==0 || ne>3) {
	    nrm_nz = Vector_3<double>::cross_product( tngs[1], tngs[0]);
	    nrm_nz.normalize();
	  }
	  
	  // Solve for vector x, to determine maximum time step
	  Vector_3<double> S[3] = { tngs[1], tngs[0], nrm_nz};
	  Vector_3<double> xs; solve( S, cs[windex], xs);

	  double temp = ws[windex];
	  ws[windex] = std::max( (6-ne)*std::max( xs[0], xs[1]), ws[windex]);
	  if(temp != ws[windex]){
	    ws2[windex] = ws[windex];
	    bcs[windex*2]= xs[0];
	    bcs[windex*2+1] = xs[1];
	    PNids[windex] = j+1;
	  }
	}

	uindex=vindex; vindex=windex;
      }
    }
  }

  _wm->reduce_on_shared_nodes( w_attr, Window_manifold_2::OP_MAX);
}

void Rocmop::evaluate_face_normals() {

  if(_verb > 1) 
    std::cout << "Entering Rocmop::evaluate_face_normals" << std::endl;

  COM_assertion_msg(_wrk_window, "Unexpected NULL pointer encountered.");

  if(_verb > 3) 
    std::cout << "Getting attribute ids" << std::endl;

  const std::string att1("facenormals");
  const std::string att2("facecenters");
  int facenormals_id = _wrk_window->attribute(att1)->id();
  int facecenters_id = _wrk_window->attribute(att2)->id();

  // Loop through the panes and its real faces
  std::vector< COM::Pane*> allpanes;
  _wrk_window->panes(allpanes);

  std::vector< COM::Pane*>::iterator it = allpanes.begin();
  for ( int i=0, local_npanes = allpanes.size(); i<local_npanes; ++i, ++it) { 
    COM::Pane *pane = *it;

    // Obtain nodal coordinates of current pane, assuming contiguous layout
    COM_assertion_msg( pane->attribute( COM_NC)->stride()==3 ||
		       pane->size_of_real_nodes()==0, 
		       "Coordinates must be stored in contiguous layout.");

    const Vector_3<double> *pnts = reinterpret_cast< Vector_3<double>* >
      (pane->attribute(COM_NC)->pointer());
    Vector_3<double> *nrms = reinterpret_cast< Vector_3<double>* >
      ( pane->attribute(facenormals_id)->pointer());
    Vector_3<double> *cnts = reinterpret_cast< Vector_3<double>* >
      ( pane->attribute(facecenters_id)->pointer());

    COM_assertion_msg(pnts, "NULL pointer to nodal coords.");
    COM_assertion_msg(nrms, "NULL pointer to face normals.");
    COM_assertion_msg(cnts, "NULL pointer to face centers.");    

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1); 
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {

      // Obtain current face normal
      Element_node_vectors_k_const<Vector_3<double> > ps; ps.set( pnts, ene, 1);
      SURF::Generic_element_2 e(ene.size_of_edges(), ene.size_of_nodes());

      Vector_2<double> nc(0.5,0.5);
      Vector_3<double> J[2];

      e.Jacobian( ps, nc, J);

      nrms[j] = Vector_3<double>::cross_product( J[0], J[1]);

      nrms[j].normalize();

      e.interpolate_to_center( ps, &cnts[j]);

    }
  }

  if(_verb > 1) 
    std::cout << "Exiting Rocmop::evaluate_face_normals" << std::endl;
}

void Rocmop::
identify_ridge_edges() 
{
  // Allocate buffer space for maximum and minedgev
  COM::Attribute *maxedgev = 
    _wrk_window->new_attribute( "maxedgev", 'n', COM_DOUBLE, 1, "");
  _wrk_window->resize_array( maxedgev, 0);

  COM::Attribute *minedgev = 
    _wrk_window->new_attribute( "minedgev", 'n', COM_DOUBLE, 1, "");
  _wrk_window->resize_array( minedgev, 0);
  _wrk_window->init_done();

  // Loop through the panes and its real faces
  std::vector< COM::Pane*> allpanes;
  _wrk_window->panes(allpanes);
  std::vector< COM::Pane*>::iterator it = allpanes.begin();
  Window_manifold_2::PM_iterator pm_it=_wm->pm_begin();

  const std::string att1("eigvecs"), att2("tangranks");

  int eigvecs_id = _wrk_window->attribute(att1)->id();
  int tangranks_id = _wrk_window->attribute(att2)->id();

  std::vector<std::map<Edge_ID, double> > edge_maps( allpanes.size());
  for ( int i=0, local_npanes = allpanes.size(); 
	i<local_npanes; ++i, ++it, ++pm_it) { 
    COM::Pane *pane = *it;
    std::map< Edge_ID, double> &edges_pn = edge_maps[i];

    // Obtain pointers to A_attr, bm_attr, and r_attr 
    Vector_3<double> *es = reinterpret_cast<Vector_3<double>*>
      ( pane->attribute(eigvecs_id)->pointer());
    int *tranks = reinterpret_cast<int*>
      ( pane->attribute(tangranks_id)->pointer());
    double *minvs = reinterpret_cast<double*>
      ( pane->attribute(minedgev->id())->pointer());
    double *maxvs = reinterpret_cast<double*>
      ( pane->attribute(maxedgev->id())->pointer());

    // Loop through real elements of the current pane
    Element_node_enumerator ene( pane, 1); 
    for ( int j=0, nj=pane->size_of_real_elements(); j<nj; ++j, ene.next()) {
      Halfedge h( &*pm_it, Edge_ID(j+1,0), SURF::REAL_PANE), h0=h;
      
      int vindex = ene[ h.id().lid()]-1;
      do {
	int windex = ene[ h.next().id().lid()]-1;
	
	// If vertex is on a ridge
	if ( tranks[ vindex]==1) {
	  int m = 1 + (tranks[windex]<2);
	  Vector_3<double> tng = h.tangent(); tng.normalize();
	  double feg = m * (es[3*vindex+2]* tng);

	  edges_pn[ h.id()] = feg;
	  if ( feg > maxvs[vindex]) maxvs[vindex] = feg;
	  if ( feg < minvs[vindex]) minvs[vindex] = feg;
	}

	// Process border edge
	if ( h.is_border() && tranks[windex]==1) {
	  int m = 1 + (tranks[vindex]<2);
	  Vector_3<double> tng = -h.tangent(); tng.normalize();
	  double feg = m * (es[3*windex+2]* tng);

	  edges_pn[ h.opposite().id()] = feg;
	  if ( feg > maxvs[windex]) maxvs[windex] = feg;
	  if ( feg < minvs[windex]) minvs[windex] = feg;
	}

	// Update vindex
	vindex = windex;
      } while ( (h=h.next()) != h0);
    }
  }

  // Reduce on shared nodes for minedgev and maxedgev
  _wm->reduce_on_shared_nodes( minedgev, Window_manifold_2::OP_MIN);
  _wm->reduce_on_shared_nodes( maxedgev, Window_manifold_2::OP_MAX);

  it = allpanes.begin(); 
  // Insert ridge edges onto _edges
  _edges.clear(); _edges.resize( allpanes.size());
  for ( int i=0, local_npanes = allpanes.size(); i<local_npanes; ++i, ++it) { 
    COM::Pane *pane = *it;
    std::map< Edge_ID, double> &emap_pn = edge_maps[i];
    std::set< Edge_ID> &eset_pn = _edges[i];

    const double *minvs = reinterpret_cast<double*>
      ( pane->attribute(minedgev->id())->pointer());
    const double *maxvs = reinterpret_cast<double*>
      ( pane->attribute(maxedgev->id())->pointer());

    // Loop throug ebuf and put ridges edges onto edges_pn
    for ( std::map< Edge_ID, double>::const_iterator eit=emap_pn.begin(), 
	    eend=emap_pn.end(); eit!=eend; ++eit) {
      Element_node_enumerator ene( pane, eit->first.eid());
      int vindex = ene[eit->first.lid()]-1;

      // Push edge onto edge
      if ( eit->second == minvs[vindex] || eit->second == maxvs[vindex])
	eset_pn.insert( eit->first);
    }
  }

  // Deallocate minedgev and maxedgev
  _wrk_window->delete_attribute( minedgev->name());
  _wrk_window->delete_attribute( maxedgev->name());
  _wrk_window->init_done();
}

int Rocmop::
eigen_analyze_vertex( Vector_3<double> A_io[3], Vector_3<double> &b_io, 
		      Vector_3<double> *nrm_nz, double ad) {

  // A_io will be replaced by its eigenvectors.
  Vector_3<double> *es = A_io;

  // Create a backup of b_io, as b_io will be replaced by eigenvalues.
  const Vector_3<double> b = b_io;
  Vector_3<double> &lambdas = b_io;

  // Compute the eigenvalues and eigenvectors of the matrix. 
  // Note that lambdas will be in decreasing order!
  compute_eigenvectors( es, lambdas);

  int nrank;
  double eps = lambdas[0]*1.e-8;
  double gs[] = { b*es[0]/lambdas[0],
		  b*es[1]/std::max(eps,lambdas[1]), 
		  b*es[2]/std::max(eps,lambdas[2]) };
  double gs_abs[] = {std::abs(gs[0]), std::abs(gs[1]), std::abs(gs[2])};
  
  // Classify offset intersection based on acos(-b.norm()/rv) and lambdas
  if ( lambdas[2] >= _saliency_crn *
       std::max(lambdas[0]-lambdas[1], lambdas[1]-lambdas[2]) || ad >= .5*pi)
    nrank = 3;
  else if ( gs_abs[1] >= gs_abs[0] || lambdas[0]*_dir_thres<lambdas[1]) {
    if ( lambdas[1] < lambdas[0]*_eig_thres) {
      lambdas[1] = lambdas[0]*_eig_thres;
      std::cerr << "Rocprop Warning: Mesh contains near cusp." 
		<< "lambdas[1]=" << lambdas[1] << " and lambdas[0]="
		<< lambdas[0] << std::endl;
    }
    nrank = 2; // ridge vertex
#if PRINT_FEATURE
    if ( gs_abs[1] < gs_abs[0]) {
      double theta = (360/pi)*std::atan(std::sqrt(lambdas[1]/lambdas[0]));
      if ( theta<fangle_min_eigen) fangle_min_eigen=theta;
      if ( theta>fangle_max_eigen) fangle_max_eigen=theta;
    }
#endif
  }
  else
    nrank = 1; // smooth vertex

  if ( !_reorthog || nrank==1)
    *nrm_nz=( es[0]*b>0)?-es[0]:es[0];
  else {
    // Solve for x within primary space
    Vector_3<double> x(0,0,0);
    for ( int k=0; k<nrank; ++k) x -= gs[k]*es[k];
    *nrm_nz = x.normalize();
  }
  return nrank;
}

#if 0
// Solve the minimization problem (x^T)Ax+2b^Tx with eigen-decomposition.
int Rocmop::
eigen_analyze_vertex( Vector_3<double> A_io[3], Vector_3<double> &b_io, 
		      Vector_3<double> *nrm_nz) {

  // A_io will be replaced by its eigenvectors.
  Vector_3<double> *es = A_io;

  // Create a backup of b_io, as b_io will be replaced by eigenvalues.
  const Vector_3<double> b = b_io;
  Vector_3<double> &lambdas = b_io;

  // Compute the eigenvalues and eigenvectors of the matrix. 
  // Note that lambdas will be in decreasing order!
  compute_eigenvectors( es, lambdas);

  int orank;
  double eps = lambdas[0]*1.e-8;
  double gs[] = { b*es[0]/lambdas[0],
		  b*es[1]/std::max(eps,lambdas[1]), 
		  b*es[2]/std::max(eps,lambdas[2]) };
  double gs_abs[] = {std::abs(gs[0]), std::abs(gs[1]), std::abs(gs[2])};

  // Classify offset intersection based on acos(-b.norm()/rv) and lambdas
  if ( lambdas[2] >= _saliency_crn*
       std::max(lambdas[0]-lambdas[1], lambdas[1]-lambdas[2]) || 
       gs_abs[2]>=gs_abs[1] && gs_abs[2] >= gs_abs[0])
    orank = 3;
  else if ( gs_abs[1] >= gs_abs[0] || lambdas[0]*_dir_thres<lambdas[1]) {
    if ( lambdas[1] < lambdas[0]*_eig_thres) { 
      lambdas[1] = lambdas[0]*_eig_thres;
      std::cerr << "Rocmop Warning: Mesh contains near cusp." << std::endl;
    }
    orank = 2; // ridge vertex
  }
  else
    orank = 1; // smooth vertex

  if ( !_reorthog || orank==1)
    *nrm_nz=( es[0]*b>0)?-es[0]:es[0];
  else {
    // Solve for x within primary space
    Vector_3<double> x(0,0,0);
    for ( int k=0; k<orank; ++k) x -= gs[k]*es[k];
    *nrm_nz = x.normalize();
  }
  return orank;
}
#endif

// Compute eigenvalues and eigenvectors of a 3x3 matrix A.
// At output, the eigenvalues are saved in lambdas in decending order.
// The columns of A are replaced by the orthonormal eigenvectors. 
void Rocmop::
compute_eigenvectors( Vector_3<double> A[3], Vector_3<double> &lambdas) {

  double abuf[3][3] = { {A[0][0], A[0][1], A[0][2]},
			{A[1][0], A[1][1], A[1][2]},
			{A[2][0], A[2][1], A[2][2]}};
  double ebuf[3][3];
  
  int info = dsyevq3( abuf, ebuf, &lambdas[0]);
  COM_assertion_msg( info==0, "Computation of eigenvectos failed");

  std::swap( ebuf[0][1], ebuf[1][0]);
  std::swap( ebuf[0][2], ebuf[2][0]);
  std::swap( ebuf[1][2], ebuf[2][1]);

  const Vector_3<double> *buf = (const Vector_3<double>*)&ebuf[0][0];

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

void Rocmop::
get_constraint_directions( int type, const Vector_3<double> &dir,
			   int &ndirs, Vector_3<double> dirs[2]) {
  if ( type == 1 || type == -1)
    COM_assertion_msg( &dir && std::abs(dir.squared_norm()-1)<=1.e-6, 
		       "Constrained direction must be normalized");

  switch ( type) {
  case 0: 
    ndirs = 0;
    break;
  case 2: 
    ndirs = 3;
    break;
  case 1:
    ndirs = 1; dirs[0] = dir;
    break;
  case -1: {
    ndirs = 2;

    double sqnrm;
    if ( std::abs( dir[2])>0.7) {
      sqnrm = dir[1]*dir[1]+dir[2]*dir[2];
      dirs[0] = Vector_3<double>( 0, -dir[2], dir[1]);
      dirs[1] = Vector_3<double>( sqnrm, -dir[0]*dir[1], -dir[0]*dir[2]);
    }
    else {
      sqnrm = dir[0]*dir[0]+dir[1]*dir[1];
      dirs[0] = Vector_3<double>( dir[1], -dir[0], 0);
      dirs[1] = Vector_3<double>( -dir[0]*dir[2], -dir[1]*dir[2], sqnrm);
    }

    double s = 1./std::sqrt( sqnrm);
    dirs[0] *= s; dirs[1] *= s;
    break;
  }
  case 'x':  
    ndirs = 1; 
    dirs[0] = Vector_3<double>(1,0,0); 
    break;
  case -'x': 
    ndirs = 2; 
    dirs[0] = Vector_3<double>(0,1,0); 
    dirs[1] = Vector_3<double>(0,0,1); 
    break;
  case 'y': 
    ndirs = 1; 
    dirs[0] = Vector_3<double>(0,1,0); 
    break;
  case -'y':
    ndirs = 2;
    dirs[0] = Vector_3<double>(1,0,0); 
    dirs[1] = Vector_3<double>(0,0,1); 
    break;
  case 'z': 
    ndirs = 1; 
    dirs[0] = Vector_3<double>(0,0,1); 
    break;
  case -'z':
    ndirs = 2;
    dirs[0] = Vector_3<double>(1,0,0); 
    dirs[1] = Vector_3<double>(0,1,0); 
    break;
  default: COM_assertion_msg( false, "Unknown type of constraint");
  }
}

MOP_END_NAMESPACE






