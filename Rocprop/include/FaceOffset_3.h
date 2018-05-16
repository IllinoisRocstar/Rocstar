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
// $Id: FaceOffset_3.h,v 1.23 2008/12/06 08:45:27 mtcampbe Exp $

/*******************************************************************
 * This file implements the Face Offsetting Method for 
 * interface propagation. It inherits some basic primitives
 * from the base class Propagation_3.
 *******************************************************************/

#ifndef _FACE_OFFSET_3_H_
#define _FACE_OFFSET_3_H_

#include "Propagation_3.h"
#include "Element_accessors.hpp"

PROP_BEGIN_NAMESPACE

class FaceOffset_3 : public Propagation_3 {
  typedef std::map<Edge_ID,std::pair<int, int> >    ObscendSet;
public:
  /// Default constructor.
  FaceOffset_3() : _is_strtd(0), _As(NULL), _boffset(NULL), _bmedial(NULL), 
		   _eigvalues(NULL), _eigvecs(NULL), _vnormals(NULL),
		   _facenormals(NULL), _facecenters(NULL), 
		   _faceareas(NULL), _vcenters(NULL), _tangranks(NULL), 
		   _ctangranks(NULL), _scales(NULL), _weights(NULL)
  { reset_default_parameters(); }

  /// Construct an object from a window manifold and a buffer window
  FaceOffset_3( Manifold *wm, COM::Window *buf);

  virtual ~FaceOffset_3() {}

  /// Main entry of the algorithm. At input, the value of smoothed specifies 
  /// whether smoothing is desired.
  virtual double time_stepping( const COM::DataItem *spds, double dt,
				COM::DataItem *disps, int *smoothed=NULL);

  /// Set the wavefrontal switch
  void set_wavefrontal_expansion( bool b) { _wf_expn = b; }

  /// Obtain the wavefrontal switch
  bool get_wavefrontal_expansion() const { return _wf_expn; }
  
  /// Set number of iterations for normal diffusion.
  void set_normal_diffusion( char i) { _nrm_dfsn = i; }

  /// Set whether or not to use feature layer
  void set_feature_layer( bool b) { _feature_layer = b; }

  /// Get number of iterations for normal diffusion
  int get_normal_diffuion() const { return _nrm_dfsn; }
  
  /// Set the threshold of eigenvalues for splitting primary and null spaces
  void set_eigen_threshold( double eps) 
  { _eig_thres = eps; COM_assertion( eps>0 && eps<1); }

  /// Get the threshold of eigenvalues for splitting primary and null spaces
  double get_eigen_threshold() const { return _eig_thres; }

  /// Set the constant in front of CFL condition (0<=c<1)
  void set_courant_constant( double c) 
  { _courant = c; _inv_sqrt_c = 1/std::sqrt(c); COM_assertion(0<=c && c<1); }

  /// Get the constant in front of CFL condition (0<=c<=1)
  double get_courant_constant() const { return _courant; }

  /// Weighting scheme
  void set_weighting_scheme( int w) { _wght_scheme = w; }
  int  get_weighting_scheme() const { return _wght_scheme; }

  /// Set the threshold for weak-feature angle
  void set_fangle_weak( double psi) {
    _tol_angle_weak = psi * pi/180;
    _dir_thres_weak = square(tan(0.5*_tol_angle_weak));
  }

  /// Set the threshold for weak-feature angle
  void set_fangle_strong( double psi) { 
    _tol_angle_strong = psi * pi/180;
  }

  /// Set the threshold for weak-feature angle
  void set_fangle_turn( double psi) { 
    _tol_cos_osta = std::cos(psi * pi/360);
  }

  /// Get the threshold for weak-feature angle
  double get_fangle_weak_radians() const 
  { return _tol_angle_weak; }

  /// Get the threshold for strong-feature angle
  double get_fangle_strong_radians() const 
  { return _tol_angle_strong; }

  /// Set the threshold for weak-feature angle
  double get_fangle_turn_radians() const 
  { return std::acos(_tol_cos_osta)*2; }

  /// Set the choice of mesh smoother.
  void set_smoother( int smoother) { _smoother = smoother; }

  /// Set mass conservation
  void set_conserve( int b) { _conserv = b; }

  // Set parameters to default values
  void reset_default_parameters( bool b=false);

protected:
  /// Introduce numerical dissipation into equation.
  void numerical_dissipation( COM::DataItem *nodal_buffer);

  /// Compute quadrics for every vertex.
  void compute_quadrics( double dt, const COM::DataItem *spds, 
			 COM::DataItem *rhs, COM::DataItem *predicted);

  /// Compute mean normals.
  void compute_directions( COM::DataItem *b_offset, bool in_principle);
  
  /// Compute volumn-vector
  void compute_volume_vector( const COM::DataItem *disps, COM::DataItem *bs);
  
  /// Decompose propagation directions based on constraints
  bool obtain_constrained_directions( COM::DataItem *disps, 
				      COM::DataItem *vcenters);
  
  /// Reduce time step based on stability limit, and rescale displacements 
  /// by the reduced time step and diffusing expansion.
  /// Returns the relative time step.
  double rescale_displacements( COM::DataItem *disps, COM::DataItem *vcenters,
				int depth=0);

  /// Filter out isolated ridge vertices and identify ridge edges
  void filter_and_identify_ridge_edges( bool filter_curves);

  // Upgrade/downgrade vertices between ridge and smooth vertices
  void reclassify_ridge_vertices( const bool upgrade_corners, 
				  const bool upgrade_ridge,
				  COM::DataItem *neighbor_feas,
				  COM::DataItem *tr_attr, 
				  bool to_filter);

  /// Build the list of ridge nighbors and obtain a list of weak-ended vertices
  /// Return the number of obscended vertices.
  int build_ridge_neighbor_list( const COM::DataItem *maxtrnangv, 
				 const COM::DataItem *mintrnangv,
				 const std::vector<std::map<Edge_ID, double> > &edge_maps,
				 const COM::DataItem *maxranglev, 
				 const COM::DataItem *minranglev,
				 const std::vector<std::map<Edge_ID, double> > &rangle_maps,
				 std::vector< ObscendSet > &obscends);

  int append_obscure_ends( const COM::DataItem *maxtrnangv, 
			   const COM::DataItem *mintrnangv,
			   const std::vector<std::map<Edge_ID, double> > &edge_maps,
			   const COM::DataItem *maxranglev, 
			   const COM::DataItem *minranglev,
			   const std::vector<std::map<Edge_ID, double> > &rangle_maps);

  /// Filter out weak-ended curves
  bool filter_obscended_ridge( const std::vector<std::map<Edge_ID, double> > &edge_maps,
			       const std::vector<std::map<Edge_ID, double> > &diangl_maps, 
			       const std::vector< ObscendSet > &obscends);

  int remove_obscure_curves( const std::vector< ObscendSet > &obscends);

  /// Mark vertices that are weak.
  void mark_weak_vertices();

  /// Redistribute smooth vertices by NuLaplacian smoothing
  void nulaplacian_smooth( const COM::DataItem *vert_normals,
			   const COM::DataItem *tangranks,
			   const COM::DataItem *disps,
			   COM::DataItem *vcenters,
			   COM::DataItem *buf,
			   COM::DataItem *vert_weights_buf,
			   const COM::DataItem *edge_weights=NULL);

  // Balance mass.
  void balance_mass();

  // Distribute volume from smooth nodes to elements
  void distribute_volume_e2n( const COM::DataItem *fvol,
			      const COM::DataItem *tranks,
			      COM::DataItem *vvol);

  void update_face_volumes( const COM::DataItem *fnormal,
			    const COM::DataItem *vdsps,
			    COM::DataItem *fvol);

  // Adjust the normal displacements for wavefrontal motion.
  void adjust_wavefrontal_motion( COM::DataItem *disps);

  // \param det: <0 indicates expansion and >0 indicates contraction
  // \param delta: "exact" motion along the normal direction
  // \param disp_nz: offset direction (unit length)
  // \param nrm_nz:  normal direction of the face (unit length).
  std::pair<double,double>
  comp_wavefrontal_motion( double det, double w, double delta, 
			   const Vector_3 &disp_nz, const Vector_3 &nrm_nz);

protected:
  // Solve the minimization problem (x^T)Ax+2b^Tx with eigen-decomposition
  // at a vertex, return the codimension of the vertex, and save mean
  // normal in x_nz.
  int eigen_analyze_vertex( Vector_3 A_io[3], Vector_3 &b_io, 
			    double angle_defect=0);

  void obtain_directions( Vector_3 es[3], const Vector_3 &lambdas, int nrank, 
			  Vector_3 &b_medial, Vector_3 *b_offset=NULL) const;

  // Obtain the normal and center of the face offset.
  // Return the displacement of the face center along normal direction.
  void obtain_face_offset( const Point_3 *pnts, const double *spd_ptr,
			   const COM::DataItem *spd, double dt, 
			   const Vector_3 *dirs, 
			   COM::Element_node_enumerator &ene,
			   Vector_3 &ns_nz, Point_3 &cnt, 
			   Vector_3 *disps, Vector_3 *ns);

  void get_tangents( const Vector_3 *es, const Vector_3 &lambdas,
		     const Vector_3 &mn, int nrank,
		     Vector_3 vs[2], double scales[2], 
		     int *is_strong=NULL) const;

  // Compute anisotropic vertex center.
  void compute_anisotropic_vertex_centers( const COM::DataItem *disps_buf);

  // Compute a normal motino to denoise surface.
  void denoise_surface( const COM::DataItem *disps_buf,
			COM::DataItem *normal_motion);

  void get_primary_component( const Vector_3 &nrm, const Vector_3 es[],
			      int trank, Vector_3 &prim) const {
    prim = Vector_3(0,0,0);

    for ( int i=0; i<3-trank; ++i) prim += nrm*es[i]*es[i];
  }

  // Append boundary edges separating different constraints into ridge edges.
  // Return the total number of inserted edges.
  int insert_boundary_edges( COM::DataItem *tr_attr);

  // Update tangential motion of ridge vertices 
  void update_vertex_centers();

protected: // Utilities.
  // Helpers for linear solvers for 2x2 and 3x3 equations.
  template <class FT>
  static void solve (const FT &a1, const FT &a2,
		     const FT &b1, const FT &b2, 
		     const FT &c1, const FT &c2, 
		     FT &x, FT &y)
  {
    FT denom = a1*b2-b1*a2;
    
    x = - (b1*c2-c1*b2)/denom;
    
    y = (a1*c2-c1*a2)/denom;
  }

  template <class FT>
  static void solve (const FT &a1, const FT &a2, const FT &a3,
		     const FT &b1, const FT &b2, const FT &b3,
		     const FT &c1, const FT &c2, const FT &c3,
		     const FT &d1, const FT &d2, const FT &d3,
		     FT &x, FT &y, FT &z)
  {
    FT denom = b2*c1*a3-b1*c2*a3+c3*b1*a2+b3*c2*a1-c1*b3*a2-b2*c3*a1;
    
    x = - (b2*c3*d1-b2*c1*d3+c1*b3*d2+b1*c2*d3-c3*b1*d2-b3*c2*d1)/denom;
    
    z = (b2*d1*a3-b2*a1*d3+b1*a2*d3-b1*d2*a3-d1*b3*a2+a1*b3*d2)/denom;
    
    y = (a2*c3*d1-a2*c1*d3-c2*d1*a3+c2*a1*d3+d2*c1*a3-d2*c3*a1)/denom;
  }

  // Wrapper for solving 3x3 equations
  static void solve ( const Vector_3 A[3],
		      const Vector_3 &q,
		      Vector_3 &x) {
    solve( A[0][0], A[0][1], A[0][2], A[1][0], A[1][1], A[1][2], 
	   A[2][0], A[2][1], A[2][2], q[0], q[1], q[2], 
	   x[0], x[1], x[2]);
  }

  // Solve for equation a*t^2+b*t+c=0 and find the roots between 0 and 1.
  std::pair<double,double>
  solve_quadratic_func( double a, double b, double c);

  double sign( double x) {
    if ( x==0) return 0;
    if ( x>0) return 1;
    return -1;
  }

  static Vector_2 proj( const Vector_3 &v, 
			const Vector_3 &d1, const Vector_3 &d2) {
    return Vector_2( v*d1, v*d2);
  }

  static double eval_angle( const Vector_3 &v1, const Vector_3 &v2) {
    double sqrnrm = v1.squared_norm()*v2.squared_norm();
    if ( sqrnrm==0) return 0;
    double s=v1*v2/std::sqrt( sqrnrm);
    if (s>1) s=1; else if (s<-1) s =-1;
    return std::acos(s);
  }

  static double eval_angle( const Vector_2 &v1, const Vector_2 &v2) {
    double sqrnrm = v1.squared_norm()*v2.squared_norm();
    if ( sqrnrm==0) return 0;
    double s=v1*v2/std::sqrt( sqrnrm);
    if (s>1) s=1; else if (s<-1) s =-1;
    return std::acos(s);
  }

  // Compute weighted one-sided normal at a vertex
  Vector_3 compute_weighted_normal( Vector_3 es[3], 
 				    const Vector_3 &lambdas, 
 				    int trank,
 				    const Vector_3 &mean_nrm, 
 				    const Vector_3 &face_nrm);
   
  // Compute eigenvalues and eigenvectors of a squared matrix A of order 3. 
  // At output, the eigenvalues are saved in lambdas in decending order,
  // and A is replaced by the orthonormal eigenvectors.
  static void compute_eigenvectors( Vector_3 A[3], 
				    Vector_3 &lambdas);

  bool is_acute_ridge( const Vector_3 *es, const Vector_3 &lambdas,
		       const Vector_3 &mn, int nrank) const {
    return  nrank!=1 && std::abs(es[0]*mn)<=std::abs(es[1]*mn);
  }
  
  bool is_acute_corner( const Vector_3 *es, const Vector_3 &lambdas, 
			const Vector_3 &mn, int nrank) const {
    
    return nrank==3 && ( std::abs(es[0]*mn)<=std::abs(es[1]*mn) ||
			 std::abs(es[0]*mn)<=std::abs(es[2]*mn));
  }

  
  bool is_strong( const Vector_3 *es, const Vector_3 &lambdas,
		  const Vector_3 &mn, int nrank) const;


  // Obtain nonuniform weight for surface fairing.
  double get_nonuniform_weight( const Vector_3 &lambdas, int trank);

protected:
  bool           _wf_expn;     //< Whether expansion is wavefrontal
  char           _wght_scheme; //< Weighting scheme
  char           _nrm_dfsn;    //< Number of iterations to perform diffusion
  char           _smoother;         //< Choice of mesh smoother.

  double         _courant;     //< Constant in front of CFL condition (0<c<1)
  double         _inv_sqrt_c;  //< reciprocal of squared root of c

  double         _dir_thres_weak;   //< Threshold for weak features
  double         _dir_thres_strong; //< Threshold for strong featuers

  double         _tol_angle_weak;   //< Threshold for weak dih-angles \theta_f
  double         _tol_angle_strong; //< Threshold for strong dih-angle \theta_F

  int            _tol_kstrong;      //< number of strong edges.
  double         _tol_kangle;       //< Strong-angle threshold
  double         _tol_eangle;       //< Semi-strong end-edge angle threshold

  double         _tol_cos_osta;     //< Threshold for turning angle
  double         _tol_turn_angle;   //< Threshold for turning angle
  double         _tol_ad;           //< Threshold for angle defect
  double         _eig_thres;        //< Threshold for range-null psace.

  double         _eig_weak_lb;      //< Lower bound for weak vertices.
  double         _eig_weak_ub;      //< Upper bound for weak vertices.

  bool           _conserv;          //< Whether or not to conserve mass
  bool           _is_strtd;         //< Whether it is for structured mesh
  bool           _feature_layer;    //< Whether to create a feature layer

  COM::DataItem *_As;         //< Stores the covariant matrix for each vertex
  COM::DataItem *_boffset;    //< Store the right-hande side of offset quadric
  COM::DataItem *_bmedial;    //< Store the right-hande side of medial quadric
  COM::DataItem *_eigvalues;  //< Stores eigenvalues
  COM::DataItem *_eigvecs;    //< Stores orthonormal eigenvectors for each
                               //< vertex in increasing order of eigenvalues

  COM::DataItem *_vnormals;      //< Stores vertex normals from eigen-decompo
                                  //< used for projection
  COM::DataItem *_facenormals;   //< Stores the normals of each offset face.
  COM::DataItem *_facecenters;   //< Stores the centers of each offset face.
  COM::DataItem *_faceareas;     //< Areas of each face
  COM::DataItem *_vcenters;      //< Vector between center of mass and vertex

  COM::DataItem *_tangranks;     //< Stores the ranks of tangent spaces
  COM::DataItem *_ctangranks;    //< Stores the ranks of tangent spaces with constraints
  COM::DataItem *_weak;          //< Marks whether a vertex is weak feature.
  COM::DataItem *_strong;        //< Marks whether a vertex is strong feature.
  COM::DataItem *_ridges;        //< List ridge edges.
  COM::DataItem *_ridgeneighbors;//< List neighbor vertices.

  COM::DataItem *_scales;    //< Stores rescaling factors
  COM::DataItem *_weights;   //< Buffer for storing weights

  std::vector<std::set< Edge_ID> > _edges; // ridge edges

  static const double pi;
};

PROP_END_NAMESPACE

#endif






