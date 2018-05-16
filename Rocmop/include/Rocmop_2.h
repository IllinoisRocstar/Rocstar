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
// $Id: Rocmop_2.h,v 1.4 2008/12/06 08:45:24 mtcampbe Exp $

/** \file Rocmop_2.h
 */

#ifndef __ROCMOP_H_
#define __ROCMOP_H_

#include "mopbasic.h"
#include "roccom_devel.h"
#include "Pane_communicator.h"
#include "commpi.h"
#include "Manifold_2.h"
#include "Rocsurf.h"
#include "Dual_connectivity.h"

MOP_BEGIN_NAMESPACE

//! A Roccom mesh optimization module.
class Rocmop : public SURF::Rocsurf {
public:
  
  //! Valid Smoothing Methods
  enum { SMOOTH_VOL_MESQ_WG = 0, SMOOTH_VOL_MESQ_NG, 
	 SMOOTH_SURF_MEDIAL, SMOOTH_NONE};

  //! Valid Wrapper Choices
  enum { WRAPPER_SHAPE = 0, WRAPPER_BOEING, WRAPPER_MAX};

  /** \name Contructors and Destructors
   * \{
   */

public:

  //! Default Constructor.
  Rocmop() : _reorthog(true), _wght_scheme(SURF::E2N_ANGLE),
    _eig_thres(1.e-12),_dir_thres(.1), _saliency_crn(.1),
    _rediter(1), _cnstr_types(NULL), _cnstr_dirs(NULL),
    _usr_window(NULL), _buf_window(NULL), _is_pmesh(false), 
    _method(SMOOTH_VOL_MESQ_WG),_wrapper(WRAPPER_SHAPE), _verb(0), _lazy(0), _tol(165.0),
    _niter(1), _ctol(.99), _ncycle(1), 
    _cookie(MOP_COOKIE),_invert_tets(0), _invert_hexes(0),
    _maxdisp(0.0), _smoothfreq(1),
    _disp_thresh(0.0)
  {;}

  //! Destructor
  virtual ~Rocmop();
  //\}

  /** \name Load and Unload
   * \{
   */

public:

  //! Loads Rocmop onto Roccom with a given module name.
  static void load( const std::string &mname);

  //! Unloads Rocmop from Roccom.
  static void unload( const std::string &mname);
  //\}

  void read_config_file(const std::string &);
  /** \name Smoothing control flow
   * \{
   */

public:

  //! Smooth the mesh in a Rocmop managed buffer.
  /** 
   * Creates a copy of the mesh, which is then optimized in place via
   * perform_smoothing().
   *
   * \param pmesh The pmesh to be smoothed.
   * \param disp The attribute where displacements to new positions
   *  should be stored.
   * \param timestep If a timestep is supplied, then disp will contain
   *  velocity vectors instead of displacement vectors
   */
  void smooth(const COM::Attribute *pmesh,
	      COM::Attribute *disp,
	      double* timestep = NULL);

protected:

  //! Perform smoother specific initialization.
  /** 
   *  Perform smoother specific initializing, for example initializing the
   *  Window_manifold_2 for a surface mesh, adding smoother specific attributes
   *  to the buffer window, etc.
   */
  void smoother_specific_init();

  //! Perform smoothing on _buf_window.
  /**
   *  Calls smoother_specific_init(), and then selects between 
   *  perform_iterative_smoothing() and perform_noniterative_smoothing()
   ** as required by the currently selected smoothing method.
   *
   */
  void perform_smoothing();

  //! Perform iterative smoothing.
  /** 
   *  Run an iterative smoothing method if needed until either
   *  the maximum number of iterations is reached or the required
   *  quality is achieved.
   */
  void perform_iterative_smoothing();

  //! Perform noniterative smoothing.
  void perform_noniterative_smoothing();

  //! Update real nodal coordinates of _buf_window from _usr_window
  /** 
   *  Use real nodal coordinates of _usr_window to update 
   *  real nodal coordinates of _buf_window
   */
  void update_buf_real_nc();

  //! Check pconn block 3 of the input mesh
  /** 
   *  Check if all ghost nodes are listed in pconn block 3
   *  of _usr_window
   */
  int check_input_pconn();

  //! Get displacement in _usr_window based on nc in _buf_window
  /** 
   *  Use nodal coordinates in _buf_window as target and 
   *  nodal coordinates in _usr_window as initial position.
   *  Determine displacement in _usr_window.
   */
  void get_usr_disp(const COM::Attribute *pmesh,
			  COM::Attribute *disp,
			  double* timestep);

  //\} 

  /** \name Miscellaneous
   * \{
   */

public:

  //! Set a Rocomp option.
  /** 
   * \param opt the option to be set.
   * \param val the new value for the option.
   * 
   * Option List:
   *
   * - name (data type, default value) 
   *
   * - method (int, SMOOTH_NONE)
   *   smoothing method to use
   *     - SMOOTH_VOL_MESQ_WG (Volume smoothing via mesquite with ghost)
   *     - SMOOTH_VOL_MESQ_NG (Volume smoothing via mesquite no ghosts) 
   *     - SMOOTH_SURF_MEDIAL      (Surface smoothing via medial quadric)
   *     - SMOOTH_NONE        (None, default).
   * - verbose (int, 0)
   *   verbose level, an integer greater than -1
   * - lazy (int, 0)      
   *   check quality before smoothing? 0 or 1
   * - tol (float, 0.0)       
   *   looping quality tolerance [0,180]
   * - niter (int, 1)     
   *   max # iterations for smoothing
   * - ctol (float, .99)    
   *   convergence subcycle tolerance [0,1]
   * - ncycle (int, 1)
   *   max # subcycles for convergence
   * - inverted (int, 0) 
   *   tets are inverted? 0 or 1
   * - monotone (int, 0)
   *   apply a non-decreasing invariant?
   */
  void set_value(const char* opt, const void* val);  
  
  //! Determine which nodes and elements are on pane borders.
  void determine_pane_border();

  //! Determine which nodes are on the physical border.
  void determine_physical_border(COM::Attribute* w_is_surface);

  //! Determine which nodes and elements are on the physical border.
  void determine_physical_border();

  //! Determine which nodes and elements are shared.
  void determine_shared_border();

  //! Mark the nodes which contain marked elems.
  void mark_elems_from_nodes(std::vector<std::vector<bool> > &marked_nodes,
			     std::vector<std::vector<bool> > &marked_elems);

protected:
  //! Obtain a reference to the manifold.
  virtual SURF::Window_manifold_2 *manifold()
  { return _wm;}
  
  //! Check that the object is valid.
  int validate_object() const {
    if ( _cookie != MOP_COOKIE) return -1;
    else return COM_Object::validate_object();
  }
  
  //! Repair inverted tets or hexes.
  void invert_elements(int conn_type);
  
  /** \name Inter-Pane Communication
   * \{
   */
protected:
  //! Perform a sum-reduction on the shared nodes for the given attribute.
  static void reduce_sum_on_shared_nodes(COM::Attribute *att);

  //! Agree on an integer value across all panes
  /**
   * \param val The local value.
   * \param op The type of reduction to perform.
   */
  void agree_int(int& val, MPI_Op op){
    if(COMMPI_Initialized()){
      int temp = val;
      MPI_Allreduce(&val, &temp, 1, MPI_INT, op, 
		    _usr_window->get_communicator());
      val = temp;
    }
  }

  //! Agree on a double value across all panes
  /**
   * \param val The local value.
   * \param op The type of reduction to perform.
   */
  void agree_double(double& val, MPI_Op op){
    if(COMMPI_Initialized()){
      double temp = val;
      MPI_Allreduce(&val, &temp, 1, MPI_DOUBLE, op, 
		    _usr_window->get_communicator());
      val = temp;
    }
  }
  //\}

  /** \name Smoothing Methods
   * \{
   */

protected:

  //! Smooth a volume via Mesquite using ghost information
  /** 
   *  This method operates in two steps. First the panes, including ghost nodes,
   *  are smoothed individually using Mesquite. Then, shared nodes are moved
   *  to the average of their positions across all panes.  These steps two steps
   *  repeat until either ncycles sub-cycles are completed or the loss in 
   *  quality between the two steps is less than ctol.
   *
   *  \param pre_quality Mesh quality prior to smoothing.
   */
  void smooth_vol_mesq_wg();

  //! Smooths a volume using Mesquite with only shared node information.
  /**
   *  This method operates in two steps.  First the shared nodes which are not
   *  on the physical boundary of the mesh are redistributed via element based
   *  laplacian smoothing.  Then, the real nodes of the panes are smoothed
   *  individually with pane boundaries fixed.
   *
   *  \param pre_quality Mesh quality prior to smoothing.
   */
  void smooth_vol_mesq_ng(){
    ;
  }

  //! Smooths a surface using the medial quadric.
  void smooth_surf_medial();
  
  void smooth_boeing(COM::Attribute* att, int* niter);

  //! Smooth the panes of the buffer window using MESQUITE.
  /**
   *  This method creates MesqPane attributes for each of the local panes, which
   *  are then smoothed serially.
   *
   * \param ghost_level whether to use ghost cells and nodes.
   */
  void smooth_mesquite(std::vector<COM::Pane*> &allpanes,
		       int ghost_level=0);
  //\}

  /** \name Quality Checking
   * \{
   */

protected:

  //! Get the largest dihedral angle of marked real elements
  /**
   * \param marked_elems the elements whose quality is to be checked.
   * \param allpanes the set of panes containing the elements
   */
  double check_marked_elem_quality(std::vector<std::vector<bool> > &marked_elems,
				   std::vector<COM::Pane*> &allpanes);

  //! Get the largest dihedral angle of all real elements
  /**
   * \param allpanes the set of panes containing the elements
   */
  double check_all_elem_quality(std::vector<const COM::Pane*> &allpanes,
				bool with_ghost = false);

  //! Get the largest dihedral angle of all real elements marked
  /**
   * \param allpanes the set of panes containing the elements
   */
  double check_all_elem_quality(std::vector<COM::Pane*> &allpanes,
				bool with_ghosts = false);

  //! Single process print message if given verbosity level is exceeded.
  void print_legible(int verb, const char* msg);

  //! Contrain displacements to _maxdisp
  void constrain_displacements(COM::Attribute * w_disp);
  
  //! Obtain the min and max dihedral angles
  void obtain_extremal_dihedrals(const COM::Attribute* att, double *min, double *max);

  //! Print the min and max dihedral angles along with their locations
  void print_extremal_dihedrals(COM::Window * window);

  //! Print the quality range of all elements, for debugging.
  void print_quality(std::string &s,const std::string &s2);

  //! Print the quality range of marked elements, for debugging.
  void print_mquality(std::string &s,std::vector<std::vector<bool> > &to_check);
  //\}

  //! Randomly perturn stationary nodes on pane boundary, not on phys. surface.
  void perturb_stationary();

  // Check and tally displacements for triggering smoothing
  bool check_displacements(COM::Attribute *disp);

  // Zero displacement array
  void zero_displacements(COM::Attribute *disp);

public:

  // Add aspect ratio measures to the users mesh.
  void add_aspect_ratios(COM::Attribute *usr_att, 
			 COM::Attribute *buf_att =NULL);

  /** \name Functions exclusive to medial quadric smoothing.
   * \{
   */

 protected:

   /// Get orthonormals of the constraints. If nnrms is 0, then not constrained.
   /// If nnrms is 3, then the point is fixed. Otherwise, the point is 
   /// constrained either in a plane or a line.
   static void get_constraint_directions( int type, const Vector_3<double> &dir,
					  int &ndirs, Vector_3<double> dirs[2]);

   /// Evaluate face normals (copied from FaceOffset_3.[hC]
   void evaluate_face_normals();

   /// Identify ridge edges
   void identify_ridge_edges();

  /// Redistribute ridge vertices within their tangent spaces
  void redistribute_vertices_ridge();

  /// Redistribute smooth vertices within their tangent spaces
  void redistribute_vertices_smooth();
  
  /// Compute medial quadric for every vertex.
  void compute_medial_quadric();
  
  // Solve the minimization problem (x^T)Ax+2b^Tx with eigen-decomposition
  // at a vertex, return the codimension of the vertex, and save mean
  // normal in x_nz.
  int eigen_analyze_vertex( Vector_3<double> A_io[3], Vector_3<double> &b_io, 
			    Vector_3<double> *nrm_nz, double angle_defect);
  
  void get_redist_safe_factor( COM::Attribute *c_attr, COM::Attribute *w_attr, 
			       int rank);

  // Compute eigenvalues and eigenvectors of a squared matrix A of order 3. 
  // At output, the eigenvalues are saved in lambdas in decending order,
  // and A is replaced by the orthonormal eigenvectors.
  static void compute_eigenvectors( Vector_3<double> A[3], 
				    Vector_3<double> &lambdas);
  
  // Compute eigenvalues and eigenvectors of a squared matrix A of order 2.
  // At output, the eigenvalues are saved in lambdas in decending order,
  // and A is replaced by the orthonormal eigenvectors.
  static void compute_eigenvectors( Vector_2<double> A[2], 
				     Vector_2<double> &lambdas);

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
  static void solve ( const Vector_3<double> A[3],
		      const Vector_3<double> &q,
		      Vector_3<double> &x) {
    solve( A[0][0], A[0][1], A[0][2], A[1][0], A[1][1], A[1][2], 
	   A[2][0], A[2][1], A[2][2], q[0], q[1], q[2], 
	   x[0], x[1], x[2]);
  }
  //\}
  
 protected: // Variables specific to medial quadric smoothing

  bool _reorthog;                           ///< Reorthogonalize?
  char _wght_scheme;                        ///< Weighting scheme
  double _eig_thres;                        ///< Eigenvalue thresholds
  double _dir_thres;                        ///< Another threshold
  double _saliency_crn;
  int _rediter;                             ///< No.iterations for vertex redistribution
   const COM::Attribute *_cnstr_types;      ///< Stores types of nodal constraints
   COM::Attribute *_cnstr_dirs;             ///< Stores directions of nodal contraints
   std::vector<std::set< Edge_ID> > _edges; ///< ridge edges

protected: // Variables not currently modifiable via set_option
   
   //std::vector<std::vector<bool> >  _elems_to_check; ///< Check quality of these.
   std::vector<std::vector<bool> > _is_shared_node;   ///< Is the node shared?
   std::vector<std::vector<bool> > _is_shared_elem;   ///< Does the element contain shared nodes?
   std::vector<std::vector<bool> > _is_phys_bnd_node; ///< Is the node on the physical boundary?
   std::vector<std::vector<bool> > _is_phys_bnd_elem; ///< Does the element contain nodes on the phys. boundary?
   std::vector<std::vector<bool> > _is_pane_bnd_node; ///< Is the node on the pane boundary?
   std::vector<std::vector<bool> > _is_pane_bnd_elem; ///< Does the element contain nodes on the pane boundary?

   const COM::Window           *_usr_window;         ///< The user's window
   COM::Window                 *_buf_window;          ///< The buffer window

   std::vector<MAP::Pane_dual_connectivity*> _dcs;   // Pane dual connectivities 
   
  bool                         _is_pmesh;            ///< pmesh or mesh ?

 protected: // Variables which can be modified via set option.

   enum { MOP_COOKIE=762667};

  int _method; ///< Choice of smoothing method.      
  // 0 == mesquite
  // 1 == mesquite, no ghosts
  // 2 == medial_quadric
  // 3 == error (default)

  int _wrapper; ///< Choice of Mesquite Smoothing Wrappers
  // 0 == ShapeImprovement Wrapper (Feasible Newton)
  // 1 == CGWrapper (Conjugate Gradient)

  int _verb; ///< Verbose level
  // 0 = no debugging info
  // 1 = routine entries
  // 2 = routine entries and exits
  // 3 = highest level in routine info
  // 4 = low level in routine info
  // 5 = 4 + high level MesqPane info
  // 6 = 5 + low level MesqPane info
  
   int _lazy; ///< Check quality before smoothing?

   float _tol; ///< Smoother iterating tolerance
   // Smoother loops until either 
   //   worst element quality > tol 
   //   OR _niter iterations reached

   int _niter; ///< Maximum number of iterations for smoothing

   float _ctol; ///< Subcycling tolerance 
   // Parallel mesquite subcycles for convergence until either 
   //      (post communication quality)/(pre communication quality) > _tol
   //      OR _ncycle cycles performed.
   // In this context, quality refers to worst quality of all real elements
   //      containing shared nodes.
   // _ctol should be set to a value between 0 and 1

   int _ncycle; ///< Max number of subcycles for convergence.

   int _cookie; ///< For Roccom.

   int _invert_tets; ///< If true (default false), then invert tets

   int _invert_hexes; ///< If true (default false), then invert hexes

   float _maxdisp; ///< Maximum displacement allowed

   int _smoothfreq;         /// Smooth every _smoothfreq'th call

  float _disp_thresh; // Threshold for smoothing trigger
};

MOP_END_NAMESPACE

#endif






