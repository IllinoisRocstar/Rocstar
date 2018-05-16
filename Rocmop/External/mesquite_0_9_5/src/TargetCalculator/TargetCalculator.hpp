/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */

/*! \file TargetCalculator.hpp

\brief The Mesquite::TargetCalculator class is the base class. Concrete classes are 
 instantiated by the user, and often implemented by the user to give 
 mesquite a measure of the perfect mesh. 

  \author Thomas Leurent
  \date   2004-09-31
 */


#ifndef TargetCalculator_hpp
#define TargetCalculator_hpp

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "Matrix3D.hpp"
#include "PatchData.hpp"
#include "PatchDataUser.hpp"

#ifdef HAVE_IEEEFP
#  include <ieeefp.h>
#endif

namespace Mesquite
{
  class PatchDataParameters;

  
  /*! \class TargetCalculator
    \brief Base class that provides the interface for computing the target corner matrices
    used in the context based smoothing.

    To implement a concrete TargetCalculator, one inherits from this class and then 
    overrides the compute_target_matrices function itself to provide corner matrices for all
    elements in the patch.

    Note that an implementation is provided in TargetCalculator::compute_default_target_matrices
    for default target matrices often used in the computation of more complex,
    reference-based target matrices.

    The target calculator is set on the QualityImprover. At runtime, it associates with
    each mesh corner (i.e. a corner of an element. Most of the time, one vertex
    corresponds to many corners.) with a TargetMatrix through an MsqTag mechanism 
    available in the MsqMeshEntity class.
   */
  class TargetCalculator : public PatchDataUser
  {
  public:

    TargetCalculator() : refMesh(0) 
      { get_all_parameters().set_global_patch_type(); }
    
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~TargetCalculator()
    {};

    static void initialize_default_target_matrices(Matrix3D &tri_M3D, Matrix3D &quad_M3D,
                                    Matrix3D &tet_M3D, Matrix3D &hex_M3D);

      //! \enum chooses whether the calculation is per element or an average
      //! for some cases of the \f$ \lambda_k \f$ coefficient.
    enum Lambda_type {
      REGULAR, //!< Each element has a lambda coefficient 
      AVERAGE  //!< The Lambda coefficient is the average on the mesh.
    };

      //! \enum chooses the type of guide matrix used in the target calculator
    enum guide_type {
      Ad, //!< Default Guide
      AK, //!< Conductivity Guide
      A0, //!< Initial Mesh Guide
      Ar, //!< Reference Mesh Guide
      As, //!< Locally-smoothed IM Guide
      Ab, //!< Bounding-box Guide
      Ac, //!< Geometric Curvature Guide
      Ap, //!< Parametric Guide
      Ae, //!< Error Indicator Guide
      Af, //!< Flux-based Guide
      Ax  //!< Solution Feature-based Guide
    };
      //! Computes the guide corner matrices A for a given element index in the reference patch.
    void compute_guide_matrices(enum guide_type type, PatchData &ref_pd, size_t elem_ind,
                                           Matrix3D A[], int num, MsqError &err);

    static double compute_Lambda(const Matrix3D &A, MsqError &err);
    Matrix3D compute_V_3D(const Matrix3D &A, MsqError &err);
    Matrix3D compute_Q_3D(const Matrix3D &A, MsqError &err);
    Matrix3D compute_Delta_3D(const Matrix3D &A, MsqError &err);
    
      //! Compute the default "isotropic" target matrices that are often used in the computation
      //! of reference-based target matrices.
      //! The resulting corner matrices are stored in tags on the elements of the PatchData.
    void compute_default_target_matrices(PatchData &pd, MsqError &err);

      //! Compute the corner matrices for the reference mesh refMesh.
      //! The refMesh data member is set by the constructors of a concrete TargetCalculator
      //! that requires a reference mesh.
    void compute_reference_corner_matrices(PatchData &pd, MsqError &err);

    //! This function wraps compute_target_matrices and checks that the determinant of each target
    //! is positive.
    void compute_target_matrices_and_check_det(PatchData& pd, MsqError& err);

    //! Reset the reference mesh so it starts from the first vertex again. 
    void reset_reference_meshset(MsqError &err);
    
    /*! \brief This function provides the corner matrices for all elements on the Patch.

         Useful functionality includes: MsqMeshEntity::set_tag, MsqTag::target_matrix,
         MsqTag::scalar .
    */
    virtual void compute_target_matrices(PatchData& pd, MsqError& err) =0;

    virtual double loop_over_mesh( MeshSet& ms, MsqError& err );
    
    virtual msq_std::string get_name();
    
    virtual AlgorithmType get_algorithm_type();

  protected:
    MeshSet* refMesh;
  };

  
  inline void TargetCalculator::initialize_default_target_matrices(Matrix3D &tri_M3D,
                                                      Matrix3D &quad_M3D,
                                                      Matrix3D &tet_M3D,
                                                      Matrix3D &hex_M3D)
  {
    static const double SIXTH_ROOT_OF_TWO = msq_stdc::pow(2., 1./6.);
    
    const double v_tri[] = {1., 0.5, 0., 0., MSQ_SQRT_THREE/2., 0., 0., 0., 1.};
    Matrix3D m1(v_tri);
    tri_M3D = m1 * MSQ_3RT_2_OVER_6RT_3;

//     const double v_tri[] = {1., 0.5, 0., 0., MSQ_SQRT_THREE/2., 0., 0., 0., sqrt(MSQ_SQRT_THREE/2.)};
//     Matrix3D m1(v_tri);
//     tri_M3D = m1 * sqrt(2./MSQ_SQRT_THREE);

    const double v_quad[] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
    Matrix3D m2(v_quad);
    quad_M3D = m2;
    
    const double v_tet[] = {1., 0.5, 0.5, 0., MSQ_SQRT_THREE/2., MSQ_SQRT_THREE/6., 0., 0., MSQ_SQRT_TWO/MSQ_SQRT_THREE};
    Matrix3D m3(v_tet);
    tet_M3D =  m3 * SIXTH_ROOT_OF_TWO;

    const double v_hex[] = {1., 0., 0.,  0., 1., 0.,  0., 0., 1.};
    Matrix3D m4(v_hex);
    hex_M3D = m4;
  }


  //! Note that this function is static, i.e. it can be used independently of an object.
  inline double TargetCalculator::compute_Lambda(const Matrix3D &A, MsqError& )
  {
    return Mesquite::cbrt(fabs(det(A)));
  }

  
  
  //!
  inline  Matrix3D TargetCalculator::compute_Q_3D(const Matrix3D &A, MsqError &err)
  {
    Vector3D a1(A[0][0], A[1][0], A[2][0]); 
    Vector3D a2(A[0][1], A[1][1], A[2][1]); 
    Vector3D a3(A[0][2], A[1][2], A[2][2]); 

    double a1_norm = A.column_length(0);
    double a2_norm = A.column_length(1);
    double a3_norm = A.column_length(2);
    
    Vector3D a1_x_a2 = a1 * a2;
    double a1_x_a2_norm = a1_x_a2.length();
    Vector3D a1_x_a3 = a1 * a3;

    double nu = Mesquite::cbrt(a1_norm*a2_norm*a3_norm);
    double det_A = det(A);
    double fac = nu / Mesquite::cbrt(fabs(det_A));
    
    Matrix3D Q;

    Q[0][0] = fac * 1.;
    Q[0][1] = fac * a1%a2 / (a1_norm*a2_norm);
    Q[0][2] = fac * a1%a3 / (a1_norm*a3_norm);
    Q[1][1] = fac * a1_x_a2_norm / (a1_norm*a2_norm);
    Q[1][2] = fac * a1_x_a2 % a1_x_a3 / (a1_x_a2_norm * a1_norm * a3_norm);
    Q[2][2] = fac * det_A / (a1_x_a2_norm*a3_norm); 
    
    if (!finite(Q[0][0]) || !finite(Q[0][1]) || !finite(Q[0][2]) ||
        !finite(Q[1][0]) || !finite(Q[1][2]) || !finite(Q[2][2])) {
      MSQ_SETERR(err)("Numerical error", MsqError::INVALID_STATE);
    }
    
    return Q;
  }

  //!
  inline  Matrix3D TargetCalculator::compute_Delta_3D(const Matrix3D &A, MsqError &err)
  {
    double a1_norm = A.column_length(0);
    double a2_norm = A.column_length(1);
    double a3_norm = A.column_length(2);
    
    double nu = Mesquite::cbrt(a1_norm*a2_norm*a3_norm);
    double fac = 1./nu ;
    if (!finite(fac)) {
      MSQ_SETERR(err)("Numerical error", MsqError::INVALID_STATE);
    }
    
    Matrix3D Delta;

    Delta[0][0] = fac * a1_norm;
    Delta[1][1] = fac * a2_norm;
    Delta[2][2] = fac * a3_norm;
    
    return Delta;
  }

  
} //namespace


#endif // TargetCalculator_hpp
