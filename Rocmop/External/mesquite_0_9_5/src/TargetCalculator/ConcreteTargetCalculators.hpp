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

/*! \file ConcreteTargetCalculators.hpp

Header file for the Mesquite::ConcreteTargetCalculator class

  \author Thomas Leurent
  \date   2004-09-31
 */


#ifndef ConcreteTargetCalculators_hpp
#define ConcreteTargetCalculators_hpp

#include "Mesquite.hpp"
#include "TargetCalculator.hpp"
#include "LVQDTargetCalculator.hpp"
#include "WTargetCalculator.hpp"

namespace Mesquite
{
  /*! \class ShapeGuides811
    \brief Shape Improvement with Unit Aspect Ratio. Use with sR-DFT
  */
  class ShapeGuides811 : public WTargetCalculator
  {
  public:
    ShapeGuides811()
    {
      guideMatrix = Ad;
    }      

    //! virtual destructor ensures use of polymorphism during destruction
    virtual ~ShapeGuides811()
      {};
  };

  /*! \class ShapeGuides812
    \brief  Shape Improvement with non-Unit Aspect Ratio. Use with sR-DFT.
  */
  class ShapeGuides812 : public LVQDTargetCalculator
  {
  public:
    ShapeGuides812(MeshSet* ref_mesh)
    {
      refMesh = ref_mesh;
      lambdaBase = REGULAR; 
      guideLambda = Ad;
      guideV = Ad;
      guideQ = Ad;
      guideDelta = A0;
    }
    
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~ShapeGuides812()
      {};
  };


  /*! \class ShapeSizeGuides821
    \brief  Shape and Size Improvement with Unit Aspect Ratio and Equidistributed Size. Use with R-DFT.
  */
  class ShapeSizeGuides821 : public LVQDTargetCalculator
  {
  public:
    ShapeSizeGuides821(MeshSet* ref_mesh)
    {
      refMesh = ref_mesh;
      lambdaBase = AVERAGE; 
      guideLambda = A0;
      guideV = Ad;
      guideQ = Ad;
      guideDelta = Ad;
    }
    
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~ShapeSizeGuides821()
      {};
  };


  /*! \class ShapeSizeGuides822
    \brief  Shape and Size Improvement with Unit Aspect Ratio and Preserved Size. Use with R-DFT.
  */
  class ShapeSizeGuides822 : public LVQDTargetCalculator
  {
  public:
    ShapeSizeGuides822(MeshSet* ref_mesh)
    {
      refMesh = ref_mesh;
      lambdaBase = REGULAR; 
      guideLambda = A0;
      guideV = Ad;
      guideQ = Ad;
      guideDelta = Ad;
    }
    
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~ShapeSizeGuides822()
      {};
  };

  /*! \class ShapeSizeGuides823
    \brief  Shape and Size Improvement with non-Unit Aspect Ratio and Equidistributed Size. Use with R-DFT.
  */
  class ShapeSizeGuides823 : public LVQDTargetCalculator
  {
  public:
    ShapeSizeGuides823(MeshSet* ref_mesh)
    {
      refMesh = ref_mesh;
      lambdaBase = AVERAGE; 
      guideLambda = A0;
      guideV = Ad;
      guideQ = Ad;
      guideDelta = A0;
    }
    
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~ShapeSizeGuides823()
      {};
  };


  /*! \class ShapeSizeGuides824
    \brief  Shape and Size Improvement with non-Unit Aspect Ratio and Preserved Size. Use with R-DFT.
  */
  class ShapeSizeGuides824 : public LVQDTargetCalculator
  {
  public:
    ShapeSizeGuides824(MeshSet* ref_mesh)
    {
      refMesh = ref_mesh;
      lambdaBase = REGULAR; 
      guideLambda = A0;
      guideV = Ad;
      guideQ = Ad;
      guideDelta = A0;
    }
    
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~ShapeSizeGuides824()
      {};
  };


  /*! \class RezoneGuides831 */

  /*! \class RezoneGuides832
    \brief  Rezone with Angle Improvement. Use with I-DFT.
  */
  class RezoneGuides832 : public LVQDTargetCalculator
  {
  public:
    RezoneGuides832(MeshSet* ref_mesh)
    {
      refMesh = ref_mesh;
      lambdaBase = REGULAR; 
      guideLambda = A0;
      guideV = A0;
      guideQ = Ad;
      guideDelta = A0;
    }
    
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~RezoneGuides832()
      {};
  };

  /*! \class RezoneGuides833 */

  /*! \class RezoneGuides834 */

  /*! \class DeformingDomainGuides841
    \brief Deforming Domain Mesh Tracking. Use with I-DFT or R-DFT
  */
  class DeformingDomainGuides841 : public WTargetCalculator
  {
  public:
    DeformingDomainGuides841(MeshSet* ref_mesh)
    {
      refMesh = ref_mesh;
      guideMatrix = Ar;
    }      

    //! virtual destructor ensures use of polymorphism during destruction
    virtual ~DeformingDomainGuides841()
      {};
  };

  /*! \class DeformingDomainGuides842 */

  /*! \class DeformingDomainGuides843
    \brief  Deforming Domain with Angle Improvement. Use with I-DFT or R-DFT.
  */
  class DeformingDomainGuides843 : public LVQDTargetCalculator
  {
  public:
    DeformingDomainGuides843(MeshSet* ref_mesh)
    {
      refMesh = ref_mesh;
      lambdaBase = REGULAR; 
      guideLambda = Ar;
      guideV = Ar;
      guideQ = Ad;
      guideDelta = Ar;
    }
    
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~DeformingDomainGuides843()
      {};
  };


  /*! \class DeformingDomainGuides844
    \brief  Deforming Domain with Angle Improvement, non-Unit AR. 
            Use with R-DFT.
  */
  class DeformingDomainGuides844 : public LVQDTargetCalculator
  {
  public:
    DeformingDomainGuides844(MeshSet* ref_mesh)
    {
      refMesh = ref_mesh;
      lambdaBase = REGULAR; 
      guideLambda = Ar;
      guideV = Ar;
      guideQ = Ad;
      guideDelta = A0;
    }
    
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~DeformingDomainGuides844()
      {};
  };
  /*! \class MorphGuides851 */

  /*! \class MorphGuides852 */

  /*! \class GeometricCurvatureGuides861 */

  /*! \class GeometricCurvatureGuides862 */

  /*! \class GeometricCurvatureGuides863 */

  /*! \class SolutionErrorGuides871 */

  /*! \class SolutionErrorGuides872 */

  /*! \class SolutionErrorGuides873 */

  /*! \class SolutionFeatureGuides881 */

  /*! \class SolutionFeatureGuides882 */

  /*! \class SolutionAlignGuides891 */



 
} //namespace


#endif // ConcreteTargetCalculator_hpp
