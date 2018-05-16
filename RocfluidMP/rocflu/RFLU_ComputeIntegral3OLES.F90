! *********************************************************************
! * Rocstar Simulation Suite                                          *
! * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
! *                                                                   *
! * Illinois Rocstar LLC                                              *
! * Champaign, IL                                                     *
! * www.illinoisrocstar.com                                           *
! * sales@illinoisrocstar.com                                         *
! *                                                                   *
! * License: See LICENSE file in top level of distribution package or *
! * http://opensource.org/licenses/NCSA                               *
! *********************************************************************
! *********************************************************************
! * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
! * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
! * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
! * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
! * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
! * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
! * Arising FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
! * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
! *********************************************************************
!******************************************************************************
!
! Purpose: Compute integral 3 of optimal LES approach.
!
! Description: None.
!
! Input: 
!   region      Current region
!
! Output: None.
!
! Notes: 
!  1. At present, assume that velocity components are stored in positions 2-4.
!  2. Note normalization of integrals - because of use of dissipation, must 
!     compute integral after computation of dissipation.
!  3. IMPORTANT: Note that since an average over prototype faces is made, need 
!     to divide by a third of nFaces, and NOT nFaces!
!
!******************************************************************************
!
! $Id: RFLU_ComputeIntegral3OLES.F90,v 1.5 2008/12/06 08:44:29 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ComputeIntegral3OLES(region)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModTools, ONLY: MakeNonZero

  USE RFLU_ModOLES

  IMPLICIT NONE

! parameters      

  TYPE(t_region) :: region

! local variables

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: a,b,c1g,c2g,c3g,d,e,g,hLoc1,hLoc2,ic1l,ic2l,ic3l,ifc,ifcp, & 
             iv1,iv2,iv3,j,k,l,m,nCells,nFaces,vLoc1,vLoc2
  INTEGER, DIMENSION(:), POINTER :: f2fpOLES
  INTEGER, DIMENSION(:,:), POINTER :: fsOLES
  REAL(RFREAL) :: avgFac,corr,term
  REAL(RFREAL), DIMENSION(:), POINTER :: vol
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: int31OLES,int32OLES             
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start, check that have primitive variables
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ComputeIntegral3OLES.F90,v $ $Revision: 1.5 $'

  global => region%global

  CALL RegisterFunction(global,'RFLU_ComputeIntegral3OLES',&
  'RFLU_ComputeIntegral3OLES.F90')

  IF ( region%mixt%cvState == CV_MIXT_STATE_CONS ) THEN 
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! region

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  fsOLES   => region%grid%fsOLES
  f2fpOLES => region%grid%f2fpOLES

  int31OLES => region%grid%int31OLES
  int32OLES => region%grid%int32OLES

  vol => region%grid%vol
  cv  => region%mixt%cv

  nCells = SIZE(region%grid%fsOLES,1)
  nFaces = region%grid%nFaces
  avgFac = 3.0_RFREAL/REAL(nFaces,KIND=RFREAL) ! NOTE factor of 3

! ******************************************************************************
! Initialize integrals
! ******************************************************************************

  int31OLES(:,:,:) = 0.0_RFREAL
  int32OLES(:,:,:) = 0.0_RFREAL  

! ******************************************************************************
! Compute integrals through loop over faces
! ******************************************************************************

  DO ifc = 1,nFaces
    ifcp = f2fpOLES(ifc) ! get prototype face

! ==============================================================================
!   Loop over cells
! ==============================================================================
 
    DO ic1l = 1,nCells
      c1g = fsOLES(ic1l,ifc)
      d = ic1l
    
      DO ic2l = 1,nCells
        c2g = fsOLES(ic2l,ifc)      
        b = ic2l
        e = ic2l
      
        DO ic3l = 1,nCells
          c3g = fsOLES(ic3l,ifc)         
          g = ic3l
          a = ic3l
          
          term = vol(c1g)*vol(c2g)*vol(c3g)       
          
! ------------------------------------------------------------------------------
!         Loop over velocity components
! ------------------------------------------------------------------------------          
          
          DO iv1 = 1,3
            l = iv1
          
            DO iv2 = 1,3
              DO iv3 = 1,3
              
! ------------- Compute correlation              
              
                corr = term*cv(iv1+1,c1g)*cv(iv2+1,c2g)*cv(iv3+1,c3g)
              
! ------------- Determine storage location for correlation              
              
                j = iv2
                k = iv3
              
                vLoc1 = RFLU_GetI1PosOLES(l,d)              
                hLoc1 = RFLU_GetQPosOLES(j,k,b,g,nCells)
                
                m = iv2
                j = iv3
                
                vLoc2 = RFLU_GetI4PosOLES(l,m,d,e,nCells)                
                hLoc2 = RFLU_GetLPosOLES(j,a)

! ------------- Store correlation              
              
                int31OLES(ifcp,vLoc1,hLoc1) = int31OLES(ifcp,vLoc1,hLoc1) & 
                                            + corr
                int32OLES(ifcp,vLoc2,hLoc2) = int32OLES(ifcp,vLoc2,hLoc2) & 
                                            + corr                                            
              END DO ! iv3
            END DO ! iv2
          END DO ! iv1
        
        END DO ! icl3
      END DO ! icl2
    END DO ! icl1

  END DO ! ifcp  

! ******************************************************************************
! Normalize and average integrals
! ******************************************************************************

  term = avgFac/(MakeNonZero(global%dissOLES)*region%grid%deltaOLES**10)

  int31OLES(:,:,:) = term*int31OLES(:,:,:)
  int32OLES(:,:,:) = term*int32OLES(:,:,:)
    
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeIntegral3OLES

!*******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeIntegral3OLES.F90,v $
! Revision 1.5  2008/12/06 08:44:29  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2003/03/15 18:31:27  haselbac
! Added footer
!
!*******************************************************************************







