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
! Purpose: Compute kinetic energy and dissipation for use in optimal LES 
!   weight computation.
!
! Description: None.
!
! Input:
!   region      Region for which energy and dissipation is to be computed
!
! Output: None.
!
! Notes: 
!   1. At present, works only for single region
!
!******************************************************************************
!
! $Id: RFLU_ComputeEnerDissOLES.F90,v 1.5 2008/12/06 08:44:29 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ComputeEnerDissOLES(region)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModParameters
  USE ModTools, ONLY: MakeNonZero

  USE RFLU_ModOLES
   
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: ic,ifc,iv1,iv2

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: a,c1g,c2g,d,errorFlag,hLoc1,ic1l,ic2l,ifcp,j,l,nCells,nFaces,vLoc1
  INTEGER, DIMENSION(:), POINTER :: f2fpOLES
  INTEGER, DIMENSION(:,:), POINTER :: fsOLES
  REAL(RFREAL) :: avgFac,corr,enerOLES,numer,sum1,sum2,sum3,term,uAvg,var, & 
                  vAvg,wAvg
  REAL(RFREAL), DIMENSION(:), POINTER :: vol  
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv
  REAL(RFREAL), DIMENSION(:,:,:), ALLOCATABLE :: int2OLES 
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: int20OLES,int21OLES   
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ComputeEnerDissOLES.F90,v $ $Revision: 1.5 $'

  global => region%global

  CALL RegisterFunction(global,'RFLU_ComputeEnerDissOLES',&
  'RFLU_ComputeEnerDissOLES.F90')

! ******************************************************************************
! Start, check that have primitive variables
! ******************************************************************************

  IF ( region%mixt%cvState == CV_MIXT_STATE_CONS ) THEN 
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! region

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  vol => region%grid%vol
  cv  => region%mixt%cv

  enerOLES = global%enerOLES

  fsOLES   => region%grid%fsOLES
  f2fpOLES => region%grid%f2fpOLES

  int20OLES => region%grid%int20OLES
  int21OLES => region%grid%int21OLES

  nCells = region%grid%nCells
  nFaces = region%grid%nFaces  
  avgFac = 1.0_RFREAL/REAL(nCells,KIND=RFREAL)

! ******************************************************************************
! Compute kinetic energy (average over cells)
! ******************************************************************************

  enerOLES = 0.0_RFREAL

  DO ic = 1,nCells
    enerOLES = enerOLES + cv(CV_MIXT_XVEL,ic)*cv(CV_MIXT_XVEL,ic) & 
                        + cv(CV_MIXT_YVEL,ic)*cv(CV_MIXT_YVEL,ic) &
                        + cv(CV_MIXT_ZVEL,ic)*cv(CV_MIXT_ZVEL,ic) 
  END DO ! ic

  global%enerOLES = avgFac*0.5_RFREAL*enerOLES  

! DEBUG
  uAvg = 0.0_RFREAL
  vAvg = 0.0_RFREAL
  wAvg = 0.0_RFREAL

  DO ic = 1,nCells
    uAvg = uAvg + cv(CV_MIXT_XVEL,ic)*cv(CV_MIXT_XVEL,ic)
    vAvg = vAvg + cv(CV_MIXT_YVEL,ic)*cv(CV_MIXT_YVEL,ic)
    wAvg = wAvg + cv(CV_MIXT_ZVEL,ic)*cv(CV_MIXT_ZVEL,ic)        
  END DO ! ic 
  
  uAvg = avgFac*uAvg
  vAvg = avgFac*vAvg
  wAvg = avgFac*wAvg
  
  global%uVarOLES = uAvg
  global%vVarOLES = vAvg
  global%wVarOLES = wAvg       
! END DEBUG

! ******************************************************************************
! Compute dissipation rate
! ******************************************************************************

  nCells = SIZE(region%grid%fsOLES,1)
  avgFac = 3.0_RFREAL/REAL(nFaces,KIND=RFREAL)

  ALLOCATE(int2OLES(3,3*nCells,3*nCells),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'int2OLES')
  END IF ! global%error

  int2OLES(:,:,:) = 0.0_RFREAL

! ==============================================================================
! Compute integral through loop over faces
! ==============================================================================

  DO ifc = 1,nFaces
    ifcp = f2fpOLES(ifc) ! get prototype face

! ------------------------------------------------------------------------------
!   Loop over cells
! ------------------------------------------------------------------------------
 
    DO ic1l = 1,nCells
      c1g = fsOLES(ic1l,ifc)
      d = ic1l
    
      DO ic2l = 1,nCells
        c2g = fsOLES(ic2l,ifc)      
        a = ic2l
                         
        term = vol(c1g)*vol(c2g)

! ----- Loop over velocity components ------------------------------------------

        DO iv1 = 1,3
          l = iv1

          DO iv2 = 1,3
            j = iv2              

! --------- Compute correlation              

            corr = term*cv(iv1+1,c1g)*cv(iv2+1,c2g)

! --------- Determine storage location for correlation              

            vLoc1 = RFLU_GetI1PosOLES(l,d)              
            hLoc1 = RFLU_GetLPosOLES(j,a)

! --------- Store correlation              

            int2OLES(ifcp,vLoc1,hLoc1) = int2OLES(ifcp,vLoc1,hLoc1) + corr                                            

          END DO ! iv2
        END DO ! iv1
        
      END DO ! icl2
    END DO ! icl1

  END DO ! ifcp  

! ==============================================================================
! Normalize and average integrals
! ==============================================================================
                                      
  int2OLES(:,:,:) = avgFac*int2OLES(:,:,:)
   
! ==============================================================================
! Attempt to extract dissipation from I20 and I21 
! ==============================================================================   
   
  sum1 = 0.0_RFREAL
  sum2 = 0.0_RFREAL
  sum3 = 0.0_RFREAL 
   
  DO ifcp = 1,3
 
    DO ic1l = 1,nCells
      d = ic1l
    
      DO ic2l = 1,nCells
        a = ic2l
      
! 1x1x2 stencil      
        IF ( ic1l /= ic2l ) THEN
! 1x1x4 stencil, innermost cells - should give the same as 1x1x2 
!        IF ( (ic1l == 2 .OR. ic1l == 3) .AND. (ic2l == 3 .OR. ic2l == 2) .AND. & 
!             (ic1l /= ic2l) ) THEN
! 1x1x4 stencil, outermost cells
!        IF ( (ic1l == 1 .OR. ic1l == 4) .AND. (ic2l == 4 .OR. ic2l == 1) .AND. & 
!             (ic1l /= ic2l) ) THEN
! 1x1x6 stencil, innermost cells  
!        IF ( (ic1l == 3 .OR. ic1l == 4) .AND. (ic2l == 4 .OR. ic2l == 3) .AND. & 
!             (ic1l /= ic2l) ) THEN 
! 1x1x6 stencil, 2nd-innermost cells  
!        IF ( (ic1l == 2 .OR. ic1l == 5) .AND. (ic2l == 5 .OR. ic2l == 2) .AND. & 
!             (ic1l /= ic2l) ) THEN                                     
! 1x1x6 stencil, outermost cells  
!        IF ( (ic1l == 1 .OR. ic1l == 6) .AND. (ic2l == 6 .OR. ic2l == 1) .AND. & 
!             (ic1l /= ic2l) ) THEN
! 1x1x8 stencil, innermost cells  
!        IF ( (ic1l == 4 .OR. ic1l == 5) .AND. (ic2l == 5 .OR. ic2l == 4) .AND. & 
!             (ic1l /= ic2l) ) THEN  
! 1x1x8 stencil, 2nd-innermost cells
!        IF ( (ic1l == 3 .OR. ic1l == 6) .AND. (ic2l == 6 .OR. ic2l == 3) .AND. & 
!             (ic1l /= ic2l) ) THEN
! 1x1x8 stencil, 3rd-innermost cells             
!        IF ( (ic1l == 2 .OR. ic1l == 7) .AND. (ic2l == 7 .OR. ic2l == 2) .AND. & 
!             (ic1l /= ic2l) ) THEN                                          
! 1x1x8 stencil, outermost cells             
!        IF ( (ic1l == 1 .OR. ic1l == 8) .AND. (ic2l == 8 .OR. ic2l == 1) .AND. & 
!             (ic1l /= ic2l) ) THEN     
          DO iv1 = 1,3
            l = iv1

            iv2 = iv1        
            j   = iv2

            vLoc1 = RFLU_GetI1PosOLES(l,d)              
            hLoc1 = RFLU_GetLPosOLES(j,a)

            sum1 = sum1 + int2OLES(ifcp,vLoc1,hLoc1)
            sum2 = sum2 + int20OLES(ifcp,vLoc1,hLoc1)
            sum3 = sum3 + int21OLES(ifcp,vLoc1,hLoc1)            
          END DO ! iv1          
        END IF ! ic1l
        
      END DO ! ic2l  
    END DO ! ic1l 

  END DO ! ifcp
   
! NOTE: Assume that rhoOLES has a constant value throughout solution region, so
! use rhoOLES(1). This will need to be modified if the grids become non-uniform    

  var   = 2.0_RFREAL*global%enerOLES/3.0_RFREAL
  numer = sum1/region%grid%deltaOLES**6 - var*sum2
  global%dissOLES = (MAX(numer/sum3,TINY(1.0_RFREAL)))**(3.0_RFREAL/2.0_RFREAL)/ & 
                    region%grid%rhoOLES(1)

  WRITE(*,'(3(2X,E15.8))') sum1/region%grid%deltaOLES**6,var*sum2,sum3
  WRITE(71,'(4(2X,E15.8))') global%currentTime,sum1/region%grid%deltaOLES**6, & 
                            var*sum2,sum3

! DEBUG
  WRITE(*,'(5(2X,E15.8))') global%enerOLES,global%dissOLES, & 
                           uAvg/(MakeNonZero(3.0_RFREAL*var)), & 
                           vAvg/(MakeNonZero(3.0_RFREAL*var)), & 
                           wAvg/(MakeNonZero(3.0_RFREAL*var))                                       
! END DEBUG


  DEALLOCATE(int2OLES,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'int2OLES')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeEnerDissOLES

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeEnerDissOLES.F90,v $
! Revision 1.5  2008/12/06 08:44:29  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/01/22 16:04:33  haselbac
! Changed declaration to eliminate warning on ALC
!
! Revision 1.2  2002/10/08 15:49:29  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.1  2002/09/09 15:39:57  haselbac
! Initial revision
!
!******************************************************************************







