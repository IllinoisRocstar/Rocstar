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
! ******************************************************************************
!
! Purpose: Collection of routines related to GENX interaction.
!
! Description: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModRocstarTools.F90,v 1.29 2008/12/06 08:44:21 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModRocstarTools

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid 
  USE ModMPI 

  IMPLICIT NONE

  INCLUDE 'comf90.h'

  PRIVATE
  PUBLIC :: RFLU_GENX_ConstrainDisp, &
            RFLU_GENX_InitBFLAG,     &
            RFLU_GENX_MoveGrid

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================  
! Private
! ==============================================================================  

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModRocstarTools.F90,v $ $Revision: 1.29 $'       
  
! ==============================================================================  
! Public
! ==============================================================================  
  

! ******************************************************************************
! Contained routines
! ******************************************************************************

  CONTAINS 
  
  
  






! ******************************************************************************
!
! Purpose: Constrain displacements to allow burn-out simulations to case.
!
! Description: 
!
! Input: 
!   pRegion 	Pointer to region data
!
! Output: None.
!
! Notes: Surface Types: Cylindrical, Ellipsoidical (Spheroid), and Spherical
!
!        Ranges: Each surface type is applied only for a particular range of 
!                the longitudinal (long axis) coordinate. The code is currently
!                configured for the RSRM where the head end is ellipsoidical,
!                the main segments are cylindrical, and the aft is spherical.
!                The code will check if a given node is in the headend, or 
!                the aft end. If neither, it assumes cylindrical constraints.
!
!        Config: To configure the constraints, one must specify the ranges for
!                which to apply each constraint type.  The relevant params are:
!
!                Constraint Surface                       Parameters        
!                ------------------                   ----------------     
!                Cylindrical                             caseRadius    
!                Spherical                         sphereCenter,caseRadius
!                Ellipsoidical                  caseTip,ellipsLong,ellipsTrans
!
!                caseRadius   = cylindrical radius of the rocket case
!                sphereCenter = longitudinal coordinate of the sphere center
!                caseTip      = the minimum longitudinal coordinate value
!                               (corresponds to the top [headend] of the case)
!                ellipsLong   = longitudinal ellipse semi axis  (a)
!                ellipsTrans  = transverse ellipse semi axis (b and c)
!              
!                From this, the head end of the rocket is assumed to have 
!                longitudinal coordinates less than (caseTip + ellipsLong) and 
!                the aft end has longitudinal coordinates greater than the 
!                sphereCenter parameter.
!
!        Method: The length of the vector R from the constraint surface origin
!                to the node is calculated.  If the node lies inside the 
!                surface, do nothing - the node is unconstrained.  
!                If the node lies outside the surface, then R is scaled such 
!                that the new node position is *on* the constraint surface.
!               
!                After all nodes have been checked/constrained, we loop over
!                each surface element and turn off burning for those which have
!                all of their vertices on the constraint surface.  These faces
!                have "burned out".
!
!
!       Changes: Added a hole for the igniter (ignHole) in the RSRM headend    
! ******************************************************************************

SUBROUTINE RFLU_GENX_ConstrainDisp(pRegion)

   USE ModTools

   IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

   TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

   INTEGER :: cntr,errorFlag,ifl,iPatch,ivg,ivl
   INTEGER :: coordTrans1,coordTrans2,coordLong,LDir
   INTEGER :: boCount, boTot, cvCount, cvTot,info
   REAL(RFREAL) :: caseRadius,caseRadius2,radius2,theta,x,y,z
   REAL(RFREAL) :: ellipsC,ellipsLong,ellipsTrans,ellipsLong2,ellipsTrans2
   REAL(RFREAL) :: caseTip, sphereCenter, longC, ignHole
   REAL(RFREAL) :: backX, backY, tolerance, dispmag
   REAL(RFREAL) :: LMinPlane, LMaxPlane, T1MinPlane
   REAL(RFREAL) :: T1MaxPlane, T2MinPlane, T2MaxPlane
   REAL(RFREAL), DIMENSION(:,:), POINTER :: pXyz,pXyzOld
   TYPE(t_global), POINTER :: global
   TYPE(t_grid), POINTER :: pGrid,pGridOld
   TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_GENX_ConstrainDisp',&
  'RFLU_ModRocstarTools.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Constraining displacements...'
  END IF ! global%myProcid
   
  pGrid    => pRegion%grid
  pXyz    => pGrid%xyz

! ==============================================================================
! Rocket orientation specification
! StarAft and RSRM are along the X axis and StarSlice is along Z
! ==============================================================================

  coordTrans1 = global%cnstrCoordT1 ! transverse-1
  coordTrans2 = global%cnstrCoordT2 ! transverse-2
  coordLong   = global%cnstrCoordL  ! longitudinal (long axis)
!  WRITE(*,*) 'coordLong = ', coordLong
  LDir = NINT(SIGN(1.0_RFREAL,global%cnstrAftEnd - global%cnstrHeadEnd))


! ==============================================================================
! Case radius setting is common to all constraint surfaces
! ==============================================================================

  caseRadius   = global%cnstrCaseRad
  caseRadius2  = caseRadius*caseRadius

! ==============================================================================
! Ellipsoidal constraints (applied for [coord_long < ellipsC])
! ==============================================================================

  ellipsLong    = global%cnstrEllipsL 
  ellipsTrans   = global%cnstrEllipsT  
  ellipsLong2   = 1.0_RFREAL/(ellipsLong * ellipsLong)
  ellipsTrans2  = 1.0_RFREAL/(ellipsTrans * ellipsTrans)
  ellipsC       = global%cnstrHeadEnd


! ==============================================================================
! Spherical constraints (applied for [coord_long > sphereC])
! Sphere radius assumed to be equal to caseRadius
! ==============================================================================

  backX         = global%cnstrAftEnd ! End of nozzle bucket at t = 0
  backY         = global%cnstrNozY ! End of nozzle bucket at t = 0, (z = 0)
  tolerance     = global%cnstrTol2 ! Wiggle room
  sphereCenter  = backX - LDir*SQRT(caseRadius2 - backY**2)
  
! ==============================================================================
! Planar constraints
! ==============================================================================
  LMinPlane     = global%cnstrLMinPlane  
  LMaxPlane     = global%cnstrLMaxPlane  
  T1MinPlane    = global%cnstrT1MinPlane 
  T1MaxPlane    = global%cnstrT1MaxPlane
  T2MinPlane    = global%cnstrT2MinPlane
  T2MaxPlane    = global%cnstrT2MaxPlane

! ==============================================================================
! Keep track of how many nodes are constrained and how many faces burn out
! ==============================================================================
  boCount = 0 
  cvCount = 0 

! *****************************************************************************=
! Loop over patches
! *****************************************************************************=
   
  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
     
! -----------------------------------------------------------------------------   
!     Loop over vertices, check and if necessary limit new position
! -----------------------------------------------------------------------------    
    
       
    DO ivl = 1,pPatch%nBVert
      ivg = pPatch%bv(ivl)            
  
      x = pXyz(coordLong,ivg)   + pPatch%dXyz(coordLong,ivl)
      y = pXyz(coordTrans1,ivg) + pPatch%dXyz(coordTrans1,ivl)
      z = pXyz(coordTrans2,ivg) + pPatch%dXyz(coordTrans2,ivl)

! -----------------------------------------------------------------------------
! Spherical constraints (aft end)
! -----------------------------------------------------------------------------           

      IF( ( pPatch%bcCoupled == BC_BURNING ) .AND.  &
          ( (LDir > 0 .AND. x > sphereCenter ) .OR.  &
            (LDir < 0 .AND. x < sphereCenter ) ) ) THEN
        info = 1
        longC = sphereCenter
        x = x - longC
        radius2 = x*x + y*y + z*z


! -----------------------------------------------------------------------------
! Ellipsoidal constraints (head end)
! -----------------------------------------------------------------------------           

      ELSE IF( ( LDir > 0 .AND. x < ellipsC) .OR.  &
               ( LDir < 0 .AND. x > ellipsC) ) THEN 
        info = 2
        longC = ellipsC
        x = x - longC
        radius2 =(x*x*ellipsLong2+y*y*ellipsTrans2+z*z*ellipsTrans2)*&
             caseRadius2


! ------------------------------------------------------------------------------
! Cylindrical constraints (everywhere else)
! ------------------------------------------------------------------------------  

      ELSE 
        info = 3
        longC = x
        x = x - longC
        radius2 = y*y+z*z
      ENDIF ! x


! ------------------------------------------------------------------------------
! Constraint section common to all surfaces
! -----------------------------------------------------------------------------           
      IF(radius2 > caseRadius2) THEN

!RAF Do not constrain nodes that have zero displacements.  Hope this
!RAF allows normal operation while preventing negative volume in BSM,
!RAF in which nozzle protrudes beyond case radius.

        dispmag = SQRT(pPatch%dXyz(coordLong,ivl)**2 +  &
                       pPatch%dXyz(coordTrans1,ivl)**2 +  &
                       pPatch%dXyz(coordTrans2,ivl)**2 )
        IF(dispmag > 1.0e-12) THEN

          cvCount = cvCount + 1
          theta = SQRT(caseRadius2/radius2)                 
!                   if(info == 1) then
!                    WRITE(*,*) 'Aft end constraints:',x,y,z,theta
!                 endif
!                   if(info == 2) then
!                     WRITE(*,*) 'Head end constraints:',x,y,z,theta
!                  endif
!                   if(info == 3) then
!                     WRITE(*,*) 'Bore constraints:',x,y,z,theta
!                  endif
         
          x = theta*x
          y = theta*y
          z = theta*z

        ENDIF ! dispmag

        x = x + longC
        pPatch%dXyz(coordLong,ivl)   = x - pXyz(coordLong,ivg)
        pPatch%dXyz(coordTrans1,ivl) = y - pXyz(coordTrans1,ivg)
        pPatch%dXyz(coordTrans2,ivl) = z - pXyz(coordTrans2,ivg)

      ELSE ! If it is not constrained in other ways, check planar constraints
!RAF         x = x + longC
!RAF    For some reason the compiler gets a NaN for longC with the above line.
        x = pXyz(coordLong,ivg)   + pPatch%dXyz(coordLong,ivl)

!RAF    Apply planar constraints only to moving nodes

        dispmag = SQRT(pPatch%dXyz(coordLong,ivl)**2 +  &
                       pPatch%dXyz(coordTrans1,ivl)**2 +  &
                       pPatch%dXyz(coordTrans2,ivl)**2 )
        IF(dispmag > 1.0e-12) THEN

          IF(x < LMinPlane)  x = LMinPlane
          IF(x > LMaxPlane)  x = LMaxPlane
          IF(y < T1MinPlane) y = T1MinPlane
          IF(y > T1MaxPlane) y = T1MaxPlane
          IF(z < T2MinPlane) z = T2MinPlane
          IF(z > T2MaxPlane) z = T2MaxPlane
!          IF(x < LMinPlane)  then
!            WRITE(*,*) 'xmin planar constraint:',x,y,z
!            WRITE(*,*) 'longC = ',longC
!            WRITE(*,*) 'pXyz(coordLong,ivg) = ',pXyz(coordLong,ivg)
!            WRITE(*,*) 'pPatch%dXyz(coordLong,ivl) = ',pPatch%dXyz(coordLong,ivl)
!            x = LMinPlane
!          endif
!          IF(x > LMaxPlane)  then
!            WRITE(*,*) 'xamx planar constraint:',x,y,z
!            x = LMaxPlane
!          endif
!          IF(y < T1MinPlane) then
!            WRITE(*,*) 'ymin planar constraint:',x,y,z
!            y = T1MinPlane
!          endif
!          IF(y > T1MaxPlane) then
!            WRITE(*,*) 'ymax planar constraint:',x,y,z
!            y = T1MaxPlane
!          endif
!          IF(z < T2MinPlane) then
!            WRITE(*,*) 'zmin planar constraint:',x,y,z
!            z = T2MinPlane
!          endif
!          IF(z > T2MaxPlane) then
!            WRITE(*,*) 'zmax planar constraint:',x,y,z
!            z = T2MaxPlane
!          endif

        ENDIF ! dispmag

        pPatch%dXyz(coordLong,ivl)   = x - pXyz(coordLong,ivg)
        pPatch%dXyz(coordTrans1,ivl) = y - pXyz(coordTrans1,ivg)
        pPatch%dXyz(coordTrans2,ivl) = z - pXyz(coordTrans2,ivg)
         
      ENDIF ! radius2
    END DO ! ivl
        
! -----------------------------------------------------------------------------   
!     Loop over burning triangles, turn off burning if all vertices on case
! -----------------------------------------------------------------------------   

    IF ( pPatch%bcCoupled == BC_BURNING ) THEN 
      DO ifl = 1,pPatch%nBTris
        cntr = 0
           
        DO ivl = 1,3
          ivg = pPatch%bTri2v(ivl,ifl)
              
          x = pXyz(coordLong,ivg)   + pPatch%dXyz(coordLong,ivl)
          y = pXyz(coordTrans1,ivg) + pPatch%dXyz(coordTrans1,ivl)
          z = pXyz(coordTrans2,ivg) + pPatch%dXyz(coordTrans2,ivl)

          IF( (LDir > 0 .AND. x > sphereCenter) .OR.  &
              (LDir < 0 .AND. x < sphereCenter) ) THEN
            x = x - sphereCenter
            radius2 = x*x + y*y + z*z
          ELSE IF ( (LDir > 0 .AND. x < ellipsC) .OR.  &
                    (LDir < 0 .AND. x > ellipsC) ) THEN
            x = x - ellipsC
            radius2 = (x*x*ellipsLong2+y*y*ellipsTrans2+z*z*ellipsTrans2) &
                 * caseRadius2
          ELSE 
            radius2 = y*y + z*z        
          ENDIF ! x
          IF ( radius2 >= (caseRadius2-global%cnstrTol1)) THEN
            cntr = cntr + 1
          ELSE
            x = pXyz(coordLong,ivg)   + pPatch%dXyz(coordLong,ivl)
            IF(x <= LMinPlane)       THEN
               cntr  = cntr + 1
            ELSE IF(x >= LMaxPlane)  THEN
               cntr  = cntr + 1
            ELSE IF(y <= T1MinPlane) THEN 
               cntr  = cntr + 1
            ELSE IF(y >= T1MaxPlane) THEN
               cntr  = cntr + 1
            ELSE IF(z <= T2MinPlane) THEN
               cntr  = cntr + 1
            ELSE IF(z >= T2MaxPlane) THEN
               cntr  = cntr + 1
            ENDIF
          END IF ! FloatEqual

        END DO ! ivl
           
        IF ( cntr == 3 ) THEN
          boCount = boCount + 1
          pPatch%bflag(ifl)   = 0 
          pPatch%mdotAlp(ifl) = 0.0_RFREAL
!        ELSE
!          pPatch%bflag(ifl) = 1 ! make sure burning is on if not burned out
        END IF ! cntr
      END DO ! ifl

! -----------------------------------------------------------------------------    
!     Loop over quadrilaterals, turn off burning if all vertices on case
! -----------------------------------------------------------------------------   

      DO ifl = 1,pPatch%nBQuads
        cntr = 0
           
        DO ivl = 1,4
          ivg = pPatch%bQuad2v(ivl,ifl)
              
          x = pXyz(coordLong,ivg)   + pPatch%dXyz(coordLong,ivl)
          y = pXyz(coordTrans1,ivg) + pPatch%dXyz(coordTrans1,ivl)
          z = pXyz(coordTrans2,ivg) + pPatch%dXyz(coordTrans2,ivl)
              
          IF( (LDir > 0 .AND. x > sphereCenter) .OR.  &
              (LDir < 0 .AND. x < sphereCenter) ) THEN
            x = x - sphereCenter
            radius2 = x*x + y*y + z*z
          ELSE IF ( (LDir > 0 .AND. x < ellipsC) .OR.  &
                    (LDir > 0 .AND. x < ellipsC) ) THEN
            x = x - ellipsC
            radius2 = (x*x*ellipsLong2+y*y*ellipsTrans2+z*z*ellipsTrans2) &
                 * caseRadius2
          ELSE 
            radius2 = y*y + z*z        
          ENDIF ! x
          IF ( radius2 >= (caseRadius2-global%cnstrTol1)) THEN
             cntr = cntr + 1        
         ELSE
            x = pXyz(coordLong,ivg)   + pPatch%dXyz(coordLong,ivl)
            IF(x <= LMinPlane)       THEN
               cntr  = cntr + 1
            ELSE IF(x >= LMaxPlane)  THEN
               cntr  = cntr + 1
            ELSE IF(y <= T1MinPlane) THEN 
               cntr  = cntr + 1
            ELSE IF(y >= T1MaxPlane) THEN
               cntr  = cntr + 1
            ELSE IF(z <= T2MinPlane) THEN
               cntr  = cntr + 1
            ELSE IF(z >= T2MaxPlane) THEN
               cntr  = cntr + 1
            ENDIF
          END IF ! FloatEqual
        END DO ! ivl
           
        IF ( cntr == 4 ) THEN ! All vertices on case, stop motion and burning
          boCount = boCount + 1
          pPatch%bflag(ifl+pPatch%nBTris)   = 0 
          pPatch%mdotAlp(ifl+pPatch%nBTris) = 0.0_RFREAL
!        ELSE
!          pPatch%bflag(ifl)   = 1
        END IF ! cntr
      END DO ! ifl 
    END IF ! pPatch%bcCoupled   
  END DO ! iPatch              
  
! *****************************************************************************
! End
! *****************************************************************************

  IF (global%checkLevel == CHECK_HIGH) THEN
    boTot = 0
    cvTot = 0
    CALL MPI_Reduce(boCount,boTot,1,MPI_INTEGER,MPI_SUM, &
         MASTERPROC,global%mpiComm,global%error)
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
    END IF ! global%errorFlag
    CALL MPI_Reduce(cvCount,cvTot,1,MPI_INTEGER,MPI_SUM, &
         MASTERPROC,global%mpiComm,global%error)
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
    END IF ! global%errorFlag
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH) THEN
      WRITE(STDOUT,'(A,3X,A,1X,I9)') SOLVER_NAME, &
           'Total number of constrained nodes:',cvTot
      WRITE(STDOUT,'(A,3X,A,1X,I9)') SOLVER_NAME, &
           'Total number of burned out faces:',boTot
    END IF ! global%myProcid
  ENDIF ! global%checkLevel
  
  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Constraining displacements done.'
  END IF ! global%myProcid
  
  CALL DeregisterFunction(global)
  
END SUBROUTINE RFLU_GENX_ConstrainDisp


SUBROUTINE RFLU_GENX_InitBFLAG(pRegion)
  
  USE ModTools
  
  IMPLICIT NONE
  
  ! ******************************************************************************
  ! Declarations and definitions
  ! ******************************************************************************
  
  ! ==============================================================================
  ! Arguments
  ! ==============================================================================
  
  TYPE(t_region), POINTER :: pRegion
  
  ! ==============================================================================
  ! Locals
  ! ==============================================================================
  
  INTEGER :: cntr,errorFlag,ifl,iPatch,ivg,ivl
  INTEGER :: coordTrans1,coordTrans2,coordLong,LDir
  INTEGER :: boCount, boTot, cvCount, cvTot,info
  REAL(RFREAL) :: caseRadius,caseRadius2,radius2,theta,x,y,z
  REAL(RFREAL) :: ellipsC,ellipsLong,ellipsTrans,ellipsLong2,ellipsTrans2
  REAL(RFREAL) :: caseTip, sphereCenter, longC, ignHole
  REAL(RFREAL) :: backX, backY, tolerance
   REAL(RFREAL) :: LMinPlane, LMaxPlane, T1MinPlane
   REAL(RFREAL) :: T1MaxPlane, T2MinPlane, T2MaxPlane
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pXyz,pXyzOld
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid,pGridOld
  TYPE(t_patch), POINTER :: pPatch
  
  ! ******************************************************************************
  ! Start
  ! ******************************************************************************
  
  global => pRegion%global
  
  CALL RegisterFunction(global,'RFLU_GENX_InitBFLAG',&
  'RFLU_ModRocstarTools.F90')
  
  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
     WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing BFLAG...'
  END IF ! global%myProcid
  
  pGrid    => pRegion%grid
  pXyz    => pGrid%xyz
  
  ! ==============================================================================
  ! Rocket orientation specification
  ! StarAft and RSRM are along the X axis and StarSlice is along Z
  ! ==============================================================================
  IF(global%cnstrCaseRad > 0) THEN
     coordTrans1 = global%cnstrCoordT1 ! transverse-1
     coordTrans2 = global%cnstrCoordT2 ! transverse-2
     coordLong   = global%cnstrCoordL  ! longitudinal (long axis)
     LDir = NINT(SIGN(1.0_RFREAL,global%cnstrAftEnd - global%cnstrHeadEnd))
     
     
     ! ==============================================================================
     ! Case radius setting is common to all constraint surfaces
     ! ==============================================================================
     
     caseRadius   = global%cnstrCaseRad
     caseRadius2  = caseRadius*caseRadius
     
     ! ==============================================================================
     ! Ellipsoidal constraints (applied for [coord_long < ellipsC])
     ! ==============================================================================
     
     ellipsLong    = global%cnstrEllipsL 
     ellipsTrans   = global%cnstrEllipsT  
     ellipsLong2   = 1.0_RFREAL/(ellipsLong * ellipsLong)
     ellipsTrans2  = 1.0_RFREAL/(ellipsTrans * ellipsTrans)
     ellipsC       = global%cnstrHeadEnd
     
     
     ! ==============================================================================
     ! Spherical constraints (applied for [coord_long > sphereC])
     ! Sphere radius assumed to be equal to caseRadius
     ! ==============================================================================
     
     backX         = global%cnstrAftEnd ! End of nozzle bucket at t = 0
     backY         = global%cnstrNozY ! End of nozzle bucket at t = 0, (z = 0)
     tolerance     = global%cnstrTol2 ! Wiggle room
     sphereCenter  = backX - LDir*SQRT(caseRadius2 - backY**2)
     
     ! ==============================================================================
     ! Planar constraints
     ! ==============================================================================
     LMinPlane     = global%cnstrLMinPlane  
     LMaxPlane     = global%cnstrLMaxPlane  
     T1MinPlane    = global%cnstrT1MinPlane 
     T1MaxPlane    = global%cnstrT1MaxPlane
     T2MinPlane    = global%cnstrT2MinPlane
     T2MaxPlane    = global%cnstrT2MaxPlane

     ! ==============================================================================
     ! Keep track of how many nodes are constrained and how many faces burn out
     ! ==============================================================================
     boCount = 0 
     cvCount = 0 
     
     ! *****************************************************************************=
     ! Loop over patches
     ! *****************************************************************************=
     
     DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)
        
        ! -----------------------------------------------------------------------------   
        !     Loop over burning triangles, turn off burning if all vertices on case
        ! -----------------------------------------------------------------------------   
        
        IF ( pPatch%bcCoupled == BC_BURNING ) THEN 
           DO ifl = 1,pPatch%nBTris
              cntr = 0
              
              DO ivl = 1,3
                 ivg = pPatch%bTri2v(ivl,ifl)
                 
                 x = pXyz(coordLong,ivg)   
                 y = pXyz(coordTrans1,ivg) 
                 z = pXyz(coordTrans2,ivg) 
                 
                 IF( (LDir > 0 .AND. x > sphereCenter) .OR.  &
                     (LDir < 0 .AND. x < sphereCenter) ) THEN
                    x = x - sphereCenter
                    radius2 = x*x + y*y + z*z
                 ELSE IF ( (LDir > 0 .AND. x < ellipsC) .OR.  &
                           (LDir < 0 .AND. x > ellipsC) ) THEN
                    x = x - ellipsC
                    radius2 = (x*x*ellipsLong2+y*y*ellipsTrans2+z*z*ellipsTrans2) &
                         * caseRadius2
                 ELSE 
                    radius2 = y*y + z*z        
                 ENDIF ! x
                 IF ( radius2 >= (caseRadius2-global%cnstrTol1)) THEN
                    cntr = cntr + 1
                 ELSE
                    x = pXyz(coordLong,ivg)
                    IF(x <= LMinPlane)       THEN
                       cntr  = cntr + 1
                    ELSE IF(x >= LMaxPlane)  THEN
                       cntr  = cntr + 1
                    ELSE IF(y <= T1MinPlane) THEN 
                       cntr  = cntr + 1
                    ELSE IF(y >= T1MaxPlane) THEN
                       cntr  = cntr + 1
                    ELSE IF(z <= T2MinPlane) THEN
                       cntr  = cntr + 1
                    ELSE IF(z >= T2MaxPlane) THEN
                       cntr  = cntr + 1
                    ENDIF
                 END IF ! FloatEqual
                 
              END DO ! ivl
              
              IF ( cntr == 3 ) THEN
                 boCount = boCount + 1
                 pPatch%bflag(ifl)   = 0 
                 pPatch%mdotAlp(ifl) = 0.0_RFREAL
              ELSE
                 pPatch%bflag(ifl) = 1 ! make sure burning is on if not burned out
              END IF ! cntr
           END DO ! ifl
           
           ! -----------------------------------------------------------------------------    
           !     Loop over quadrilaterals, turn off burning if all vertices on case
           ! -----------------------------------------------------------------------------   
           
           DO ifl = 1,pPatch%nBQuads
              cntr = 0
              
              DO ivl = 1,4
                 ivg = pPatch%bQuad2v(ivl,ifl)
                 
                 x = pXyz(coordLong,ivg)
                 y = pXyz(coordTrans1,ivg)
                 z = pXyz(coordTrans2,ivg)
                 
                 IF( (LDir > 0 .AND. x > sphereCenter) .OR.  &
                     (LDir < 0 .AND. x < sphereCenter) ) THEN
                    x = x - sphereCenter
                    radius2 = x*x + y*y + z*z
                 ELSE IF ( (LDir > 0 .AND. x < ellipsC) .OR.  &
                           (LDir < 0 .AND. x < ellipsC) ) THEN
                    x = x - ellipsC
                    radius2 = (x*x*ellipsLong2+y*y*ellipsTrans2+z*z*ellipsTrans2) &
                         * caseRadius2
                 ELSE 
                    radius2 = y*y + z*z        
                 ENDIF ! x
                 IF ( radius2 >= (caseRadius2-global%cnstrTol1)) THEN
                    cntr = cntr + 1        
                 ELSE
                    x = pXyz(coordLong,ivg)
                    IF(x <= LMinPlane)       THEN
                       cntr  = cntr + 1
                    ELSE IF(x >= LMaxPlane)  THEN
                       cntr  = cntr + 1
                    ELSE IF(y <= T1MinPlane) THEN 
                       cntr  = cntr + 1
                    ELSE IF(y >= T1MaxPlane) THEN
                       cntr  = cntr + 1
                    ELSE IF(z <= T2MinPlane) THEN
                       cntr  = cntr + 1
                    ELSE IF(z >= T2MaxPlane) THEN
                       cntr  = cntr + 1
                    ENDIF
                 END IF ! FloatEqual
              END DO ! ivl
              
              IF ( cntr == 4 ) THEN ! All vertices on case, stop motion and burning
                 boCount = boCount + 1
                 pPatch%bflag(ifl+pPatch%nBTris)   = 0 
                 pPatch%mdotAlp(ifl+pPatch%nBTris) = 0.0_RFREAL
              ELSE
                 pPatch%bflag(ifl)   = 1
              END IF ! cntr
           END DO ! ifl 
        END IF ! pPatch%bcCoupled   
     END DO ! iPatch              
     
     ! *****************************************************************************
     ! End
     ! *****************************************************************************
     
     IF (global%checkLevel == CHECK_HIGH) THEN
        boTot = 0
        CALL MPI_Reduce(boCount,boTot,1,MPI_INTEGER,MPI_SUM, &
             MASTERPROC,global%mpiComm,global%error)
        IF ( global%error /= ERR_NONE ) THEN
           CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
        END IF ! global%errorFlag
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel >= VERBOSE_HIGH) THEN
           WRITE(STDOUT,'(A,3X,A,1X,I9)') SOLVER_NAME, &
                'Total number of burned out faces:',boTot
        END IF ! global%myProcid
     ENDIF ! global%checkLevel

  ENDIF ! (global%cnstrCaseRad > 0)

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
     WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing BFLAG done.'
  END IF ! global%myProcid
  
  CALL DeregisterFunction(global)
  
END SUBROUTINE RFLU_GENX_InitBFLAG








! ******************************************************************************
!
! Purpose: Move grid through Roccom.
!
! Description: 
!
! Input: 
!   regions     Region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_GENX_MoveGrid(regions)

  IMPLICIT NONE
  
! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: winName
  INTEGER :: errorFlag,icg,iPatch,iReg,ivg,ivl
  INTEGER :: handleDisp,handleOption,handlePMesh,handleSmooth
  REAL(RFREAL) :: dispTot,dispMax
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pDisp,pRhs,pXyz,pXyzOld
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid,pGridOld
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! *****************************************************************************
! Start
! *****************************************************************************

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_GENX_MoveGrid',&
  'RFLU_ModRocstarTools.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Moving grid...'
  END IF ! global%myProcid
   
! *****************************************************************************
! Get function and dataitem handles
! *****************************************************************************
   
  winName = TRIM(global%volWinName)   
       
  handleSmooth = COM_get_function_handle('Rocflu-MOP.smooth')    
  handlePMesh  = COM_get_dataitem_handle(TRIM(winName)//'.pmesh')
  handleDisp   = COM_get_dataitem_handle(TRIM(winName)//'.disp')  
    
  handleOption = COM_get_function_handle('Rocflu-MOP.set_value')
  CALL COM_call_function(handleOption,2,'inverted',1) 
    
! *****************************************************************************
! Copy grid data to old grid
!
! This will only be done for non-Newton-Krylov solvers since that implicit
! solver will use displacement ramping and hense will have to control when
! the grid is stored.
! *****************************************************************************

  IF ( global%solverType .NE. SOLV_IMPLICIT_NK ) THEN
    DO iReg = 1,global%nRegionsLocal
    pGrid    => regions(iReg)%grid
    pGridOld => regions(iReg)%gridOld

    DO icg = 1,pGrid%nCellsTot ! Explicit copy to avoid ASCI White problem
      pGridOld%vol(icg) = pGrid%vol(icg)
    END DO ! icg

    DO ivg = 1,pGrid%nVertTot ! Explicit copy to avoid ASCI White problem
      pGridOld%xyz(XCOORD,ivg) = pGrid%xyz(XCOORD,ivg)
      pGridOld%xyz(YCOORD,ivg) = pGrid%xyz(YCOORD,ivg)
      pGridOld%xyz(ZCOORD,ivg) = pGrid%xyz(ZCOORD,ivg)
    END DO ! iv
  END DO ! iReg
  END IF ! global%solverType
      
! ******************************************************************************
! Initialize displacements to zero and impose patch displacements on vertex
! coordinates. NOTE this is different from other mesh motion routines because 
! Rocmop ignores displacement array on entry and requires that coordinates 
! already have boundary displacement added. NOTE pass in non-zero displacements
! to Rocmop for thresholding.
! ******************************************************************************    
  dispMax = 0.0_RFREAL
  DO iReg = 1,global%nRegionsLocal 
    pGrid    => regions(iReg)%grid
    pGridOld => regions(iReg)%gridOld 
   
    pDisp   => pGrid%disp  
    pXyz    => pGrid%xyz
    pXyzOld => pGridOld%xyz   

    DO ivg = 1,pGrid%nVertTot
      pDisp(XCOORD,ivg) = 0.0_RFREAL 
      pDisp(YCOORD,ivg) = 0.0_RFREAL       
      pDisp(ZCOORD,ivg) = 0.0_RFREAL          
    END DO ! ivg
    
    DO iPatch = 1,pGrid%nPatches
      pPatch => regions(iReg)%patches(iPatch)
 
      DO ivl = 1,pPatch%nBVert
        ivg = pPatch%bv(ivl)            

        pXyz(XCOORD,ivg) = pXyzOld(XCOORD,ivg) + pPatch%dXyz(XCOORD,ivl)
        pXyz(YCOORD,ivg) = pXyzOld(YCOORD,ivg) + pPatch%dXyz(YCOORD,ivl)
        pXyz(ZCOORD,ivg) = pXyzOld(ZCOORD,ivg) + pPatch%dXyz(ZCOORD,ivl)


!       Pass in displacements for displacement thresholding
        pDisp(XCOORD,ivg) = pPatch%dXyz(XCOORD,ivl)
        pDisp(YCOORD,ivg) = pPatch%dXyz(YCOORD,ivl)
        pDisp(ZCOORD,ivg) = pPatch%dXyz(ZCOORD,ivl)
!
      END DO ! ivl
    END DO ! iPatch         
    IF (global%checkLevel == CHECK_HIGH) THEN
       DO iPatch = 1,pGrid%nPatches
          pPatch => regions(iReg)%patches(iPatch)

          DO ivl = 1,pPatch%nBVert
             ivg = pPatch%bv(ivl)

             dispTot = pDisp(XCOORD,ivg)*pDisp(XCOORD,ivg)  +            &
                  pDisp(YCOORD,ivg)*pDisp(YCOORD,ivg) +                  &
                  pDisp(ZCOORD,ivg)*pDisp(ZCOORD,ivg)
             IF(dispTot > dispMax) THEN
                dispMax = dispTot
             ENDIF

          END DO ! ivl
       END DO ! iPatch
    ENDIF

  END DO ! iReg        
     
  IF (global%checkLevel == CHECK_HIGH) THEN
     dispTot = dispMax
     CALL MPI_Reduce(dispTot,dispMax,1,MPI_RFREAL,MPI_MAX, &
          MASTERPROC,global%mpiComm,global%error)
     IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
     END IF ! global%errorFlag
     IF ( global%myProcid == MASTERPROC .AND. &
          global%verbLevel >= VERBOSE_HIGH) THEN
        dispMax = SQRT(dispMax)
        WRITE(STDOUT,'(A,3X,A,1X,E13.6)') SOLVER_NAME, &
             'Maximum nodal displacement from Rocstar:',dispMax
     END IF ! global%myProcid
  ENDIF

! ******************************************************************************
! Smooth displacements through Roccom
! ******************************************************************************     

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Smoothing grid with Rocmop...'
  END IF ! global%myProcid
   
  CALL COM_call_function(handleSmooth,2,handlePMesh,handleDisp)     
     
  CALL MPI_Barrier(global%mpiComm,global%error)
  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Done smoothing grid with Rocmop.'
  END IF ! global%myProcid

  dispTot = 0.0
  dispMax = 0.0
  IF (global%checkLevel == CHECK_HIGH) THEN
     DO iReg = 1,global%nRegionsLocal
        pGrid    => regions(iReg)%grid
        pDisp   => pGrid%disp
        DO iPatch = 1,pGrid%nPatches
           pPatch => regions(iReg)%patches(iPatch)
           DO ivl = 1,pPatch%nBVert
              ivg = pPatch%bv(ivl)
              pDisp(XCOORD,ivg) = 0.0_RFREAL
              pDisp(YCOORD,ivg) = 0.0_RFREAL
              pDisp(ZCOORD,ivg) = 0.0_RFREAL
              dispTot = pDisp(XCOORD,ivg)*pDisp(XCOORD,ivg) +            &
                   pDisp(YCOORD,ivg)*pDisp(YCOORD,ivg)      +            &
                   pDisp(ZCOORD,ivg)*pDisp(ZCOORD,ivg)
              IF(dispTot > dispMax) THEN
                 dispMax = dispTot
              ENDIF
           END DO ! ivl
        END DO ! iPatch
     END DO ! iReg
     dispTot = dispMax
     CALL MPI_Reduce(dispTot,dispMax,1,MPI_RFREAL,MPI_MAX, &
          MASTERPROC,global%mpiComm,global%error)
     IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
     END IF ! global%errorFlag
     IF (dispMax /= 0.0) THEN
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel >= VERBOSE_NONE) THEN
           dispMax = SQRT(dispMax)
           WRITE(STDOUT,'(A,3X,A,1X,E13.6)') SOLVER_NAME, &
                '***WARNING - NONZERO SURFACE DISP FROM ROCMOP:',&
                dispMax
        END IF ! global%myProcid
     ENDIF
  ENDIF
! ******************************************************************************
! Update coordinates
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal 
    pGrid => regions(iReg)%grid
    pDisp => pGrid%disp
    pXyz  => pGrid%xyz

    DO ivg = 1,pGrid%nVertTot
       pXyz(XCOORD,ivg) = pXyz(XCOORD,ivg) + pDisp(XCOORD,ivg)
       pXyz(YCOORD,ivg) = pXyz(YCOORD,ivg) + pDisp(YCOORD,ivg)
       pXyz(ZCOORD,ivg) = pXyz(ZCOORD,ivg) + pDisp(ZCOORD,ivg) 
    END DO ! ivg
  END DO ! iReg
  
  IF (global%checkLevel ==  CHECK_HIGH) THEN
    dispTot = 0.0
    dispMax = 0.0

    DO iReg = 1,global%nRegionsLocal
       pGrid => regions(iReg)%grid
       pDisp => pGrid%disp

       DO ivg = 1,pGrid%nVertTot
          dispTot = pDisp(XCOORD,ivg)*pDisp(XCOORD,ivg)+ &
               pDisp(YCOORD,ivg)*pDisp(YCOORD,ivg)+ &
               pDisp(ZCOORD,ivg)*pDisp(ZCOORD,ivg)
          IF(dispTot > dispMax) THEN
             dispMax = dispTot
          ENDIF
       END DO ! ivg

    END DO ! iReg

    dispTot = dispMax
    CALL MPI_Reduce(dispTot,dispMax,1,MPI_RFREAL,MPI_MAX, &
         MASTERPROC,global%mpiComm,global%error)
    IF ( global%error /= ERR_NONE ) THEN
       CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
    END IF ! global%errorFlag
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH) THEN
       dispMax = SQRT(dispMax)
       WRITE(STDOUT,'(A,3X,A,1X,E13.6)') SOLVER_NAME, &
            'Maximum nodal displacement from Rocmop:',dispMax
    END IF ! global%myProcid

 END IF ! global%checkLevel

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Moving grid done.'
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GENX_MoveGrid








END MODULE RFLU_ModRocstarTools

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModGENXTools.F90,v $
! Revision 1.29  2008/12/06 08:44:21  mtcampbe
! Updated license.
!
! Revision 1.28  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.27  2008/09/26 15:53:02  rfiedler
! Do not apply planar constraints to nodes with vanishing displacements (for BSM).
!
! Revision 1.26  2008/09/17 20:14:41  rfiedler
! Reduce tolerance for small displacements.
!
! Revision 1.25  2008/09/11 22:42:04  rfiedler
! Do not constrain nodes that have 0 displacements -- for BSM.
!
! Revision 1.24  2008/06/20 19:31:32  rfiedler
! Uses sign of HeadEnd - AftEnd to determine direction of rocket.
!
! Revision 1.23  2008/01/31 22:47:28  mtcampbe
! added default coordinate shift of 0.0 which should be used when planar
! constraints are considered
!
! Revision 1.22  2008/01/14 22:08:30  mtcampbe
! added planar constr
!
! Revision 1.21  2007/04/20 16:07:49  mtcampbe
! Updating for burnout support function RFLU_GENX_InitBFLAG
!
! Revision 1.20  2007/04/14 14:08:50  mtcampbe
! Cleaned up, bug fixes to constraints for mixed meshes, updated for ROCKET inp section
!
! Revision 1.19  2007/02/06 18:12:10  mtcampbe
! Attempt to fix surface node mismatch problem by zeroing tiny surface displacements
! received from Rocmop.  Currently, this will happen only if checking is high.
!
! Revision 1.18  2006/11/30 16:50:37  mtcampbe
! Fixed 2 constraint code bugs: XyzOld changed to Xyz, and now looping over
! all patches to limit motion instead of just burning patches.
!
! Revision 1.17  2006/11/13 19:26:15  mtcampbe
! New constraint codes and displacement checking
!
! Revision 1.16  2006/09/10 16:57:51  mtcampbe
! Added the constraint displacement routine.
!
! Revision 1.15  2006/02/08 21:13:02  hdewey2
! Deactivating grid copy for implicit solver bcos of disp ramping
!
! Revision 1.14  2005/09/22 16:16:17  mtcampbe
! Fixed checking bug
!
! Revision 1.13  2005/09/22 16:02:43  mtcampbe
! Re-enabled displacement passing, added displacement checking for grid motion
!
! Revision 1.12  2005/09/16 00:44:20  mtcampbe
! Emergency bugfix - disabled displacement passing to Rocmop
!
! Revision 1.11  2005/09/14 21:32:22  haselbac
! Bug fix: Disps were applied several times for bverts shared betw patches
!
! Revision 1.10  2005/09/04 16:32:25  mtcampbe
! Changed to pass in mesh displacements to Rocmop
!
! Revision 1.9  2005/09/01 16:04:42  mtcampbe
! Removed control calls to Rocmop which have been replaced by RocmopControl.txt
!
! Revision 1.8  2005/06/20 21:24:50  haselbac
! Bug fix: Removed IF on pPatch%movePatchDir
!
! Revision 1.7  2005/06/09 20:20:49  haselbac
! Now use movePatchDir instead of movePatch and smoothGrid
!
! Revision 1.6  2005/06/06 22:44:38  rfiedler
! Reduce default Rocmop verbosity to 0.
!
! Revision 1.5  2005/04/15 15:06:57  haselbac
! Cosmetics only
!
! Revision 1.4  2005/04/02 22:14:29  haselbac
! Added setting of method option for mesh smoother
!
! Revision 1.3  2004/12/08 14:43:01  haselbac
! Bug fixes: Pass flag for inverted orientation and patch smoothing
!
! Revision 1.2  2004/10/22 14:01:57  haselbac
! Adapted code to work with Rocmop
!
! Revision 1.1  2004/10/19 19:27:32  haselbac
! Initial revision
!
! ******************************************************************************









