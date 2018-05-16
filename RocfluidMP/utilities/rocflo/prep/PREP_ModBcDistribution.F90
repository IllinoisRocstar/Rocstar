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
! Purpose: Suite of bc distribution routines.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: PREP_ModBcDistribution.F90,v 1.9 2008/12/06 08:44:50 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE PREP_ModBcDistribution

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region  
  USE ModBndPatch, ONLY  : t_patch
  USE ModGrid, ONLY: t_grid
  USE ModMPI
  USE PREP_ModParameters
  
  IMPLICIT NONE

  PRIVATE 
  PUBLIC :: BcDistributionFiles
 
! private : BcNoslipDistrib, &
!           BcInflowDistrib, &
!           BcInjectDistrib, &
!           BcOutflowDistrib, &
!           BcFarfDistrib, &
!           BcSlipWallDistrib, &
!           ReadPatchSection, &
!           WriteBcToFile
 
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: PREP_ModBcDistribution.F90,v $ $Revision: 1.9 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS

!******************************************************************************
!
! Purpose: make bc distribution files
!
! Description: none.
!
! Input: boundary condition file.
!
! Output: regions = BC files containing data distribution.
!
! Notes: 
!
!******************************************************************************

SUBROUTINE BcDistributionFiles( regions )

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  CHARACTER(2*CHRLEN+9) :: fname
  CHARACTER(256)        :: line

  INTEGER :: distrib, errorFlag

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'BcDistributionFiles',&
  'PREP_ModBcDistribution.F90' )

! open file

  WRITE(fname,'(A)') TRIM(global%inDir)//TRIM(global%casename)//'.bc'
  OPEN(IF_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! allocate and initialize bc-plane edges

  ALLOCATE( global%infloPlanEdges(XCOORD:ZCOORD,XCOORD:ZCOORD,2), &
            stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( global%xyzMinmax(XCOORD:ZCOORD,2),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  global%infloPlanEdges = 0._RFREAL
  global%xyzMinmax(:,1) = REAL_LARGE
  global%xyzMinmax(:,2) = REAL_SMALL

! read file looking for keywords

  CALL BcCaseLoop( 1 )
  IF (ABS( SUM( global%infloPlanEdges ) ) > 1.E-10_RFREAL) THEN
    REWIND(IF_INPUT)
    CALL BcCaseLoop( 2)
  ENDIF

! close file ------------------------------------------------------------------

  CLOSE(IF_INPUT,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )

! finalization & error handling -----------------------------------------------

  CALL DeregisterFunction( global )

! -----------------------------------------------------------------------------
  CONTAINS

  SUBROUTINE BcCaseLoop( n )

    INTEGER :: n

    DO
      READ(IF_INPUT,'(A256)',err=10,end=86) line
      SELECT CASE(TRIM(line))

      CASE ('# BC_SLIPWALL')
        distrib=0
        CALL BcSlipWallDistrib( regions,n,BC_SLIPWALL,distrib )
        IF (distrib/=0) CALL ErrorStop( global,ERR_UNKNOWN_OPTION,__LINE__, &
                        ' generator for slip wall distribution not ready yet' )

      CASE ('# BC_NOSLIPWALL')
        distrib=0
        CALL BcNoslipDistrib( regions,n,BC_NOSLIPWALL,distrib )
        IF (distrib/=0) CALL ErrorStop( global,ERR_UNKNOWN_OPTION,__LINE__, &
                      ' generator for noslip wall distribution not ready yet' )

! TEMPORARY - Keep this for backward compatibility
      CASE ('# BC_INFLOW')
        distrib=0
        CALL BcInflowTotAngDistrib( regions,n,BC_INFLOW,distrib )
        IF (distrib/=0) CALL ErrorStop( global,ERR_UNKNOWN_OPTION,__LINE__, &
                    ' generator for TotAng-inflow distribution not ready yet' )
! END TEMPORARY          
    
      CASE ('# BC_INFLOW_TOTANG')
        distrib=0
        CALL BcInflowTotAngDistrib( regions,n,BC_INFLOW_TOTANG,distrib )
        IF (distrib/=0) CALL ErrorStop( global,ERR_UNKNOWN_OPTION,__LINE__, &
                    ' generator for TotAng-inflow distribution not ready yet' )

      CASE ('# BC_INFLOW_VELTEMP')
        distrib=0
        CALL BcInflowVelDistrib( regions,n,BC_INFLOW_VELTEMP,distrib )

      CASE ('# BC_INFLOW_VELPRESS')
        distrib=0
        CALL BcInflowVelDistrib( regions,n,BC_INFLOW_VELPRESS,distrib )

      CASE ('# BC_OUTFLOW')
        distrib=0
        CALL BcOutflowDistrib( regions,n,BC_OUTFLOW,distrib )
        IF (distrib/=0) CALL ErrorStop( global,ERR_UNKNOWN_OPTION,__LINE__, &
                        ' generator for outflow distribution not ready yet' )
    
      CASE ('# BC_FARFIELD')
        distrib=0
        CALL BcFarfDistrib( regions,n,BC_FARFIELD,distrib )
        IF (distrib/=0) CALL ErrorStop( global,ERR_UNKNOWN_OPTION,__LINE__, &
                        ' generator for farfield distribution not ready yet' )
    
      CASE ('# BC_INJECTION')
        distrib=0
        CALL BcInjectDistrib( regions,n,BC_INJECTION,distrib )
        IF (distrib/=0) CALL ErrorStop( global,ERR_UNKNOWN_OPTION,__LINE__, &
                        ' generator for injection distribution not ready yet' )

      END SELECT
    ENDDO

86   CONTINUE

  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) )

999  CONTINUE

  END SUBROUTINE BcCaseLoop

END SUBROUTINE BcDistributionFiles

!******************************************************************************
!
! Purpose: read in user input related to slip-wall boundary condition
!          and create bc distribution file if applicable
!
! Description: none.
!
! Input: boundary condition file.
!
! Output: regions = BC data.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE BcSlipWallDistrib( regions,ithRead,bcTitle,distrib )

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  INTEGER :: ithRead, bcTitle, distrib

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  CHARACTER(10)  :: keys(2)
  CHARACTER(256) :: fname

  INTEGER :: brbeg, brend, prbeg, prend, profType, errorFlag

  LOGICAL :: defined(2)

  REAL(RFREAL) :: vals(2)

  TYPE(t_patch), POINTER :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'BcSlipWallDistrib',&
  'PREP_ModBcDistribution.F90' )

! specify keywords and search for them ----------------------------------------

  keys(1) = 'EXTRAPOL'
  keys(2) = 'MAXCHANGE'

  distrib  = 0
  profType = 0
  CALL ReadPatchSection( global,IF_INPUT,2,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,profType,fname,defined )
  
! check if all values defined -------------------------------------------------

  IF (.NOT. defined(1) .OR. &
      .NOT. defined(2)) CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )

! copy values/distribution to variables ---------------------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_SLIPWALL .AND. &
           patch%bcType<=BC_SLIPWALL+BC_RANGE) .AND. &    ! my boundary type
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        IF (patch%mixt%bcSet) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,__LINE__,'Slip-wall boundary.' )

        patch%mixt%nData     = 0
        patch%mixt%nSwitches = 1
        patch%mixt%distrib   = BCDAT_CONSTANT
        IF (ithRead==2) patch%mixt%bcSet = .true.

! ----- get value of switch

        ALLOCATE( patch%mixt%switches(patch%mixt%nSwitches), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

          patch%mixt%switches(BCSWI_SLIPW_EXTRAP) = EXTRAPOL_CONST
        IF (vals(1) > 0.1) &
          patch%mixt%switches(BCSWI_SLIPW_EXTRAP) = EXTRAPOL_LINEAR

        patch%mixt%maxChange = vals(2)

      ENDIF    ! bcType, active region on my processor

    ENDDO      ! iPatch
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE BcSlipWallDistrib

!******************************************************************************
!
! Purpose: read in user input related to noslip wall boundary condition.
!          and create bc distribution file if applicable
!
! Description: none.
!
! Input: boundary condition file.
!
! Output: regions = BC data.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE BcNoslipDistrib( regions,ithRead,bcTitle,distrib )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  INTEGER :: ithRead, bcTitle, distrib

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  CHARACTER(10)  :: keys(2)
  CHARACTER(256) :: fname

  INTEGER :: brbeg, brend, prbeg, prend, profType, switch
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag

  LOGICAL :: defined(2)

  REAL(RFREAL) :: vals(2)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'BcNoslipDistrib',&
  'PREP_ModBcDistribution.F90' )

! specify keywords and search for them ----------------------------------------

  keys(1) = 'ADIABAT'
  keys(2) = 'TWALL'

  distrib  = 0
  profType = 0
  CALL ReadPatchSection( global,IF_INPUT,2,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,profType,fname,defined )

! get switches & check if all necessary values defined ------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_NOSLIPWALL .AND. &
           patch%bcType<=BC_NOSLIPWALL+BC_RANGE) .AND. &   ! my boundary type,
          regions(iReg)%procid==global%myProcid  .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN               ! on my processor

        IF (patch%mixt%bcSet) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,__LINE__,'Noslip boundary.' )

        patch%mixt%nData     = 0
        patch%mixt%nSwitches = 1
        IF (patch%bcCoupled == BC_EXTERNAL) THEN   ! data from outside
          patch%mixt%distrib = BCDAT_DISTRIB    ! => always distribution
        ELSE
          patch%mixt%distrib = distrib
        ENDIF

        ALLOCATE( patch%mixt%switches(patch%mixt%nSwitches), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! ----- check if switch defined
        IF (defined(1)) THEN
          patch%mixt%switches(BCSWI_NOSLIP_ADIABAT) = BCOPT_ADIABAT
          IF (vals(1) < 0.1) &
            patch%mixt%switches(BCSWI_NOSLIP_ADIABAT) = BCOPT_NON_ADIABAT
        ELSE
          CALL ErrorStop( global,ERR_NO_BCSWITCH,__LINE__,'(adiabatic wall yes/no).' )
        ENDIF

! ----- check if Twall specified (value or file with distribution)
        IF (patch%mixt%switches(BCSWI_NOSLIP_ADIABAT) == &
            BCOPT_NON_ADIABAT) THEN
          IF (patch%mixt%distrib==BCDAT_CONSTANT .AND. &
              (.NOT. defined(2))) CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )
        ENDIF

! ----- set flag to BC specified
        IF (ithRead==2) patch%mixt%bcSet = .true.

      ENDIF   ! my BC & processor, active
    ENDDO
  ENDDO

! copy values/distribution to variables ---------------------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_NOSLIPWALL .AND. &
           patch%bcType<=BC_NOSLIPWALL+BC_RANGE) .AND. &   ! my boundary type,
          regions(iReg)%procid==global%myProcid  .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN               ! on my processor
        switch = patch%mixt%switches(BCSWI_NOSLIP_ADIABAT)
      ELSE
        switch = BCOPT_ADIABAT
      ENDIF

      IF ((patch%bcType>=BC_NOSLIPWALL .AND. &
           patch%bcType<=BC_NOSLIPWALL+BC_RANGE) .AND. &   ! my boundary type,
          switch==BCOPT_NON_ADIABAT              .AND. &   ! Twall required,
          regions(iReg)%procid==global%myProcid  .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN               ! on my processor

! ----- allocate memory for the values

        patch%mixt%nData = 1

        IF (patch%mixt%distrib == BCDAT_DISTRIB) THEN
          n1    = ABS(patch%l1end-patch%l1beg)
          n2    = ABS(patch%l2end-patch%l2beg)
          iOff  = n1 + 1
          ijBeg = IndIJ( 0, 0,iOff)
          ijEnd = IndIJ(n1,n2,iOff)
        ELSE
          ijBeg = 0
          ijEnd = 1
        ENDIF
        ALLOCATE( patch%mixt%vals(patch%mixt%nData,ijBeg:ijEnd), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! ----- distribution from file

        IF (patch%mixt%distrib==BCDAT_DISTRIB .AND. &
            patch%bcCoupled      /=BC_EXTERNAL  ) THEN
!          CALL WriteBcToFile( global,fname,patch )

! ----- distribution from external source / constant value

        ELSE
          patch%mixt%vals(BCDAT_NOSLIP_TWALL,:) = vals(2)
        ENDIF  ! distribution?

      ENDIF    ! bcType, Twall req., active region on my processor

    ENDDO      ! iPatch
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE BcNoslipDistrib

!******************************************************************************
!
! Purpose: read in user input related to inflow boundary condition.
!          and create bc distribution file if applicable
!
! Description: present inflow bc is based on total pressure, total temperature
!              and flow angle.
!
! Input: boundary condition file.
!
! Output: regions = BC data.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE BcInflowTotAngDistrib( regions,ithRead,bcTitle,distrib )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  INTEGER :: ithRead, bcTitle, distrib

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  CHARACTER(10)  :: keys(7)
  CHARACTER(256) :: fname

  INTEGER :: brbeg, brend, prbeg, prend, profType, switch
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag

  LOGICAL :: defined(7)

  REAL(RFREAL) :: vals(7)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'BcInflowTotAngDistrib',&
  'PREP_ModBcDistribution.F90' )

! specify keywords and search for them ----------------------------------------

  keys(1) = 'TYPE'
  keys(2) = 'FIXED'
  keys(3) = 'PTOT'
  keys(4) = 'TTOT'
  keys(5) = 'BETAH'
  keys(6) = 'BETAV'
  keys(7) = 'MACH'

  distrib  = 0
  profType = 0
  CALL ReadPatchSection( global,IF_INPUT,7,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,profType,fname,defined )

! get switches & check if all necessary values defined ------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_INFLOW .AND. &
           patch%bcType<=BC_INFLOW+BC_RANGE)    .AND. &   ! my boundary type,
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        patch%bcType = BC_INFLOW_TOTANG

        IF (patch%mixt%bcSet) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,__LINE__,'Inflow boundary.' )

        patch%mixt%nSwitches = 2
        IF (patch%bcCoupled == BC_EXTERNAL) THEN   ! data from outside
          patch%mixt%distrib = BCDAT_DISTRIB    ! => always distribution
        ELSE
          patch%mixt%distrib = distrib
        ENDIF

        ALLOCATE( patch%mixt%switches(patch%mixt%nSwitches), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! ----- check if switch defined
        IF (defined(1)) THEN
          patch%mixt%switches(BCSWI_INFLOW_TYPE)   = BCOPT_SUBSONIC
          IF (vals(1) < 0.1) &
            patch%mixt%switches(BCSWI_INFLOW_TYPE) = BCOPT_SUPERSONIC
          IF (vals(1) > 1.9) &
            patch%mixt%switches(BCSWI_INFLOW_TYPE) = BCOPT_MIXED
        ELSE
          CALL ErrorStop( global,ERR_NO_BCSWITCH,__LINE__,'(inflow type).' )
        ENDIF

        IF (defined(2)) THEN
          IF (vals(2) < 0.1) &
            patch%mixt%switches(BCSWI_INFLOW_FIXED) = BCOPT_FIXED_NO
          IF (vals(2) > 0.9) &
            patch%mixt%switches(BCSWI_INFLOW_FIXED) = BCOPT_FIXED_YES
        ELSE
          patch%mixt%switches(BCSWI_INFLOW_FIXED) = BCOPT_FIXED_NO
        ENDIF

! ----- check if appropriate values specified
        IF (patch%mixt%switches(BCSWI_INFLOW_TYPE) == BCOPT_SUBSONIC) THEN
          IF (patch%mixt%distrib==BCDAT_CONSTANT .AND. &
              (.NOT. defined(3) .OR. &
               .NOT. defined(4) .OR. &
               .NOT. defined(5) .OR. &
               .NOT. defined(6))) CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )
        ENDIF
        IF (patch%mixt%switches(BCSWI_INFLOW_TYPE) == BCOPT_SUPERSONIC .OR. &
            patch%mixt%switches(BCSWI_INFLOW_TYPE) == BCOPT_MIXED) THEN
          IF (patch%mixt%distrib==BCDAT_CONSTANT .AND. &
              (.NOT. defined(3) .OR. &
               .NOT. defined(4) .OR. &
               .NOT. defined(5) .OR. &
               .NOT. defined(6) .OR. &
               .NOT. defined(7))) CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )
        ENDIF

! ----- set flag to BC specified
        IF (ithRead==2) patch%mixt%bcSet = .true.

      ENDIF   ! my BC & processor, active
    ENDDO
  ENDDO

! copy values/distribution to variables ---------------------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_INFLOW .AND. &
           patch%bcType<=BC_INFLOW+BC_RANGE)    .AND. &   ! my boundary type,
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        switch = patch%mixt%switches(BCSWI_INFLOW_TYPE)
        IF (switch == BCOPT_SUBSONIC) THEN
          patch%mixt%nData = 4
        ELSE
          patch%mixt%nData = 5
        ENDIF

! ----- allocate memory for the values

        IF (patch%mixt%distrib == BCDAT_DISTRIB) THEN
          n1    = ABS(patch%l1end-patch%l1beg)
          n2    = ABS(patch%l2end-patch%l2beg)
          iOff  = n1 + 1
          ijBeg = IndIJ( 0, 0,iOff)
          ijEnd = IndIJ(n1,n2,iOff)
        ELSE
          ijBeg = 0
          ijEnd = 1
        ENDIF
        ALLOCATE( patch%mixt%vals(patch%mixt%nData,ijBeg:ijEnd), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! ----- distribution from file

        IF (patch%mixt%distrib==BCDAT_DISTRIB .AND. &
            patch%bcCoupled      /=BC_EXTERNAL  ) THEN
!          CALL WriteBcToFile( global,fname,patch )

          IF (switch == BCOPT_SUBSONIC) THEN
            patch%mixt%vals(BCDAT_INFLOW_BETAH,:) = &
              patch%mixt%vals(BCDAT_INFLOW_BETAH,:)*global%rad
            patch%mixt%vals(BCDAT_INFLOW_BETAV,:) = &
              patch%mixt%vals(BCDAT_INFLOW_BETAV,:)*global%rad
          ENDIF

! ----- distribution from external source / constant value

        ELSE
          patch%mixt%vals(BCDAT_INFLOW_PTOT ,:) = vals(3)
          patch%mixt%vals(BCDAT_INFLOW_TTOT ,:) = vals(4)
          patch%mixt%vals(BCDAT_INFLOW_BETAH,:) = vals(5)*global%rad
          patch%mixt%vals(BCDAT_INFLOW_BETAV,:) = vals(6)*global%rad
          IF (switch /= BCOPT_SUBSONIC) THEN
            patch%mixt%vals(BCDAT_INFLOW_MACH,:) = vals(7)
          ENDIF
        ENDIF  ! distribution?

      ENDIF    ! bcType, active region on my processor

    ENDDO      ! iPatch
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE BcInflowTotAngDistrib

!******************************************************************************
!
! Purpose: read in user input related to inflow boundary condition.
!          and create bc distribution file if applicable
!
! Description: present inflow bc is based on prescribed velocities and 
!              either temperature or pressure.
!
! Input: boundary condition file.
!
! Output: regions = BC data.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE BcInflowVelDistrib( regions,ithRead,bcTitle,distrib )

  USE PREP_ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  INTEGER :: ithRead, bcTitle, distrib

! ... loop variables
  INTEGER :: iReg, iPatch, i, j, k, l

! ... local variables
  INTEGER, PARAMETER :: NVALS_MAX = 8

  CHARACTER(10)  :: keys(NVALS_MAX)
  CHARACTER(256) :: fname

  INTEGER :: nvals, brbeg, brend, prbeg, prend, profType, switch
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iLev, iNOff, ijNOff
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend, ijkN

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL) :: vals(NVALS_MAX)
  REAL(RFREAL), POINTER :: xyz(:,:)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'BcInflowVelDistrib',&
  'PREP_ModBcDistribution.F90' )

! specify keywords and search for them ----------------------------------------

  keys(1) = 'TYPE'
  keys(2) = 'VELX'
  keys(3) = 'VELY'
  keys(4) = 'VELZ'
  IF (bcTitle==BC_INFLOW_VELTEMP) THEN
    keys(5) = 'TEMP'
    keys(6) = 'PRESS'
  ELSEIF (bcTitle==BC_INFLOW_VELPRESS) THEN
    keys(5) = 'PRESS'
    keys(6) = 'TEMP'
  ELSE
    CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__ )
  ENDIF
  keys(7) = 'AXIALPOWER'
  keys(8) = 'NORMALFACT'

  nvals = NVALS_MAX

  distrib  = 0
  profType = 0
  CALL ReadPatchSection( global,IF_INPUT,nvals,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,profType,fname,defined )

! get switches & check if all necessary values defined ------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_INFLOW .AND. &
           patch%bcType<=BC_INFLOW+BC_RANGE)    .AND. &   ! my boundary type,
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        patch%bcType = bcTitle

        IF (patch%mixt%bcSet) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,__LINE__,'Inflow boundary.' )

        patch%mixt%nSwitches = 2
        IF (patch%bcCoupled == BC_EXTERNAL) THEN   ! data from outside
          patch%mixt%distrib = BCDAT_DISTRIB    ! => always distribution
        ELSE
          patch%mixt%distrib = distrib
        ENDIF

        ALLOCATE( patch%mixt%switches(patch%mixt%nSwitches), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! ----- check if switch defined
        IF (defined(1)) THEN
          patch%mixt%switches(BCSWI_INFLOW_TYPE)   = BCOPT_SUBSONIC
          IF (vals(1) < 0.1) &
            patch%mixt%switches(BCSWI_INFLOW_TYPE) = BCOPT_SUPERSONIC
          IF (vals(1) > 1.9) &
            patch%mixt%switches(BCSWI_INFLOW_TYPE) = BCOPT_MIXED
        ELSE
          CALL ErrorStop( global,ERR_NO_BCSWITCH,__LINE__,'(inflow type).' )
        ENDIF

! ----- check if appropriate values specified
        IF (patch%mixt%switches(BCSWI_INFLOW_TYPE) == BCOPT_SUBSONIC) THEN
          IF (patch%mixt%distrib==BCDAT_CONSTANT .AND. &
              (.NOT. defined(2) .OR. &
               .NOT. defined(3) .OR. &
               .NOT. defined(4) .OR. &
               .NOT. defined(5))) CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )
        ENDIF
        IF (patch%mixt%switches(BCSWI_INFLOW_TYPE) == BCOPT_SUPERSONIC .OR. &
            patch%mixt%switches(BCSWI_INFLOW_TYPE) == BCOPT_MIXED) THEN
          IF (patch%mixt%distrib==BCDAT_CONSTANT .AND. &
              (.NOT. defined(2) .OR. &
               .NOT. defined(3) .OR. &
               .NOT. defined(4) .OR. &
               .NOT. defined(5) .OR. &
               .NOT. defined(6))) CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )
        ENDIF

! ----- set flag to BC specified
        IF (ithRead==2) patch%mixt%bcSet = .true.

      ENDIF   ! my BC & processor, active
    ENDDO
  ENDDO

! obtain geometrical parameters -----------------------------------------------

! allocate memory for the geometrical edges

  DO iReg=brbeg,brend

    iLev =  regions(iReg)%currLevel
    xyz  => regions(iReg)%levels(iLev)%grid%xyz

    CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )

    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_INFLOW .AND. &
           patch%bcType<=BC_INFLOW+BC_RANGE)    .AND. &   ! my boundary type,
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        IF (patch%mixt%distrib == BCDAT_DISTRIB .AND. &
            patch%bcCoupled       /= BC_EXTERNAL  ) THEN

          IF (.NOT. defined(6) .OR. profType==0) &
              CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__, &
                             ' PROFILE or TEMP or PRESS should be defined' )
          IF (.NOT. defined(7)) vals(7) = 2._RFREAL
          IF (.NOT. defined(8)) vals(8) = 1._RFREAL

! ------- search for global geometrical edges on the bc plane

          CALL RFLO_GetPatchIndicesNodes( regions(iReg),patch,iLev, &
                                          ibeg,iend,jbeg,jend,kbeg,kend )
          DO k=kbeg,kend
            DO j=jbeg,jend
              DO i=ibeg,iend
                ijkN = IndIJK(i,j,k,iNOff,ijNOff)
                IF (xyz(XCOORD,ijkN) < global%xyzMinmax(XCOORD,1)) THEN
                  global%xyzMinmax(XCOORD,1) = xyz(XCOORD,ijkN)
                  DO l = XCOORD,ZCOORD
                    global%infloPlanEdges(l,XCOORD,1) = xyz(l,ijkN)
                  ENDDO
                ENDIF
                IF (xyz(XCOORD,ijkN) > global%xyzMinmax(XCOORD,2)) THEN
                  global%xyzMinmax(XCOORD,2) = xyz(XCOORD,ijkN)
                  DO l = XCOORD,ZCOORD
                    global%infloPlanEdges(l,XCOORD,2) = xyz(l,ijkN)
                  ENDDO
                ENDIF
                IF (xyz(YCOORD,ijkN) < global%xyzMinmax(YCOORD,1)) THEN
                  global%xyzMinmax(YCOORD,1) = xyz(YCOORD,ijkN)
                  DO l = XCOORD,ZCOORD
                    global%infloPlanEdges(l,YCOORD,1) = xyz(l,ijkN)
                  ENDDO
                ENDIF
                IF (xyz(YCOORD,ijkN) > global%xyzMinmax(YCOORD,2)) THEN
                  global%xyzMinmax(YCOORD,2) = xyz(YCOORD,ijkN)
                  DO l = XCOORD,ZCOORD
                    global%infloPlanEdges(l,YCOORD,2) = xyz(l,ijkN)
                  ENDDO
                ENDIF
                IF (xyz(ZCOORD,ijkN) < global%xyzMinmax(ZCOORD,1)) THEN
                  global%xyzMinmax(ZCOORD,1) = xyz(ZCOORD,ijkN)
                  DO l = XCOORD,ZCOORD
                    global%infloPlanEdges(l,ZCOORD,1) = xyz(l,ijkN)
                  ENDDO
                ENDIF
                IF (xyz(ZCOORD,ijkN) > global%xyzMinmax(ZCOORD,2)) THEN
                  global%xyzMinmax(ZCOORD,2) = xyz(ZCOORD,ijkN)
                  DO l = XCOORD,ZCOORD
                    global%infloPlanEdges(l,ZCOORD,2) = xyz(l,ijkN)
                  ENDDO
                ENDIF

              ENDDO  ! i
            ENDDO    ! j
          ENDDO      ! k

        ENDIF  ! distrib and not-coupled

      ENDIF    ! bcType, active region on my processor

    ENDDO      ! iPatch
  ENDDO        ! iReg

! copy values/distribution to variables ---------------------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_INFLOW .AND. &
           patch%bcType<=BC_INFLOW+BC_RANGE)    .AND. &   ! my boundary type,
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        switch = patch%mixt%switches(BCSWI_INFLOW_TYPE)
        IF (switch == BCOPT_SUBSONIC) THEN
          patch%mixt%nData = 4
        ELSE
          patch%mixt%nData = 5
        ENDIF

! ----- allocate memory for the values

        IF (patch%mixt%distrib == BCDAT_DISTRIB) THEN
          n1    = ABS(patch%l1end-patch%l1beg)
          n2    = ABS(patch%l2end-patch%l2beg)
          iOff  = n1 + 1
          ijBeg = IndIJ( 0, 0,iOff)
          ijEnd = IndIJ(n1,n2,iOff)
        ELSE
          ijBeg = 0
          ijEnd = 1
        ENDIF
        ALLOCATE( patch%mixt%vals(patch%mixt%nData,ijBeg:ijEnd), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! ----- distribution from file

        IF (patch%mixt%distrib==BCDAT_DISTRIB .AND. &
            patch%bcCoupled      /=BC_EXTERNAL  ) THEN

          IF (patch%mixt%bcSet) &
            WRITE(STDOUT,*)'inlet region and patch #:',iReg,iPatch

          IF (profType==INFLO_TAYLOR_CYL) THEN
            CALL ProfInflowVTTaylorCyl(  global,patch,nvals,switch,keys, &
                                         defined,vals )
          ELSEIF (profType==INFLO_TAYLOR_PLAN) THEN
!            CALL ProfInflowVTTaylorPlan( global,patch,nvals,switch,keys, &
!                                         defined,vals )
            CALL ErrorStop( global,ERR_UNKNOWN_OPTION,__LINE__, &
                           'INFLO_TAYLOR_PLAN' )
          ELSEIF (profType==INFLO_BLLAM_CYL) THEN
!            CALL ProfInflowVTBlLamCyl(   global,patch,nvals,switch,keys, &
!                                         defined,vals )
            CALL ErrorStop( global,ERR_UNKNOWN_OPTION,__LINE__, &
                           'INFLO_BLLAM_CYL' )
          ELSEIF (profType==INFLO_BLLAM_PLAN) THEN
!            CALL ProfInflowVTBlLamPlan(  global,patch,nvals,switch,keys, &
!                                         defined,vals )
            CALL ErrorStop( global,ERR_UNKNOWN_OPTION,__LINE__, &
                           'INFLO_BLLAM_PLAM' )
          ELSEIF (profType==INFLO_BLTURB_CYL) THEN
!            CALL ProfInflowVTBlTurbCyl(  global,patch,nvals,switch,keys, &
!                                         defined,vals )
            CALL ErrorStop( global,ERR_UNKNOWN_OPTION,__LINE__, &
                           'INFLO_BLTURB_CYL' )
          ELSEIF (profType==INFLO_BLTURB_PLAN) THEN
!            CALL ProfInflowVTBlTurbPlan( global,patch,nvals,switch,keys, &
!                                         defined,vals )
            CALL ErrorStop( global,ERR_UNKNOWN_OPTION,__LINE__, &
                           'INFLO_BLTURB_PLAN' )
          ENDIF

          CALL WriteBcToFile( global,fname,patch )

! ----- distribution from external source / constant value

        ELSE
          patch%mixt%vals(BCDAT_INFLOW_U,:) = vals(2)
          patch%mixt%vals(BCDAT_INFLOW_V,:) = vals(3)
          patch%mixt%vals(BCDAT_INFLOW_W,:) = vals(4)
          IF (bcTitle==BC_INFLOW_VELTEMP) THEN
            patch%mixt%vals(BCDAT_INFLOW_T,:) = vals(5)
          ELSEIF (bcTitle==BC_INFLOW_VELPRESS) THEN
            patch%mixt%vals(BCDAT_INFLOW_P,:) = vals(5)
          ENDIF
          IF (switch /= BCOPT_SUBSONIC) THEN
            IF (bcTitle==BC_INFLOW_VELTEMP) THEN
              patch%mixt%vals(BCDAT_INFLOW_P,:) = vals(6)
            ELSEIF (bcTitle==BC_INFLOW_VELPRESS) THEN
              patch%mixt%vals(BCDAT_INFLOW_T,:) = vals(6)
            ENDIF
          ENDIF
        ENDIF  ! distribution?

      ENDIF    ! bcType, active region on my processor

    ENDDO      ! iPatch
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE BcInflowVelDistrib

!******************************************************************************
!
! Purpose: read in user input related to outflow boundary condition.
!          and create bc distribution file if applicable
!
! Description: none.
!
! Input: boundary condition file.
!
! Output: regions = BC data.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE BcOutflowDistrib( regions,ithRead,bcTitle,distrib )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  INTEGER :: ithRead, bcTitle, distrib

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  CHARACTER(10)  :: keys(2)
  CHARACTER(256) :: fname

  INTEGER :: brbeg, brend, prbeg, prend, profType, switch
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag

  LOGICAL :: defined(2)

  REAL(RFREAL) :: vals(2)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'BcOutflowDistrib',&
  'PREP_ModBcDistribution.F90' )

! specify keywords and search for them ----------------------------------------

  keys(1) = 'TYPE'
  keys(2) = 'PRESS'

  distrib  = 0
  profType = 0
  CALL ReadPatchSection( global,IF_INPUT,2,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,profType,fname,defined )

! get switches & check if all necessary values defined ------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_OUTFLOW .AND. &
           patch%bcType<=BC_OUTFLOW+BC_RANGE)    .AND. &  ! my boundary type,
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        IF (patch%mixt%bcSet) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,__LINE__,'Outflow boundary.' )

        patch%mixt%nData     = 0
        patch%mixt%nSwitches = 1
        IF (patch%bcCoupled == BC_EXTERNAL) THEN   ! data from outside
          patch%mixt%distrib = BCDAT_DISTRIB    ! => always distribution
        ELSE
          patch%mixt%distrib = distrib
        ENDIF

        ALLOCATE( patch%mixt%switches(patch%mixt%nSwitches), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! ----- check if switch defined
        IF (defined(1)) THEN
          patch%mixt%switches(BCSWI_OUTFLOW_TYPE)   = BCOPT_SUBSONIC
          IF (vals(1) < 0.1) &
            patch%mixt%switches(BCSWI_OUTFLOW_TYPE) = BCOPT_SUPERSONIC
          IF (vals(1) > 1.9) &
            patch%mixt%switches(BCSWI_OUTFLOW_TYPE) = BCOPT_MIXED
        ELSE
          CALL ErrorStop( global,ERR_NO_BCSWITCH,__LINE__,'(outflow type).' )
        ENDIF

! ----- check if appropriate values specified
        IF (patch%mixt%switches(BCSWI_OUTFLOW_TYPE) /= BCOPT_SUPERSONIC) THEN
          IF (patch%mixt%distrib==BCDAT_CONSTANT .AND. &
              .NOT. defined(2)) CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )
        ENDIF

! ----- set flag to BC specified
        IF (ithRead==2) patch%mixt%bcSet = .true.

      ENDIF   ! my BC & processor, active
    ENDDO
  ENDDO

! copy values/distribution to variables ---------------------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_OUTFLOW .AND. &
           patch%bcType<=BC_OUTFLOW+BC_RANGE)   .AND. &   ! my boundary type,
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor
        switch = patch%mixt%switches(BCSWI_OUTFLOW_TYPE)
      ELSE
        switch = BCOPT_SUPERSONIC
      ENDIF

      IF ((patch%bcType>=BC_OUTFLOW .AND. &
           patch%bcType<=BC_OUTFLOW+BC_RANGE)   .AND. &   ! my boundary type,
           switch/=BCOPT_SUPERSONIC             .AND. &   ! p required
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

! ----- allocate memory for the values

        patch%mixt%nData = 1

        IF (patch%mixt%distrib == BCDAT_DISTRIB) THEN
          n1    = ABS(patch%l1end-patch%l1beg)
          n2    = ABS(patch%l2end-patch%l2beg)
          iOff  = n1 + 1
          ijBeg = IndIJ( 0, 0,iOff)
          ijEnd = IndIJ(n1,n2,iOff)
        ELSE
          ijBeg = 0
          ijEnd = 1
        ENDIF
        ALLOCATE( patch%mixt%vals(patch%mixt%nData,ijBeg:ijEnd), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! ----- distribution from file

        IF (patch%mixt%distrib==BCDAT_DISTRIB .AND. &
            patch%bcCoupled      /=BC_EXTERNAL  ) THEN
!          CALL WriteBcToFile( global,fname,patch )

! ----- distribution from external source / constant value

        ELSE
          patch%mixt%vals(BCDAT_OUTFLOW_PRESS,:) = vals(2)
        ENDIF  ! distribution?

      ENDIF    ! bcType, p req., active region on my processor

    ENDDO      ! iPatch
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE BcOutflowDistrib

!******************************************************************************
!
! Purpose: read in user input related to far field boundary condition.
!          and create bc distribution file if applicable
!
! Description: none.
!
! Input: boundary condition file.
!
! Output: regions = BC data
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE BcFarfDistrib( regions,ithRead,bcTitle,distrib )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  INTEGER :: ithRead, bcTitle, distrib

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  CHARACTER(10)  :: keys(5)
  CHARACTER(256) :: fname

  INTEGER :: brbeg, brend, prbeg, prend, profType
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag

  LOGICAL :: defined(5)

  REAL(RFREAL) :: vals(5)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'BcFarfDistrib',&
  'PREP_ModBcDistribution.F90' )

! specify keywords and search for them ----------------------------------------

  keys(1) = 'MACH'
  keys(2) = 'ATTACK'
  keys(3) = 'SLIP'
  keys(4) = 'PRESS'
  keys(5) = 'TEMP'

  distrib  = 0
  profType = 0
  CALL ReadPatchSection( global,IF_INPUT,5,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,profType,fname,defined )

! check if all values defined -------------------------------------------------

  IF (distrib==BCDAT_CONSTANT .AND. &
      (.NOT. defined(1) .OR. &
       .NOT. defined(2) .OR. &
       .NOT. defined(3) .OR. &
       .NOT. defined(4) .OR. &
       .NOT. defined(5))) CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )

! copy values/distribution to variables ---------------------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_FARFIELD .AND. &
           patch%bcType<=BC_FARFIELD+BC_RANGE)  .AND. &   ! my boundary type
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        IF (patch%mixt%bcSet) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,__LINE__,'Farfield boundary.' )

        patch%mixt%nData     = 5
        patch%mixt%nSwitches = 0
        IF (ithRead==2) patch%mixt%bcSet = .true.

        IF (patch%bcCoupled == BC_EXTERNAL) THEN   ! data from outside
          patch%mixt%distrib = BCDAT_DISTRIB    ! => always distribution
        ELSE
          patch%mixt%distrib = distrib
        ENDIF

! ----- allocate memory for the values

        IF (patch%mixt%distrib == BCDAT_DISTRIB) THEN
          n1    = ABS(patch%l1end-patch%l1beg)
          n2    = ABS(patch%l2end-patch%l2beg)
          iOff  = n1 + 1
          ijBeg = IndIJ( 0, 0,iOff)
          ijEnd = IndIJ(n1,n2,iOff)
        ELSE
          ijBeg = 0
          ijEnd = 1
        ENDIF
        ALLOCATE( patch%mixt%vals(patch%mixt%nData,ijBeg:ijEnd), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! ----- distribution from file

        IF (patch%mixt%distrib==BCDAT_DISTRIB .AND. &
            patch%bcCoupled      /=BC_EXTERNAL  ) THEN
!          CALL WriteBcToFile( global,fname,patch )

          patch%mixt%vals(BCDAT_FARF_ATTACK,:) = &
            patch%mixt%vals(BCDAT_FARF_ATTACK,:)*global%rad
          patch%mixt%vals(BCDAT_FARF_SLIP  ,:) = &
            patch%mixt%vals(BCDAT_FARF_SLIP  ,:)*global%rad

! ----- distribution from external source / constant value

        ELSE
          patch%mixt%vals(BCDAT_FARF_MACH  ,:) = vals(1)
          patch%mixt%vals(BCDAT_FARF_ATTACK,:) = vals(2)*global%rad
          patch%mixt%vals(BCDAT_FARF_SLIP  ,:) = vals(3)*global%rad
          patch%mixt%vals(BCDAT_FARF_PRESS ,:) = vals(4)
          patch%mixt%vals(BCDAT_FARF_TEMP  ,:) = vals(5)
        ENDIF  ! distribution?

      ENDIF    ! bcType, active region on my processor

    ENDDO      ! iPatch
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE BcFarfDistrib

!******************************************************************************
!
! Purpose: read in user input related to injection boundary condition.
!          and create bc distribution file if applicable
!
! Description: none.
!
! Input: boundary condition file.
!
! Output: regions = BC data.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE BcInjectDistrib( regions,ithRead,bcTitle,distrib )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  INTEGER :: ithRead, bcTitle, distrib

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  CHARACTER(10)  :: keys(4)
  CHARACTER(256) :: fname

  INTEGER :: brbeg, brend, prbeg, prend, profType
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag

  LOGICAL :: defined(4)

  REAL(RFREAL) :: vals(4)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'BcInjectDistrib',&
  'PREP_ModBcDistribution.F90' )

! specify keywords and search for them ----------------------------------------

  keys(1) = 'MFRATE'
  keys(2) = 'TEMP'
  keys(3) = 'EXTRAPOL'
  keys(4) = 'MAXCHANGE'

  distrib  = 0
  profType = 0
  CALL ReadPatchSection( global,IF_INPUT,4,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,profType,fname,defined )
  
! check if all values defined -------------------------------------------------

  IF (distrib==BCDAT_CONSTANT .AND. &
      (.NOT. defined(1) .OR. &
       .NOT. defined(2))) CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )

  IF (.NOT. defined(3) .OR. &
      .NOT. defined(4))   CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )

! copy values/distribution to variables ---------------------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_INJECTION .AND. &
           patch%bcType<=BC_INJECTION+BC_RANGE) .AND. &   ! my boundary type
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        IF (patch%mixt%bcSet) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,__LINE__,'Injection boundary.' )

        patch%mixt%nData     = 5
        patch%mixt%nSwitches = 1
        IF (ithRead==2) patch%mixt%bcSet = .true.

        IF (patch%bcCoupled == BC_EXTERNAL) THEN   ! data from outside
          patch%mixt%distrib = BCDAT_DISTRIB    ! => always distribution
        ELSE
          patch%mixt%distrib = distrib
        ENDIF

! ----- get value of switch

        ALLOCATE( patch%mixt%switches(patch%mixt%nSwitches), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

          patch%mixt%switches(BCSWI_INJECT_EXTRAP) = EXTRAPOL_CONST
        IF (vals(3) > 0.1) &
          patch%mixt%switches(BCSWI_INJECT_EXTRAP) = EXTRAPOL_LINEAR

        patch%mixt%maxChange = vals(4)

! ----- allocate memory for the values

        IF (patch%mixt%distrib == BCDAT_DISTRIB) THEN
          n1    = ABS(patch%l1end-patch%l1beg)
          n2    = ABS(patch%l2end-patch%l2beg)
          iOff  = n1 + 1
          ijBeg = IndIJ( 0, 0,iOff)
          ijEnd = IndIJ(n1,n2,iOff)
        ELSE
          ijBeg = 0
          ijEnd = 1
        ENDIF
        ALLOCATE( patch%mixt%vals(patch%mixt%nData,ijBeg:ijEnd), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! ----- distribution from file

        IF (patch%mixt%distrib==BCDAT_DISTRIB .AND. &
            patch%bcCoupled      /=BC_EXTERNAL  ) THEN
          patch%mixt%nData = 2
!          CALL WriteBcToFile( global,fname,patch )
          patch%mixt%nData = 5

! ----- distribution from external source / constant value

        ELSE
          patch%mixt%vals(BCDAT_INJECT_MFRATE,:) = vals(1)
          patch%mixt%vals(BCDAT_INJECT_TEMP  ,:) = vals(2)
        ENDIF  ! distribution?

      ENDIF    ! bcType, active region on my processor

    ENDDO      ! iPatch
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE BcInjectDistrib

!******************************************************************************
!
! Purpose: read in a section of a file (until # is encountered), read
!          keywords and store the associated numerical values.
!
! Description:
!  - ReadPatchSection = section applies to a range of patches (prbeg:prend)
!                       within a range of regions (brbeg:brend)
!
! Input: fileID = file number
!        nvals  = number of values to search for and to store
!        keys   = keywords to search for
!
! Output: vals    = values associated with keywords (reals only)
!         defined = flag if for certain keyword a value was read in
!         brbeg   = begin of region range (values set for these regions)
!         brend   = end of region range
!         prbeg   = begin of patch range (values set for these patches)
!         prend   = end of patch range
!         distrib = single value for a patch (=0) or distribution (>0)
!         profType= profile type of vals data
!         fname   = file with distribution for a patch
!
!******************************************************************************

SUBROUTINE ReadPatchSection( global,fileID,nvals,keys,vals,brbeg,brend, &
                             prbeg,prend,distrib,profType,fname,defined )

  IMPLICIT NONE

! ... parameters
  INTEGER      :: brbeg, brend
  INTEGER      :: fileID, nvals, prbeg, prend, distrib, profType
  CHARACTER(*) :: keys(nvals), fname

  LOGICAL      :: defined(nvals)
  REAL(RFREAL) :: vals(nvals)
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: ival

! ... local variables
  CHARACTER(256)    :: line

  INTEGER :: errorFlag, nc

!******************************************************************************

  CALL RegisterFunction( global,'ReadPatchSection',&
  'PREP_ModBcDistribution.F90' )

! read lines from file until # or EOF found

  brbeg = 1               ! region range: input applies to all regions (default)
  brend = global%nRegions

  prbeg = 1               ! patch range: input applies to all patches (default)
  prend = 999999          ! can have different # of patches in each region

  distrib    = 0          ! no distribution as a default
  fname      = ''         ! no file name

  IF ( nvals /= 0 ) THEN
    defined(:) = .false.  ! keeps track of values being provided by the user
  END IF ! nvals

  DO
    READ(fileID,'(A256)',iostat=errorFlag) line
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_READ,__LINE__ )
    IF (line(1:1) == '#') EXIT

    IF (line(1:5) == 'BLOCK') THEN
      READ(line(6:256),*) brbeg,brend
      brend = MIN(brend,global%nRegions)
      IF (brbeg <= 0    ) brbeg = 1
      IF (brend <= 0    ) brend = global%nRegions
      IF (brend <  brbeg) brend = brbeg
    ELSE IF (line(1:5) == 'PATCH') THEN
      READ(line(6:256),*) prbeg,prend
      IF (prbeg <= 0    ) prbeg = 1
      IF (prend <= 0    ) prend = 999999
      IF (prend <  prbeg) prend = prbeg
    ELSE IF (line(1:7) == 'DISTRIB') THEN
      READ(line(8:256),*) distrib
      distrib = MAX(distrib,0)
      distrib = MIN(distrib,1)
    ELSE IF (line(1:7) == 'PROFILE') THEN
      READ(line(8:256),*) profType
    ELSE IF (line(1:4) == 'FILE') THEN
      READ(line(5:256),*) fname
    ELSE
      DO ival=1,nvals
        nc = LEN_TRIM(keys(ival))
        IF (line(1:nc) == TRIM(keys(ival))) THEN    ! found matching keyword
          READ(line(nc+1:256),*) vals(ival)
          defined(ival) = .true.
        ENDIF
      ENDDO
    ENDIF
  ENDDO

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadPatchSection

!******************************************************************************
!
! Purpose: write boundary condition data to a file.
!
! Description: none.
!
! Input: global   = global variables (needed for error function)
!        fname    = file name
!        patch    = BC patch for which the data is to be read in.
!
! Output: patch%mixt%vals = BC data for the mixture.
!
! Notes: currently only the mixture BC data are written out.
!
!******************************************************************************

SUBROUTINE WriteBcToFile( global,fname,patch )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  CHARACTER(*) :: fname

  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch

! ... loop variables
  INTEGER :: iReg, iPatch, n, i, j, ij

! ... local variables
  INTEGER :: n1, n2, iOff, errorFlag

!******************************************************************************

  CALL RegisterFunction( global,'WriteBcToFile',&
  'PREP_ModBcDistribution.F90' )

! dimensions

  n1   = ABS(patch%l1end-patch%l1beg)
  n2   = ABS(patch%l2end-patch%l2beg)
  iOff = n1 + 1

! write to file

  OPEN(IF_DISTR,file=fname,form='formatted',status='unknown',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

  WRITE(IF_DISTR,*,err=10) n1+1,n2+1

  DO n=1,patch%mixt%nData
    DO j=0,n2
      DO i=0,n1
        ij = IndIJ(i,j,iOff)
        WRITE(IF_DISTR,*,err=10) patch%mixt%vals(n,ij)
      ENDDO
    ENDDO
  ENDDO

  CLOSE(IF_DISTR,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )

GOTO 999

! error handling

10  CONTINUE
  CALL ErrorStop( global,ERR_FILE_WRITE,__LINE__,'File: '//TRIM(fname) )

999 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE WriteBcToFile

!******************************************************************************
!
! Purpose: write boundary condition data to a file.
!
! Description: none.
!
! Input: global   = global variables (needed for error function)
!        patch    = BC patch for which the data is to be read in.
!        nvals    = number of values to search for and to store
!        keys     = keywords to search for
!        defined  = flag if for certain keyword a value was read in
!
! Output: patch%mixt%vals = BC data for the mixture.
!
! Notes: currently only the mixture BC data are written out.
!
!******************************************************************************

SUBROUTINE ProfInflowVTTaylorCyl( global,patch,nvals,switch,keys,defined,vals )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch
  INTEGER      :: nvals, switch
  CHARACTER(*) :: keys(nvals)
  LOGICAL      :: defined(nvals)

  REAL(RFREAL) :: vals(nvals)

! ... loop variables
  INTEGER :: iReg, iPatch, n, i, j, ij

! ... local variables
  INTEGER :: n1, n2, iOff, nf1, nf2, lbound, errorFlag
  INTEGER :: ivelx, ively, ivelz, itemp, ipres, ipow, ifact
  INTEGER :: in, jn, ip, jp, ng1, ng2
  REAL(RFREAL) :: vaxi, vinj, xc, yc, zc, xs, ys, zs, radius, rd
  REAL(RFREAL) :: sign, costerm, sinterm, cosa, sina
  REAL(RFREAL) :: ymin(XCOORD:ZCOORD), ymax(XCOORD:ZCOORD)
  REAL(RFREAL) :: zmin(XCOORD:ZCOORD), zmax(XCOORD:ZCOORD)

!******************************************************************************

  CALL RegisterFunction( global,'ProfInflowVTTaylorCyl',&
  'PREP_ModBcDistribution.F90' )

! remember keys

  DO i=1,nvals
    IF (keys(i)=='VELX')        ivelx =i
    IF (keys(i)=='VELY')        ively =i
    IF (keys(i)=='VELZ')        ivelz =i
    IF (keys(i)=='TEMP')        itemp =i
    IF (keys(i)=='PRESS')       ipres =i
    IF (keys(i)=='AXIALPOWER')  ipow  =i
    IF (keys(i)=='NORMALFACT')  ifact =i
  ENDDO

  lbound = patch%lbound

  IF (lbound==1 .OR. lbound==2) THEN
    vaxi = vals(ivelx)
    vinj = vals(ively)
  ELSEIF (lbound==3 .OR. lbound==4) THEN
    vaxi = vals(ively)
    vinj = vals(ivelz)
  ELSEIF (lbound==5 .OR. lbound==6) THEN
    vaxi = vals(ivelz)
    vinj = vals(ivelx)
  ENDIF

! dimensions

  n1   = ABS(patch%l1end-patch%l1beg)
  n2   = ABS(patch%l2end-patch%l2beg)
  iOff = n1 + 1

! define cirkel center and radius

  IF (lbound/=3 .AND. lbound/=4) THEN
    ymin(XCOORD:ZCOORD) = global%infloPlanEdges(XCOORD:ZCOORD,YCOORD,1)
    ymax(XCOORD:ZCOORD) = global%infloPlanEdges(XCOORD:ZCOORD,YCOORD,2)
    xc     = 0.5_RFREAL*(ymin(XCOORD) + ymax(XCOORD))
    yc     = 0.5_RFREAL*(ymin(YCOORD) + ymax(YCOORD))
    zc     = 0.5_RFREAL*(ymin(ZCOORD) + ymax(ZCOORD))
    radius = 0.5_RFREAL*SQRT( (ymax(XCOORD)-ymin(XCOORD))**2 + &
                              (ymax(YCOORD)-ymin(YCOORD))**2 + &
                              (ymax(ZCOORD)-ymin(ZCOORD))**2 )
  ELSEIF (lbound==3 .OR. lbound==4) THEN
    zmin(XCOORD:ZCOORD) = global%infloPlanEdges(XCOORD:ZCOORD,ZCOORD,1)
    zmax(XCOORD:ZCOORD) = global%infloPlanEdges(XCOORD:ZCOORD,ZCOORD,2)
    xc     = 0.5_RFREAL*(zmin(XCOORD) + zmax(XCOORD))
    yc     = 0.5_RFREAL*(zmin(YCOORD) + zmax(YCOORD))
    zc     = 0.5_RFREAL*(zmin(ZCOORD) + zmax(ZCOORD))
    radius = 0.5_RFREAL*SQRT( (zmax(XCOORD)-zmin(XCOORD))**2 + &
                              (zmax(YCOORD)-zmin(YCOORD))**2 + &
                              (zmax(ZCOORD)-zmin(ZCOORD))**2 )
  ENDIF

  IF (patch%mixt%bcSet) &
    WRITE(STDOUT,100)'(x,y,z)_cirkel-center, radius:', xc,yc,zc,radius

! compute inflow profiles

  DO j=0,n2
    DO i=0,n1
      ij = IndIJ(i,j,iOff)

      IF      (lbound==1 .OR. lbound==2) THEN
        IF (lbound == 2) THEN
          ng1 = i - 0 + 1
        ELSE
          ng1 = n1 - i + 1
        ENDIF
        ng2 = j - 0 + 1
      ELSE IF (lbound==3 .OR. lbound==4) THEN
        ng1 = i - 0 + 1
        IF (lbound == 4) THEN
          ng2 = j - 0 + 1
        ELSE
          ng2 = n2 - j + 1
        ENDIF
      ELSE IF (lbound==5 .OR. lbound==6) THEN
        IF (lbound == 6) THEN
          ng1 = i - 0 + 1
        ELSE
          ng1 = n1 - i + 1
        ENDIF
        ng2 = j - 0 + 1
      ENDIF

      in = ng1
      jn = ng2
      ip = ng1+1
      jp = ng2+1
      xs = 0.25_RFREAL*(patch%surfCoord(XCOORD,in,jn) + &
                        patch%surfCoord(XCOORD,ip,jn) + &
                        patch%surfCoord(XCOORD,in,jp) + &
                        patch%surfCoord(XCOORD,ip,jp))
      ys = 0.25_RFREAL*(patch%surfCoord(YCOORD,in,jn) + &
                        patch%surfCoord(YCOORD,ip,jn) + &
                        patch%surfCoord(YCOORD,in,jp) + &
                        patch%surfCoord(YCOORD,ip,jp))
      zs = 0.25_RFREAL*(patch%surfCoord(ZCOORD,in,jn) + &
                        patch%surfCoord(ZCOORD,ip,jn) + &
                        patch%surfCoord(ZCOORD,in,jp) + &
                        patch%surfCoord(ZCOORD,ip,jp))
      rd = SQRT((xs-xc)**2+(ys-yc)**2+(zs-zc)**2)/radius

      costerm = COS( 0.5_RFREAL*global%pi*rd**vals(ipow) )
      sinterm = SIN( 0.5_RFREAL/vals(ifact)*global%pi*rd*rd )/rd/ &
                SIN( 0.5_RFREAL/vals(ifact)*global%pi )

      sign    = 1._RFREAL
      IF (lbound==1 .OR. lbound==2) THEN
        cosa = (ys-yc)/rd/radius
        sina = (zs-zc)/rd/radius
        IF (lbound==2) sign=-1._RFREAL
        patch%mixt%vals(BCDAT_INFLOW_U,ij) =  vaxi*costerm*sign
        patch%mixt%vals(BCDAT_INFLOW_V,ij) = -vinj*sinterm*cosa
        patch%mixt%vals(BCDAT_INFLOW_W,ij) = -vinj*sinterm*sina
      ELSEIF (lbound==3 .OR. lbound==4) THEN
        cosa = (zs-zc)/rd/radius
        sina = (xs-xc)/rd/radius
        IF (lbound==4) sign=-1._RFREAL
        patch%mixt%vals(BCDAT_INFLOW_V,ij) =  vaxi*costerm*sign
        patch%mixt%vals(BCDAT_INFLOW_W,ij) = -vinj*sinterm*cosa
        patch%mixt%vals(BCDAT_INFLOW_U,ij) = -vinj*sinterm*sina
      ELSEIF (lbound==5 .OR. lbound==6) THEN
        cosa = (xs-xc)/rd/radius
        sina = (ys-yc)/rd/radius
        IF (lbound==6) sign=-1._RFREAL
        patch%mixt%vals(BCDAT_INFLOW_W,ij) =  vaxi*costerm*sign
        patch%mixt%vals(BCDAT_INFLOW_U,ij) = -vinj*sinterm*cosa
        patch%mixt%vals(BCDAT_INFLOW_V,ij) = -vinj*sinterm*sina
      ENDIF
      patch%mixt%vals(BCDAT_INFLOW_T,ij) = vals(itemp)
      IF (switch/=BCOPT_SUBSONIC) &
        patch%mixt%vals(BCDAT_INFLOW_P,ij) = vals(ipres)
    ENDDO
  ENDDO

100 FORMAT( A,4E17.10 )

  CALL DeregisterFunction( global )

END SUBROUTINE ProfInflowVTTaylorCyl

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE PREP_ModBcDistribution

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PREP_ModBcDistribution.F90,v $
! Revision 1.9  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2006/08/19 15:41:08  mparmar
! Renamed patch variables
!
! Revision 1.6  2005/05/06 00:32:52  wasistho
! fixed orientation of patch%surfCoord TaylorVTCyl routine
!
! Revision 1.5  2005/05/04 19:04:35  wasistho
! fixed bug, integer headers n1,n2 in WriteBcToFile
!
! Revision 1.4  2005/05/03 08:18:22  wasistho
! make more efficient
!
! Revision 1.3  2005/05/03 03:20:09  wasistho
! enabled modified cyl.Taylor inflow profile
!
! Revision 1.2  2005/05/02 18:08:48  wasistho
! added cylindrical Taylor inflow profile capability
!
! Revision 1.1  2005/04/29 03:32:33  wasistho
! added distribution bc file generator
!
!
! ******************************************************************************

















