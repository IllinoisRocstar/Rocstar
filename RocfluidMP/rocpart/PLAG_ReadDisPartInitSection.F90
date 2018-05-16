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
! Purpose: read in user input related to discrete particle module 
!          for initialization phase.
!
! Description: none.
!
! Input: user input file.
!
! Output: regions = number of initial particles, their positions
!                  diameters and temperatures.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_ReadDisPartInitSection.F90,v 1.9 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_ReadDisPartInitSection( regions )

  USE ModDataTypes 
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag_input
  USE ModInterfaces, ONLY : ReadListSection, ReadSection
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
! Local variables
! ==============================================================================

  INTEGER, PARAMETER :: NVALS_MAX = 2

  CHARACTER(CHRLEN) :: RCSIdentString
  CHARACTER(15)     :: keysInit(NVALS_MAX),keysInitMat

  INTEGER :: brbeg,brend,errorFlag,iReg,iRow1,iRow2,iVal,nCols,&
             nPclsIni,nRegions,nRows,nVals

  LOGICAL :: definedInit(NVALS_MAX),definedInitMat
  
  REAL(RFREAL) :: valsInit(NVALS_MAX)
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: valsInitMat
  
  TYPE(t_global),   POINTER :: global 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ReadDisPartInitSection.F90,v $ $Revision: 1.9 $'

  global => regions(1)%global
  
  CALL RegisterFunction( global, 'PLAG_ReadDisPartInitSection',&
  'PLAG_ReadDisPartInitSection.F90' )
 
! ******************************************************************************
! Initialize
! ******************************************************************************
 
  nRegions = global%nRegions
   
#ifdef RFLO    
  brbeg = 1
  brend = nRegions
#endif
#ifdef RFLU
  brbeg = LBOUND(regions,1)
  brend = UBOUND(regions,1)
#endif

! ******************************************************************************
! Read section pertinent to initial solution flag
! ******************************************************************************

  nVals          = NVALS_MAX
  definedInit(:) = .FALSE.
  keysInit(1)    = 'FLAG'
  keysInit(2)    = 'NPCLSRAND'
  valsInit(:)    = 0.0_RFREAL

#ifdef RFLO
  CALL ReadregionSection( global,IF_INPUT,nVals,keysInit,valsInit, &
                          brbeg,brend,definedInit )

  IF (definedInit(1) .EQV. .TRUE.) &
    global%initPlagFlag = NINT(valsInit(1))
  
  IF (definedInit(2) .EQV. .TRUE.) &
    regions(brbeg:brend)%plagInput%nPclsIni = NINT(valsInit(2))
#endif

#ifdef RFLU
  CALL ReadSection(global,IF_INPUT,nVals,keysInit,valsInit,definedInit )

  IF ( definedInit(1) .EQV. .FALSE. ) THEN
    CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__,'DISPARTINIT-FLAG')
  ELSE
    IF ( NINT(valsInit(1)) == PLAG_INIT_FROMSCRATCH ) THEN 
      global%initPlagFlag = PLAG_INIT_FROMSCRATCH
    ELSE IF ( NINT(valsInit(1)) == PLAG_INIT_FROMFILE ) THEN 
      global%initPlagFlag = PLAG_INIT_FROMFILE
    ELSE IF ( NINT(valsInit(1)) == PLAG_INIT_FROMHARDCODE ) THEN 
      global%initPlagFlag = PLAG_INIT_FROMHARDCODE
    ELSE IF ( NINT(valsInit(1)) == PLAG_INIT_FROMRANDOMSTATE ) THEN 
      global%initPlagFlag = PLAG_INIT_FROMRANDOMSTATE
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! valsInit     
  END IF ! defined

  IF ( ( NINT(valsInit(1)) == PLAG_INIT_FROMRANDOMSTATE ) .AND. &
       ( definedInit(2) .EQV. .FALSE. )                         ) THEN
    CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__,'DISPARTINIT-NPCLSRAND')
  ELSE
    IF ( NINT(valsInit(2)) >= 0 ) THEN 
      regions(brbeg:brend)%plagInput%nPclsIni =  NINT(valsInit(2))
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! valsInit     
  END IF ! iValInit 
#endif

! ******************************************************************************
! Read section pertinent to particle positions, diameters and temperature
! ******************************************************************************

  SELECT CASE ( NINT(valsInit(1)) )

! ==============================================================================
!   PLAG_INIT_FROMSCRATCH:
!   Read section pertinent to scratch initialization
! ==============================================================================

    CASE ( PLAG_INIT_FROMSCRATCH ) 
      definedInitMat = .FALSE.
      keysInitMat    = 'NUMBER'
      nCols          = 9

      CALL ReadListSection( global, IF_INPUT,keysInitMat,nCols,nRows, &
                            valsInitMat,definedInitMat )

! ==============================================================================
!     Always load number of particles for scratch field
! ==============================================================================

      regions(brbeg:brend)%plagInput%nPclsIni = nRows

      IF (definedInitMat .EQV. .TRUE.) THEN

! ==============================================================================
!       Allocate arrays
! ==============================================================================

        DO iReg = brbeg,brend
          ALLOCATE( regions(iReg)%plagInput%iniPosX(nRows),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ , &
                            'regions%plagInput%iniPosX' ) 
          END IF ! global%error

          ALLOCATE( regions(iReg)%plagInput%iniPosY(nRows),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ , &
                            'regions%plagInput%iniPosY' ) 
          END IF ! global%error

          ALLOCATE( regions(iReg)%plagInput%iniPosZ(nRows),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ , &
                            'regions%plagInput%iniPosZ' ) 
          END IF ! global%error

          ALLOCATE( regions(iReg)%plagInput%iniDiam(nRows),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,&
                            'regions%plagInput%iniDiam' ) 
          END IF ! global%error

          ALLOCATE( regions(iReg)%plagInput%iniTemp(nRows),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ , &
                          'regions%plagInput%iniTemp' )
          END IF ! global%error
 
          ALLOCATE( regions(iReg)%plagInput%iniSpLoad(nRows),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ , &
                            'regions%plagInput%iniSpLoad' )
          END IF ! global%error

          ALLOCATE( regions(iReg)%plagInput%iniVelX(nRows),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ , &
                            'regions%plagInput%iniVelX' ) 
          END IF ! global%error

          ALLOCATE( regions(iReg)%plagInput%iniVelY(nRows),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ , &
                            'regions%plagInput%iniVelY' ) 
          END IF ! global%error

          ALLOCATE( regions(iReg)%plagInput%iniVelZ(nRows),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) THEN
            CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ , &
                            'regions%plagInput%iniVelZ' ) 
          END IF ! global%error
        END DO !iReg

! ==============================================================================
!       Initialize arrays
! ==============================================================================

        DO iReg = brbeg,brend
          regions(iReg)%plagInput%iniPosX   = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniPosY   = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniPosZ   = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniDiam   = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniTemp   = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniSpLoad = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniVelX   = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniVelY   = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniVelZ  = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        END DO ! iReg

! ==============================================================================
!       Load arrays
! ==============================================================================
       
        DO iReg= brbeg,brend
          nPclsIni = regions(iReg)%plagInput%nPclsIni  

          DO iVal=1, nPclsIni 
            regions(iReg)%plagInput%iniPosX(iVal) = valsInitMat(iVal,1)
            regions(iReg)%plagInput%iniPosY(iVal) = valsInitMat(iVal,2)
            regions(iReg)%plagInput%iniPosZ(iVal) = valsInitMat(iVal,3)
            regions(iReg)%plagInput%iniDiam(iVal) = ABS(valsInitMat(iVal,4))
            regions(iReg)%plagInput%iniTemp(iVal) = ABS(valsInitMat(iVal,5))
            regions(iReg)%plagInput%iniSpLoad(iVal) = ABS(valsInitMat(iVal,6))
            regions(iReg)%plagInput%iniVelX(iVal) = valsInitMat(iVal,7)
            regions(iReg)%plagInput%iniVelY(iVal) = valsInitMat(iVal,8) 
            regions(iReg)%plagInput%iniVelZ(iVal) = valsInitMat(iVal,9)
          ENDDO ! iVal
        ENDDO ! iReg
      ENDIF ! definedInitMat 

! ******************************************************************************
!   PLAG_INIT_FROMRANDOM:
!   Read section pertinent to random state initialization
! ******************************************************************************

    CASE ( PLAG_INIT_FROMRANDOMSTATE ) 
      definedInitMat = .FALSE.
      keysInitMat    = 'NUMBER'
      nCols          = 9

      CALL ReadListSection( global, IF_INPUT,keysInitMat,nCols,nRows, &
                            valsInitMat,definedInitMat )

      IF (definedInitMat .EQV. .TRUE.) THEN

        IF ( nRows /= 2 ) THEN
          CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__,'DISPARTINIT-Random nRows')
        ENDIF ! nRows

! ==============================================================================
!       Initialize arrays
! ==============================================================================

        DO iReg = brbeg,brend
          regions(iReg)%plagInput%iniRandDiamMax   = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandDiamMin   = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandTempMax   = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandTempMin   = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandSploadMax = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandSploadMin = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandXMax      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandXMin      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandYMax      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandYMin      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandZMax      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandZMin      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandUMax      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandUMin      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandVMax      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandVMin      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandWMax      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
          regions(iReg)%plagInput%iniRandWMin      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        END DO ! iReg

! ==============================================================================
!       Load arrays
!       First row contains minimum values
!       Second row contains maximum values
! ==============================================================================

        iRow1 = 1
        iRow2 = 2       

        DO iReg= brbeg,brend
          regions(iReg)%plagInput%iniRandXMin      =     valsInitMat(iRow1,1) 
          regions(iReg)%plagInput%iniRandYMin      =     valsInitMat(iRow1,2) 
          regions(iReg)%plagInput%iniRandZMin      =     valsInitMat(iRow1,3) 
          regions(iReg)%plagInput%iniRandDiamMin   = ABS(valsInitMat(iRow1,4))
          regions(iReg)%plagInput%iniRandTempMin   = ABS(valsInitMat(iRow1,5))
          regions(iReg)%plagInput%iniRandSploadMin = ABS(valsInitMat(iRow1,6))
          regions(iReg)%plagInput%iniRandUMin      =     valsInitMat(iRow1,7)
          regions(iReg)%plagInput%iniRandVMin      =     valsInitMat(iRow1,8)
          regions(iReg)%plagInput%iniRandWMin      =     valsInitMat(iRow1,9)
          regions(iReg)%plagInput%iniRandXMax      =     valsInitMat(iRow2,1) 
          regions(iReg)%plagInput%iniRandYMax      =     valsInitMat(iRow2,2) 
          regions(iReg)%plagInput%iniRandZMax      =     valsInitMat(iRow2,3) 
          regions(iReg)%plagInput%iniRandDiamMax   = ABS(valsInitMat(iRow2,4))
          regions(iReg)%plagInput%iniRandTempMax   = ABS(valsInitMat(iRow2,5))
          regions(iReg)%plagInput%iniRandSploadMax = ABS(valsInitMat(iRow2,6))   
          regions(iReg)%plagInput%iniRandUMax      =     valsInitMat(iRow2,7)
          regions(iReg)%plagInput%iniRandVMax      =     valsInitMat(iRow2,8)
          regions(iReg)%plagInput%iniRandWMax      =     valsInitMat(iRow2,9)
        ENDDO ! iReg

      ELSE 
        CALL ErrorStop(global,ERR_MISSING_VALUE,__LINE__)

      ENDIF ! definedInitMat

  END SELECT ! valsInit

! ******************************************************************************
!  Deallocate pointer array for non-null nRows
! ******************************************************************************

   IF (nRows > 0) THEN 
     DEALLOCATE( valsInitMat, stat=errorFlag ) 
     global%error = errorFlag
     IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__ ,'valsInitMat' ) 
     END IF ! global%error    
   ENDIF ! nRows

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_ReadDisPartInitSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ReadDisPartInitSection.F90,v $
! Revision 1.9  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2006/10/26 14:58:51  fnajjar
! Added initial random min-max velocities with proper initialization
!
! Revision 1.6  2006/07/17 15:50:17  fnajjar
! Removed ABS from the x, y and z coordinate extents
!
! Revision 1.5  2006/07/15 21:53:54  fnajjar
! Added min-max extents for x-y-z coordinates for random field initialization
!
! Revision 1.4  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.3  2005/12/01 18:39:54  fnajjar
! Added error trap for missing input value and set nPclsIni to nRows when initializing from scratch
!
! Revision 1.2  2005/03/31 20:25:59  fnajjar
! Added initial particle velocities for scratch solution
!
! Revision 1.1  2004/12/01 20:58:03  fnajjar
! Initial revision after changing case
!
! Revision 1.7  2004/11/06 18:47:48  fnajjar
! Bug fix: deallocated valsInitMat array for non-zero nRows
!
! Revision 1.6  2004/10/11 22:12:23  haselbac
! Bug fix
!
! Revision 1.5  2004/10/11 19:38:48  fnajjar
! Renamed ininPlag to nPclsIni to follow naming convention
!
! Revision 1.4  2004/10/09 16:39:05  fnajjar
! Streamlined routine and included reading variables for random state initialization
!
! Revision 1.3  2004/08/23 20:09:47  fnajjar
! Removed ABS from iniPos to allow negative positions
!
! Revision 1.2  2004/08/20 23:35:02  fnajjar
! Moved deallocation inside definedInitMat IF statement for null ininPlag
!
! Revision 1.1  2004/08/20 23:29:24  fnajjar
! Initial import of Plag prep tool
!
!******************************************************************************







