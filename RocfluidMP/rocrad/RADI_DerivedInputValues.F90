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
! Purpose: Set values derived from user input.
!
! Description: Derived variables are set based on user input parameters.
!              Derived variables are components of radiInput data type 
!              (region%radiInput).
!
! Input: regions = input parameters for all regions.
!
! Output: regions = derived variables stored as part of radiInput data.
!
! Notes: Unlike mixture, derived parameters/variables are not necessarily dv.
!
!******************************************************************************
!
! $Id: RADI_DerivedInputValues.F90,v 1.4 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

#ifdef RFLO
SUBROUTINE RADI_DerivedInputValues( regions )
#endif
#ifdef RFLU
SUBROUTINE RADI_DerivedInputValues
#endif

  USE ModDataTypes
#ifdef RFLO  
  USE ModDataStruct, ONLY : t_region
#endif
  USE ModGlobal, ONLY     : t_global
  USE ModRadiation, ONLY  : t_radi_input
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

#ifdef RFLO
! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, l, m, n
#endif

! ... local variables
  TYPE(t_global), POINTER     :: global
  TYPE(t_radi_input), POINTER :: input

  INTEGER      :: errorFlag, angles(3,2)
  REAL(RFREAL) :: pi, twopi

!******************************************************************************

  global => regions(1)%global
  CALL RegisterFunction( global,'RADI_DerivedInputValues',&
  'RADI_DerivedInputValues.F90' )

! set local constants ---------------------------------------------------------

  pi    = global%pi
  twopi = 2._RFREAL*pi

#ifdef RFLO
! global values ---------------------------------------------------------------

! region related data (all levels) --------------------------------------------

  DO iReg=1,global%nRegions

    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

      input => regions(iReg)%radiInput

! --- Stefan-Boltzman constant
      input%stBoltz = 5.67E-8_RFREAL

! --- set logical value of radiUsed

      IF (input%radiModel /= RADI_MODEL_NONE) THEN
        regions(iReg)%mixtInput%radiUsed = .TRUE.
      ENDIF

! --- set number of cv, dv and grad components

      input%nCv   = 0
      input%nGrad = 0 
      IF (input%radiModel == RADI_MODEL_FLDTRAN) THEN
        input%nCv   = 1
        input%nGrad = 3
      ENDIF

      input%nDv   = 0
      IF ((input%radiModel /= RADI_MODEL_ROSS)    .AND. &
          (input%radiModel /= RADI_MODEL_FLDSRC)  .AND. &
          (input%radiModel /= RADI_MODEL_FLDTRAN)) THEN
        input%nDv   = 1
      ENDIF

! --- discrete ordinates and/or intensity angles of RTE models

      IF ((input%radiModel == RADI_MODEL_RTEGRAY)  .OR. &
          (input%radiModel == RADI_MODEL_RTEBAND)) THEN

        IF (input%solMethod == RADI_NUM_DOM4) THEN
          input%nOrdin = 4
          input%nAng   = input%nOrdin
        ELSEIF (input%solMethod == RADI_NUM_DOM8) THEN
          input%nOrdin = 8          
          input%nAng   = input%nOrdin
        ELSEIF (input%solMethod == RADI_NUM_DOM16) THEN
          input%nOrdin = 16          
          input%nAng   = input%nOrdin
        ELSEIF (input%solMethod == RADI_NUM_FVM) THEN
          input%nAng = (input%nPol+1)*(input%nAzi+1)
        ENDIF ! solMethod
      ENDIF   ! radiModel

! --- assign angles read from input file to data structure

      ALLOCATE( input%angles(input%nAng,RADI_ANGLE_NCOMP),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      IF ((input%radiModel == RADI_MODEL_ROSS)  .OR. &
          (input%radiModel == RADI_MODEL_FLDSRC) .OR. &
          (input%radiModel == RADI_MODEL_FLDTRAN)) THEN
        READ(input%line(1),*,err=10,end=20) (input%angles(l,RADI_ANGLE_POLAR), &
                                             l=1,input%nAng )  
        READ(input%line(2),*,err=10,end=20) (input%angles(l,RADI_ANGLE_AZIMU), &
                                             l=1,input%nAng )  
        input%angles = input%angles*global%rad 

      ELSEIF ((input%radiModel == RADI_MODEL_RTEGRAY)  .OR. &
              (input%radiModel == RADI_MODEL_RTEBAND)) THEN
        IF (input%solMethod == RADI_NUM_FVM) THEN
          DO m = 1,input%nPol+1
            DO l = 1,input%nAzi+1
              n = (m-1)*(input%nAzi+1) + l
              input%angles(n,RADI_ANGLE_AZIMU) =   pi*DBLE(l-1)/DBLE(input%nAzi)
              input%angles(n,RADI_ANGLE_POLAR) =twopi*DBLE(m-1)/DBLE(input%nPol)
            ENDDO
          ENDDO
        ELSE  ! DOM

        ENDIF ! solMethod
      ENDIF   ! radiModel

    ENDIF ! region active and my processor
  ENDDO   ! iReg
  
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__, &
                  'Error in reading real numbers from string' )
20   CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__, &
                  'Number of intensity angles is inconsistent' )

999  CONTINUE

END SUBROUTINE RADI_DerivedInputValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_DerivedInputValues.F90,v $
! Revision 1.4  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/09/30 17:10:04  wasistho
! prepared for full FLD radiation model
!
! Revision 1.1  2004/09/22 02:35:49  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.6  2004/09/22 01:31:23  wasistho
! switch LFD to FLD for flux limited diffusion
!
! Revision 1.5  2004/09/18 17:41:18  wasistho
! install Limited Flux Diffusion radiation
!
! Revision 1.4  2003/08/01 22:16:10  wasistho
! prepared rocrad for Genx
!
! Revision 1.3  2003/07/23 03:13:49  wasistho
! cured baby illness
!
! Revision 1.2  2003/07/22 03:05:41  wasistho
! include logical write-parameter
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!
!******************************************************************************







