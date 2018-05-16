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
! Purpose: Read in user input within RADI section (done on all processors).
!
! Description: none.
!
! Input: regions = user input file of all regions.
!
! Output: regions = RADI input parameters.
!
! Notes: Mother routine = ReadInputFile.
!
!******************************************************************************
!
! $Id: RADI_ReadRadiSection.F90,v 1.4 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

#ifdef RFLO
SUBROUTINE RADI_ReadRadiSection( regions )
#endif
#ifdef RFLU
SUBROUTINE RADI_ReadRadiSection
#endif

  USE ModDataTypes
#ifdef RFLO  
  USE ModDataStruct, ONLY : t_region
#endif
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO 
  USE ModInterfaces, ONLY : ReadRegionSection, ReadStringSection
#endif
#ifdef RFLU
  USE ModInterfaces, ONLY : ReadSection, ReadStringSection
#endif
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  
  IMPLICIT NONE

! ... parameters
#ifdef RFLO
  TYPE(t_region), POINTER :: regions(:)
#endif

! ... loop variables
  INTEGER :: iReg

! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: nVals
#ifdef RFLO
  INTEGER :: brbeg, brend
#endif
  INTEGER, PARAMETER :: NVALS_MAX = 20

  REAL(RFREAL)      :: vals(NVALS_MAX)
  LOGICAL           :: defined(NVALS_MAX)
  CHARACTER(20)     :: keys(NVALS_MAX)
  CHARACTER(256)    :: sectionline, line(2)

!******************************************************************************

  global => regions(1)%global
  CALL RegisterFunction( global,'RADI_ReadRadiSection',&
  'RADI_ReadRadiSection.F90' )

! specify keywords

  nVals = NVALS_MAX
  keys( 1) = 'RADIMODEL'
  keys( 2) = 'MEDIA'
  keys( 3) = 'FLUXLIMITER'
  keys( 4) = 'SMOOCF'
  keys( 5) = 'DISCR'
  keys( 6) = 'K2'
  keys( 7) = '1/K4'
  keys( 8) = 'ORDER'
  keys( 9) = 'CONPARTVFRAC'
  keys(10) = 'CONPARTDIAM'
  keys(11) = 'CONPARTQE'
  keys(12) = 'DISPARTVFRAC'
  keys(13) = 'DISPARTDIAM'
  keys(14) = 'DISPARTQE'
  keys(15) = 'SOLMETHOD'
  keys(16) = 'NPOLAR'
  keys(17) = 'NAZIMUTHAL'
  keys(18) = 'NINTANGLES'
  keys(19) = 'POLANGLES'
  keys(20) = 'AZIANGLES'

! 4  -  8    only relevant for FLD involving transport eq. for Er
! 9,10,12,13 only relevant for artificial media
! 15         only relevant for RTE radiation model
! 16 - 17    only relevant for RTE radiation model/ FVM sol. method
! 18         only relevant for diffusion approximation methods
! 19 - 20    only relevant for diffusion approximation methods

! safety check

  IF (nvals-1 /= 14) THEN
    CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
           'number of input parameters being read is inconsistent' )
  ENDIF

! search for keywords

#ifdef RFLO
  CALL ReadRegionSection( global,IF_INPUT,nVals-2,keys(1:nVals-2), &
                          vals(1:nVals-2),brbeg,brend,defined(1:nVals-2) )

  IF (defined(1)) THEN
    regions(brbeg:brend)%radiInput%radiModel = INT(vals(1)+0.5_RFREAL)
  ENDIF
  IF (defined(2)) THEN
    IF (INT(vals(2)+0.5_RFREAL) <= 1) THEN
      regions(brbeg:brend)%radiInput%media = RADI_MEDIA_ARTIF
    ELSE
      regions(brbeg:brend)%radiInput%media = RADI_MEDIA_REAL
    ENDIF
  ENDIF
  IF (defined(3)) THEN
    IF (INT(vals(3)+0.5_RFREAL) <= 0) THEN
      regions(brbeg:brend)%radiInput%fluxLim = FLD_LIM_NONE
    ELSE
      regions(brbeg:brend)%radiInput%fluxLim = FLD_LIM_LP
    ENDIF
  ENDIF

  IF (defined(4)) regions(brbeg:brend)%radiInput%smoocf = vals(4)

  IF (defined(5)) THEN
    regions(brbeg:brend)%radiInput%spaceDiscr  = INT(vals(5)+0.5_RFREAL)
  ENDIF

  IF (defined(6)) regions(brbeg:brend)%radiInput%vis2   = ABS(vals(6))

  IF (defined(7)) THEN
    IF (vals(7) > 1.E-10_RFREAL) THEN
      regions(brbeg:brend)%radiInput%vis4 = 1._RFREAL/vals(7)
    ELSEIF (vals(7) > 0._RFREAL .AND. vals(7) <= 1.E-10_RFREAL) THEN
      regions(brbeg:brend)%radiInput%vis4 = 1.E+10_RFREAL
    ELSEIF (vals(7) <= 0._RFREAL ) THEN
      regions(brbeg:brend)%radiInput%vis4 = 0.0_RFREAL
    ENDIF
  ENDIF

  IF (defined(8)) THEN
    regions(brbeg:brend)%radiInput%spaceOrder  = INT(vals(8)+0.5_RFREAL)
  ENDIF

  IF (defined(9)) THEN
    DO iReg = brbeg,brend
      regions(iReg)%radiInput%optConst(PHASE_PROP_V,RADI_PHASE_CONPART)= &
                                       ABS(vals(9))
    ENDDO
  ENDIF
  IF (defined(10)) THEN
    DO iReg = brbeg,brend
      regions(iReg)%radiInput%optConst(PHASE_PROP_D,RADI_PHASE_CONPART)= &
                                       ABS(vals(10))
    ENDDO
  ENDIF
  IF (defined(11)) THEN
    DO iReg = brbeg,brend
      regions(iReg)%radiInput%optConst(PHASE_PROP_Q,RADI_PHASE_CONPART)= &
                                       ABS(vals(11))
    ENDDO
  ENDIF
  IF (defined(12)) THEN
    DO iReg = brbeg,brend
      regions(iReg)%radiInput%optConst(PHASE_PROP_V,RADI_PHASE_DISPART)= &
                                       ABS(vals(12))
    ENDDO
  ENDIF
  IF (defined(13)) THEN
    DO iReg = brbeg,brend
      regions(iReg)%radiInput%optConst(PHASE_PROP_D,RADI_PHASE_DISPART)= &
                                       ABS(vals(13))
    ENDDO
  ENDIF
  IF (defined(14)) THEN
    DO iReg = brbeg,brend
      regions(iReg)%radiInput%optConst(PHASE_PROP_Q,RADI_PHASE_DISPART)= &
                                       ABS(vals(14))
    ENDDO
  ENDIF

  IF (defined(15)) THEN
    regions(brbeg:brend)%radiInput%solMethod = INT(vals(15)+0.5_RFREAL)
  ENDIF
  IF (defined(16)) THEN
    regions(brbeg:brend)%radiInput%nPol      = INT(vals(16)+0.5_RFREAL)
  ENDIF
  IF (defined(17)) THEN
    regions(brbeg:brend)%radiInput%nAzi      = INT(vals(17)+0.5_RFREAL)
  ENDIF
  IF (defined(18)) THEN
    regions(brbeg:brend)%radiInput%nAng      = INT(vals(18)+0.5_RFREAL)
  ENDIF

  REWIND(IF_INPUT, err=10)
  DO
    READ(IF_INPUT,'(A256)',err=10,end=79) sectionline

    SELECT CASE(TRIM(sectionline))
    CASE ('# RADIATION')
      CALL ReadStringSection( global,IF_INPUT,2,keys(nvals-1:nvals), &
                              line(1:2),defined(nvals-1:nvals) )
    END SELECT
  ENDDO

79 CONTINUE

  IF (defined(19) .AND. defined(nvals-1)) THEN
    regions(brbeg:brend)%radiInput%line(1)  = line(1)
  ELSE
    regions(brbeg:brend)%radiInput%nAng     = 1
    regions(brbeg:brend)%radiInput%line(1)  = '45'
  ENDIF
  IF (defined(20) .AND. defined(nvals)) THEN
    regions(brbeg:brend)%radiInput%line(2)  = line(2)
  ELSE
    regions(brbeg:brend)%radiInput%nAng     = 1
    regions(brbeg:brend)%radiInput%line(2)  = '45'
  ENDIF
#endif

#ifdef RFLU
  CALL ReadSection( global,IF_INPUT,nVals-2,keys(1:nVals-2),vals(1:nVals-2), & 
                    defined(1:nVals-2) ) 
  
  IF (defined(1)) THEN 
    radiInput%radiModel = NINT(vals(1))
  END IF
  IF (defined(2)) THEN
    IF (NINT(vals(2)) <= 1) THEN
      radiInput%media = RADI_MEDIA_ARTIF
    ELSE
      radiInput%media = RADI_MEDIA_REAL
    ENDIF
  ENDIF
  IF (defined(3)) THEN
    IF (NINT(vals(3)) <= 0) THEN
      radiInput%fluxLim = FLD_LIM_NONE
    ELSE
      radiInput%fluxLim = FLD_LIM_LP
    ENDIF
  ENDIF

  IF (defined(4)) THEN
    radiInput%smoocf = vals(4)
  ENDIF

  IF (defined(5)) THEN
    radiInput%spaceDiscr = NINT(vals(5))
  ENDIF

  IF (defined(6)) THEN
    radiInput%vis2 = ABS(vals(6))
  ENDIF

  IF (defined(7)) THEN
    IF (vals(7) > 1.E-10_RFREAL) THEN
      radiInput%vis4 = 1._RFREAL/vals(7)
    ELSEIF (vals(7) > 0._RFREAL .AND. vals(7) <= 1.E-10_RFREAL) THEN
      radiInput%vis4 = 1.E+10_RFREAL
    ELSEIF (vals(7) <= 0._RFREAL ) THEN
      radiInput%vis4 = 0.0_RFREAL
    ENDIF
  ENDIF

  IF (defined(8)) THEN
    radiInput%spaceOrder  = NINT(vals(8))
  ENDIF

  IF (defined(9)) &
    radiInput%optConst(PHASE_PROP_V,RADI_PHASE_CONPART)= ABS(vals(9))
  IF (defined(10)) &
    radiInput%optConst(PHASE_PROP_D,RADI_PHASE_CONPART)= ABS(vals(10))
  IF (defined(11)) &
    radiInput%optConst(PHASE_PROP_Q,RADI_PHASE_CONPART)= ABS(vals(11))
  IF (defined(12)) &
    radiInput%optConst(PHASE_PROP_V,RADI_PHASE_DISPART)= ABS(vals(12))
  IF (defined(13)) &
    radiInput%optConst(PHASE_PROP_D,RADI_PHASE_DISPART)= ABS(vals(13))
  IF (defined(14)) &
    radiInput%optConst(PHASE_PROP_Q,RADI_PHASE_DISPART)= ABS(vals(14))

  IF (defined(15)) THEN 
    radiInput%solMethod = NINT(vals(15))
  END IF

  IF (defined(16)) THEN 
    radiInput%nPol      = NINT(vals(16))
  END IF

  IF (defined(17)) THEN 
    radiInput%nAzi      = NINT(vals(17))
  END IF

  IF (defined(18)) THEN 
    radiInput%nAng      = NINT(vals(18))
  END IF

  REWIND(IF_INPUT, err=10)
  DO
    READ(IF_INPUT,'(A256)',err=10,end=89) sectionline

    SELECT CASE(TRIM(sectionline))
    CASE ('# RADIATION')
      CALL ReadStringSection( global,IF_INPUT,2,keys(nvals-1:nvals), &
                              line(1:2),defined(nvals-1:nvals) )
    END SELECT
  ENDDO

89 CONTINUE

  IF (defined(19) .AND. defined(nvals-1)) THEN
    radiInput%line(1)   = line(1)
  ELSE
    radiInput%nAng      = 1
    radiInput%line(1)   = '45'
  ENDIF
  IF (defined(20) .AND. defined(nvals)) THEN
    radiInput%line(2)   = line(2)
  ELSE
    radiInput%nAng      = 1
    radiInput%line(2)   = '45'
  ENDIF
#endif 

! finalize -------------------------------------------------------

  GOTO 999

! error handling

10  CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__, &
       'reading Radiation section in Input File' )

999 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_ReadRadiSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_ReadRadiSection.F90,v $
! Revision 1.4  2008/12/06 08:44:38  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/09/30 17:10:42  wasistho
! prepared for full FLD radiation model
!
! Revision 1.1  2004/09/22 02:35:50  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.6  2004/09/22 01:30:58  wasistho
! switch LFD to FLD for flux limited diffusion
!
! Revision 1.5  2004/09/18 17:41:04  wasistho
! install Limited Flux Diffusion radiation
!
! Revision 1.4  2003/07/30 22:23:55  wasistho
! enter part and smoke data into radiation
!
! Revision 1.3  2003/07/23 03:13:36  wasistho
! cured baby illness
!
! Revision 1.2  2003/07/18 01:39:31  wasistho
! removed bcModel from input data structure
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!
!******************************************************************************







