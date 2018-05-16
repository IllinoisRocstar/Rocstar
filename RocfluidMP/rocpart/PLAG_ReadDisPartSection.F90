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
! Purpose: read in user input related to discrete particle module.
!
! Description: none.
!
! Input: user input file.
!
! Output: regions = total number of particles, drag model,injection model.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_ReadDisPartSection.F90,v 1.7 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_ReadDisPartSection( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt_input
  USE ModPartLag, ONLY    : t_plag_input
  USE ModError
  USE ModParameters

  USE PLAG_ModParameters
  USE PLAG_ModInterfaces, ONLY : PLAG_ReadPdfFromFile

#ifdef RFLO
  USE ModInterfaces, ONLY : ReadRegionSection
#endif
#ifdef RFLU
  USE ModInterfaces, ONLY: ReadSection
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop  variables

! ... local variables
  INTEGER, PARAMETER :: NVALS_MAX=18

  CHARACTER(CHRLEN) :: RCSIdentString
  CHARACTER(20)     :: keys(NVALS_MAX)

  INTEGER           :: brbeg, brend, nVals

  LOGICAL           :: defined(NVALS_MAX)

  REAL(RFREAL)      :: vals(NVALS_MAX)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ReadDisPartSection.F90,v $ $Revision: 1.7 $'

  global => regions(1)%global

  CALL RegisterFunction( global, 'PLAG_ReadDisPartSection',&
  'PLAG_ReadDisPartSection.F90' )

! specify keywords and search for them ----------------------------------------

  nVals   = NVALS_MAX

  defined(:) = .FALSE.

  keys( 1)   = 'USED'
  keys( 2)   = 'NPCLSMAX'
  keys( 3)   = 'EJECMODEL'
  keys( 4)   = 'INJCVELRATIO'
  keys( 5)   = 'SPLOAD'
  keys( 6)   = 'INJCBETA'
  keys( 7)   = 'INJCDIAMDIST'
  keys( 8)   = 'INJCDIAMMEAN'
  keys( 9)   = 'INJCDIAMMIN'
  keys(10)   = 'INJCDIAMMAX'
  keys(11)   = 'INJCSTDDEV'  
  keys(12)   = 'INTRPLMIXTMODEL'
  keys(13)   = 'NPCLSBUFFTOT'
  keys(14)   = 'NPCLSBUFFCECELLSMAX'
  keys(15)   = 'BREAKUPMODEL'
  keys(16)   = 'BREAKUPFAC'
  keys(17)   = 'BREAKUPWEBSWI'
  keys(18)   = 'FINDPCLMETHOD'

#ifdef RFLO
  CALL ReadRegionSection( global, IF_INPUT,nVals,keys,vals, &
       brbeg,brend,defined )
#endif
#ifdef RFLU
  CALL ReadSection(global,IF_INPUT,nVals,keys(1:nVals),vals(1:nVals), &
       defined(1:nVals))

  brbeg = LBOUND(regions,1)
  brend = UBOUND(regions,1)
#endif

  IF (defined(1)) THEN
     IF (vals(1) < 0.1) THEN
        regions(brbeg:brend)%plagInput%readStatus = 0
     ELSE
        regions(brbeg:brend)%plagInput%readStatus = 1
     ENDIF ! vals
  ENDIF ! defined

  IF (defined(2)) &
       regions(brbeg:brend)%plagInput%nPclsMax = MAX(1,INT(ABS(vals(2))+0.5_RFREAL))

  IF (defined(3)) THEN
     regions(brbeg:brend)%plagInput%ejecModel     = NINT(ABS(vals(3)))
     IF ( vals(3) > 0.9 .AND. vals(3) < 1.1 ) &
        regions(brbeg:brend)%plagInput%ejecModel  = PLAG_EJEC_MODEL1
     IF ( vals(3) > 1.9 .AND. vals(3) < 2.1 ) &
          regions(brbeg:brend)%plagInput%ejecModel  = PLAG_EJEC_CRE
  ENDIF ! defined

  IF (defined(4)) &
       regions(brbeg:brend)%plagInput%injcVelRatio = ABS(vals(4))

  IF (defined(5)) &
       regions(brbeg:brend)%plagInput%spLoad       = ABS(vals(5))

  IF (defined(6)) &
       regions(brbeg:brend)%plagInput%injcBeta     = ABS(vals(6))

  IF (defined(7)) THEN
     regions(brbeg:brend)%plagInput%injcDiamDist    = NINT(ABS(vals(7)))
     IF ( vals(7) > 0.9 .AND. vals(7) < 1.1 ) &
       regions(brbeg:brend)%plagInput%injcDiamDist  = PLAG_INJC_LOGNORM   
     IF ( vals(7) > 1.9 .AND. vals(7) < 2.1 ) &
       regions(brbeg:brend)%plagInput%injcDiamDist  = PLAG_INJC_LOGSKWD
     IF ( vals(7) > 2.9 .AND. vals(7) < 3.1 ) THEN
       regions(brbeg:brend)%plagInput%injcDiamDist  = PLAG_INJC_PDF
       CALL PLAG_ReadPdfFromFile( regions,brbeg,brend )
     ENDIF

  ENDIF ! defined

  IF (defined(8)) &
    regions(brbeg:brend)%plagInput%injcDiamMean = ABS(vals(8))

  IF (defined(9)) &
    regions(brbeg:brend)%plagInput%injcDiamMin  = ABS(vals(9))

  IF (defined(10)) &
    regions(brbeg:brend)%plagInput%injcDiamMax  = ABS(vals(10))

  IF (defined(11)) &
    regions(brbeg:brend)%plagInput%injcStdDev   = ABS(vals(11))

  IF (defined(12)) THEN
    regions(brbeg:brend)%plagInput%intrplMixtModel    = NINT(ABS(vals(12)))
    IF ( vals(12) > -0.1 .AND. vals(12) < 0.1 ) &
      regions(brbeg:brend)%plagInput%intrplMixtModel  = ZEROTH_ORDER
    IF ( vals(12) > 0.9 .AND. vals(12) < 1.1 ) &
      regions(brbeg:brend)%plagInput%intrplMixtModel  = FIRST_ORDER
  ENDIF ! defined

  IF (defined(13)) &
    regions(brbeg:brend)%plagInput%nPclsBuffTot = MAX(1,INT(ABS(vals(13))+0.5_RFREAL))

  IF (defined(14)) &
    regions(brbeg:brend)%plagInput%nPclsBuffCECellsMax = MAX(1,INT(ABS(vals(14))+0.5_RFREAL))  

  IF (defined(15)) THEN
    regions(brbeg:brend)%plagInput%breakupModel    = NINT(ABS(vals(15))) 
    IF ( vals(15) > -0.1 .AND. vals(15) < 0.1 ) &
      regions(brbeg:brend)%plagInput%breakupModel  = PLAG_BREAKUP_NOMODEL 
    IF ( vals(15) > 0.9 .AND. vals(15) < 1.1 ) &
      regions(brbeg:brend)%plagInput%breakupModel  = PLAG_BREAKUP_MODEL1
  ENDIF ! defined

  IF (defined(16)) &
    regions(brbeg:brend)%plagInput%breakupFac = ABS(vals(16))
              
  IF (defined(17)) THEN
    regions(brbeg:brend)%plagInput%breakupWebSwi   = NINT(ABS(vals(17))) 
    IF ( vals(17) > -0.1 .AND. vals(17) < 0.1 ) &
      regions(brbeg:brend)%plagInput%breakupWebSwi = PLAG_BREAKUP_NOWEBSWI 
    IF ( vals(17) > 0.9 .AND. vals(17) < 1.1 ) &
      regions(brbeg:brend)%plagInput%breakupWebSwi = PLAG_BREAKUP_WEBSWI1
  ENDIF ! defined

  IF ( defined(18) .EQV. .TRUE. ) THEN 
    regions(brbeg:brend)%plagInput%findPclMethod = vals(18)  
  END IF ! defined

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_ReadDisPartSection
!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ReadDisPartSection.F90,v $
! Revision 1.7  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2007/03/06 23:15:32  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.4  2005/05/17 18:35:41  fnajjar
! Further IO improvement
!
! Revision 1.3  2005/05/17 17:07:37  fnajjar
! Improved loading of values for proper error trapping
!
! Revision 1.2  2005/04/25 18:39:10  luca1
! Imposed PDF from file option for random particle ejection
!
! Revision 1.1  2004/12/01 20:58:05  fnajjar
! Initial revision after changing case
!
! Revision 1.13  2004/10/08 22:13:33  haselbac
! Added reading of findPclMethod
!
! Revision 1.12  2004/10/08 20:03:14  fnajjar
! Bug fix for nVals size and defined NVALS_MAX
!
! Revision 1.11  2004/07/29 16:58:07  fnajjar
! Inconsistent numbers in array vals leading to io clobbering
!
! Revision 1.10  2004/06/17 15:19:48  fnajjar
! Added infrastructure for ejection model and revamped DISPART input section
!
! Revision 1.9  2004/06/16 23:07:17  fnajjar
! Renamed variabled for CRE kernel
!
! Revision 1.8  2004/03/10 23:09:50  fnajjar
! Added maximum buffer size for corner-edge cells
!
! Revision 1.7  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2004/02/26 21:02:19  haselbac
! Added RFLU support
!
! Revision 1.5  2003/09/17 21:05:13  fnajjar
! Added infrastructure for skewed Log distribution in injection model
!
! Revision 1.4  2003/09/13 20:14:22  fnajjar
! Added infrastructure for Breakup model
!
! Revision 1.3  2003/09/10 23:35:50  fnajjar
! Removed flags that are subsumed with Rocinteract
!
! Revision 1.2  2003/01/24 19:41:55  f-najjar
! Bug fix for nPclsBuffTot to read vals(14)
!
! Revision 1.1  2002/10/25 14:19:16  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************







