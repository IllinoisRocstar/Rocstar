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
! Purpose: Read in user input within TURB section (done on all processors).
!
! Description: none.
!
! Input: regions = user input file of all regions.
!
! Output: regions = turbulence input parameters.
!
! Notes: Mother routine = ReadInputFile.
!
!******************************************************************************
!
! $Id: TURB_ReadTurbSection.F90,v 1.7 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_ReadTurbSection( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO 
  USE ModInterfaces, ONLY : ReadRegionSection
#endif
#ifdef RFLU
  USE ModInterfaces, ONLY : ReadSection
#endif
  USE ModTurbulence
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global


  INTEGER :: nVals
  INTEGER :: brbeg, brend
  INTEGER, PARAMETER :: NVALS_MAX = 25

  LOGICAL           :: defined(NVALS_MAX)
  REAL(RFREAL)      :: vals(NVALS_MAX)
  CHARACTER(20)     :: keys(NVALS_MAX)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_ReadTurbSection.F90,v $'

  global => regions(1)%global
  CALL RegisterFunction( global,'TURB_ReadTurbSection',&
  'TURB_ReadTurbSection.F90' )

! specify keywords and search for them

  nVals = NVALS_MAX

  defined(:) = .FALSE.

  keys( 1) = 'TURBMODEL'
  keys( 2) = 'OUTPUTNUMBER'
  keys( 3) = 'CSMAGORINSKY'
  keys( 4) = 'XTRIPLOC'
  keys( 5) = 'YTRIPLOC'
  keys( 6) = 'ZTRIPLOC'
  keys( 7) = 'FILTERTYPE'
  keys( 8) = 'DELTATYPE'
  keys( 9) = 'IFILTERWIDTH'
  keys(10) = 'JFILTERWIDTH'
  keys(11) = 'KFILTERWIDTH'
  keys(12) = 'IHOMOGENDIR'
  keys(13) = 'JHOMOGENDIR'
  keys(14) = 'KHOMOGENDIR'
  keys(15) = 'ENERGYMODEL'
  keys(16) = 'CALCVORTIC'
  keys(17) = 'WALLDISTMETHOD'
  keys(18) = 'WALLDISTFREQ'
  keys(19) = 'VISCFUNCTION'
  keys(20) = 'CDES'
  keys(21) = 'SMOOCF'
  keys(22) = 'DISCR'
  keys(23) = 'K2'
  keys(24) = '1/K4'
  keys(25) = 'ORDER'

#ifdef RFLO
  CALL ReadRegionSection( global,IF_INPUT,nVals,keys(1:nVals),vals(1:nVals), & 
                          brbeg,brend,defined(1:nVals) )
#endif
#ifdef RFLU  
  CALL ReadSection( global,IF_INPUT,nVals,keys(1:nVals),vals(1:nVals), & 
                    defined(1:nVals) ) 

  brbeg = LBOUND(regions,1)   ! temporary for now before zonal modeling apply
  brend = UBOUND(regions,1)
#endif
  IF (defined(1) .eqv. .true.) THEN
    regions(brbeg:brend)%mixtInput%turbModel  = INT(vals(1)+0.5_RFREAL)
  ENDIF

  IF (defined(2)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%nOutField  = INT(vals(2)+0.5_RFREAL)
  ENDIF

  IF (defined(3)  .eqv. .true.) regions(brbeg:brend)%turbInput%cSmag  = ABS(vals(3))

  IF (defined(4)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%xyzSmag(XCOORD) = ABS(vals(4))
  ENDIF

  IF (defined(5)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%xyzSmag(YCOORD) = ABS(vals(5))
  ENDIF

  IF (defined(6)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%xyzSmag(ZCOORD) = ABS(vals(6))
  ENDIF

#ifdef RFLO
  IF (defined(7)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%filterType = INT(vals(7)+0.5_RFREAL)
  ENDIF
  IF (defined(8)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%deltaType  = INT(vals(8)+0.5_RFREAL)
  ENDIF
#endif
  IF (defined(9)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%filterWidth(DIRI) = INT(vals(9)+0.5_RFREAL)
  ENDIF
#ifdef RFLO
  IF (defined(10)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%filterWidth(DIRJ) = INT(vals(10)+0.5_RFREAL)
  ENDIF
  IF (defined(11)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%filterWidth(DIRK) = INT(vals(11)+0.5_RFREAL)
  ENDIF
  IF (defined(12)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%homDir(DIRI) = INT(vals(12)+0.5_RFREAL)
  ENDIF
  IF (defined(13)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%homDir(DIRJ) = INT(vals(13)+0.5_RFREAL)
  ENDIF
  IF (defined(14)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%homDir(DIRK) = INT(vals(14)+0.5_RFREAL)
  ENDIF
#endif
  IF (defined(15)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%engModel     = INT(vals(15)+0.5_RFREAL)
  ENDIF
  IF (defined(16)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%calcVort     = INT(vals(16)+0.5_RFREAL)
  ENDIF
  IF (defined(17)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%wDistMethod  = INT(vals(17)+0.5_RFREAL)
  ENDIF
  IF (defined(18)  .eqv. .true.) THEN
    global%turbCalcWDistFreq = MAX( global%turbCalcWDistFreq, &
                                                  INT(vals(18)+0.5_RFREAL) )
  ENDIF
  IF (defined(19)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%functV1      = INT(vals(19)+0.5_RFREAL)
  ENDIF

  IF (defined(20)  .eqv. .true.) regions(brbeg:brend)%turbInput%cDes = ABS(vals(20))
  IF (defined(21)  .eqv. .true.) regions(brbeg:brend)%turbInput%smoocf = vals(21)

  IF (defined(22)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%spaceDiscr  = INT(vals(22)+0.5_RFREAL)
  ENDIF

  IF (defined(23)  .eqv. .true.) regions(brbeg:brend)%turbInput%vis2   = ABS(vals(23))

  IF (defined(24)  .eqv. .true.) THEN
    IF (vals(24) > 1.E-10_RFREAL) THEN
      regions(brbeg:brend)%turbInput%vis4 = 1._RFREAL/vals(24)
    ELSEIF (vals(24) > 0._RFREAL .AND. vals(24) <= 1.E-10_RFREAL) THEN
      regions(brbeg:brend)%turbInput%vis4 = 1.E+10_RFREAL
    ELSEIF (vals(24) <= 0._RFREAL ) THEN
      regions(brbeg:brend)%turbInput%vis4 = 0.0_RFREAL
    ENDIF
  ENDIF
  IF (defined(25)  .eqv. .true.) THEN
    regions(brbeg:brend)%turbInput%spaceOrder  = INT(vals(25)+0.5_RFREAL)
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_ReadTurbSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_ReadTurbSection.F90,v $
! Revision 1.7  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2008/10/23 18:20:57  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.4  2006/01/12 09:48:15  wasistho
! enabled tripping fixed Smagorinsky
!
! Revision 1.3  2004/04/20 20:49:57  wasistho
! added user option for frequency in computing wall distance
!
! Revision 1.2  2004/03/19 02:48:11  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.11  2004/02/19 04:03:34  wasistho
! added new rans/SA parameter VISCFUNCTION
!
! Revision 1.10  2004/02/11 03:23:50  wasistho
! added feature: variable number of turbulence output fields
!
! Revision 1.9  2003/10/27 23:12:12  wasistho
! bug fixed in reading vis2
!
! Revision 1.8  2003/10/26 00:08:26  wasistho
! added multiple discr.types and order
!
! Revision 1.7  2003/10/15 03:40:59  wasistho
! added 2nd order dissipation coeff. k2
!
! Revision 1.6  2003/10/09 20:49:38  wasistho
! added DES lengthscale coefficient CDES
!
! Revision 1.5  2003/10/07 02:07:06  wasistho
! initial installation of RaNS-SA and DES
!
! Revision 1.4  2003/08/06 15:55:50  wasistho
! added vorticities computation
!
! Revision 1.3  2003/08/01 22:17:52  wasistho
! prepared rocturb for Genx
!
! Revision 1.2  2003/07/22 02:59:55  wasistho
! prepare more accurate rocturb restart
!
! Revision 1.1  2002/10/14 23:55:30  wasistho
! Install Rocturb
!
!
!******************************************************************************







