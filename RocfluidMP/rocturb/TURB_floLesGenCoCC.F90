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
! Purpose: Get filter coefficients of cell-cell filtering.
!
! Description: The coefficients are obtained in i, j and k direction.
!              For each direction, distinction is made between the low delta
!              (1 or 2 grid-spacing) and high delta (2 or 4 grid-spacing).
!
! Input: region  = data of current region
!
! Output: turb%ccCofi1,2,4, turb%ccCofj1,2,4, turb%ccCofk1,2,4. 
!
! Notes: This routine is only relevant if non-uniform filter is selected.
!
!******************************************************************************
!
! $Id: TURB_floLesGenCoCC.F90,v 1.5 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloLesGenCoCC( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy,RFLO_GetDimensDummyNodes, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset
  USE TURB_ModInterfaces, ONLY : TURB_FloLesGenCoCCLo, TURB_FloLesGenCoCCHi
  USE ModTurbulence, ONLY : t_turb
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), TARGET :: region

! ... loop variables

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER           :: ibeg,iend,jbeg,jend,kbeg,kend
  INTEGER           :: idcbeg,idcend,jdcbeg,jdcend,kdcbeg,kdcend
  INTEGER           :: idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend
  INTEGER           :: iLev,iCOff,ijCOff,iNOff,ijNOff,minIdx,maxIdx,segId
  INTEGER           :: homDir(DIRI:DIRK),filterWidth(DIRI:DIRK)
  TYPE(t_turb), POINTER     :: turb
  REAL(RFREAL), POINTER     :: segm(:,:),ccCofA(:,:),ccCofB(:,:)
  REAL(RFREAL), ALLOCATABLE :: ds(:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floLesGenCoCC.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloLesGenCoCC',&
  'TURB_floLesGenCoCC.F90' )

! get indices, parameters and pointers ----------------------------------------

  iLev           =  region%currLevel
  homDir(:)      =  region%turbInput%homDir(:)
  filterWidth(:) =  region%turbInput%filterWidth(:)
  turb           => region%levels(iLev)%turb
  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  minIdx = MIN(idnbeg,jdnbeg,kdnbeg)-2   
  maxIdx = MAX(idnend,jdnend,kdnend)+2   
  ALLOCATE( ds(minIdx:maxIdx) )

  ibeg = idcbeg
  iend = idcend
  jbeg = jdcbeg
  jend = jdcend
  kbeg = kdcbeg
  kend = kdcend

! filter coefficients in I direction

  IF (homDir(DIRI) == OFF) THEN
    IF (filterWidth(DIRI)==FILWIDTH_ONE) THEN

! --- filter coefficients at cell centers
      segm   => region%levels(iLev)%turb%workK
      ccCofA => turb%ccCofi1
      ccCofB => turb%ccCofi2
      segId  =  1
      CALL TURB_FloLesGenCoCCLo( global,DIRI,ibeg,iend,jbeg,jend,kbeg,kend, &
                              minIdx,maxIdx,segId,iCOff,ijCOff,iNOff,ijNOff,ds, &
                              segm,ccCofA,ccCofB )  

    ELSEIF ((filterWidth(DIRI)==FILWIDTH_TWO) .OR. &
            (filterWidth(DIRI)==FILWIDTH_ZERO)) THEN

! --- filter coefficients at cell centers
      segm   => region%levels(iLev)%turb%workK
      ccCofA => turb%ccCofi2
      ccCofB => turb%ccCofi4
      segId  =  1
      CALL TURB_FloLesGenCoCCHi( global,DIRI,ibeg,iend,jbeg,jend,kbeg,kend, &
                              minIdx,maxIdx,segId,iCOff,ijCOff,iNOff,ijNOff,ds, &
                              segm,ccCofA,ccCofB )  

    ENDIF
  ENDIF

! filter coefficients in J direction

  IF (homDir(DIRJ) == OFF) THEN
    IF (filterWidth(DIRJ)==FILWIDTH_ONE) THEN

! --- filter coefficients at cell centers
      segm   => region%levels(iLev)%turb%workI
      ccCofA => turb%ccCofj1
      ccCofB => turb%ccCofj2
      segId  =  1
      CALL TURB_FloLesGenCoCCLo( global,DIRJ,ibeg,iend,jbeg,jend,kbeg,kend, &
                              minIdx,maxIdx,segId,iCOff,ijCOff,iNOff,ijNOff,ds, &
                              segm,ccCofA,ccCofB )  

    ELSEIF ((filterWidth(DIRJ)==FILWIDTH_TWO) .OR. &
            (filterWidth(DIRJ)==FILWIDTH_ZERO)) THEN

! --- filter coefficients at cell centers
      segm   => region%levels(iLev)%turb%workI
      ccCofA => turb%ccCofj2
      ccCofB => turb%ccCofj4
      segId  =  1
      CALL TURB_FloLesGenCoCCHi( global,DIRJ,ibeg,iend,jbeg,jend,kbeg,kend, &
                              minIdx,maxIdx,segId,iCOff,ijCOff,iNOff,ijNOff,ds, &
                              segm,ccCofA,ccCofB )  

    ENDIF
  ENDIF

! filter coefficients in K direction

  IF (homDir(DIRK) == OFF) THEN
    IF (filterWidth(DIRK)==FILWIDTH_ONE) THEN

! --- filter coefficients at cell centers
      segm   => region%levels(iLev)%turb%workJ
      ccCofA => turb%ccCofk1
      ccCofB => turb%ccCofk2
      segId  =  1
      CALL TURB_FloLesGenCoCCLo( global,DIRK,ibeg,iend,jbeg,jend,kbeg,kend, &
                              minIdx,maxIdx,segId,iCOff,ijCOff,iNOff,ijNOff,ds, &
                              segm,ccCofA,ccCofB )  

    ELSEIF ((filterWidth(DIRK)==FILWIDTH_TWO) .OR. &
            (filterWidth(DIRK)==FILWIDTH_ZERO)) THEN

! --- filter coefficients at cell centers
      segm   => region%levels(iLev)%turb%workJ
      ccCofA => turb%ccCofk2
      ccCofB => turb%ccCofk4
      segId  =  1
      CALL TURB_FloLesGenCoCCHi( global,DIRK,ibeg,iend,jbeg,jend,kbeg,kend, &
                              minIdx,maxIdx,segId,iCOff,ijCOff,iNOff,ijNOff,ds, &
                              segm,ccCofA,ccCofB )  

    ENDIF
  ENDIF

! deallocate temporary arrays

  DEALLOCATE( ds )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloLesGenCoCC

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floLesGenCoCC.F90,v $
! Revision 1.5  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/08/04 02:50:34  wasistho
! removed turb%avgCoI,J,K replaced with turb%workI,J,K
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.2  2003/08/29 01:41:24  wasistho
! Added TARGET attribute to region variable, since pointers are cached into it
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************







