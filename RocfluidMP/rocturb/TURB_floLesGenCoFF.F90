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
! Purpose: Get filter coefficients of face-face filtering.
!
! Description: The coefficients are obtained in i, j and k direction.
!              For each direction, distinction is made between the low delta
!              (1 or 2 grid-spacing) and high delta (2 or 4 grid-spacing).
!
! Input: region  = data of current region
!
! Output: turb%ffCofi1,2,4, turb%ffCofj1,2,4, turb%ffCofk1,2,4.
!
! Notes: This routine is only relevant if non-uniform filter is selected.
!
!******************************************************************************
!
! $Id: TURB_floLesGenCoFF.F90,v 1.5 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloLesGenCoFF( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset
  USE TURB_ModInterfaces, ONLY : TURB_FloLesGenCoFFLo, TURB_FloLesGenCoFCLo, &
                                 TURB_FloLesGenCoFFHi, TURB_FloLesGenCoFCHi
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
  INTEGER           :: idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend
  INTEGER           :: iLev,iNOff,ijNOff,minIdx,maxIdx,segId
  INTEGER           :: homDir(DIRI:DIRK),filterWidth(DIRI:DIRK)
  TYPE(t_turb), POINTER     :: turb
  REAL(RFREAL), POINTER     :: segm(:,:),ffCofA(:,:),ffCofB(:,:)
  REAL(RFREAL), ALLOCATABLE :: ds(:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floLesGenCoFF.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloLesGenCoFF',&
  'TURB_floLesGenCoFF.F90' )

! get indices, parameters and pointers ----------------------------------------

  iLev           =  region%currLevel
  homDir(:)      =  region%turbInput%homDir(:)
  filterWidth(:) =  region%turbInput%filterWidth(:)
  turb           => region%levels(iLev)%turb
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  minIdx = MIN(idnbeg,jdnbeg,kdnbeg)-2   
  maxIdx = MAX(idnend,jdnend,kdnend)+2   
  ALLOCATE( ds(minIdx:maxIdx) )

  ibeg = idnbeg
  iend = idnend
  jbeg = jdnbeg
  jend = jdnend
  kbeg = kdnbeg
  kend = kdnend

! filter coefficients in I direction

  IF (homDir(DIRI) == OFF) THEN
    IF (filterWidth(DIRI)==FILWIDTH_ONE) THEN

! --- filter coefficients at I-faces
      segm   => region%levels(iLev)%turb%workK
      ffCofA => turb%ffCofi1I
      ffCofB => turb%ffCofi2I
      segId  =  1
      CALL TURB_FloLesGenCoFFLo( global,DIRI,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

! --- filter coefficients at J-faces
      segm   => region%levels(iLev)%turb%workJ
      ffCofA => turb%ffCofi1J
      ffCofB => turb%ffCofi2J
      segId  =  2
      CALL TURB_FloLesGenCoFCLo( global,DIRI,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

! --- filter coefficients at K-faces
      segm   => region%levels(iLev)%turb%workK
      ffCofA => turb%ffCofi1K
      ffCofB => turb%ffCofi2K
      segId  =  1
      CALL TURB_FloLesGenCoFCLo( global,DIRI,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

    ELSEIF ((filterWidth(DIRI)==FILWIDTH_TWO) .OR. &
            (filterWidth(DIRI)==FILWIDTH_ZERO)) THEN

! --- filter coefficients at I-faces
      segm   => region%levels(iLev)%turb%workK
      ffCofA => turb%ffCofi2I
      ffCofB => turb%ffCofi4I
      segId  =  1
      CALL TURB_FloLesGenCoFFHi( global,DIRI,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

! --- filter coefficients at J-faces
      segm   => region%levels(iLev)%turb%workJ
      ffCofA => turb%ffCofi2J
      ffCofB => turb%ffCofi4J
      segId  =  2
      CALL TURB_FloLesGenCoFCHi( global,DIRI,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

! --- filter coefficients at K-faces
      segm   => region%levels(iLev)%turb%workK
      ffCofA => turb%ffCofi2K
      ffCofB => turb%ffCofi4K
      segId  =  1
      CALL TURB_FloLesGenCoFCHi( global,DIRI,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

    ENDIF
  ENDIF

! filter coefficients in J direction

  IF (homDir(DIRJ) == OFF) THEN
    IF (filterWidth(DIRJ)==FILWIDTH_ONE) THEN

! --- filter coefficients at J-faces
      segm   => region%levels(iLev)%turb%workI
      ffCofA => turb%ffCofj1J
      ffCofB => turb%ffCofj2J
      segId  =  1
      CALL TURB_FloLesGenCoFFLo( global,DIRJ,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

! --- filter coefficients at I-faces
      segm   => region%levels(iLev)%turb%workI
      ffCofA => turb%ffCofj1I
      ffCofB => turb%ffCofj2I
      segId  =  1
      CALL TURB_FloLesGenCoFCLo( global,DIRJ,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

! --- filter coefficients at K-faces
      segm   => region%levels(iLev)%turb%workK
      ffCofA => turb%ffCofj1K
      ffCofB => turb%ffCofj2K
      segId  =  2
      CALL TURB_FloLesGenCoFCLo( global,DIRJ,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

    ELSEIF ((filterWidth(DIRJ)==FILWIDTH_TWO) .OR. &
            (filterWidth(DIRJ)==FILWIDTH_ZERO)) THEN

! --- filter coefficients at J-faces
      segm   => region%levels(iLev)%turb%workI
      ffCofA => turb%ffCofj2J
      ffCofB => turb%ffCofj4J
      segId  =  1
      CALL TURB_FloLesGenCoFFHi( global,DIRJ,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

! --- filter coefficients at I-faces
      segm   => region%levels(iLev)%turb%workI
      ffCofA => turb%ffCofj2I
      ffCofB => turb%ffCofj4I
      segId  =  1
      CALL TURB_FloLesGenCoFCHi( global,DIRJ,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

! --- filter coefficients at K-faces
      segm   => region%levels(iLev)%turb%workK
      ffCofA => turb%ffCofj2K
      ffCofB => turb%ffCofj4K
      segId  =  2
      CALL TURB_FloLesGenCoFCHi( global,DIRJ,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

    ENDIF
  ENDIF

! filter coefficients in K direction

  IF (homDir(DIRK) == OFF) THEN
    IF (filterWidth(DIRK)==FILWIDTH_ONE) THEN

! --- filter coefficients at K-faces
      segm   => region%levels(iLev)%turb%workJ
      ffCofA => turb%ffCofk1K
      ffCofB => turb%ffCofk2K
      segId  =  1
      CALL TURB_FloLesGenCoFFLo( global,DIRK,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

! --- filter coefficients at I-faces
      segm   => region%levels(iLev)%turb%workI
      ffCofA => turb%ffCofk1I
      ffCofB => turb%ffCofk2I
      segId  =  2
      CALL TURB_FloLesGenCoFCLo( global,DIRK,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

! --- filter coefficients at J-faces
      segm   => region%levels(iLev)%turb%workJ
      ffCofA => turb%ffCofk1J
      ffCofB => turb%ffCofk2J
      segId  =  1
      CALL TURB_FloLesGenCoFCLo( global,DIRK,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

    ELSEIF ((filterWidth(DIRK)==FILWIDTH_TWO) .OR. &
            (filterWidth(DIRK)==FILWIDTH_ZERO)) THEN

! --- filter coefficients at K-faces
      segm   => region%levels(iLev)%turb%workJ
      ffCofA => turb%ffCofk2K
      ffCofB => turb%ffCofk4K
      segId  =  1
      CALL TURB_FloLesGenCoFFHi( global,DIRK,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

! --- filter coefficients at I-faces
      segm   => region%levels(iLev)%turb%workI
      ffCofA => turb%ffCofk2I
      ffCofB => turb%ffCofk4I
      segId  =  2
      CALL TURB_FloLesGenCoFCHi( global,DIRK,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

! --- filter coefficients at J-faces
      segm   => region%levels(iLev)%turb%workJ
      ffCofA => turb%ffCofk2J
      ffCofB => turb%ffCofk4J
      segId  =  1
      CALL TURB_FloLesGenCoFCHi( global,DIRK,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )  

    ENDIF
  ENDIF

! deallocate temporary arrays

  DEALLOCATE( ds )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloLesGenCoFF

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floLesGenCoFF.F90,v $
! Revision 1.5  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/08/04 02:50:49  wasistho
! removed turb%avgCoI,J,K replaced with turb%workI,J,K
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.2  2003/08/29 01:41:29  wasistho
! Added TARGET attribute to region variable, since pointers are cached into it
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************







