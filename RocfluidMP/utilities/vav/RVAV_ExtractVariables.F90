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
! Purpose: extract that primitive and derived variables for a thermally and 
!          calorically perfect gas mixture. the extracted variables (ev) are
!          stored in 3D arrays ev(ibeg:iend,jbeg:jend,kbeg:kend). for the 
!          moment we have not incorporated the ijump, jjunp and kjump options
!
! Description: this subroutine work with Stream1 & Stream2
!
! Input: region            ! the region from which we are extracting
!        ibeg, iend, ijump ! the dimensions of the subregion of the region
!        jbeg, jend, jjump !                       
!        kbeg, kend, kjump !                       
!
! Output: ev = extracted variables (e.g.: rho, u, v, w, p, T, s, omega )
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RVAV_ExtractVariables.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RVAV_ExtractVariables( global, region,            &
                                          ibeg, iend, ijump, &
                                          jbeg, jend, jjump, &
                                          kbeg, kend, kjump, &
                                          iCOff, ijCOff, &
                                          variableIndex,     &
                                          fileType,          &
                                          indCp, indMol, ev )

  USE ModMPI
  USE ModParameters
  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE ModMixture, ONLY    : t_mixt
  USE ModError
  USE ModInterfaces, ONLY : RFLO_GetCellOffset 
  USE RVAV_ModParameters
  USE RVAV_ModGlobal
  IMPLICIT NONE

#include "Indexing.h"

! ... variables in the parameter list of the subroutine
  TYPE(t_region) :: region
  TYPE(t_global), POINTER :: global
  INTEGER :: ibeg, iend, ijump, jbeg, jend, jjump, kbeg, kend, kjump
  INTEGER :: iCOff, ijCOff
  INTEGER :: indCp, indMol
  INTEGER :: variableIndex
  INTEGER :: fileType
  REAL(RFREAL), POINTER :: ev(:,:,:) ! this is an array for a single extracted
                                     ! variable with i,j,k dimensions

! ... local variables
  CHARACTER(CHRLEN) :: msg

  INTEGER :: ipc, jpc, kpc ! the indices of the physical cells
  INTEGER :: ic            ! the cell index which is obtained from (ipc,jpc,kpc)
                           ! by using the IndIJK function in Indexing.h
  INTEGER :: is, js, ks    ! the indices of the shifted cells
  INTEGER :: iLev, iOffset, ijOffset
  TYPE(t_grid),   POINTER :: grid
  TYPE(t_mixt),   POINTER :: mixt
  REAL(RFREAL), POINTER :: cv(:,:), gv(:,:)

  REAL(RFREAL) :: rgas, gamma, rhoq, pressure, temperature, sound_speed, &
                  Mach_number

!******************************************************************************

  CALL RegisterFunction( global, 'RVAV_ExtractVariables',&
  'RVAV_ExtractVariables.F90' )

! ... set local pointers to the global pointers

  iLev = global%startLevel
  grid => region%levels(iLev)%grid
  mixt => region%levels(iLev)%mixt
  cv   => mixt%cv
  gv   => mixt%gv

! ... get the offsets for the region
  
  IF ( global%verbLevel/=VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(/,A)')        'RVAV Extract'
    WRITE(STDOUT,'(A,3(I5,3X))') '  ibeg,iend,ijump',ibeg,iend,ijump
    WRITE(STDOUT,'(A,3(I5,3X))') '  jbeg,jend,jjump',jbeg,jend,jjump
    WRITE(STDOUT,'(A,3(I5,3X))') '  kbeg,kend,kjump',kbeg,kend,kjump
    WRITE(STDOUT,'(A,3(I5,3X))') '  iLev,iCOff,ijCOff',iLev,iCOff,ijCOff
  END IF !verbLevel
  
  SELECT CASE ( fileType )

    CASE (FILE_COMPUTED)
        iOffset  = iCOff
        ijOffset = ijCOff

    CASE (FILE_ANALYTICAL)
        iOffset  = iend-ibeg+1
        ijOffset = iOffset* (jend-jbeg+1)
        IF (globalRVAV%SimilarityTypeS2==RVAV_BLASIUS) THEN
          iOffset = iCOff
          ijOffset = iCOff*(jend-jbeg+1) 
        ENDIF
    CASE (FILE_EXPERIMENTAL)
        iOffset  = iend-ibeg+1
        ijOffset = iOffset* (jend-jbeg+1)

  END SELECT

  IF ( global%verbLevel/=VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3(I5,3X))') '  iLev,iOffset,ijOffset',iLev,iOffset,ijOffset
  END IF !verbLevel
  
! ... loop over the ipc,jpc,kpc indices to get the extracted variable (ev)

  DO kpc = kbeg, kend, kjump
    DO jpc = jbeg, jend, jjump
      DO ipc = ibeg, iend, ijump

! ... convert the ipc, jpc, kpc indices to the cell index ic
        ic = IndIJK( ipc, jpc, kpc, iOffset, ijOffset )

        is = INT(REAL(ipc-ibeg,KIND=RFREAL)/REAL(ijump,KIND=RFREAL))+1
        js = INT(REAL(jpc-jbeg,KIND=RFREAL)/REAL(jjump,KIND=RFREAL))+1
        ks = INT(REAL(kpc-kbeg,KIND=RFREAL)/REAL(kjump,KIND=RFREAL))+1

! ... calculate some of the gas constants needed for variable extraction

        IF ( fileType == FILE_COMPUTED  ) THEN

          rgas = 8314.3_RFREAL/gv(GV_MIXT_MOL, ic*indMol)
          gamma = gv(GV_MIXT_CP,ic*indCp)/( gv(GV_MIXT_CP,ic*indCp)-rgas )
  
! ... select appropriate cases 

        SELECT CASE( VariableIndex )
          CASE( RVAV_RHO )
            ev( is, js, ks ) =  cv(CV_MIXT_DENS, ic)

          CASE( RVAV_RHOU )
            ev( is, js, ks ) =  cv(CV_MIXT_XMOM, ic)

          CASE( RVAV_RHOV )
            ev( is, js, ks ) =  cv(CV_MIXT_YMOM, ic)

          CASE( RVAV_RHOW )
            ev( is, js, ks ) =  cv(CV_MIXT_ZMOM, ic)

          CASE( RVAV_RHOET )
            ev( is, js, ks ) =  cv(CV_MIXT_ENER, ic)

          CASE( RVAV_U )
            ev( is, js, ks ) =  cv(CV_MIXT_XMOM, ic)/cv(CV_MIXT_DENS, ic)

          CASE( RVAV_V )
            ev( is, js, ks ) =  cv(CV_MIXT_YMOM, ic)/cv(CV_MIXT_DENS, ic)

          CASE( RVAV_W )
            ev( is, js, ks ) =  cv(CV_MIXT_ZMOM, ic)/cv(CV_MIXT_DENS, ic)

          CASE( RVAV_P )
            rhoq = 0.5_RFREAL* &
                   ( 1.0_RFREAL/cv(CV_MIXT_DENS,ic) )* &
                   ( cv(CV_MIXT_XMOM,ic)*cv(CV_MIXT_XMOM,ic) + &
                     cv(CV_MIXT_YMOM,ic)*cv(CV_MIXT_YMOM,ic) + &
                     cv(CV_MIXT_ZMOM,ic)*cv(CV_MIXT_ZMOM,ic) )
            ev( is, js, ks ) =  ( gamma-1._RFREAL )* &
                                   ( cv(CV_MIXT_ENER, ic) - rhoq )

          CASE( RVAV_T )
            rhoq = 0.5_RFREAL* &
                   ( 1.0_RFREAL/cv(CV_MIXT_DENS,ic) )* &
                   ( cv(CV_MIXT_XMOM,ic)*cv(CV_MIXT_XMOM,ic) + &
                     cv(CV_MIXT_YMOM,ic)*cv(CV_MIXT_YMOM,ic) + &
                     cv(CV_MIXT_ZMOM,ic)*cv(CV_MIXT_ZMOM,ic) )
            pressure =  ( gamma-1._RFREAL )* &
                        ( cv(CV_MIXT_ENER, ic) - rhoq )
            ev( is, js, ks ) = pressure/( cv(CV_MIXT_DENS, ic)*rgas )

          CASE( RVAV_TSTAG )
            rhoq = 0.5_RFREAL* &
                   ( 1.0_RFREAL/cv(CV_MIXT_DENS,ic) )* &
                   ( cv(CV_MIXT_XMOM,ic)*cv(CV_MIXT_XMOM,ic) + &
                     cv(CV_MIXT_YMOM,ic)*cv(CV_MIXT_YMOM,ic) + &
                     cv(CV_MIXT_ZMOM,ic)*cv(CV_MIXT_ZMOM,ic) )
            pressure =  ( gamma-1._RFREAL )* &
                        ( cv(CV_MIXT_ENER, ic) - rhoq )
            temperature = pressure/( cv(CV_MIXT_DENS, ic)*rgas )
            sound_speed = SQRT(gamma*rgas*temperature)
            Mach_number = SQRT(2.0_RFREAL*rhoq/cv(CV_MIXT_DENS,ic))/sound_speed
            ev( is, js, ks ) = (1.0_RFREAL + &
                                  0.5_RFREAL*(gamma-1.0_RFREAL)* &
                                  Mach_number*Mach_number)*&
                                  temperature

          CASE( RVAV_PSTAG )
            rhoq = 0.5_RFREAL* &
                   ( 1.0_RFREAL/cv(CV_MIXT_DENS,ic) )* &
                   ( cv(CV_MIXT_XMOM,ic)*cv(CV_MIXT_XMOM,ic) + &
                     cv(CV_MIXT_YMOM,ic)*cv(CV_MIXT_YMOM,ic) + &
                     cv(CV_MIXT_ZMOM,ic)*cv(CV_MIXT_ZMOM,ic) )
            pressure =  ( gamma-1._RFREAL )* &
                        ( cv(CV_MIXT_ENER, ic) - rhoq )
            temperature = pressure/( cv(CV_MIXT_DENS, ic)*rgas )
            sound_speed = SQRT(gamma*rgas*temperature)
            Mach_number = SQRT(2.0_RFREAL*rhoq/cv(CV_MIXT_DENS,ic))/sound_speed
            ev( is, js, ks ) = 1.0_RFREAL + &
                                  0.5_RFREAL*(gamma-1.0_RFREAL)* &
                                  Mach_number*Mach_number
            ev( is, js, ks ) = (ev( is, js, ks)**(gamma/(gamma-1.0_RFREAL)))*&
                                  pressure
        END SELECT
        END IF ! fileType

! - For Analytical and Experimental Extract as is while shifting cv by 3.

        IF ( fileType == FILE_ANALYTICAL .OR. &
             fileType == FILE_EXPERIMENTAL     ) THEN

          ev( is, js, ks ) = cv(variableIndex-ZCOORD,ic)

       END IF ! fileType


! ... DO NOT DEALLOCATE ev. It is a pointer that is required outside this code
! ... No derived variables  such as entropy and vorticity.

      ENDDO
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RVAV_ExtractVariables

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_ExtractVariables.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:43:21  fnajjar
! Initial revision after changing case
!
! Revision 1.11  2003/11/21 23:40:45  fnajjar
! Remove quotes from comment lines as it confuses the compiler
!
! Revision 1.10  2003/11/20 16:40:41  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/05/15 02:57:08  jblazek
! Inlined index function.
!
! Revision 1.5  2002/09/10 00:01:45  f-najjar
! Variable global moved into regions()
!
! Revision 1.4  2002/06/24 15:50:10  f-najjar
! Fixed Offset for Blasius
!
! Revision 1.3  2002/06/19 20:22:40  f-najjar
! Cleaned unnecessary PRINT statement
!
! Revision 1.2  2002/06/15 17:28:52  f-najjar
! Bug Fix of iOffSet and ijOffset
!
! Revision 1.1.1.1  2002/06/03 21:41:29  f-najjar
! Initial Import of RocVaV
!







