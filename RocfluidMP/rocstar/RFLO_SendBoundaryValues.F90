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
! Purpose: store outgoing data into buffers for GenX.
!
! Description: none.
!
! Input: region     = dimensions of patches, types of BC`s, flow variables
!        initialize = initial data for GenX (true/false)..
!
! Output: regions%levels%patches = data in buffers for GenX.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_SendBoundaryValues.F90,v 1.11 2010/02/18 21:47:38 juzhang Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_SendBoundaryValues( region,initialize )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
        RFLO_GetPatchIndicesNodes, RFLO_GetCellOffset, RFLO_GetNodeOffset, &
        MixtPerf_R_M, MixtPerf_D_PRT
#ifdef PLAG
  USE ModInterfaces, ONLY : PLAG_SetSizeGenx
#endif
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

  LOGICAL :: initialize

! ... loop variables
  INTEGER :: iPatch, i, j, k, ijkC, ijkD, ijkN, n1, n2, ng1, ng2

! ... local variables
  INTEGER :: iLev, bcType, lbnd, iCOff, ijCOff, iNOff, ijNOff, ijkN1
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: inode, jnode, knode, i2d, distrib, indMol, nOff

  REAL(RFREAL)          :: sgn, dS, mRate, tBurn, rgas, sxn, syn, szn, &
                           edgevector(3), deltn
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), gv(:,:), xyz(:,:), &
                           sFace(:,:), vals(:,:)

  TYPE(t_global), POINTER :: global
  TYPE(t_patch) , POINTER :: patch

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_SendBoundaryValues',&
  'RFLO_SendBoundaryValues.F90' )

! store pointers to variables -------------------------------------------------

  iLev = region%currLevel
  !IF(initialize .eqv. .true.) THEN
  !   WRITE(*,*) 'RFLO: WARNING!! INITIALIZING in SENDBOUNDARYVALS!!'
  !ENDIF
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  indMol = region%levels(iLev)%mixt%indMol

  cv  => region%levels(iLev)%mixt%cv
  dv  => region%levels(iLev)%mixt%dv
  gv  => region%levels(iLev)%mixt%gv
  xyz => region%levels(iLev)%grid%xyz

! loop over all cells of the patch (if an interface)

  DO iPatch=1,region%nPatches

    patch   => region%levels(iLev)%patches(iPatch)
    bcType  =  patch%bcType
    lbnd  =  patch%lbound
    distrib =  patch%mixt%distrib
    nOff    =  ABS(patch%l1end-patch%l1beg) + 1
    vals    => patch%mixt%vals

    IF ( patch%bcCoupled == BC_EXTERNAL .OR. &   ! data from outside
        (patch%bcCoupled == BC_INTERNAL .AND. &  ! data from internal APN
         bcType == BC_INJECTION_APN)) THEN        

      CALL RFLO_GetPatchIndices( region,patch,iLev, &
                                 ibeg,iend,jbeg,jend,kbeg,kend )
      CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )

! --- to take the right face vector and make it point outwards

      sgn   = +1._RFREAL
      inode = 0
      jnode = 0
      knode = 0
      IF (lbnd==2 .OR. lbnd==4 .OR. lbnd==6) THEN
        sgn   = -1._RFREAL
        inode = -idir
        jnode = -jdir
        knode = -kdir
      ENDIF

! --- get the appropriate face vector

      IF (lbnd==1 .OR. lbnd==2) sFace => region%levels(iLev)%grid%si
      IF (lbnd==3 .OR. lbnd==4) sFace => region%levels(iLev)%grid%sj
      IF (lbnd==5 .OR. lbnd==6) sFace => region%levels(iLev)%grid%sk

! --- outgoing data

      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkC = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)
            ijkD = IndIJK(i-idir ,j-jdir ,k-kdir ,iCOff,ijCOff)
            ijkN = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
            dS   = SQRT(sFace(XCOORD,ijkN)*sFace(XCOORD,ijkN)+ &
                        sFace(YCOORD,ijkN)*sFace(YCOORD,ijkN)+ &
                        sFace(ZCOORD,ijkN)*sFace(ZCOORD,ijkN))
            sxn   = sgn*sFace(XCOORD,ijkN)/dS
            syn   = sgn*sFace(YCOORD,ijkN)/dS
            szn   = sgn*sFace(ZCOORD,ijkN)/dS
  
            IF      (lbnd==1 .OR. lbnd==2) THEN
              n1 = j - jbeg
              n2 = k - kbeg
              ijkN1 = IndIJK(i+inode+1*idir,j+jnode,k+knode  ,iNOff,ijNOff)
              IF (lbnd == 2) THEN
                ng1 = j - jbeg + 1
              ELSE
                ng1 = jend - j + 1
              ENDIF
              ng2 = k - kbeg + 1
            ELSE IF (lbnd==3 .OR. lbnd==4) THEN
              n1 = k - kbeg
              n2 = i - ibeg
              ijkN1 = IndIJK(i+inode,j+jnode+1*jdir,k+knode  ,iNOff,ijNOff)
             IF (lbnd == 4) THEN
                ng2 = i - ibeg + 1
              ELSE
                ng2 = iend - i + 1
              ENDIF
              ng1 = k - kbeg + 1
            ELSE IF (lbnd==5 .OR. lbnd==6) THEN
              n1 = i - ibeg
              n2 = j - jbeg
              ijkN1 = IndIJK(i+inode,j+jnode,k+knode+1*kdir  ,iNOff,ijNOff)
              IF (lbnd == 6) THEN
                ng1 = i - ibeg + 1
              ELSE
                ng1 = iend - i + 1
              ENDIF
              ng2 = j - jbeg + 1
            ENDIF
            IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) THEN
              i2d   = distrib * IndIJ(n1,n2,nOff)
              mRate = vals(BCDAT_INJECT_MFRATE,i2d)
              tBurn = vals(BCDAT_INJECT_TEMP  ,i2d)

!  ----   cal distance from the centerroid to the center of the wall face (assuming body-fitted meshes)
              edgevector(1) = xyz(XCOORD,ijkN1) - xyz(XCOORD,ijkN)                
              edgevector(2) = xyz(YCOORD,ijkN1) - xyz(YCOORD,ijkN)                
              edgevector(3) = xyz(ZCOORD,ijkN1) - xyz(ZCOORD,ijkN)                
              deltn = 0.5*abs(sxn*edgevector(1)+syn*edgevector(2)+szn*edgevector(3))

              IF (bcType==BC_INJECTION_APN) THEN
                patch%mdotAlp(    ng1,ng2) = mRate
                patch%tflmAlp(    ng1,ng2) = vals(BCDAT_INJECT_TEMP  ,i2d)
                patch%rhofvfAlp(1,ng1,ng2) = vals(BCDAT_INJECT_RFVFU ,i2d)
                patch%rhofvfAlp(2,ng1,ng2) = vals(BCDAT_INJECT_RFVFV ,i2d)
                patch%rhofvfAlp(3,ng1,ng2) = vals(BCDAT_INJECT_RFVFW ,i2d)
!                IF (mRate > 0._RFREAL) patch%bFlag(ng1,ng2) = 1
              ENDIF
              IF (initialize .eqv. .true.) THEN
                IF (mRate > 0._RFREAL) THEN
                  patch%bFlag(ng1,ng2) = 1   ! face ignited
                ELSE
                  patch%bFlag(ng1,ng2) = 0   ! face not ignited yet
                ENDIF
              ENDIF

              IF (patch%bFlag(ng1,ng2) == 1) THEN    ! burning
                rgas                   = MixtPerf_R_M( gv(GV_MIXT_MOL,ijkC*indMol) )
                patch%pf     (ng1,ng2) = dv(DV_MIXT_PRES,ijkC)
                patch%rhofAlp(ng1,ng2) = MixtPerf_D_PRT( dv(DV_MIXT_PRES,ijkC), &
                                                         rgas,tBurn )
                patch%tempf  (ng1,ng2) = tBurn
              ELSE                                   ! not burning
                patch%pf     (ng1,ng2) = 0.5_RFREAL*(dv(DV_MIXT_PRES,ijkC)+ &
                                                     dv(DV_MIXT_PRES,ijkD))
                patch%rhofAlp(ng1,ng2) = 0.5_RFREAL*(cv(CV_MIXT_DENS,ijkC)+ &
                                                     cv(CV_MIXT_DENS,ijkD))
                patch%tempf  (ng1,ng2) = 0.5_RFREAL*(dv(DV_MIXT_TEMP,ijkC)+ &
                                                     dv(DV_MIXT_TEMP,ijkD))
              ENDIF
                patch%tempv  (ng1,ng2) = dv(DV_MIXT_TEMP,ijkC)
              IF (bcType==BC_INJECTION_HT) THEN
                patch%dnml   (ng1,ng2) = -deltn
             ELSE
                patch%dnml   (ng1,ng2) = deltn
              ENDIF
            ELSE                           ! not an injection boundary
              IF (initialize) THEN
                patch%bFlag(ng1,ng2) = 0   ! can never burn
              ENDIF
              patch%pf     (ng1,ng2) = 0.5_RFREAL*(dv(DV_MIXT_PRES,ijkC)+ &
                                                   dv(DV_MIXT_PRES,ijkD))
              patch%rhofAlp(ng1,ng2) = 0.5_RFREAL*(cv(CV_MIXT_DENS,ijkC)+ &
                                                   cv(CV_MIXT_DENS,ijkD))
              patch%tempf  (ng1,ng2) = 0.5_RFREAL*(dv(DV_MIXT_TEMP,ijkC)+ &
                                                   dv(DV_MIXT_TEMP,ijkD))
            ENDIF
            patch%qc   (  ng1,ng2) = 0._RFREAL
            patch%nfAlp(1,ng1,ng2) = sgn*sFace(XCOORD,ijkN)/dS
            patch%nfAlp(2,ng1,ng2) = sgn*sFace(YCOORD,ijkN)/dS
            patch%nfAlp(3,ng1,ng2) = sgn*sFace(ZCOORD,ijkN)/dS
            patch%tracf(1,ng1,ng2) = patch%pf(ng1,ng2)*patch%nfAlp(1,ng1,ng2)
            patch%tracf(2,ng1,ng2) = patch%pf(ng1,ng2)*patch%nfAlp(2,ng1,ng2)
            patch%tracf(3,ng1,ng2) = patch%pf(ng1,ng2)*patch%nfAlp(3,ng1,ng2)
          ENDDO
        ENDDO
      ENDDO

    ENDIF  ! external BC

    IF (patch%bcMotion == BC_EXTERNAL) THEN        ! data from outside

! ----surface grid

      CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                      ibeg,iend,jbeg,jend,kbeg,kend )

      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkN = IndIJK(i,j,k,iNOff,ijNOff)
            IF      (lbnd==1 .OR. lbnd==2) THEN
              IF (lbnd == 2) THEN
                ng1 = j - jbeg + 1
              ELSE
                ng1 = jend - j + 1
              ENDIF
              ng2 = k - kbeg + 1
            ELSE IF (lbnd==3 .OR. lbnd==4) THEN
              ng1 = k - kbeg + 1
              IF (lbnd == 4) THEN
                ng2 = i - ibeg + 1
              ELSE
                ng2 = iend - i + 1
              ENDIF
            ELSE IF (lbnd==5 .OR. lbnd==6) THEN
              IF (lbnd == 6) THEN
                ng1 = i - ibeg + 1
              ELSE
                ng1 = iend - i + 1
              ENDIF
              ng2 = j - jbeg + 1
            ENDIF
            patch%surfCoord(1,ng1,ng2) = xyz(XCOORD,ijkN)
            patch%surfCoord(2,ng1,ng2) = xyz(YCOORD,ijkN)
            patch%surfCoord(3,ng1,ng2) = xyz(ZCOORD,ijkN)
          ENDDO
        ENDDO
      ENDDO

    ENDIF  ! external Bc and motion
  ENDDO    ! iPatch

#ifdef PLAG
  IF (global%plagUsed)  CALL PLAG_SetSizeGenx( region )
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_SendBoundaryValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_SendBoundaryValues.F90,v $
! Revision 1.11  2010/02/18 21:47:38  juzhang
! Heat transfer bc for non-propellant surface documented in Rocburn_PY_HT.pdf in Rocburn_PY directory is implemented within Rocburn_PY. Major changes were made to Rocburn, Rocman3, RocfluidMP/genx, RocfluidMP/modflo directories.
!
! Revision 1.10  2008/12/06 08:44:01  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2006/08/19 15:37:50  mparmar
! Renamed patch variables
!
! Revision 1.7  2006/01/25 06:04:03  wasistho
! set bflag for injectionAPN
!
! Revision 1.6  2006/01/25 04:32:11  wasistho
! added sending rhofvfalp and tflmalp for injectionAPN
!
! Revision 1.5  2006/01/23 22:44:21  wasistho
! added condition for internal injectionAPN
!
! Revision 1.4  2005/07/03 03:42:52  wasistho
! bcMotion external for surfCoord
!
! Revision 1.3  2005/02/01 21:18:00  wasistho
! reactivated bflag
!
! Revision 1.2  2005/01/20 20:05:49  wasistho
! temporary outcommented initial bflag condition
!
! Revision 1.1  2004/12/01 21:23:54  haselbac
! Initial revision after changing case
!
! Revision 1.18  2004/11/13 22:37:55  wasistho
! invert orientation of genx-surface-variables
!
! Revision 1.17  2004/07/02 22:07:19  fnajjar
! Added call to PLAG_SetSizeGenx for Roccom3 import
!
! Revision 1.16  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.13  2003/08/14 20:06:58  jblazek
! Corrected bug associated with radiation flux qr.
!
! Revision 1.12  2003/08/09 02:14:14  wasistho
! restore version 1.10
!
! Revision 1.11  2003/08/09 02:09:26  wasistho
! added TURB and RADI_initGenxInterface
!
! Revision 1.10  2003/05/15 02:57:00  jblazek
! Inlined index function.
!
! Revision 1.9  2002/10/15 23:22:59  jblazek
! dded Rocturb to GenX compilation path.
!
! Revision 1.8  2002/10/15 17:51:35  jblazek
! And flippt the normals again ...
!
! Revision 1.7  2002/10/15 00:38:18  jblazek
! Added routine to send newest fluids density to GenX.
!
! Revision 1.6  2002/10/10 23:49:48  jblazek
! Changed orientation of surface grid.
!
! Revision 1.5  2002/10/03 21:33:48  jblazek
! Init. of bcflag moved from SendBoundaryValues to InitGenxInterface.
!
! Revision 1.4  2002/10/01 21:53:27  jiao
! Modified surface temperature evaluation for burning surfaces.
!
! Revision 1.3  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.2  2002/09/24 23:18:01  jblazek
! Changed bcflag to a pointer.
!
! Revision 1.1  2002/09/20 22:22:34  jblazek
! Finalized integration into GenX.
!
!******************************************************************************







