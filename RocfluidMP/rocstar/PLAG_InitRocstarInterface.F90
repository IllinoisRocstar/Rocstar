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
! Purpose: register PLAG windows with GenX.
!
! Description: none.
!
! Input: regions    = data of all regions,
!        wins, winp = GenX window registrations.
!
! Output: to Roccom.
!
! Notes: Surface registration for Tiles works only for External coupled bc.
!        Need to activate for both Internal and External bc.
!
!******************************************************************************
!
! $Id: PLAG_InitGenxInterface.F90,v 1.5 2009/10/26 00:19:31 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_InitGenxInterface( regions, wins, inPlag, obtain_dataitem )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  USE ModPartLag, ONLY    : t_plag, t_tile_plag
  IMPLICIT NONE
  INCLUDE 'comf90.h'

! ... parameters

  TYPE(t_region), POINTER :: regions(:)

  CHARACTER(CHRLEN) :: wins, inPlag

  INTEGER :: obtain_dataitem

! ... loop variables

  INTEGER :: iCont, iReg, iPatch

! ... local variables
  
  CHARACTER(CHRLEN) :: winp

  INTEGER, PARAMETER :: ASCII_ZERO = 48   ! char representation of zero
  INTEGER :: bcType, errorFlag, ilb, icount, iLev,pid
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: nAiv, nArv, nCont, nContMax, nCv, nCvTile, nDv, nDvTile
  INTEGER :: nDimPlag
  INTEGER :: plagStatus
  INTEGER, POINTER, DIMENSION(:) :: pCvTileMass
  INTEGER, POINTER :: pInt

  TYPE(t_global),     POINTER :: pGlobal
  TYPE(t_patch) ,     POINTER :: pPatch
  TYPE(t_plag)  ,     POINTER :: pPlag
  TYPE(t_tile_plag),  POINTER :: pTilePlag

!******************************************************************************

  pGlobal => regions(1)%global

  CALL RegisterFunction( pGlobal,'PLAG_InitGenxInterface',&
  'PLAG_InitRocstarInterface.F90' )

! set PLAG window -------------------------------------------------------------

  winp = TRIM(pGlobal%winName)//'_plag'
  pGlobal%winp = winp
  
! input data (currently none) -------------------------------------------------

! output surface data ---------------------------------------------------------

  CALL COM_new_dataitem( TRIM(wins)//'.plag_momnrm', 'e',  &
                          COM_DOUBLE, 1, 'kg/(m^2 s)' )
  CALL COM_new_dataitem( TRIM(wins)//'.plag_ener',   'e',  &
                          COM_DOUBLE, 1, 'J/kg'       )

  nContMax = 0
  IF ( pGlobal%plagUsed ) THEN
    DO iReg=1,pGlobal%nRegions
      nContMax = MAX(nContMax,regions(iReg)%plagInput%nCont)
    ENDDO ! iReg
  ENDIF ! plagUsed

  DO iCont=1,nContMax
    CALL COM_new_dataitem( TRIM(wins)//'.plag_mass'//CHAR(iCont+ASCII_ZERO), 'e', &
                            COM_DOUBLE, 1, 'kg/m^3' )
  ENDDO ! iCont

  CALL COM_new_dataitem( TRIM(wins)//'.plag_dv_timefctr', 'e',   &
                          COM_DOUBLE, 1, ''  )

  CALL COM_new_dataitem( TRIM(wins)//'.plag_dv_diam'    , 'e',   &
                          COM_DOUBLE, 1, 'm' )

  CALL COM_new_dataitem( TRIM(wins)//'.plag_dv_spload'  , 'e',   &
                          COM_DOUBLE, 1, ''  )

! restart data (av, cv ) and additional plot data (dv= diam) ------------------
  CALL COM_new_window( TRIM(winp))

  CALL COM_new_dataitem( TRIM(winp)//'.plag_aiv_pidini', 'n' ,  &
                          COM_INTEGER, 1, '' )
  CALL COM_new_dataitem( TRIM(winp)//'.plag_aiv_regini', 'n' ,  &
                          COM_INTEGER, 1, '' )
  CALL COM_new_dataitem( TRIM(winp)//'.plag_aiv_regcrt', 'n' ,  &
                          COM_INTEGER, 1, '' )
  CALL COM_new_dataitem( TRIM(winp)//'.plag_aiv_icells', 'n' ,  &
                          COM_INTEGER, 1, '' )
  CALL COM_new_dataitem( TRIM(winp)//'.plag_aiv_indexi', 'n' ,  &
                          COM_INTEGER, 1, '' )
  CALL COM_new_dataitem( TRIM(winp)//'.plag_aiv_indexj', 'n' ,  &
                          COM_INTEGER, 1, '' )
  CALL COM_new_dataitem( TRIM(winp)//'.plag_aiv_indexk', 'n' ,  &
                          COM_INTEGER, 1, '' )
  CALL COM_new_dataitem( TRIM(winp)//'.plag_aiv_burnstat', 'n' ,&
                          COM_INTEGER, 1, '' )
  CALL COM_new_dataitem( TRIM(winp)//'.plag_arv_spload', 'n' ,  &
                          COM_DOUBLE, 1, '' )

  CALL COM_new_dataitem( TRIM(winp)//'.plag_xmom', 'n',         &
                          COM_DOUBLE, 1, 'kg/(m^2 s)'  )
  CALL COM_new_dataitem( TRIM(winp)//'.plag_ymom', 'n',         &
                          COM_DOUBLE, 1, 'kg/(m^2 s)'  )
  CALL COM_new_dataitem( TRIM(winp)//'.plag_zmom', 'n',         &
                          COM_DOUBLE, 1, 'kg/(m^2 s)'  )
  CALL COM_new_dataitem( TRIM(winp)//'.plag_ener', 'n',         &
                          COM_DOUBLE, 1, 'J/kg'        )
  CALL COM_new_dataitem( TRIM(winp)//'.plag_enervapor', 'n',    &
                          COM_DOUBLE, 1, 'J/kg'        )
  DO iCont = 1, nContMax
    CALL COM_new_dataitem( TRIM(winp)//'.plag_mass'//CHAR(iCont+ASCII_ZERO), 'n', &
                            COM_DOUBLE,1, 'kg/m^3'     )
  ENDDO ! iCont

  CALL COM_new_dataitem( TRIM(winp)//'.plag_diam', 'n',         &
                          COM_DOUBLE, 1,'m'            )

  CALL COM_new_dataitem( TRIM(winp)//'.plag_nextid', 'p',COM_INTEGER, 1,'')

! loop over all regions -------------------------------------------------------

  DO iReg=1,pGlobal%nRegions
    IF ( regions(iReg)%procid==pGlobal%myProcid .AND. &   ! region active and
         regions(iReg)%active==ACTIVE           .AND. &   ! on my processor
         pGlobal%plagUsed ) THEN                          ! and particles used
 
      iLev   = regions(iReg)%currLevel
      icount = 0

! - set pointer ---------------------------------------------------------------

      pPlag => regions(iReg)%levels(iLev)%plag



! - get dimensions ------------------------------------------------------------

      nDimPlag = regions(iReg)%plagInput%nPclsMax
      nCont    = regions(iReg)%plagInput%nCont
      nAiv     = pPlag%nAiv
      nArv     = pPlag%nArv
      nCv      = pPlag%nCv
      nDv      = pPlag%nDv
      nCvTile  = CV_TILE_LAST+nCont
      nDvTile  = DV_TILE_LAST

! - surface data for tile infrastructure --------------------------------------

      DO iPatch=1,regions(iReg)%nPatches
        pPatch => regions(iReg)%levels(iLev)%patches(iPatch)
        bcType =  pPatch%bcType

        IF ( pPatch%bcCoupled == BC_EXTERNAL   .AND.                      &
             ( bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE ) ) THEN
          icount  = icount + 1
          pid     = iReg*REGOFF + icount

          pTilePlag   => pPatch%tilePlag
          pCvTileMass => pTilePlag%cvTileMass

! -- output tile data ---------------------------------------------------------

          ilb = LBOUND(pTilePlag%cv,2)
          CALL COM_set_array( TRIM(wins)//'.plag_momnrm',pid, &
                              pTilePlag%cv(CV_TILE_MOMNRM,ilb),nCvTile)
          CALL COM_set_array( TRIM(wins)//'.plag_ener',  pid, &
                              pTilePlag%cv(CV_TILE_ENER,ilb),nCvTile)

          DO iCont = 1, nCont
            CALL COM_set_array( TRIM(wins)//'.plag_mass'//CHAR(iCont+ASCII_ZERO),&
                                pid, pTilePlag%cv(pCvTileMass(iCont),ilb),nCvTile)
          ENDDO ! iCont

          ilb = LBOUND(pTilePlag%dv,2)
          CALL COM_set_array( TRIM(wins)//'.plag_dv_timefctr',  &
                              pid, pTilePlag%dv(DV_TILE_COUNTDOWN,ilb),nDvTile)
          CALL COM_set_array( TRIM(wins)//'.plag_dv_diam',  &
                              pid, pTilePlag%dv(DV_TILE_DIAM    ,ilb),nDvTile)
          CALL COM_set_array( TRIM(wins)//'.plag_dv_spload',    &
                              pid, pTilePlag%dv(DV_TILE_SPLOAD  ,ilb),nDvTile)

        ENDIF   ! external BC
      ENDDO     ! iPatch

! -- volume data --------------------------------------------------------------

! -- create PLAG pane data ----------------------------------------------------

      pid = iReg*REGOFF+1

! Obtain status
      plagStatus = COM_get_status( TRIM(inPlag),pid )
      
! Obtain the number of particles
! check if status is inexistent first.

      IF ( plagStatus == -1 ) THEN
#ifndef NATIVE_MP_IO
        pPlag%nPcls = 0
#endif
      ELSE 
        CALL COM_get_size( TRIM(inPlag)//'.nc',pid,pPlag%nPcls)
      END IF ! plagStatus

! COM_set_size procedure must be called in RFLO_sendBoundaryValues
! if nPcls has changed.

      CALL COM_set_size( TRIM(winp)//'.nc',pid,pPlag%nPcls)

! COM_set_array on '.nc' processes all three components at once
      ilb = LBOUND(pPlag%cv,2)
      CALL COM_set_array( TRIM(winp)//'.nc',pid, &
                          pPlag%cv(CV_PLAG_XPOS,ilb),nCv,nDimPlag )

      pInt => pPlag%nextIdNumber
      CALL COM_set_size( TRIM(winp)//'.plag_nextid',pid,1)
      CALL COM_set_array( TRIM(winp)//'.plag_nextid', pid, pInt)

! -- aiv data -----------------------------------------------------------------

      ilb = LBOUND(pPlag%aiv,2)
      CALL COM_set_array( TRIM(winp)//'.plag_aiv_pidini',pid,          &
                          pPlag%aiv(AIV_PLAG_PIDINI,ilb),nAiv,nDimPlag )
      CALL COM_set_array( TRIM(winp)//'.plag_aiv_regini',pid,          &
                          pPlag%aiv(AIV_PLAG_REGINI,ilb),nAiv,nDimPlag )
      CALL COM_set_array( TRIM(winp)//'.plag_aiv_regcrt',pid,          &
                          pPlag%aiv(AIV_PLAG_REGCRT,ilb),nAiv,nDimPlag )
      CALL COM_set_array( TRIM(winp)//'.plag_aiv_icells',pid,          &
                          pPlag%aiv(AIV_PLAG_ICELLS,ilb),nAiv,nDimPlag )
      CALL COM_set_array( TRIM(winp)//'.plag_aiv_indexi',pid,          &
                          pPlag%aiv(AIV_PLAG_INDEXI,ilb),nAiv,nDimPlag )
      CALL COM_set_array( TRIM(winp)//'.plag_aiv_indexj',pid,          &
                          pPlag%aiv(AIV_PLAG_INDEXJ,ilb),nAiv,nDimPlag )
      CALL COM_set_array( TRIM(winp)//'.plag_aiv_indexk',pid,          &
                          pPlag%aiv(AIV_PLAG_INDEXK,ilb),nAiv,nDimPlag )
      CALL COM_set_array( TRIM(winp)//'.plag_aiv_burnstat',pid,          &
                          pPlag%aiv(AIV_PLAG_BURNSTAT,ilb),nAiv,nDimPlag )

! -- arv data -----------------------------------------------------------------

      ilb = LBOUND(pPlag%arv,2)
      CALL COM_set_array( TRIM(winp)//'.plag_arv_spload',pid,          &
                          pPlag%arv(ARV_PLAG_SPLOAD,ilb),nArv,nDimPlag )

! -- cv data ------------------------------------------------------------------

      ilb = LBOUND(pPlag%cv,2)
      CALL COM_set_array( TRIM(winp)//'.plag_xmom',pid,           &
                          pPlag%cv(CV_PLAG_XMOM,ilb),nCv,nDimPlag )
      CALL COM_set_array( TRIM(winp)//'.plag_ymom',pid,           &
                          pPlag%cv(CV_PLAG_YMOM,ilb),nCv,nDimPlag )
      CALL COM_set_array( TRIM(winp)//'.plag_zmom',pid,           &
                          pPlag%cv(CV_PLAG_ZMOM,ilb),nCv,nDimPlag )
      CALL COM_set_array( TRIM(winp)//'.plag_ener',pid,           &
                          pPlag%cv(CV_PLAG_ENER,ilb),nCv,nDimPlag )
      CALL COM_set_array( TRIM(winp)//'.plag_enervapor',pid,           &
                          pPlag%cv(CV_PLAG_ENERVAPOR,ilb),nCv,nDimPlag )
      DO iCont = 1, nCont
        CALL COM_set_array( TRIM(winp)//'.plag_mass'//CHAR(iCont+ASCII_ZERO),  &
                            pid, pPlag%cv(CV_PLAG_LAST+iCont,ilb),nCv,nDimPlag )
      ENDDO ! iCont

! -- dv data (for plotting) ---------------------------------------------------

      ilb = LBOUND(pPlag%dv,2)
      CALL COM_set_array( TRIM(winp)//'.plag_diam',pid,           &
                          pPlag%dv(DV_PLAG_DIAM,ilb),nDv,nDimPlag )

    ENDIF      ! region on this processor and active

  ENDDO        ! iReg

  CALL COM_window_init_done( TRIM(winp))

#ifndef NATIVE_MP_IO
  CALL COM_call_function( obtain_dataitem,2, &
                          COM_get_dataitem_handle_const(TRIM(inPlag)//".all"), &
                          COM_get_dataitem_handle(TRIM(winp)//".all") )
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( pGlobal )

END SUBROUTINE PLAG_InitGenxInterface

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InitGenxInterface.F90,v $
! Revision 1.5  2009/10/26 00:19:31  mtcampbe
! Updates for completion of NATIVE_MP_IO
!
! Revision 1.4  2008/12/06 08:44:00  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/03/13 18:54:19  fnajjar
! Fixed bug to change nPclsTot to nPclsMax
!
! Revision 1.1  2004/12/01 21:23:44  haselbac
! Initial revision after changing case
!
! Revision 1.14  2004/11/05 22:29:31  fnajjar
! Included a COM_get_status call before calling COM_get_size
!
! Revision 1.13  2004/11/05 19:12:06  fnajjar
! Deleted REGOFF definition since it now resides in ModParameters
!
! Revision 1.12  2004/10/03 03:04:46  fnajjar
! Bug fix for dataitem surface variable changed to plag_dv_diam
!
! Revision 1.11  2004/07/02 22:05:56  fnajjar
! Modified routine for Roccom3 import
!
! Revision 1.10  2004/06/16 23:08:48  fnajjar
! Renamed DV_TILE_TIMEFCTR to DV_TILE_COUNTDOWN for CRE kernel
!
! Revision 1.9  2004/05/12 14:13:29  fnajjar
! Added missing array entries for aiv and cv to Genx IO
!
! Revision 1.8  2004/03/16 21:24:51  fnajjar
! Bug fix for plagUsed as it needs to use pGlobal pointer rather than global pointer
!
! Revision 1.7  2004/03/05 22:08:58  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2004/03/05 16:27:18  fnajjar
! Added dv(diam) and dv(spload) from tile datastructure to insure proper restart
!
! Revision 1.5  2004/01/21 16:37:20  fnajjar
! Fixed semantics of iCont in Roccom window to bypass Frost compilation error
!
! Revision 1.4  2004/01/20 23:22:48  fnajjar
! Fixed CHAR operation in RocCom window naming to correct Frost compilation error
!
! Revision 1.3  2003/12/03 03:02:07  jiao
! Removed all calls involving COM_NULL.
!
! Revision 1.2  2003/11/26 22:04:12  fnajjar
! Defined nDvTile correctly for Roccom restart
!
! Revision 1.1  2003/11/21 22:21:42  fnajjar
! Initial import of Rocpart GenX interfaces
!
!******************************************************************************







