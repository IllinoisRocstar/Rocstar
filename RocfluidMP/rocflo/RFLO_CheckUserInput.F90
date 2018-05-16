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
! Purpose: check parameters specified by the user.
!
! Description: none.
!
! Input: regions = input parameters for all grid regions.
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_CheckUserInput.F90,v 1.16 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CheckUserInput( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMPI
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  INTEGER :: istop
  INTEGER :: viscModel
  LOGICAL :: computeTv

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_CheckUserInput',&
  'RFLO_CheckUserInput.F90' )

! fixed parameters check ------------------------------------------------------

! check consistency of mixture gradient and dv indices

  istop = 0

  IF (((DV_MIXT_VVEL - DV_MIXT_UVEL)/=1).OR. &
      ((DV_MIXT_WVEL - DV_MIXT_VVEL)/=1).OR. &
      ((DV_MIXT_TEMP - DV_MIXT_WVEL)/=1)) istop = 1
  IF (((GR_MIXT_VX - GR_MIXT_UX)/=1).OR.((GR_MIXT_WX - GR_MIXT_VX)/=1).OR. &
      ((GR_MIXT_TX - GR_MIXT_WX)/=1).OR.((GR_MIXT_UY - GR_MIXT_TX)/=1).OR. &
      ((GR_MIXT_VY - GR_MIXT_UY)/=1).OR.((GR_MIXT_WY - GR_MIXT_VY)/=1).OR. &
      ((GR_MIXT_TY - GR_MIXT_WY)/=1).OR.((GR_MIXT_UZ - GR_MIXT_TY)/=1).OR. &
      ((GR_MIXT_VZ - GR_MIXT_UZ)/=1).OR.((GR_MIXT_WZ - GR_MIXT_VZ)/=1).OR. &
      ((GR_MIXT_TZ - GR_MIXT_WZ)/=1) ) istop = 1

  IF (istop == 1) CALL ErrorStop( global,ERR_GRAD_INDEX,__LINE__ )

! user input check ------------------------------------------------------------

! check mapping to processors

  IF (global%nRegionsProc < 0) THEN         ! manual region to processor mapping
    CALL ErrorStop( global,ERR_NO_PROCMAP,__LINE__ )
  ELSE IF (global%nRegionsProc > 0) THEN    ! number regions per processor given
    IF (global%nProcAlloc*global%nRegionsProc /= global%nRegions) &
      CALL ErrorStop( global,ERR_NO_PROCMATCH,__LINE__ )
  ENDIF

! dual time-stepping

#ifdef GENX
  IF (global%flowType == FLOW_UNSTEADY .AND. &
      global%solverType == SOLV_IMPLICIT) THEN
    IF (global%predictSol.eqv..true.) THEN
      CALL ErrorStop( global,ERR_ILLEGAL_VALUE,__LINE__, &
                      'for Rocstar, safer to use PREDICTSOL=0' )
    ENDIF
  ENDIF
#endif

! grid motion

  IF (global%moveGridScheme == MOVEGRID_FRAME .OR. &
      global%moveGridScheme == MOVEGRID_FOMS) THEN
    IF (global%nRegions < 8) THEN
      CALL ErrorStop( global,ERR_ILLEGAL_VALUE,__LINE__, &
      'Selected grid motion type only applicable for cases with >= 8 regions' )
    ENDIF
    IF (global%moveGridNbour > (26 + (global%nRegions-8)*2)) THEN
      CALL ErrorStop( global,ERR_ILLEGAL_VALUE,__LINE__, &
      'Maximum number of NEIGHBOR in type 2 grid motion is 26+2*(nReg-8)' )
    ENDIF
    IF (global%moveGridNbour < 2) THEN
      CALL ErrorStop( global,ERR_ILLEGAL_VALUE,__LINE__, &
      'Minimum number of NEIGHBOR in type 2 grid motion is 2' )
    ENDIF
    IF (global%moveGridNsmatch < 2) THEN
      CALL ErrorStop( global,ERR_ILLEGAL_VALUE,__LINE__, &
      'Minimum number NSURFMATCH in type 2 grid motion is 2' )
    ENDIF
  ENDIF

! random number generator
  
  IF ( global%randSeedType < RAND_SEED_TYPE_FIXED .OR. &
       global%randSeedType > RAND_SEED_TYPE_CLOCK      ) THEN
    CALL ErrorStop(global,ERR_RAND_SEED_TYPE_INVALID,__LINE__)
  END IF ! global%randSeedType

! aerodynamic coefficients

  IF ( global%aeroCoeffs==ACTIVE .AND. global%nProbes==0) THEN
    CALL ErrorStop(global,ERR_DEPENDENT_INPUT,__LINE__, &
      'activate Probe section (NUMBER>0) when AEROCOEFFS/=0 in Force section')
  ENDIF

! region dependent ------------------------------------------------------------

  DO iReg=1,global%nRegions

! - check numeric

    IF (regions(iReg)%mixtInput%spaceDiscr/=DISCR_CEN_SCAL .AND. &
        regions(iReg)%mixtInput%spaceDiscr/= DISCR_UPW_ROE) THEN
      CALL ErrorStop(global,ERR_UNKNOWN_OPTION,__LINE__, &
                     'only central and Roe-upwind implemented')
    ENDIF

    IF (regions(iReg)%mixtInput%spaceOrder < DISCR_ORDER_1 .OR. &
        regions(iReg)%mixtInput%spaceOrder > DISCR_ORDER_2) THEN
      CALL ErrorStop(global,ERR_UNKNOWN_OPTION,__LINE__, &
                     'only 1st and 2nd order space discr. implemented')
    ENDIF

! - check dummy cells

    IF (regions(iReg)%mixtInput%spaceOrder <= DISCR_ORDER_2 .AND. &
        regions(iReg)%nDumCells < 2) THEN          ! <= 2nd-order
      CALL ErrorStop( global,ERR_GRID_DUMCELLS,__LINE__, &
                      'number of dummy layers should be at least 2'  )
    ENDIF

! - check consistency of viscosity models

    computeTv = regions(iReg)%mixtInput%computeTv
    viscModel = regions(iReg)%mixtInput%viscModel
    IF ((computeTv.eqv..true.) .AND. &
       (viscModel <  VISC_SUTHR .OR. viscModel > VISC_ANTIB )) &
      CALL ErrorStop( global,ERR_UNKNOWN_VISCMODEL,__LINE__ )

! - check consistency of gridmotion

    IF (regions(iReg)%mixtInput%moveGrid.eqv..true.) THEN
      IF (MOD( global%nRegions,global%nProcAlloc )/=0) THEN
        CALL ErrorStop( global,ERR_NO_PROCMATCH,__LINE__, &
        'Sorry, currently gridmotion only allows nRegs multiple of nProcs' )
      ENDIF
    ENDIF

! - check if mpi-tag-shifts are sufficiently large

    IF (regions(iReg)%localNumber > BASE_TAG_SHIFT) THEN
      CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__, &
           'increase TAG_SHIFTs parameters or reduce maximum nRegions/proc.' )
    ENDIF
  ENDDO

! tag shift of physical modules

  IF (BASE_TAG_SHIFT < 0 .OR. &
      MOD( PEUL_TAG_SHIFT,BASE_TAG_SHIFT ) /= 1 .OR. &
      MOD( TURB_TAG_SHIFT,BASE_TAG_SHIFT ) /= 1 .OR. &
      MOD( RADI_TAG_SHIFT,BASE_TAG_SHIFT ) /= 1 .OR. &
      MOD( PLAG_TAG_SHIFT,BASE_TAG_SHIFT ) /= 1) THEN
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__, &
         'TAG_SHIFTs should be n*BASE_TAG_SHIFT+1, n>0, BASE_TAG_SHIFT>0' )
  ENDIF

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_CheckUserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CheckUserInput.F90,v $
! Revision 1.16  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.15  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.14  2008/10/23 18:20:57  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.13  2006/03/25 01:06:05  wasistho
! changed errcode to ERR_DEPENDENT_INPUT
!
! Revision 1.12  2006/03/25 01:02:04  wasistho
! modified test to aeroCoeffs==ACTIVE
!
! Revision 1.11  2006/03/24 23:31:48  wasistho
! added consistency check between aeroCoeffs and probes
!
! Revision 1.10  2005/12/01 17:12:49  fnajjar
! Added check for randSeedType
!
! Revision 1.9  2005/11/08 21:34:24  wasistho
! forgot line continuation sign &
!
! Revision 1.8  2005/11/08 21:32:13  wasistho
! added check for movegrid_foms
!
! Revision 1.7  2005/08/19 00:02:03  wasistho
! made gridmotion max neighbour dependent on nReg
!
! Revision 1.6  2005/08/18 19:48:52  wasistho
! added check moveGridNsmatch
!
! Revision 1.5  2005/06/25 03:57:10  wasistho
! check nRegions multiple of nProcs for moving grid
!
! Revision 1.4  2005/06/23 01:45:02  wasistho
! added checks specific for type 2 grid motion
!
! Revision 1.3  2005/06/20 16:55:42  wasistho
! added check nprocs vs nregions for type 2 grid motion
!
! Revision 1.2  2004/12/15 09:03:44  wasistho
! added trap for predictsol (dual tst) within genx
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.21  2004/10/01 18:08:18  wasistho
! trapped unimplemented space discretization options
!
! Revision 1.20  2004/04/08 03:16:42  wasistho
! nDummyCells in Rocflo read from INITFLOW section
!
! Revision 1.19  2004/03/06 02:34:06  wasistho
! added PLAG tag shift
!
! Revision 1.18  2004/03/05 22:09:02  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.17  2004/03/05 21:09:03  wasistho
! added check consistency of mpi-tag-shifts
!
! Revision 1.16  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.15  2003/04/10 23:30:51  fnajjar
! Checking consistency of viscosity models
!
! Revision 1.14  2003/02/11 22:49:53  jferry
! Re-worked BC and TBC input routines to add multi-physics capability
!
! Revision 1.1  2003/02/11 22:30:21  jferry
! Re-worked BC and TBC input routines to add multi-physics capability
!
! Revision 1.13  2002/09/26 19:29:40  jferry
! removed TBC check
!
! Revision 1.11  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.10  2002/09/02 23:35:22  wasistho
! Add grad index check and removed TURB compilation check
!
! Revision 1.9  2002/08/30 01:47:58  jblazek
! Added support for moving grids.
!
! Revision 1.8  2002/08/24 03:15:30  wasistho
! modify TURB error msg
!
! Revision 1.7  2002/08/18 02:30:47  wasistho
! Added check TURB module activation
!
! Revision 1.6  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.5  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.3  2002/01/31 20:56:30  jblazek
! Added basic boundary conditions.
!
! Revision 1.2  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.1  2001/12/08 00:18:41  jblazek
! Added routines to read BC input file.
!
!******************************************************************************







