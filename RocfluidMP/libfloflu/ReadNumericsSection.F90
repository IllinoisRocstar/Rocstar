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
! Purpose: read in user input related to numerical procedure.
!
! Description: none.
!
! Input: user input file.
!
! Output: regions = numerical switches and parameters.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ReadNumericsSection.F90,v 1.13 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadNumericsSection( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO  
  USE ModInterfaces, ONLY : ReadRegionSection
#endif
#ifdef RFLU
  USE ModInterfaces, ONLY : ReadSection
#endif
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  INTEGER :: nVals

#ifdef RFLO
  INTEGER :: brbeg, brend
  INTEGER, PARAMETER :: NVALS_MAX = 11
#endif
#ifdef RFLU
  INTEGER :: iReg
  INTEGER, PARAMETER :: NVALS_MAX = 16
#endif

  CHARACTER(10) :: keys(NVALS_MAX)

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL) :: vals(NVALS_MAX)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'ReadNumericsSection',&
  'ReadNumericsSection.F90' )

! specify keywords and search for them

  nVals = NVALS_MAX

#ifdef RFLO
  keys( 1) = 'CFL'
  keys( 2) = 'SMOOCF'
  keys( 3) = 'DISCR'
  keys( 4) = 'ORDER'
  keys( 5) = 'LIMFAC'
  keys( 6) = 'ENTROPY'
  keys( 7) = 'K2'
  keys( 8) = '1/K4'
  keys( 9) = 'PSWTYPE'
  keys(10) = 'PSWOMEGA'
  keys(11) = 'FEAVERAG'

  CALL ReadRegionSection( global,IF_INPUT,nVals,keys(1:nVals),vals(1:nVals), & 
                          brbeg,brend,defined(1:nVals) )

  IF (defined(3).eqv..true.) THEN
    SELECT CASE(INT(vals(3)+0.5_RFREAL))
    CASE (0)
      regions(brbeg:brend)%mixtInput%spaceDiscr = DISCR_CEN_SCAL
    CASE (1)
      regions(brbeg:brend)%mixtInput%spaceDiscr = DISCR_UPW_ROE
    CASE (2)
      regions(brbeg:brend)%mixtInput%spaceDiscr = DISCR_UPW_MAPS
    END SELECT
  ENDIF
  IF (defined(4).eqv..true.) THEN
    SELECT CASE(INT(vals(4)+0.5_RFREAL))
    CASE (1)
      regions(brbeg:brend)%mixtInput%spaceOrder = DISCR_ORDER_1
    CASE (2)
      regions(brbeg:brend)%mixtInput%spaceOrder = DISCR_ORDER_2
    CASE (4)
      regions(brbeg:brend)%mixtInput%spaceOrder = DISCR_ORDER_4
    END SELECT
  ENDIF
  IF (defined(8).eqv..true.) THEN
    IF (vals(8) > 1.E-10_RFREAL) THEN
      regions(brbeg:brend)%mixtInput%vis4 = 1._RFREAL/vals(8)
    ELSEIF (vals(8) > 0._RFREAL .AND. vals(8) <= 1.E-10_RFREAL) THEN
      regions(brbeg:brend)%mixtInput%vis4 = 1.E+10_RFREAL
    ELSEIF (vals(8) <= 0._RFREAL ) THEN
      regions(brbeg:brend)%mixtInput%vis4 = 0.0_RFREAL
    ENDIF
  ENDIF
  IF (defined(9).eqv..true.) THEN
    IF (vals(9) < 0.5_RFREAL) THEN
      regions(brbeg:brend)%mixtInput%pSwitchType = PSWITCH_STD
    ELSE
      regions(brbeg:brend)%mixtInput%pSwitchType = PSWITCH_TVD
    ENDIF
  ENDIF
  IF (defined( 1).eqv..true.) regions(brbeg:brend)%mixtInput%cfl          = ABS(vals(1))
  IF (defined( 2).eqv..true.) regions(brbeg:brend)%mixtInput%smoocf       = vals(2)
  IF (defined( 5).eqv..true.) regions(brbeg:brend)%mixtInput%limfac       = ABS(vals(5))
  IF (defined( 6).eqv..true.) regions(brbeg:brend)%mixtInput%epsentr      = ABS(vals(6))
  IF (defined( 7).eqv..true.) regions(brbeg:brend)%mixtInput%vis2         = ABS(vals(7))
  IF (defined(10).eqv..true.) regions(brbeg:brend)%mixtInput%pSwitchOmega = ABS(vals(10))
  IF (defined(11).eqv..true.) THEN
    IF (vals(11) < 0.5_RFREAL) THEN
      regions(brbeg:brend)%mixtInput%faceEdgeAvg = FE_AVG_UNIFORM
    ELSE
      regions(brbeg:brend)%mixtInput%faceEdgeAvg = FE_AVG_LINEAR
    ENDIF
  ENDIF
#endif

#ifdef RFLU
  keys( 1) = 'CFL'
  keys( 2) = 'DISCR'
  keys( 3) = 'ORDER'
  keys( 4) = 'ENTROPY'
  keys( 5) = 'DISSFACT'
  keys( 6) = 'RECONST'
  keys( 7) = 'DIMENS'
  keys( 8) = 'CRECONSTC'
  keys( 9) = 'CRECONSTF'
  keys(10) = 'CRECONSTCW'
  keys(11) = 'CRECONSTFW'
  keys(12) = 'TOLERICT' 
  keys(13) = 'SDIMENSC'
  keys(14) = 'SDIMENSF'
  keys(15) = 'SDIMENSBF'
  keys(16) = 'ORDERBF'
  
  CALL ReadSection( global,IF_INPUT,nVals,keys(1:nVals),vals(1:nVals), & 
                    defined(1:nVals) ) 
  
  IF ( defined(1) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%cfl = ABS(vals(1))
    END DO ! iReg
  END IF ! defined
  
  IF ( defined(2) .EQV. .TRUE. ) THEN 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%spaceDiscr = NINT(vals(2))
    END DO ! iReg
  END IF ! defined 
  
  IF ( defined(3) .EQV. .TRUE. ) THEN 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%spaceOrder = NINT(vals(3)) 
    END DO ! iReg
  END IF ! defined 
  
  IF ( defined(4) .EQV. .TRUE. ) THEN 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%epsentr = ABS(vals(4))
    END DO ! iReg
  END IF ! defined 
  
  IF ( defined(5) .EQV. .TRUE. ) THEN 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%dissFact = ABS(vals(5))
    END DO ! iReg
  END IF ! defined

  IF ( defined(6) .EQV. .TRUE. ) THEN 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%reconst = NINT(vals(6))
    END DO ! iReg
  END IF ! defined
  
  IF ( defined(7) .EQV. .TRUE. ) THEN 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%dimens = NINT(vals(7))
    END DO ! iReg    
  END IF ! defined 
  
  IF ( defined(8) .EQV. .TRUE. ) THEN 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%cReconstCells = NINT(vals(8))
    END DO ! iReg    
  END IF ! defined          

  IF ( defined(9) .EQV. .TRUE. ) THEN 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%cReconstFaces = NINT(vals(9))
    END DO ! iReg    
  END IF ! defined          
  
  IF ( defined(10) .EQV. .TRUE. ) THEN 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%cReconstCellsWeight = ABS(vals(10))
    END DO ! iReg    
  END IF ! defined          

  IF ( defined(11) .EQV. .TRUE. ) THEN 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%cReconstFacesWeight = ABS(vals(11))
    END DO ! iReg    
  END IF ! defined
  
  IF ( defined(12) .EQV. .TRUE. ) THEN 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%tolerICT = ABS(vals(12))
    END DO ! iReg
  END IF ! defined 
  
  IF ( defined(13) .EQV. .TRUE. ) THEN 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%stencilDimensCells = ABS(vals(13))
    END DO ! iReg
  END IF ! defined 
  
  IF ( defined(14) .EQV. .TRUE. ) THEN 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%stencilDimensFaces = ABS(vals(14))
    END DO ! iReg
  END IF ! defined
  
  IF ( defined(15) .EQV. .TRUE. ) THEN 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%stencilDimensBFaces = ABS(vals(15))
    END DO ! iReg
  END IF ! defined                 
  
  IF ( defined(16) .EQV. .TRUE. ) THEN 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%spaceOrderBFaces = ABS(vals(16))
    END DO ! iReg
  END IF ! defined                 
#endif 

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadNumericsSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadNumericsSection.F90,v $
! Revision 1.13  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.12  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.11  2008/10/23 18:20:55  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.10  2006/10/20 21:20:48  mparmar
! Renamed keyword SORDERBF to ORDERBF
!
! Revision 1.9  2006/08/19 15:38:29  mparmar
! Added reading of mixtInput%spaceOrderBFaces
!
! Revision 1.8  2006/04/07 14:37:44  haselbac
! Added new stencilDimens parameters
!
! Revision 1.7  2006/01/06 22:04:19  haselbac
! Added SDIMENS
!
! Revision 1.6  2005/12/25 15:20:21  haselbac
! Added constrained reconstuction input variables, minor bug fixes (NINT instead of ABS)
!
! Revision 1.5  2005/12/24 21:23:35  haselbac
! Added reading of ICT tolerance
!
! Revision 1.4  2005/10/27 18:53:24  haselbac
! Added input parameter for constrained reconstruction
!
! Revision 1.3  2005/07/11 19:22:20  mparmar
! Added reading of RECONST keyword
!
! Revision 1.2  2005/03/09 14:51:41  haselbac
! Added DIMENS variable
!
! Revision 1.1  2004/12/01 16:50:37  haselbac
! Initial revision after changing case
!
! Revision 1.19  2004/09/02 02:34:11  wasistho
! added face-edge averaging input-option parameter in Rocflo
!
! Revision 1.18  2004/07/08 02:18:37  haselbac
! Added dissFact, cosmetics
!
! Revision 1.17  2003/12/05 16:54:34  haselbac
! Fixed bug in DO loop ranges
!
! Revision 1.16  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.13  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.12  2003/04/09 20:42:54  wasistho
! enable k4 = 0.0
!
! Revision 1.11  2002/09/09 14:03:46  haselbac
! mixtInput now under regions
!
! Revision 1.10  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.9  2002/07/25 00:38:00  jblazek
! Option for TVD type pressure switch.
!
! Revision 1.8  2002/05/04 16:21:58  haselbac
! Bug fix: spaceOrder was set twice
!
! Revision 1.7  2002/03/26 19:06:47  haselbac
! Added ROCFLU functionality
!
! Revision 1.6  2002/02/27 18:38:19  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
! Revision 1.5  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.3  2002/01/11 17:18:31  jblazek
! Updated description of I/O variables.
!
! Revision 1.2  2001/12/22 00:09:38  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1  2001/12/07 16:54:31  jblazek
! Added files to read user input.
!
!******************************************************************************







