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
! Purpose: LES inter processor communication. 
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = current region number
!
! Output: global% = new global parameters after global communication
!
! Notes: This procedure is always done regardless turbulence model selected.
!        The input/output global variables are always defined. 
!        If not needed, the result is simply not used.
!
!******************************************************************************
!
! $Id: TURB_LesCommunication.F90,v 1.6 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_LesCommunication( regions ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  
! ... local variables
  TYPE(t_global), POINTER :: global
  REAL(RFREAL) :: esg1Glo, esg4Glo

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'TURB_LesCommunication',&
  'TURB_LesCommunication.F90' )
  
! global sum of LES energy sgs model 1 and 4 ----------------------------------

#ifdef MPI
  CALL MPI_ALLREDUCE( global%esg1sum, esg1Glo, 1, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, global%mpiComm, global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
  CALL MPI_ALLREDUCE( global%esg4sum, esg4Glo, 1, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, global%mpiComm, global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
  global%esg1sum = esg1Glo
  global%esg4sum = esg4Glo
#endif

!#ifdef RFLU
!#ifdef CHARM
!  CALL FEM_REDUCE( global%fieldFlagTurbReal,global%esg1sum,esg1Glo,FEM_SUM )
!  CALL FEM_REDUCE( global%fieldFlagTurbReal,global%esg4sum,esg4Glo,FEM_SUM )
!  global%esg1sum = esg1Glo
!  global%esg4sum = esg4Glo
!#endif
!#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_LesCommunication

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_LesCommunication.F90,v $
! Revision 1.6  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/12/29 19:45:51  wasistho
! removed CHARM stuff
!
! Revision 1.3  2005/04/15 15:07:35  haselbac
! Removed Charm/FEM stuff
!
! Revision 1.2  2004/03/19 02:49:28  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.1  2004/02/26 21:32:51  wasistho
! install TURB_lesCommunication
!
!
!******************************************************************************







