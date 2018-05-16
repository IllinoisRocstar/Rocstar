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
! Purpose: define the properties of materials used throughout the code
!
! Description: none.
!
! Notes:
!
!   this module will grow as the modeling for materials becomes more
!   sophisticated
!
!******************************************************************************
!
! $Id: ModMaterials.F90,v 1.7 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

MODULE ModMaterials

  USE ModDataTypes
  IMPLICIT NONE

! data types ------------------------------------------------------------------

  TYPE t_material
    CHARACTER(CHRLEN) :: name  ! Name of material
    INTEGER           :: phase ! 1 = Gas, 2 = Liquid, 3 = Solid
    INTEGER           :: index ! Index of material in global%materials(:)
    REAL(RFREAL)      :: molw  ! Molecular weight
    REAL(RFREAL)      :: dens  ! Density
    REAL(RFREAL)      :: spht  ! Specific heat
    REAL(RFREAL)      :: surftens  ! Surface tension
    REAL(RFREAL)      :: Tboil ! Boiling (or condensation) temperature
    REAL(RFREAL)      :: Tmelt ! Melting (or freezing) temperature
  END TYPE t_material

END MODULE ModMaterials

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModMaterials.F90,v $
! Revision 1.7  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.4  2004/03/02 21:49:45  jferry
! Added melting and boiling point to material definitions
!
! Revision 1.3  2003/09/13 20:15:55  fnajjar
! Added surface tension to Materials datastructure
!
! Revision 1.2  2003/03/12 18:06:13  jferry
! added index component to t_material
!
! Revision 1.1  2003/03/11 15:53:25  jferry
! Created data type for material properties
!
!******************************************************************************






