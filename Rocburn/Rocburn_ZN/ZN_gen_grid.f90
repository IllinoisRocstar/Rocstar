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
!
! ---------------------------------------------------------------------------
!
!   SUBROUTINE  : ZN_gen_grid
!
!   purpose     : grid generation  
!
!   Authors          :  P. Loner
!
!   Creation Date    :  Sep. 3, 2002
!
!   Modifications    :
!
!    No.     Date         Authors       Description
!     1    09/03/02       K. Tang       Global variables G_ZN
!
! ---------------------------------------------------------------------------
!
!
!   arguments   :
!
!      G_ZN     : Global variables for Rocburn_1D_ZN
!      gridtype : = 1  expoential grid 
!                 = 2  boundary layer grid
!      numx     : maximum number of grid points
!      beta     : stretch parameter 
!      xmax     : the maximum dimension in x-coordinate
!      x        : the x-coordinate array
!      z        : the z-coordinate arary
!      zx       : the first derivative of z
!      zxx      : the second derivative of z 
!      delz     : the z spacing
!                                 
! ---------------------------------------------------------------------------
!
  SUBROUTINE ZN_gen_grid(G_ZN, gridtype, numx, x, z, zx, zxx)

    USE M_Rocburn_ZN_Global_Data

    IMPLICIT NONE
    INCLUDE 'mpif.h'

    TYPE (G_BURN_1D), POINTER :: G_ZN
    INTEGER, INTENT(IN) :: gridtype, numx
    REAL(DBL), INTENT (OUT) :: x(:), z(:), zx(:), zxx(:)
!
!
! ---------------------------------------------------------------------------
!   local variables

    REAL(DBL) :: bpm,bp1,bm1,czx,czxx,term,term2
    INTEGER :: i

! ---------------------------------------------------------------------------
    IF(G_ZN%xmax >= 0.0) then
         PRINT *,' ROCBURN_ZN: rank=',G_ZN%rank,' xmax >= 0'
         PRINT *,' ROCBURN_ZN: job aborted!'
         STOP
    END IF

    IF(numx.gt.(G_ZN%nxmax)) then
      PRINT *,' ROCBURN_ZN: rankk=',G_ZN%rank,' numx=',numx, &
              ' nxmax=',G_ZN%nxmax
      PRINT *,' ROCBURN_ZN: nx > maximum number of', &
              ' allowed spatial nodes'
      PRINT *,' ROCBURN_ZN: job aborted!'
      STOP
    END IF

    G_ZN%delz = 1.0/dble(numx-1)

    IF(G_ZN%rank ==0) THEN
      PRINT *,' ROCBURN_ZN: delz=',G_ZN%delz
      PRINT *,' ROCBURN_ZN: nxmax=',G_ZN%nxmax
    END IF

    IF (gridtype == 1) THEN
!
!     exponential grid
!

      DO i=1,numx
        z(i) = dble(i-1)*G_ZN%delz
        x(i) = log(1.-z(i))/G_ZN%beta
        IF (i.eq.numx) THEN
            x(i) = -10.
        END IF
        zx(i)  = -G_ZN%beta*exp(G_ZN%beta*x(i))
        zxx(i) = -G_ZN%beta**2*exp(G_ZN%beta*x(i))
      END DO ! i

    ELSE
!
!     Anderson, Tanehill and Pletcher boundary layer grid control
!

      bp1=G_ZN%beta+1.0
      bm1=G_ZN%beta-1.0
      bpm=bp1/bm1
      czx=2.0*G_ZN%beta/log(bpm)/G_ZN%xmax
      czxx=-2.0/G_ZN%xmax*czx

      DO i=1,numx
        z(i)=dble(i-1)*G_ZN%delz
        term=bpm**(1.0-z(i))
        x(i)=G_ZN%xmax*(bp1-bm1*term)/(term+1.0)
        term2=G_ZN%beta**2-(1.0-x(i)/G_ZN%xmax)**2
        zx(i)=czx/term2
        zxx(i)= czxx*(1.0-x(i)/G_ZN%xmax)/term2**2
      END DO ! i

    END IF

    IF(G_ZN%rank == 0) THEN
       WRITE (*,*)
       WRITE (*,*)' ROCBURN_ZN: Smallest delta x = ',x(1)-x(2)
       WRITE (*,*)' ROCBURN_ZN: Largest delta x  = ',x(numx-1)-x(numx)
       WRITE (*,*)
    END IF

    RETURN

  END SUBROUTINE ZN_gen_grid






