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
!     ------------------------------------------------------------------------

  PROGRAM BURN_INITIALIZE

    USE M_ROCBURN_1D_PY
    IMPLICIT NONE

    INTEGER :: nxmax,i
    REAL(DBL) :: To_read,P, PATM, Tflame,dzero,time
    REAL(DBL) :: Toa, fr,rhoc,rhoc_cgs,rb, dt,p_coor(3),fact
    REAL(DBL),ALLOCATABLE :: Tn(:),Tnp1(:)
    INTEGER :: bflag,ntime_steps
    CHARACTER (LEN=80) PASS_DIR
    TYPE(parameter_structure),POINTER :: pass_in_out

    
     PASS_DIR = './'

!----READ INPUT FILE AND GENERATE GRID
     CALL BURN_INIT_0D(pass_in_out,1,TRIM(PASS_DIR), nxmax, To_read)

!----ALLOCATION
     ALLOCATE(Tn(nxmax))
     ALLOCATE(Tnp1(nxmax))

!----USE a burning interface
     bflag =1   
     dzero = zero

!----SET Patm
     Patm = 30.0
     P = Patm /pa2atm

!----SET rhoc, in genx is passed by Rocfrac
     rhoc_cgs = 1.7026
     rhoc = rhoc*gcc2kgmc

!----INITIALIZE temperature profile
     CALL BURN_INIT_1D(pass_in_out, bflag, &
                    P, To_read, rhoc, p_coor, rb, &
                    Toa, fr, Tn(:), Tflame)
  
     time = 0 
     ntime_steps = 100000
     dt = 1.0e-6

!----RESET P if starting from steady conditions
     if(bflag ==1) then
        fact = 2.0d0
        P = P*fact
     endif

!----MAIN TIME LOOP
     do i = 1,ntime_steps
        time = time+dt
        CALL burn_get_burning_rate1d(pass_in_out, dt, P, &
             To_read, Tn(:),   &
             dzero, dzero,     &
             dzero, dzero,     &
             rhoc, dzero,      &
             rb, dzero,        &
             bflag, Tnp1(:), Tflame)
        Tn = Tnp1
        print*,i,time,Tn(1)
     enddo
 
  STOP
END PROGRAM BURN_INITIALIZE
!*****************************************************************************






