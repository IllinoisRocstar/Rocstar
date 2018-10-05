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
!   SUBROUTINE  : APN_input_0d                                                  
!                                                                              
!   This subroutine read inputs for burning rate model Rocburn_1D_APN
!
!   Authors          : 
!
!   Creation Date    :  Sep. 10, 2002
!   
!   Modifications    :
!
!    No.     Date         Authors       Description
!
!                                                                              
! ---------------------------------------------------------------------------  
!                                                                              
!                                                                              
!   arguments   :                                                              
!                                                                              
!      G_APN     : Global variables for Rocburn_1D_APN
!      Indir     : directory for input data file
!                                                                              
! ---------------------------------------------------------------------------  
!
  SUBROUTINE APN_input_0d(G_APN, Indir)

    USE M_Rocburn_APN_Global_Data

    IMPLICIT NONE
    INCLUDE 'mpif.h'

!
!   Global data for Rocburn_1D_APN passed as a pointer
!

    TYPE (G_BURN_1D), POINTER :: G_APN
!
!   arguments
!
    CHARACTER(*), INTENT(IN)  :: Indir
    CHARACTER(*), PARAMETER :: ControlFile = "RocburnAPNControl.txt"

!
!   
!
!
! ----------------------------------------------------------------
!   local variables
    CHARACTER(LEN=80) :: Infile

    INTEGER   :: ir, ioerr
    INTEGER   :: mat

!
!       read propellant thermophysical properties
!

    ir = 10

    if (Indir(LEN_TRIM(Indir):LEN_TRIM(Indir)) == '/') then
       Infile= TRIM(Indir) // ControlFile
    else
       Infile= TRIM(Indir) // '/' // ControlFile
    endif

    OPEN (unit=ir,file=Infile,status='old')

!RAF
!RAF Extend for multiple propellants, but keep it backward compatible
!RAF Specify same value for nxmax and To for all materials.
!RAF

    DO mat = 1, MATMAX
      READ(ir,*,IOSTAT=ioerr)   G_APN%a_p(mat)
      IF (ioerr /= 0) THEN
        G_APN%nmat = mat - 1
        EXIT
      ENDIF
      G_APN%nmat = mat
      READ(ir,*)   G_APN%n_p(mat)
      READ(ir,*)   G_APN%nxmax
      READ(ir,*)   G_APN%Tf_adiabatic(mat)
      READ(ir,*)   G_APN%To
      READ(ir,*,IOSTAT=ioerr) G_APN%xmax(mat)
      IF (ioerr /= 0) THEN
        G_APN%xmax(mat) = 1.0e+15
      ENDIF
      IF (ioerr /= 0) EXIT
    END DO

    READ(ir,*,IOSTAT=ioerr) G_APN%verbosity
    IF(ioerr /= 0) THEN
      G_APN%verbosity = 1
    ENDIF

    CLOSE(ir)

    IF(G_APN%rank .eq. 0 .AND. G_APN%verbosity .gt. 0) THEN
       WRITE(6,'(A)') 'RocburnAPN: *********** Using APN Burn Model **************'
    ENDIF
    IF(G_APN%rank .eq. 0 .AND. G_APN%verbosity .gt. 1) THEN
       WRITE(6,'(A,i3,A)') 'RocburnAPN:  Found a total of ',G_APN%nmat,' materials'
       DO mat = 1, G_APN%nmat
          WRITE(6,'(A,i3,A,f12.4)') 'RocburnAPN:  a_p(',mat,')  =',G_APN%a_p(mat)
          WRITE(6,'(A,i3,A,f12.4)') 'RocburnAPN:  n_p(',mat,')  =',G_APN%n_p(mat)
          WRITE(6,'(A,i3)') 'RocburnAPN:  nxmax= ',G_APN%nxmax
          WRITE(6,'(A,i3,A,f12.4)') 'RocburnAPN:  Tf_adiabatic(',mat,') = ',&
                                     G_APN%Tf_adiabatic(mat)
          WRITE(6,'(A,f12.4)') 'RocburnAPN:  To= ',G_APN%To
          WRITE(6,'(A,i3,A,f12.4)') 'RocburnAPN:  xmax(',mat,')  =',&
                                    G_APN%xmax(mat)
       END DO
    END IF

    IF (G_APN%verbosity.gt.2) PRINT *,'RocburnAPN: rank=',G_APN%rank,       &
            ' done input ', ControlFile

    RETURN 

  END SUBROUTINE APN_input_0d






