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
      program summary
!
!  Summarizes timing data in GENXTimingDatannnn.txt files
!
!  Finds the longest wall clock time over all processes
!  in a given run.  If the longest time is shorter in the
!  2nd step than in the 3rd, it uses the 2nd step.
!
!  Use executable summary.x with timing_script.
!
!  Written by Robert Fiedler, revised 2/7/02.
!
      integer m,n,np,nmax,mmax
      parameter(nmax=10000,mmax=3)
      double precision times(nmax,mmax),tmax2,tmax3,tmax
!      character*55 junk
!      character*56 junk
      character*33 junk
      character*33 junk1
!
      times = 0.
!
!  Make a system call to grep to get list of timings
!  If "system" is not available on the machine, comment this
!  out and do it by hand.
!
!      call system('grep "PrecCorr loop time" *.txt > summary')
!

!
!  This list has mmax entries per processor.  Read in the times.
!
      open(10,file='summary')
      do n=1,nmax  ! Loop over timing files (processors)
!        if (n == 1) then
!
! See if it's in the format for 1 processor
!
!          do m=1,mmax
!            read(10,'(a,d25.16)',iostat=ios) junk1,times(n,m)
!            if (ios /= 0) exit
!          end do
!          if (ios /= 0) rewind(10)
!        endif
!
! Try format for more than 1 processor
!
        do m=1,mmax
          read(10,'(a,d25.16)',iostat=ios) junk,times(n,m)
          if (ios /= 0) exit
        end do
        if (ios /= 0) exit
      end do
      np = max(1,n-1)
!
!  Find the longest time for the 2nd and 3rd PreCorr loops 
!
      tmax2 = 0.
      tmax3 = 0.
      do n=1,np
        if (times(n,2) > tmax2) tmax2 = times(n,2)
        if (times(n,3) > tmax3) tmax3 = times(n,3)
      end do
      tmax = tmax2
      if (tmax3 .ne. 0.) tmax = min(tmax2,tmax3)
      write(*,10) np, tmax
 10   format(i4, 2x, 1pd16.5)
!      print *,np,tmax
!
      end






