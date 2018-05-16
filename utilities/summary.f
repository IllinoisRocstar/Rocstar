      program summary
!
!  Summarizes timing data in GEN1TimingData*.txt files
!  Written by Robert Fiedler, 5/14/01.
!
      integer m,n,np,nmax
      parameter(nmax=512)
      double precision times(nmax,6),tmax4,tmax6,tmax
!      character*55 junk
      character*56 junk
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
!  This list has 6 entries per processor.  Read in the times.
!
      open(10,file='summary')
      do n=1,nmax  ! Loop over timing files (processors)
        if (n == 1) then
!
! See if it's in the format for 1 processor
!
          do m=1,6
            read(10,'(a,d25.16)',iostat=ios) junk1,times(n,m)
            if (ios /= 0) exit
          end do
          if (ios /= 0) rewind(10)
        endif
!
! Try format for more than 1 processor
!
        do m=1,6
          read(10,'(a,d25.16)',iostat=ios) junk,times(n,m)
          if (ios /= 0) exit
        end do
        if (ios /= 0) exit
      end do
      np = max(1,n-1)
!
!  Find the longest time for the 4th and 6th PreCorr loops
!
      tmax4 = 0.
      tmax6 = 0.
      do n=1,np
        if (times(n,4) > tmax4) tmax4 = times(n,4)
        if (times(n,6) > tmax6) tmax6 = times(n,6)
      end do
      tmax = tmax4
      if (tmax6 .ne. 0.) tmax = min(tmax4,tmax6)
      write(*,10) np, tmax
 10   format(i4, 2x, 1pd16.5)
!      print *,np,tmax
!
      end
