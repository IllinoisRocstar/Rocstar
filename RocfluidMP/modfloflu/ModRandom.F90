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
! Purpose: Portable random number generator
!
! Description: This module provides the following routines:
!
!   RandomSeed: seeds the random number generator
!
!   RandUniform:    generates an array of uniformly distributed random numbers
!                   on the interval [0,1)
!
!   RandNormal:     generates an array of normally distributed random numbers
!
!   RandLogNormal:  generates an array of log-normally distributed numbers
!
!   RandImposedPDF: generates an array of number distributed according to a dscrete PDF
!
!   Rand1Uniform:    function form of RandUniform -- returns a single number
!
!   Rand1Normal:     function form of RandNormal
!
!   Rand1LogNormal:  function form of RandLogNormal
!
!   Rand1ImposedPDF: function form of RandImposedPDF
!
! Notes:
!
!   This is an implementation of the Mersenne Twister random number
!   generator, MT19937, based on the Mersenne Prime 2^19937 - 1.
!
!   It is (somewhat painfully) coded such that it gives the same results
!   whether integers are 32- or 64-bit, and these results are identical
!   to the standard, C version of the Mersenne Twister.
!
!   It outputs random numbers with 53-bit precision.
!
!******************************************************************************
!
! $Id: ModRandom.F90,v 1.6 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

MODULE ModRandom

  USE ModDataTypes
  IMPLICIT NONE

! Period parameter
  INTEGER, PARAMETER :: RAND_BUFFER_SIZE = 624
  INTEGER, PARAMETER :: RAND_MTI_INDEX   = RAND_BUFFER_SIZE
  INTEGER, PARAMETER :: RAND_SEED_INDEX  = RAND_BUFFER_SIZE + 1
  INTEGER, PARAMETER :: RAND_CALLS_INDEX = RAND_BUFFER_SIZE + 2
  INTEGER, PARAMETER :: RAND_EXTRA_INTS  = 3
  INTEGER, PARAMETER :: RAND_TOTAL_SIZE  = RAND_BUFFER_SIZE + RAND_EXTRA_INTS

  TYPE t_rand_data
    INTEGER :: mt(0:RAND_TOTAL_SIZE-1)
  END TYPE t_rand_data

CONTAINS

! ******************************************************************************
! Initialization subroutine
! ******************************************************************************
!   sets initial seeds to rdata%mt using the generator Line 25 of Table 1 in
!   [KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102]

  SUBROUTINE RandomSeed(seed,rdata)

! ... parameters
    INTEGER,           INTENT(IN)  :: seed
    TYPE(t_rand_data), INTENT(OUT) :: rdata

! ... loop variables
    INTEGER :: i

! ... local variables
    INTEGER, PARAMETER :: ONE = 1
    INTEGER, PARAMETER :: MULT = 1812433253

    INTEGER :: MASK32

    MASK32 = ISHFT(ONE,32) - ONE

    rdata%mt(0) = IAND(seed,MASK32)

    DO i=1,RAND_BUFFER_SIZE-1
! --- See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier
! --- In the previous versions, MSBs of the seed affect
! --- only MSBs of the array mt[].
! --- 2002/01/09 modified by Makoto Matsumoto

      rdata%mt(i) = MULT * (IEOR(rdata%mt(i-1), &
                                 ISHFT(rdata%mt(i-1),-30))) + i

! --- for >32 bit machines
      rdata%mt(i) = IAND(rdata%mt(i),MASK32)

    END DO ! i

    rdata%mt(RAND_MTI_INDEX)   = RAND_BUFFER_SIZE
    rdata%mt(RAND_SEED_INDEX)  = IAND(seed,MASK32)
    rdata%mt(RAND_CALLS_INDEX) = 0

  END SUBROUTINE RandomSeed

! ******************************************************************************
! Uniform random number generator on [0,1) with 53-bit precision
! ******************************************************************************

  SUBROUTINE RandUniform(a,rdata)

! ... parameters
    REAL(RFREAL),      INTENT(OUT)   :: a(:)
    TYPE(t_rand_data), INTENT(INOUT) :: rdata

! ... loop variables
    INTEGER :: i,j,kk

! ... local variables

! - Period parameters
    INTEGER, PARAMETER :: ONE = 1
    INTEGER, PARAMETER :: M   = 397
    INTEGER, PARAMETER :: LMASK      =  2147483647 ! least significant r bits
    INTEGER, PARAMETER :: PRE_MATA   = -1727483681
    INTEGER, PARAMETER :: PRE_TMASKB = -1658038656
    INTEGER, PARAMETER :: PRE_TMASKC = -272236544

! - Normalization parameters
    REAL(RFREAL), PARAMETER :: TWO32 = 4294967296.0_RFREAL
    REAL(RFREAL), PARAMETER :: NFAC1 = 67108864.0_RFREAL

! - These should be parameters, but cannot be initialized on all platforms
    INTEGER      :: MASK32, MATA, UMASK, TMASKB, TMASKC, mag01(0:1)
    REAL(RFREAL) :: NFAC2

    INTEGER      :: y,mti,iamin,iamax
    REAL(RFREAL) :: yr(2)

! - Set values of things that should, ideally, be parameters

    UMASK    = ISHFT(ONE,31)           ! most significant w-r bits
    MASK32   = ISHFT(ONE,32) - ONE     ! mask to force 32-bit behavior
    MATA     = IAND(PRE_MATA,  MASK32) ! constant vector
    TMASKB   = IAND(PRE_TMASKB,MASK32) ! tempering parameter
    TMASKC   = IAND(PRE_TMASKC,MASK32) ! tempering parameter
    mag01(0) = 0
    mag01(1) = MATA
    NFAC2    = 1.0_RFREAL/9007199254740992.0_RFREAL

! - Begin

    mti = rdata%mt(RAND_MTI_INDEX)

    iamin = LBOUND(a,1)
    iamax = UBOUND(a,1)

    DO i = iamin,iamax

      DO j = 1,2

        IF (mti >= RAND_BUFFER_SIZE) THEN ! make RAND_BUFFER_SIZE integers

          DO kk=0,RAND_BUFFER_SIZE-M-1
            y = IOR(IAND(rdata%mt(kk),UMASK),IAND(rdata%mt(kk+1),LMASK))
            rdata%mt(kk) = IEOR(IEOR(rdata%mt(kk+M),ISHFT(y,-1)), &
                                mag01(IAND(y,ONE)))
          END DO ! kk

          DO kk=RAND_BUFFER_SIZE-M,RAND_BUFFER_SIZE-2
            y = IOR(IAND(rdata%mt(kk),UMASK),IAND(rdata%mt(kk+1),LMASK))
            rdata%mt(kk) = IEOR(IEOR(rdata%mt(kk+(M-RAND_BUFFER_SIZE)), &
                                     ISHFT(y,-1)), mag01(IAND(y,ONE)))
          END DO ! kk

          y = IOR(IAND(rdata%mt(RAND_BUFFER_SIZE-1),UMASK), &
                  IAND(rdata%mt(0),LMASK))
          rdata%mt(RAND_BUFFER_SIZE-1) = &
            IEOR(IEOR(rdata%mt(M-1), ISHFT(y,-1)), mag01(IAND(y,ONE)))

          mti = 0

        END IF ! mti

        y = rdata%mt(mti)

! ----- Tempering

        y = IEOR(y,     ISHFT(y,-11))
        y = IEOR(y,IAND(ISHFT(y,  7),TMASKB))
        y = IEOR(y,IAND(ISHFT(y, 15),TMASKC))
        y = IEOR(y,     ISHFT(y,-18))

        y = ISHFT(y,-4-j)

        IF (y < 0) THEN
          yr(j) = REAL(y,KIND=RFREAL) + TWO32
        ELSE
          yr(j) = REAL(y,KIND=RFREAL)
        END IF ! y

        mti = mti + 1

      END DO ! j

! --- Form random number from two 32-bit integers

      a(i) = (yr(1)*NFAC1 + yr(2)) * NFAC2

    END DO ! i

    rdata%mt(RAND_MTI_INDEX)   = mti
    rdata%mt(RAND_CALLS_INDEX) = rdata%mt(RAND_CALLS_INDEX) + iamax - iamin + 1

  END SUBROUTINE RandUniform

! ******************************************************************************
! Normal distribution
! ******************************************************************************

  SUBROUTINE RandNormal(a,med,dev,rdata)
! med = median
! dev = standard deviation

    REAL(RFREAL),      INTENT(OUT)   :: a(:)
    REAL(RFREAL),      INTENT(IN)    :: med,dev
    TYPE(t_rand_data), INTENT(INOUT) :: rdata

    INTEGER :: i,iamin,iamax
    REAL(RFREAL) :: r,th
    REAL(RFREAL) :: twopi

    iamin = LBOUND(a,1)
    iamax = UBOUND(a,1)
    twopi = 8.0_RFREAL*ATAN(1.0_RFREAL)
    CALL RandUniform(a,rdata)

    DO i=iamin,iamax-1,2
      th = twopi*a(i)
      r  = dev*SQRT(-2.0_RFREAL*LOG(1.0_RFREAL-a(i+1)))
      a(i)   = r*COS(th) + med
      a(i+1) = r*SIN(th) + med
    END DO ! i

    IF (MODULO(iamax-iamin+1,2)==1) THEN
      th = twopi*a(iamax)
      r  = Rand1Uniform(rdata)
      r  = dev*SQRT(-2.0_RFREAL*LOG(1.0_RFREAL-r))
      a(iamax) = r*COS(th) + med
    END IF

  END SUBROUTINE RandNormal

! ******************************************************************************
! Log normal distribution
! ******************************************************************************

  SUBROUTINE RandLogNormal(a,logmed,gdev,rdata)
! logmed = log of median
! gdev = standard deviation of corresponding Gaussian

    REAL(RFREAL),      INTENT(OUT)   :: a(:)
    REAL(RFREAL),      INTENT(IN)    :: logmed,gdev
    TYPE(t_rand_data), INTENT(INOUT) :: rdata

    CALL RandNormal(a,logmed,gdev,rdata)
    a = EXP(a)

  END SUBROUTINE RandLogNormal

! ******************************************************************************
! Discrete PDF distribution from file
! ******************************************************************************

  SUBROUTINE RandImposedPDF(a,pdfvalues,locMax,valmax,rdata)

    REAL(RFREAL),      INTENT(IN)    :: pdfvalues(:,:)
    REAL(RFREAL),      INTENT(OUT)   :: a(:)
    TYPE(t_rand_data), INTENT(INOUT) :: rdata
    INTEGER,           INTENT(IN)    :: locMax
    REAL(RFREAL),      INTENT(IN)    :: valmax

    INTEGER :: i,iamin,iamax,k1,k2

    iamin = LBOUND(a,1)
    iamax = UBOUND(a,1)
    CALL RandUniform(a,rdata)

    DO i=iamin,iamax
       call FINDPDF(a(i),pdfvalues(:,2),k1,k2)
       a(i) = (a(i) - pdfvalues(k1,2))*(pdfvalues(k2,3)-pdfvalues(k1,3))&
             /(pdfvalues(k2,2) - pdfvalues(k1,2)) + pdfvalues(k1,3)
    END DO ! i

    RETURN
    
  CONTAINS
    SUBROUTINE FINDPDF(datum,vec,k1,k2)
      REAL(RFREAL),  INTENT(IN)    :: datum
      REAL(RFREAL),  INTENT(IN)    :: vec(:)
      INTEGER,       INTENT(OUT)   :: k1,k2  
      integer :: isgn

      isgn = sign(1.1_RFREAL,datum-valmax)
      k1 = locMax
      k2 = k1+isgn
      do while (int(sign(1.1_RFREAL,datum-vec(k2)))*isgn > 0)
         k1 = k2
         k2 = k2+isgn
      enddo
    END SUBROUTINE FINDPDF
      

  END SUBROUTINE RandImposedPDF

! ******************************************************************************
! Function form of RandUniform 
! ******************************************************************************

  REAL(RFREAL) FUNCTION Rand1Uniform(rdata)

! ... parameters
    TYPE(t_rand_data), INTENT(INOUT) :: rdata

! ... local variables
    REAL(RFREAL) :: a(1)

    CALL RandUniform(a,rdata)
    Rand1Uniform = a(1)

  END FUNCTION Rand1Uniform

! ******************************************************************************
! Function form of RandNormal
! ******************************************************************************

  REAL(RFREAL) FUNCTION Rand1Normal(med,dev,rdata)

! ... parameters
    REAL(RFREAL),      INTENT(IN)    :: med,dev
    TYPE(t_rand_data), INTENT(INOUT) :: rdata

! ... local variables
    REAL(RFREAL) :: a(1)

    CALL RandNormal(a,med,dev,rdata)
    Rand1Normal = a(1)

  END FUNCTION Rand1Normal

! ******************************************************************************
! Function form of RandLogNormal
! ******************************************************************************

  REAL(RFREAL) FUNCTION Rand1LogNormal(logmed,gdev,rdata)

! ... parameters
    REAL(RFREAL),      INTENT(IN)    :: logmed,gdev
    TYPE(t_rand_data), INTENT(INOUT) :: rdata

! ... local variables
    REAL(RFREAL) :: a(1)

    CALL RandLogNormal(a,logmed,gdev,rdata)
    Rand1LogNormal = a(1)

  END FUNCTION Rand1LogNormal

! ******************************************************************************
! Function form of RandUniform 
! ******************************************************************************

  REAL(RFREAL) FUNCTION Rand1ImposedPDF(rdata,pdfvalues,locMax,valmax)

! ... parameters
    TYPE(t_rand_data), INTENT(INOUT) :: rdata
    REAL(RFREAL),      INTENT(IN)    :: pdfvalues(:,:)
    INTEGER,           INTENT(IN)    :: locMax
    REAL(RFREAL),      INTENT(IN)    :: valmax

! ... local variables
    REAL(RFREAL) :: a(1)

    CALL RandImposedPDF(a,pdfvalues,locMax,valmax,rdata)
    Rand1ImposedPDF = a(1)

  END FUNCTION Rand1ImposedPDF
END MODULE ModRandom

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModRandom.F90,v $
! Revision 1.6  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/04/25 18:39:09  luca1
! Imposed PDF from file option for random particle ejection
!
! Revision 1.3  2003/12/08 16:39:47  jferry
! Removed non-standard KIND intrinsic
!
! Revision 1.2  2003/11/21 22:35:51  fnajjar
! Update Random Number Generator
!
! Revision 1.1  2003/02/17 19:31:11  jferry
! Implemented portable random number generator ModRandom
!
!******************************************************************************






