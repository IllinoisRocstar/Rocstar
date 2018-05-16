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
! Purpose: read in user defined porbability density function (PDF) of
!          the mass injection process.
!
! Description: none
!
! Input: user input file.
!
! Output: structure with information relative to the imposed pdf
!
! Notes:  none
!
!******************************************************************************
!
! $Id: PLAG_ReadPdfFromFile.F90,v 1.4 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_ReadPdfFromFile(regions, brbeg, brend )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters

  USE PLAG_ModParameters

  TYPE(t_region), POINTER :: regions(:)
  INTEGER,INTENT(IN)      :: brbeg,brend

! ... local variables
  TYPE(t_global), POINTER :: global
  CHARACTER(CHRLEN)       :: fname
  REAL(RFREAL),POINTER    :: tmpvec(:,:)
  REAL(RFREAL)            :: tmpscal
  INTEGER                 :: k,tmpint,errorFlag,nbins,nrow
!******************************************************************************


  global => regions(1)%global

  CALL RegisterFunction( global, 'PLAG_ReadPdfFromFile',&
  'PLAG_ReadPdfFromFile.F90' )

!OPEN PDF file
  WRITE(fname,'(A,A,A)') TRIM(global%inDir),TRIM(global%casename),'.plag_injcpdf'
  OPEN(IF_PLAG_INJCPDF,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
       CALL ErrorStop( global, ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )
  READ(IF_PLAG_INJCPDF,*,err=10,end=10) nbins

!Allocation
  do k = brbeg,brend
     ALLOCATE(regions(k)%plagInput%PDF%pdfvalues(nbins+1,3), stat=errorFlag )
  enddo
  global%error = errorFlag
  ALLOCATE(tmpvec(nbins+1,5), stat=errorFlag)
  global%error = abs(global%error) + abs(errorFlag)
  IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )

!Read in the discrete PDF
  do k = 1,nbins
     READ(IF_PLAG_INJCPDF,*,err=10,end=10) nrow,tmpvec(k,1:2)
  enddo
        
  tmpvec(1,3) = (tmpvec(2,1)-tmpvec(1,1))        
  do k = 2,nbins
     tmpvec(k,3) = 2_RFREAL*(tmpvec(k,1)-tmpvec(k-1,1)) - tmpvec(k-1,3)
  enddo

!evaluate the two following cumulative sums 
  tmpvec(1,4:5) = 0_RFREAL
  do k = 1,nbins
     tmpvec(k+1,4) = tmpvec(k,4) + tmpvec(k,2)*tmpvec(k,3)
     tmpvec(k+1,5) = tmpvec(k,5) + tmpvec(k,3)
  enddo
  tmpvec(:,4) = tmpvec(:,4) / tmpvec(nbins+1,4)
  tmpvec(:,5) = tmpvec(:,5) + tmpvec(1,1) - tmpvec(1,3)/2_RFREAL


!assign the tmp vector entries to the region pointer
  regions(brbeg:brend)%plagInput%PDF%nbins = nbins        
  do k = brbeg,brend
     regions(k)%plagInput%PDF%pdfvalues(:,1:3) = tmpvec(:,3:5)
  enddo

!--find the maximum in the PDF curve. This is the most probable diameter value,
!--the one from which the search will start;                                                        
  tmpscal = 0_RFREAL
  do k = 1,nbins
     if(tmpvec(k+1,4) - tmpvec(k,4) > tmpscal) then
        tmpint = k+1
        tmpscal = tmpvec(k+1,4) - tmpvec(k,4)
     endif
  enddo
  tmpscal = tmpvec(tmpint,4)
  regions(brbeg:brend)%plagInput%PDF%locmax = tmpint
  regions(brbeg:brend)%plagInput%PDF%valmax = tmpscal

  DEALLOCATE(tmpvec)
  close(IF_PLAG_INJCPDF)

GOTO 999

10  CONTINUE
  CALL ErrorStop( global, ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) )

! finalize --------------------------------------------------------------------

999 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_ReadPdfFromFile

!******************************************************************************







