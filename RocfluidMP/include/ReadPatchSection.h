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
INTERFACE
#ifdef RFLO
  SUBROUTINE ReadPatchSection( global,fileID,nvals,keys,vals,brbeg,brend, &
                               prbeg,prend,distrib,fname,defined )
#endif
#ifdef RFLU
  SUBROUTINE ReadPatchSection( global,fileID,nvals,keys,vals, &
                               prbeg,prend,distrib,fname,bcName,defined )
#endif
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
#ifdef RFLO
    INTEGER      :: brbeg, brend
#endif    
    INTEGER      :: fileID, nvals, prbeg, prend, distrib
    CHARACTER(*) :: keys(nvals), fname
#ifdef RFLU
    CHARACTER(*) :: bcName
#endif    
    LOGICAL      :: defined(nvals)
    REAL(RFREAL) :: vals(nvals)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadPatchSection
END INTERFACE






