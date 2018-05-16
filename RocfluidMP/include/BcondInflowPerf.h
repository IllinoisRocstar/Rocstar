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
  SUBROUTINE BcondInflowPerf( bcOpt,ptot,ttot,betah,betav,mach,sxn,syn,szn, & 
                              cpgas,mm,rl,rul,rvl,rwl,rr,rur,rvr,rwr, & 
                              rer,pr )
    USE ModDataTypes
    INTEGER, INTENT(IN) :: bcOpt
    REAL(RFREAL), INTENT(IN)  :: betah, betav, cpgas, mach, mm, sxn, syn, szn, &
                                 ptot, rl, rul, rvl, rwl, ttot
    REAL(RFREAL), INTENT(OUT) :: rer, rr, rur, rvr, rwr, pr
  END SUBROUTINE BcondInflowPerf
END INTERFACE






