!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
SUBROUTINE huang_const_model(dfgrad,matrix,nmatrix, &
     particle,nparticle,nparticletype,interfac,ninterfac,stress, strain, ilarge, ismall)

!!****f* Rocfrac/Source/huang_const_model.f90
!!
!!  NAME
!!     huang_const_model
!!
!!  FUNCTION
!!
!!    this program is for the constitutive modeling of solid propellants
!!     accounting for the effect of particle/matrix interface debonding.  the
!!     program can handle one or two types of particles that have the same
!!     elastic moduli but different radii.
!!
!!  INPUTS
!!	
!!-----  all variables above except the last few, i.e., 
!!              dfgrad,matrix,nmatrix,particle,nparticle,
!!              interfac,ninterfac
!!	dfgrad: deformation gradient
!!	nmatrix: number of matrix properties
!!     matrix: properties of the matrix
!!	   matrix(1): young's modulus of the matrix
!!	   matrix(2): poisson's ratio of the matrix
!!	   matrix(3): volume fraction of the matrix
!!     nparticletype: number of the particle types
!!                    nparticletype is no more than 2 in this program!
!!	nparticle: number particle properties for each type 
!!	particle: properties of particles
!!	   particle(1,i): young's modulus of type-i particles
!!	   particle(2,i): poisson's ratio of type-i particles
!!	   particle(3,i): volume fraction of type-i particles
!!	   particle(4,i): radius of type-i particles
!!	ninterfac: number of interface properties
!!	interfac: properties of the interface
!!	   interfac(1): strength of the interface
!!	   interfac(2): linear modulus of the interface
!!	   interfac(3): softening modulus of the interface
!!	stress: stress tensor
!!   strain: strain tensor
!!	ilarge, ismall: damage information for large and small particles
!!	ilarge=1: large particles not failed
!!	ilarge=2: large particles is failing
!!	ilarge=3: large particles completely failed
!!	ismall=1: small particles not failed
!!	ismall=2: small particles is failing
!!	ismall=3: small particles completely failed
!!
!!  OUTPUT
!!
!!    stress -- Stress
!!    Strain -- Strain
!!    ilarge, ismall -- failing state of large and small particles.
!!****

!-----  arrays for input/output
  REAL*8 dfgrad(3,3),matrix(nmatrix), &
       particle(nparticle,nparticletype),interfac(ninterfac), &
       stress(3,3)

  REAL*8 :: ilarge, ismall

!-----  arrays for the internal use in the subroutine
  REAL*8 alpha(nparticletype),alphaprime(nparticletype), &
       bulkp(nparticletype),shearp(nparticletype), &
       strain(3,3),stressprime(3,3),strainprime(3,3), &
       mstrainpt(5),mstresspt(5)

!-----  variables for internal use in the subroutine
  REAL*8 em,posm,bulkm,shearm,bulkcomposite,shearcomposite,strength, &
       kinthard,kintsoft,stressmean,strainmean

!-----  the green strain 
  DO i=1,3
     DO  j=1,3
        strain(i,j)= dfgrad(1,i)*dfgrad(1,j) &
             +dfgrad(2,i)*dfgrad(2,j) &
             +dfgrad(3,i)*dfgrad(3,j) 
        IF (i.EQ.j) THEN
           strain(i,j)=(strain(i,j)-1)*0.5d0
        ELSE
           strain(i,j)=strain(i,j)*0.5d0
        END IF
     ENDDO
  ENDDO

!-----  the mean and deviatoric strains
  strainmean=(strain(1,1)+strain(2,2)+strain(3,3))/3.d0
  DO  i=1,3
     DO  j=1,3
        IF (i.EQ.j) THEN 
           strainprime(i,j)=strain(i,j)-strainmean
        ELSE
           strainprime(i,j)=strain(i,j)
        END IF
     ENDDO
  ENDDO
				
!-----  parameters for particles: alpha and alphaprime
  em   = matrix(1)
  posm = matrix(2)
  bulkm= em/(3.d0*(1.d0-2.d0*posm))
  shearm=em/(2.d0*(1.d0+posm))
  strength=interfac(1)
  kinthard=interfac(2)
  kintsoft=interfac(3)
  
  DO i=1,nparticletype
     bulkp(i)=particle(1,i)/(3*(1-2*particle(2,i)))
     shearp(i)=particle(1,i)/(2*(1+particle(2,i)))
     alpha(i)=3.d0*(1-posm)/(2*em) &
          /( 1.d0/(kinthard*particle(4,i)) &
          +1.d0/(3.d0*bulkp(i)) &
          +1.d0/(4.d0*shearm))
     alphaprime(i)=-3.d0*(1.d0-posm)/(2.d0*em) &
          /(-1.d0/(kintsoft*particle(4,i)) &
          +1.d0/(3.d0*bulkp(i)) &
          +1.d0/(4.d0*shearm))
     
  ENDDO


  c1=particle(3,1)
  c2=particle(3,2)
  a1=particle(4,1)
  a2=particle(4,2)
  c=c1+c2

  alpha1=alpha(1)
  alpha2=alpha(2)
  alphaprime1=alphaprime(1)
  alphaprime2=alphaprime(2)

!-----  critical particle radius
  aprime=1.d0/(1.d0/(4.d0*shearm)+1.d0/(3.d0*bulkp(1)))/kintsoft        

!-----  the mean stress and strain of the composite at the transition point 
!     between different stages.
!     path2:(i,i)->(ii,i)->(iii,i)->(iii,ii)->(iii,iii)
  mstrainpt(1)=0.d0
  mstresspt(1)=0.d0
!     (i,i)->(ii,i)
  mstrainpt(2)=strength/(2.d0*em*alpha1) &
       *(2.d0*(1.d0-2.d0*posm)-(1.d0+posm)*(-c+c2*alpha2+c1*alpha1))
  mstresspt(2)=strength*(1.d0-c+c2*alpha2+c1*alpha1)/alpha1
!     (ii,i)->(iii,i)
  mstrainpt(3)=strength*(1.d0/alpha1+1.d0/alphaprime1)  &
       *(2.d0*(1.d0-2.d0*posm)+(1.d0+posm)*(c-c2*alpha2))/(2.d0*em)	 
  mstresspt(3)=strength*(1.d0/alpha1+1.d0/alphaprime1)*(1.d0-c+c2*alpha2)
!     (iii,i)->(iii,ii)
  mstrainpt(4)=strength*(2.d0*(1.d0-2.d0*posm)+c*(1.d0+posm)-(1.d0+posm)*c2*alpha2)  &
       /(2.d0*em*alpha2)
  mstresspt(4)=strength*(1.d0-c+c2*alpha2)/alpha2
!     (iii,ii)->(iii,iii)
  mstrainpt(5)=strength/(2.d0*em) &
       *(1/alpha2+1/alphaprime2)*(2*(1-2*posm)+c*(1+posm))
  mstresspt(5)=(1.d0-c)*strength*(1.d0/alpha2+1.d0/alphaprime2)

!     path1:(i,i)->(ii,i)->(ii,ii)->(iii,ii)->(iii,iii)
  IF (a1.LT.aprime) THEN
!     (ii,i)->(ii,ii) 
     mstrainpt(3)=strength/(2.d0*em)  &
          *(2.d0*(1.d0-2.d0*posm)/alpha2  &
          -(1.d0+posm)*(c+c1*alphaprime1*(1.d0/alpha1-1.d0/alpha2)  &
          -c/alpha2))        
     mstresspt(3)=strength  &
          *((1.d0-c)/alpha2+c+c1*alphaprime1 &
          *(1.d0/alpha1-1.d0/alpha2))     
!     (ii,ii)->(iii,ii)
     mstrainpt(4)=strength/(2.d0*em) &
          *((2.d0*(1.d0-2.d0*posm)+c*(1.d0+posm)) &
          *(1.d0/alpha1+1.d0/alphaprime1) &
          -(1.d0+posm)*c2*(1.d0-alphaprime2/alpha1 &
          -alphaprime2/alphaprime1 &
          +alphaprime2/alpha2))
     mstresspt(4)=strength*((1.d0-c)*(1.d0/alpha1+1.d0/alphaprime1) &
          +c2*(1.d0-alphaprime2/alpha1 &
          -alphaprime2/alphaprime1 &
          +alphaprime2/alpha2))
  END IF

!----  mean stress in the composite
  DO i=1,4
     IF (strainmean.GE.mstrainpt(i).AND.  &
          strainmean.LT.mstrainpt(i+1)) THEN
        stressmean=mstresspt(i)+(strainmean-mstrainpt(i)) &
             *(mstresspt(i+1)-mstresspt(i)) &
             /(mstrainpt(i+1)-mstrainpt(i))

		IF(i.EQ.1) THEN
			ilarge=1
			ismall=1
		ELSE IF(i.EQ.2) THEN
			ilarge=2
			ismall=1
		ELSE IF(i.EQ.3) THEN
			IF(a1.LT.aprime) THEN
				ilarge=2
				ismall=2
			ELSE
				ilarge=3
				ismall=1
			ENDIF
		ELSE IF(i.EQ.4) THEN
			ilarge=3
			ismall=2
		ENDIF
     END IF
  ENDDO
  IF (strainmean.GE.mstrainpt(5)) THEN
     stressmean=mstresspt(5)*strainmean/mstrainpt(5)
		ilarge=3
		ismall=3
  END IF
  IF (strainmean.LT.0) THEN
     bulkcomposite=bulkm+c*(bulkp(1)-bulkm)*3.d0*(1-posm) &
          /(3.d0*(1.d0-posm)+(1.d0+posm) &
          *(bulkp(1)/bulkm-1.d0)*(1.d0-c))
     stressmean=3.d0*bulkcomposite*strainmean
		ilarge=1
		ismall=1
  END IF

!-----  deviatoric stress
  shearcomposite=shearm+c*(shearp(1)-shearm)*15.d0*(1.d0-posm) &
       /(15.d0*(1.d0-posm)+2.d0*(4.d0-5.d0*posm) &
       *(shearp(1)/shearm-1.d0)*(1.d0-c))
  DO i=1,3
     DO j=1,3
        stressprime(i,j)=2.d0*shearcomposite*strainprime(i,j)
     ENDDO
  ENDDO

!-----  update the stress
  DO i=1,3
     DO j=1,3
        IF (i.EQ.j) THEN 
           stress(i,j)=stressprime(i,j)+stressmean
        ELSE
           stress(i,j)=stressprime(i,j)
        END IF
     ENDDO
  ENDDO

END SUBROUTINE huang_const_model

