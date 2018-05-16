      program hep
c
c Compute the steady-state head-end pressure of a rocket with
c a given burning area.
c
c Also compute optimum thrust.
c
c Usage:
c
c % hep_area.x < parameters
c
c Robert Fiedler, revised 4/16/05
c
      implicit NONE
      real pi,n,a,p_ref,rho_s,gamma,mw,t_flame,r_t
      real r_gas,gm1,sb,b,p_c,area_t,area_c
      real a_c,cp,mf_ig
      real gp1, thrust, p_a, r_e
      integer ios
c
c Read from standard input instead of parameters file
c
c      open(8,file='parameters_BATES',iostat=ios)
c      if (ios /= 0) then
c        print *,'Need to read from file called: parameters_rsrm'
c        stop
c      endif
      pi = 2.*asin(1.)
c
c DEFINITION OF PARAMETERS
C
c n        exponent in burn rate a*(p/p_ref)^n
c a        burn rate in m/s at pressure p_ref
c p_ref    pressure for normalization in burn rate (Pa)
c rho_s    propellant density (kg/m^3)
c gamma    ratio of specific heats
c mw       molecular weight (g/mol)
c cp       heat capacity at constant pressure (MKS units)
c t_flame  flame temperature (K)
c a_c      chamber area (m^2)
c r_t      nozzle throat radius (m)
c r_e      nozzle exit radius (m)
c P_a      ambient pressure (Pa)
c
      read(*,*) n
      read(*,*) a
      read(*,*) p_ref
      read(*,*) rho_s
      read(*,*) gamma
c      read(*,*) mw
      read(*,*) cp
      read(*,*) t_flame
      read(*,*) a_c
      read(*,*) r_t
      read(*,*,iostat=ios) r_e
      read(*,*,iostat=ios) p_a

c Hard code typical igniter mass flux

c      mf_ig = 200.
      mf_ig = 0.
c
c Compute gas constant
c
      gm1 = gamma - 1.

c      r_gas = 1000.*8.3143/mw
      if (cp < 100.) then
c
c The molecular weight was input as cp
c
        r_gas = 1000.*8.3143/cp
      else
c
c The heat capacity was input
c
        r_gas = cp * gm1 / gamma
      endif
c
c Compute throat area and chamber area
c
      area_t = pi * r_t**2
      area_c = a_c

c Use 0-D theory to get steady chamber pressure step by step
c
      sb = 1. + 0.5*gm1
c
      b = sqrt(gamma/(r_gas*(t_flame/sb)))
      b = b /sb**(gamma/gm1)
      b = b*area_t
      b = b*p_ref**n/(area_c*rho_s*a + mf_ig)
c
c b is p_c**(n-1.)
c
      p_c = b**(1./(n-1.))
c
      print *,'Steady pressure is ',p_c,' Pa, or ',p_c/6895.,' psi'
c
      if (ios == 0) then
        gp1 = gamma + 1.0
        thrust = 1.0 - (p_a/p_c)**(gm1/gamma)
        thrust = thrust * (2.0/gp1)**(gp1/gm1)
        thrust = thrust * 2.0 * gamma**2 / gm1
        thrust = area_t * p_c * sqrt (thrust)
c
        print *,'Optimum thrust is ',thrust,' N'
c
      endif
c
      end
