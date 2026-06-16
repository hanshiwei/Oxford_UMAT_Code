!     May 03rd, 2022
!     Eralp Demir
!
      module backstress
      implicit none
      contains
!
!     Armstrong-Frederic backstress model
!     Two input parameters are required
      subroutine backstressmodel1(backstressparam,nslip,X,gdot,dt,dX)
      use userinputs, only : maxnparam
      implicit none
!     Inputs
!     Backstress parameters
      real(8), intent(in) :: backstressparam(maxnparam)
!     Number of slip systems
      integer, intent(in) :: nslip
!     Current value of slip rate
      real(8), intent(in) :: gdot(nslip)
!     Current value of backstress
      real(8), intent(in) :: X(nslip)
!     time increment
      real(8), intent(in) :: dt
!
!     Output
!     Backtress increment
      real(8), intent(out) :: dX(nslip)
!
!     Variables used in this subroutine
!     Backstress evolution parameter
      real(8) :: h
!     Backstress relief parameter
      real(8) :: hD
      integer :: is
!
!     Backstress evolution
      h = backstressparam(1)
!
!     Backstress annihiliation
      hD = backstressparam(2)
!
      dX = 0.
      do is=1,nslip
!
          dX(is) = (h*gdot(is) - hD*X(is)*abs(gdot(is)))*dt
!
      end do
!
!
!
      return
      end subroutine backstressmodel1
!
!
!
!
!
!    NON-LOCAL backstress calculation based on GND density
      subroutine backstressmodel2
      use userinputs, only: maxnslip
      use globalvariables, only: numel, numpt,
     + materialid, numslip_all, numscrew_all, screw_all,
     + burgerv_all, gf_all, G12_all, v12_all,
     + backstressparam_all, statev_backstress,
     + statev_gnd, statev_gammasum
      use utilities, only: matvec6
      implicit none
      integer :: matid, nslip, nscrew
      real(8) :: burgerv(maxnslip)
      integer :: screw(maxnslip)
      real(8) :: X(maxnslip)
      real(8) :: gf, G12, v12, xi
      real(8) :: rhoGNDe, rhoGNDs, gsum
      real(8) :: rhoGND(maxnslip)
      integer :: ie, ip, is, i, k, l
!
!
!
!
!
!
      do ie=1, numel
!
!         Reset arrays
          burgerv=0.;screw=0
!
!
!
!         Assume the same material for al the Gaussian points of an element
          matid = materialid(ie,1)
!
!         Number of slip systems
          nslip = numslip_all(matid)
!
!         Number of screw systems
          nscrew = numscrew_all(matid)
!
!         Screw systems
          screw = screw_all(matid,1:nscrew)
!
!         Geometric factor
          gf = gf_all(matid)
!
!         Shear Modulus
          G12 = G12_all(matid)
!
!         Poisson's ratio
          v12 = v12_all(matid)
!
!         Burgers vector
          burgerv(1:nslip) = burgerv_all(matid,1:nslip)
!
!         Backstress parameter
          xi = backstressparam_all(matid,1)
!
!
          do ip=1,numpt
!
!
!
!
!
!
!             Reset backstress
              X = 0.
!
!
              rhoGND=0.
!             Edges
              do is = 1, nslip
!
!                 Sum of shear rates
                  gsum = statev_gammasum(ie,ip,is)

!                 Edge dislocation density
                  rhoGNDe = statev_gnd(ie,ip,is)
!
!!
!                  rhoGND(is) = abs(statev_gnd(ie,ip,is))
!
!
!                 Backstress due to edges
                  X(is) = xi*G12/(1.-v12)*burgerv(is)*
     + sqrt(abs(rhoGNDe))*sign(1.0,gsum)
!
!
!
              end do
!
!
!
!             Screws
              do i = 1, nscrew
!
                  is = screw(i)
!
!                 Sum of shear rates
                  gsum = statev_gammasum(ie,ip,is)

!                 Screw dislocation density
                  rhoGNDs = statev_gnd(ie,ip,nslip+i)
!
!                  rhoGND(is) = rhoGND(is) + 
!     + abs(statev_gnd(ie,ip,nslip+i))
!
!                 Backstress due to screws
                  X(is) = X(is) + xi*G12*burgerv(is)*
     + sqrt(abs(rhoGNDs))*sign(1.0,gsum)
!
!
!
              end do
!
!
!!             Backstress
!              do is = 1, nslip
!!
!!                 Sum of shear rates
!                  gsum = statev_gammasum(ie,ip,is)
!!
!                  X(is) = xi*G12*burgerv(is)*sqrt(rhoGND(is))
!     + *sign(1.0,gsum)
!!
!              end do
!
!
!             Assign to the global vector
              statev_backstress(ie,ip,1:nslip) = X(1:nslip)
!
!
!
!
!
!
          end do
!
!
!
!
!
!
!
      end do
!
!
!
!
!
!
!
      return
!
      end subroutine backstressmodel2
!

!<<<<<by Shiwei 2026/06/09
!     dislocation well model
!     activate this backstress model by setting variable backstressmodel = 11 in userinputs.f
!     
!     Two input parameters are required
      subroutine dislocationwellmodel(backstressparam, nslip, X, gdot,
     + dt, dX, GamImp_t, GamImp, dGamImp, CapImp)
      use userinputs, only : maxnparam,nimpede,bpredef
      use errors, only : error

      implicit none
!     Inputs
!     Backstress parameters
      real(8), intent(in) :: backstressparam(maxnparam)
!     Number of slip systems
      integer, intent(in) :: nslip
!     Current value of slip rate
      real(8), intent(in) :: gdot(nslip)
!     Current value of backstress
      real(8), intent(in) :: X(nslip)
!     time increment
      real(8), intent(in) :: dt
!     impeded shear strain on each slip system converged at time t
      real(8), intent(in) :: GamImp_t(nslip,nimpede)
!     (trial) impeded shear strain on each slip system from outer iteration
      real(8), intent(in) :: GamImp(nslip,nimpede)
!     Capcity of total impeded strain on each slip system 
      real(8), intent(in) :: CapImp(nslip,nimpede)    
!     Output
!     Backtress increment
      real(8), intent(out) :: dX(nslip)
!     increment of the impeded shear strain on each slip system
      real(8), intent(out) :: dGamImp(nslip,nimpede)

      
!     Variables used in this subroutine
!     hardening moduli
      real(8) :: c
      
!     percentage 1st term
      real(8) :: f
!     impeded strain capacity
      real(8) :: h
!     tunneling index
      real(8) :: q

      
      integer :: is,k
      
      real(8) :: signHeav
      real(8) :: gamImpRatAbs,gamImpCoe,GamImpTri
      

!     Hardening moduli
      c = backstressparam(1)
      
!
      dX = 0.
      
      
      do k=1,nimpede
          f=   backstressparam(1+(K-1)*3+1)          
          !h=   backstressparam(1+(K-1)*3+2)
          q=1./backstressparam(1+(K-1)*3+3)
          do is=1,nslip
              
              h=CapImp(is,k)    
              !  check if h is less than or equal to 0, if so, throw an error    
              if (h <= 0.) then
                  call error(101)
              end if    
              
              !   Heaviside symbol H(gammadot*f*gammaImp)
              signHeav = 0.
              if ((gdot(is)*GamImp(is,k))>0.) then
                  signHeav=1.
              end if 
              
              !   |gammaImp/h|
              gamImpRatAbs=abs(GamImp(is,k)/h)
              
              !   
              if (gamImpRatAbs>=1.) then
                  gamImpCoe= signHeav       
              else 
                  gamImpCoe= gamImpRatAbs**q*signHeav
              end if
              !   [1 - |gammaImp/h|^q*Heaviside(gammadot*GamImp)]*GammaDot*fraction
              dGamImp(is,K)=(1.- gamImpCoe)*gdot(is)*dt*f
              
              !   avoid the overshoot of the impeded strain
              !   if Gimp + dGimp > h or <-h, then dGimp should be adjusted accordingly
              GamImpTri = GamImp_t(is,K) + dGamImp(is,K)
              
              if (GamImpTri > +h) GamImpTri = +h
              if (GamImpTri < -h) GamImpTri = -h
              
              dGamImp(is,K) = GamImpTri - GamImp_t(is,K)
              
              !   form the backstress by linear hardening moduli
              dX(is) = dX(is) + C*dGamImp(is,K)
    !         
              
          end do
      end do
!
!
      return
      end subroutine dislocationwellmodel
!<<<<<by Shiwei 2026/06/09
!
      end module backstress