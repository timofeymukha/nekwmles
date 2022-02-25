!> @brief Rough log law for temperature 
      subroutine set_heat_flux(h, ix, iy, iz, ie)
      implicit none

      include 'SIZE'
      include 'FRAMELP'
      include 'WMLES'

      ! sampling height
      real h 

      ! local sampled and surface temperature
      real th, ts

      ! the indices of the gll point
      integer ix, iy, iz, ie

      real utau

      ! log law parameters
      real kappa, logb, z1

      ! dummy variables for retrieving parameters
      integer itmp
      logical ltmp
      character*20 ctmp
!-----------------------------------------------------------------------

      ! assign kappa and B and z1
      call rprm_rp_get(itmp,kappa,ltmp,ctmp,wmles_logkappa_id,rpar_real)
      call rprm_rp_get(itmp,logb,ltmp,ctmp,wmles_logb_id,rpar_real)
      call rprm_rp_get(itmp,z1,ltmp,ctmp,wmles_z1_id,rpar_real)

      ! Get utau (should be previously computed by set_momentum_flux)
      utau = tau(1, ix, iy, iz, ie)**2 +
     $       tau(2, ix, iy, iz, ie)**2 +
     $       tau(3, ix, iy, iz, ie)**2
      utau = sqrt(sqrt(utau))
     
      th = temph(ix, iy, iz, ie)
      
      if (ifviscosity) then
        ts = temps(ix, iy, iz, ie)
      else
        ts = wmles_surface_temp
      endif

      heat_flux(ix, iy, iz, ie) = kappa*utau*(ts - th)/log(h/z1)

      end

