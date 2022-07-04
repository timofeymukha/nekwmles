!> @brief Rough log law for temperature 
      subroutine wmles_set_heat_flux(i)
      implicit none

      include 'SIZE'
      include 'FRAMELP'
      include 'WMLES'

      ! sampling height
      real h 

      ! local sampled and surface temperature
      real th, ts

      ! the indices of the gll point
      integer i, ix, iy, iz, ie

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

      h = wmles_sampling_h(i)

      ! Get utau (should be previously computed by set_momentum_flux)
      utau = wmles_tau(i, 1)**2 +
     $       wmles_tau(i, 2)**2 +
     $       wmles_tau(i, 3)**2
      utau = sqrt(sqrt(utau))
     
      th = wmles_solh(i, 4)
      
      if (wmles_ifviscosity) then
        ts = wmles_solh(i, 5)
      else
        ts = wmles_surface_temp
      endif


      wmles_q(i) = kappa*utau*(ts - th)/log(h/z1)

      end

