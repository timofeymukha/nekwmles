!> @brief Rough log law for temperature in implicit form
      subroutine set_heat_flux(h, ix, iy, iz, ie)
      implicit none

      include 'SIZE'
      include 'FRAMELP'
      include 'WMLES'

      ! sampling height
      real h 

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

      ! assign kappa and B and z0
      call rprm_rp_get(itmp,kappa,ltmp,ctmp,wmles_logkappa_id,rpar_real)
      call rprm_rp_get(itmp,logb,ltmp,ctmp,wmles_logb_id,rpar_real)
      call rprm_rp_get(itmp,z1,ltmp,ctmp,wmles_z1_id,rpar_real)

      ! Magnitude of the sampled velocity
      utau = tau(1, ix, iy, iz, ie)**2 +
     $       tau(2, ix, iy, iz, ie)**2 +
     $       tau(3, ix, iy, iz, ie)**2
      utau = sqrt(sqrt(utau))
     
      utau = (magvh - logb)*kappa/log(h/z0)

      heat_flux(ix, iy, iz, ie) = -utau**2*th(1, ix, iy, iz, ie)/magvh

      end

