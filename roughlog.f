!> @brief Rough log law in implicit form
      subroutine wmles_set_momentum_flux(i)
      implicit none

      include 'SIZE'
      include 'FRAMELP'
      include 'WMLES'

      ! sampling height
      real h 

      ! the indices of the gll point
      integer i, ix, iy, iz, ie

      ! sampled velocity
      real magvh

      real utau

      ! log law parameters
      real kappa, logb, z0

      ! dummy variables for retrieving parameters
      integer itmp
      logical ltmp
      character*20 ctmp
!-----------------------------------------------------------------------

      ! assign kappa and B and z0
      call rprm_rp_get(itmp,kappa,ltmp,ctmp,wmles_logkappa_id,rpar_real)
      call rprm_rp_get(itmp,logb,ltmp,ctmp,wmles_logb_id,rpar_real)
      call rprm_rp_get(itmp,z0,ltmp,ctmp,wmles_z0_id,rpar_real)
      
      ix = wmles_indices(i, 1)
      iy = wmles_indices(i, 2)
      iz = wmles_indices(i, 3)
      ie = wmles_indices(i, 4)

      h = wmles_sampling_h(i)

      ! Magnitude of the sampled velocity
      magvh = wmles_solh(i, 1)**2 +
     $        wmles_solh(i, 2)**2 +
     $        wmles_solh(i, 3)**2
      magvh = sqrt(magvh)
     
      utau = (magvh - logb)*kappa/log(h/z0)

      ! Assign proportional to the velocity magnitudes at
      ! the sampling point
      wmles_tau(i, 1) = -utau**2*wmles_solh(i, 1)/magvh
      wmles_tau(i, 2) = -utau**2*wmles_solh(i, 2)/magvh
      wmles_tau(i, 3) = -utau**2*wmles_solh(i, 3)/magvh

      end

