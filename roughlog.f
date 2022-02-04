!> @brief Rough log law in implicit form
      real function set_wall_quantities(h, ix, iy, iz, ie)
      implicit none

      include 'SIZE'
      include 'FRAMELP'
      include 'WMLES'

      ! sampling height
      real h 

      ! the indices of the gll point
      integer ix, iy, iz, ie

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

      ! Magnitude of the sampled velocity
      magvh = vh(1, ix, iy, iz, ie)**2 +
     $        vh(2, ix, iy, iz, ie)**2 +
     $        vh(3, ix, iy, iz, ie)**2
      magvh = sqrt(magvh)
     
      utau = (magvh - logb)*kappa/log(h/z0)

      ! Assign proportional to the velocity magnitudes at
      ! the sampling point
      tau(1, ix, iy, iz, ie) = -utau**2*vh(1, ix, iy, iz, ie)/magvh
      tau(2, ix, iy, iz, ie) = -utau**2*vh(2, ix, iy, iz, ie)/magvh
      tau(3, ix, iy, iz, ie) = -utau**2*vh(3, ix, iy, iz, ie)/magvh

      end

