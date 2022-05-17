!> @brief Spaldings law in implicit form
      real function set_momentum_flux(h, i, ix, iy, iz, ie)
      implicit none

      include 'SIZE'
      include 'WMLES'

      ! sampling height
      real h 

      ! linear index for the sampled solution values
      integer i

      ! the indices of the gll point
      integer ix, iy, iz, ie

      ! solver for the law of the wall
      real newton

      ! Laws of the wall
      external spalding_value, spalding_derivative

      ! sampled velocity
      real magvh

      ! utau and a guess for it
      real guess, utau
!-----------------------------------------------------------------------

      ! Magnitude of the sampled velocity
      magvh = solh(i, 1)**2 +
     $        solh(i, 2)**2 +
     $        solh(i, 3)**2
      magvh = sqrt(magvh)
     
      ! take the stress mag at the previous step as a guess
      ! inital value at simulation start is the GUESS param
      guess = tau(1, ix, iy, iz, ie)**2 +
     $        tau(2, ix, iy, iz, ie)**2 +
     $        tau(3, ix, iy, iz, ie)**2

      ! Double sqrt to get utau.
      guess = sqrt(sqrt(guess))
      utau = newton(spalding_value, spalding_derivative,
     $              magvh, h, guess, 1e-4, 50)

      ! Assign proportional to the velocity magnitudes at
      ! the sampling point
      tau(1, ix, iy, iz, ie) = -utau**2*solh(i, 1)/magvh
      tau(2, ix, iy, iz, ie) = -utau**2*solh(i, 2)/magvh
      tau(3, ix, iy, iz, ie) = -utau**2*solh(i, 3)/magvh

      end

!-----------------------------------------------------------------------

!> @brief Spaldings law in implicit form
!! @param[in]   u                 velocity magnitude
!! @param[in]   y                 wall-normal distance
!! @param[in]   utau              friction velocity
      real function spalding_value(u, y, utau)
      implicit none
      include 'SIZE'
      include 'FRAMELP'
      include 'WMLES'
      include 'INPUT'

      ! input
      real u, y, utau

      ! model constants
      real kappa, B

      ! scaled quantities
      real uplus, yplus 

      ! dummy variables for retrieving parameters
      integer itmp
      logical ltmp
      character*20 ctmp
!-----------------------------------------------------------------------

      ! assign kappa and B
      call rprm_rp_get(itmp,kappa,ltmp,ctmp,wmles_logkappa_id,rpar_real)
      call rprm_rp_get(itmp,B,ltmp,ctmp,wmles_logb_id,rpar_real)

      uplus = u/utau
      yplus = y*utau/param(2)

      spalding_value = 
     $  (uplus + exp(-kappa*B)*(exp(kappa*uplus) - 1.0 -
     $   kappa*uplus - 0.5*(kappa*uPlus)**2 -
     $   1./6*(kappa*uPlus)**3) -  yplus)

      end

!> @brief Derivative of Spaldings law in implicit form
!! @param[in]   u                 velocity magnitude
!! @param[in]   y                 wall-normal distance
!! @param[in]   utau              friction velocity
      real function spalding_derivative(u, y, utau)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'WMLES'

      ! input
      real u, y, utau

      ! model constants
      real kappa, B

      ! scaled quantities
      real uplus, yplus 

      ! dummy variables for retrieving parameters
      integer itmp
      logical ltmp
      character*20 ctmp
!-----------------------------------------------------------------------

      ! Assign kappa and B
      call rprm_rp_get(itmp,kappa,ltmp,ctmp,wmles_logkappa_id,rpar_real)
      call rprm_rp_get(itmp,B,ltmp,ctmp,wmles_logb_id,rpar_real)

      uplus = u/utau

      spalding_derivative = 
     $  (-y/param(2) - u/utau**2 - kappa*uplus/utau*exp(-kappa*B) *
     $  (exp(kappa*uplus) - 1 - kappa*uplus - 0.5*(kappa*uPlus)**2))
      end
