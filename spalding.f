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