      real function newton(u, y, nu, utau, tol, maxiter, debug)

      implicit none

      ! Stuff to pass to the wall model
      real u, y , nu
      

      real spalding_value
      real spalding_derivative

      ! The friction velocty to be found
      real utau
      
      ! Misc parameters
      real tol
      integer maxiter
      logical debug

      real deltax, fx, fxprime
      integer k

      if (debug) then
        print 11, utau
11      format('Initial guess: utau = ', e22.15)
      endif

      newton = utau

      do k=1,maxiter

        ! evaluate function and its derivative:
        fx = spalding_value(u, y, nu, newton)
        fxprime = spalding_derivative(u, y, nu, newton)

        if (abs(fx) < tol) then
          exit  ! jump out of do loop
        endif

        ! compute Newton increment x:
        deltax = fx/fxprime

        ! update x:
        newton = newton - deltax

        if (debug) then
          print 12, k, newton 
 12       format('After', i3, ' iterations, x = ', e22.15)
        endif

      enddo

      end


      real function spalding_value(u, y, nu, utau)
      ! Return the value of the implicit function defined by Spalding's
      ! law. To be used with Newton's law.

      implicit none

      ! Sampled velocity value
      real u

      ! Wall-normal distance to the sampling point
      real y

      ! Kinematic viscosity
      real nu

      ! The friction velocity
      real utau 

      ! model constants
      real kappa, B

      ! scaled quantities
      real uplus, yplus 

      kappa  = 0.41
      B = 5.5

      uplus = u/utau
      yplus = y*utau/nu
      spalding_value = 
     $  (uplus + exp(-kappa*B)*(exp(kappa*uplus) - 1.0 -
     $   kappa*uplus - 0.5*(kappa*uPlus)**2 -
     $   1./6*(kappa*uPlus)**3) -  yPlus)

      end

      real function spalding_derivative(u, y, nu, utau)
      !Return the value of the derivative of the implicit function
      ! defined by Spalding's law. To be used with Newton's method.

      implicit none

      ! Sampled velocity value
      real u

      ! Wall-normal distance to the sampling point
      real y

      ! Kinematic viscosity
      real nu

      ! The friction velocity
      real utau 

      ! model constants
      real kappa, B

      ! scaled quantities
      real uplus, yplus 

      kappa  = 0.41
      B = 5.5

      uplus = u/utau

      spalding_derivative = 
     $  (-y/nu - u/utau**2 - kappa*uplus/utau*exp(-kappa*B) *
     $  (exp(kappa*uplus) - 1 - kappa*uplus - 0.5*(kappa*uPlus)**2))

      end