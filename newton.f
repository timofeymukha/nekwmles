!> @brief Newton-Raphson solver for non-linear algeraic equations
!! @param[in]   u              velocity magnitude
!! @param[in]   y              wall-normal distance
!! @param[in]   utau           initial guess for the friction velocity
!! @param[in]   tol            error tolerance
!! @param[in]   maxitter       max number of iterations
      real function newton(f, d, u, y, utau, tol, maxiter)

      implicit none

      ! input
      real f, d
      real u, y, utau, tol
      integer maxiter 

      real deltax, fx, fxprime
      integer k
      
      newton = utau

      do k=1,maxiter

        ! evaluate function and its derivative
        fx = f(u, y, newton)
        fxprime = d(u, y, newton)

        if (abs(fx) < tol) then
          exit
        endif

        ! compute Newton increment
        deltax = fx/fxprime

        ! update solution:
        newton = newton - deltax

      enddo

      end
