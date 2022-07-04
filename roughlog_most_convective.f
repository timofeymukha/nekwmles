!> @brief Rough log law in implicit form
      subroutine wmles_set_momentum_flux(i)
      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'FRAMELP'
      include 'WMLES'

      ! sampling height
      real h
      
      ! obukhov length
      real l_obukhov, l_upper, l_lower, l_backup, l_old

      ! the indices of the gll point
      integer i, ix, iy, iz, ie, count

      ! sampled velocity
      real magvh
      
      ! richardson number
      real rib

      real utau, utau_old, g, th, ts, q
      
      ! newton iteration stuff
      real f, dfdl, fd_h

      ! log law parameters
      real kappa, z0, z1
      
      ! similarity law for velocity
      real similarity_law

      ! dummy variables for retrieving parameters
      integer itmp
      logical ltmp
      character*20 ctmp
!-----------------------------------------------------------------------
      

      ! assign kappa and B and z0
      call rprm_rp_get(itmp,kappa,ltmp,ctmp,wmles_logkappa_id,rpar_real)
      call rprm_rp_get(itmp,z0,ltmp,ctmp,wmles_z0_id,rpar_real)
      call rprm_rp_get(itmp,z1,ltmp,ctmp,wmles_z1_id,rpar_real)
      
      g = 9.80665

      ix = wmles_indices(i, 1)
      iy = wmles_indices(i, 2)
      iz = wmles_indices(i, 3)
      ie = wmles_indices(i, 4)

      h = wmles_sampling_h(i)


      ! Sample the values at the sampling point
      ! Velocity
      magvh = wmles_solh(i, 1)**2 +
     $        wmles_solh(i, 2)**2 +
     $        wmles_solh(i, 3)**2
      magvh = sqrt(magvh)

      ! Temperature
      th = wmles_solh(i, 4)
      
      ! Surface temperature
      if (wmles_ifviscosity) then
        ts = wmles_solh(i, 5) !not actually supported yet!
      else
        ts = wmles_surface_temp
      endif

      ! Get uncorrected utau for a first guess other wise from last
      ! timestep
      if (ISTEP .lt. 3) then
        utau = wmles_tau(i, 1)**2 +
     $         wmles_tau(i, 2)**2 +
     $         wmles_tau(i, 3)**2
        utau = sqrt(sqrt(utau))
      else
        utau = magvh*kappa/log(h/z0)
      end if
      
      
      if (ISTEP .gt. 3) then
        ! q is known
        q = wmles_q(i)
        rib = -g*h/th*q/(magvh**3*kappa**2)
        
        ! Obukhov l based on the previous-step utau
        l_obukhov = -(wmles_theta0*utau**3)/(kappa*g*q)
        
        ! In case the iteration diverges we will just use this
        l_backup = l_obukhov
        
        l_old = 0
        count = 0
c        !write(*,*) "ERR",  abs(l_old - l_obukhov)/l_obukhov 
        do while ((abs(l_old - l_obukhov)/abs(l_obukhov) .gt. 1e-3)
     $             .and. (count .lt. 20))

          l_old = l_obukhov
          count = count + 1

          ! for the central diff for evaluating dfdl
          fd_h = 1e-3*l_obukhov
          l_upper = l_obukhov + fd_h
          l_lower = l_obukhov - fd_h

          f = (rib - h/l_obukhov/similarity_law(l_obukhov, h, z0)**3)

          dfdl = (-h/l_upper/similarity_law(l_upper, h, z0)**3)
          dfdl = dfdl + (h/l_lower/similarity_law(l_lower, h, z0)**3)
          dfdl = dfdl/(2*fd_h)

          l_obukhov = l_obukhov - f/dfdl
c          write(*,*) l_backup, l_obukhov, count, rib

          ! This is an adhoc upper bound for L, at which point we
          ! consider N-R to be diverged
          if (abs(l_obukhov) .gt. 20000 .or.
     $        abs(l_obukhov) .lt. 1e-5) then
            count = 20
          end if
        enddo
        
        
        ! if we did not converge
        if (count .eq. 20) then
          write(*,*) "Unconverged :("
          l_obukhov = l_backup
        endif

        ! compute u* with the new obukhov length
        utau = kappa*magvh/similarity_law(l_obukhov, h, z0) 
      
      endif

      ! Assign tau proportional to the velocity magnitudes at
      ! the sampling point
      wmles_tau(i, 1) = -utau**2*wmles_solh(i, 1)/magvh
      wmles_tau(i, 2) = -utau**2*wmles_solh(i, 2)/magvh
      wmles_tau(i, 3) = -utau**2*wmles_solh(i, 3)/magvh
      end
      

!> @brief Compute correction for the u log law in the convective case
      real function correction(z, l)
      implicit none

      real z, l, xi, pi
!-----------------------------------------------------------------------

      pi = 4*atan(1.0)
      xi = (1.0 - 16.0*z/l)**0.25
      correction = 2*log(0.5*(1 + xi)) + log(0.5*(1 + xi**2)) -
     $             2*atan(xi) + pi/2

      end
      
!> @brief Compute the similarity law for velocity
      real function similarity_law(l_obukhov, h, z0)
      implicit none
      
      real l_obukhov, h, z0
      real correction
!-----------------------------------------------------------------------
      
      similarity_law = log(h/z0) - correction(h, l_obukhov)
     $                           + correction(z0, l_obukhov)
      
      end