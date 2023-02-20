!> @brief Sets both utau and q for the stable ABL case
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

      real utau, g, th, ts, q
      
      ! newton iteration stuff
      real f, dfdl

      ! log law parameters
      real kappa, z0, z1
      
      ! parameters in front of the correction functions for heat and 
      ! momentum.
      real a, b
      
      ! similarity law and coorection functions for it
      real similarity_law_u, similarity_law_q
      
      ! check wether the Newton iteration diverged
      logical diverged

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
c      magvh = wmles_solh(i, 1)**2 +
c     $        wmles_solh(i, 2)**2 +
c     $        wmles_solh(i, 3)**2

      magvh = wmles_uh_average(1)**2 +
     $        wmles_uh_average(2)**2 +
     $        wmles_uh_average(3)**2

      magvh = sqrt(magvh)

c      write(*,*) i
c      if (i .eq. 1 .and. nid .eq. 64) then
c        write(*,*) "Sampled u", magvh, nid
c      endif

      ! Temperature
c      th = wmles_solh(i, 4)
      th = wmles_th_average
      
      ! Surface temperature
      if (wmles_ifviscosity) then
        ts = wmles_solh(i, 5)
        ts = wmles_ts_average
      else
        ts = wmles_surface_temp
      endif

      ! Get uncorrected utau for a first guess
      utau = magvh*kappa/log(h/z0)
      q = kappa*utau*(ts - th)/log(h/z1)
      
c      write(*,*) utau, q

      a = 5.0
      b = 5.0
      
      diverged = .false.
      
      if (ISTEP .gt. 3) then

        rib = g*h/th*(th - ts)/magvh**2

        ! Obukhov l based on the previous-step utau and q
        l_obukhov = -(wmles_theta0*utau**3)/(kappa*g*q)
        
        wmles_lobukhov(i) = l_obukhov
        
        if (l_obukhov .ge. 20000) then
            diverged = .true.
        endif

c        write(*,*) "l", l_obukhov
        
        l_old = 0
        count = 0
c       !write(*,*) "ERR",  abs(l_old - l_obukhov)/l_obukhov 
        do while ((abs(l_old - l_obukhov)/abs(l_obukhov) .gt. 1e-3)
     $           .and. (count .lt. 20) .and. (.not. diverged))

          l_old = l_obukhov
          count = count + 1

          ! for the central diff for evaluating dfdl
c          l_upper = l_obukhov + 1e-3*l_obukhov
c          l_lower = l_obukhov - 1e-3*l_obukhov

          f = rib - h/l_obukhov*
     $              similarity_law_q(l_obukhov, h, z0)/
     $              similarity_law_u(l_obukhov, h, z0)**2

          dfdl = ((h*log(h/z0)*(2*a*h - b*h + l_obukhov*log(h/z0)))/
     $            (b*h + l_obukhov*log(h/z0))**3)
      
c          dfdl = -h/l_upper*
c     $             similarity_law_u(l_upper, h, z0)/
c     $             similarity_law_q(l_upper, h, z0)**2
c          dfdl = dfdl + h/l_lower*
c     $            similarity_law_u(l_lower, h, z0)/
c     $            similarity_law_q(l_lower, h, z0)**2
c          dfdl = dfdl/(l_upper - l_lower)/2

          l_obukhov = l_obukhov - f/dfdl
c          write(*,*) count, l_obukhov
          
          ! This is an adhoc upper bound for L, at which point we
          ! consider N-R to be diverged
          if (abs(l_obukhov) > 20000 ) then
            diverged = .true.
          end if
        enddo
        
        if (count .eq. 20) then
            diverged = .true.
        endif
        
        !write(*,*) l_backup, l_obukhov, count, rib
        
        ! if we did not converge
        if (diverged) then
c          write(*,*) "Unconverged :("
        else
              
          if (l_obukhov < 5) then
c            write(*,*) l_obukhov, magvh, th, count, wmles_uh_average
          endif
          ! udate stored values
          wmles_lobukhov(i) = l_obukhov

          ! compute u* with the new obukhov length
          utau = kappa*magvh/similarity_law_u(l_obukhov, h, z0) 

          ! compute q with the new obukhov length
          q = kappa*utau*(ts - th)/similarity_law_q(l_obukhov, h, z0) 
        endif
      endif

      ! Assign tau proportional to the velocity magnitudes at
      ! the sampling point
      wmles_tau(i, 1) = -utau**2*wmles_solh(i, 1)/magvh
      wmles_tau(i, 2) = -utau**2*wmles_solh(i, 2)/magvh
      wmles_tau(i, 3) = -utau**2*wmles_solh(i, 3)/magvh
      wmles_q(i) = q
      end subroutine
      

!> @brief Compute correction for the u log law in the stable case
      real function correction_u(z, l)
      implicit none

      real z, l
!-----------------------------------------------------------------------

      correction_u = -5*z/l
      end function

!> @brief Compute correction for the u log law in the convective case
      real function correction_q(z, l)
      implicit none

      real z, l
!-----------------------------------------------------------------------

      correction_q = -5*z/l
      end function
      
!> @brief Compute the similarity law for velocity
      real function similarity_law_u(l_obukhov, h, z0)
      implicit none
      
      real l_obukhov, h, z0
      real correction_u
!-----------------------------------------------------------------------
      
      similarity_law_u = log(h/z0) - correction_u(h, l_obukhov)
c     $                             + correction_u(z0, l_obukhov)
      
      end function

!> @brief Compute the similarity law for heat
      real function similarity_law_q(l_obukhov, h, z0)
      implicit none
      
      real l_obukhov, h, z0
      real correction_q
!-----------------------------------------------------------------------
      
      similarity_law_q = log(h/z0) - correction_q(h, l_obukhov)
c     $                             + correction_q(z0, l_obukhov)
      

      end function
      
      
!> @brief Dummy for the function to compute q
      subroutine wmles_set_heat_flux(h, ix, iy, iz, ie)
      end subroutine