!> @brief Algebraic wall-model based on a law of the wall
      subroutine wmles_algebraic()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'WMLES'
      include 'FRAMELP'

      ! linear index for the solution at the boundary
      integer i_linear

      ! Timer function and current time place holder
      real dnekclock, ltim
      
      ! Start the timer 
      ltim = dnekclock()

      ! Sample the solution values at h
      call wmles_sample_distance_based

      i_linear = 0

      do i_linear = 1, wmles_nbpoints
         call wmles_set_momentum_flux(i_linear)

         if (ifheat) then
           call wmles_set_heat_flux(i_linear)
         end if
      enddo

      ! If we use viscosity to impose the shear stress, then
      ! we expect a no-slip velocity bc and we call the set_viscosity
      ! function. Otherwise, we expect a 'sh ' bc and nothing needs
      ! be done here, instead traction is to be assigned in userbc.
      if (wmles_ifviscosity) then
        call wmles_set_viscosity()
      end if

      ! Stop the timer and add to total
      ltim = dnekclock() - ltim
      call mntr_tmr_add(wmles_tmr_tot_id, 1, ltim)
      end subroutine
!=======================================================================
      subroutine wmles_print_averaged_quantities()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'WMLES'
      
      ! alias for wmles_nbpoints to shorten stuff
      integer nbp

      ! inner product and sum-reduction
      real glsc2, glsum, vlsum

      ! just a work array
      real total(10)

      real totalarea, utau
      
      nbp = wmles_nbpoints

      totalarea = glsum(vlsum(wmles_areas(1:nbp), nbp), 1)
      total(1) = glsc2(wmles_solh(:, 1), wmles_areas, nbp)
      total(2) = glsc2(wmles_solh(:, 2), wmles_areas, nbp)
      total(3) = glsc2(wmles_solh(:, 3), wmles_areas, nbp)
      if (ifheat) then
        total(4) = glsc2(wmles_solh(:, 4), wmles_areas, nbp)
        total(5) = glsc2(wmles_solh(:, 5), wmles_areas, nbp)
        total(9) = glsc2(wmles_q, wmles_areas, nbp)
      endif
      total(6) = glsc2(wmles_tau(:, 1), wmles_areas, nbp)
      total(7) = glsc2(wmles_tau(:, 2), wmles_areas, nbp)
      total(8) = glsc2(wmles_tau(:, 3), wmles_areas, nbp)

      utau = sqrt(sqrt(total(6)**2 + total(7)**2 + total(8)**2))

      if (nid .eq. 0) then
        write(*,*) "[WMLES] average tau =", total(6)/totalarea,
     $    total(7)/totalarea, total(8)/totalarea, utau/param(2)
        write(*,*) "[WMLES] average v =", total(1)/totalarea,
     $    total(2)/totalarea, total(3)/totalarea
        if (ifheat) then
        write(*,*) "[WMLES] average q =", total(9)/totalarea
        write(*,*) "[WMLES] average t, ts =", total(4)/totalarea,
     $    total(5)/totalarea
        endif
      
      endif

c      if (ifheat) then
c        totalt = glsum(vlsum(wmles_solh(1:wmles_nbpoints, 4),
c     $  wmles_nbpoints), 1)
c        totalts = glsum(vlsum(wmles_solh(1:wmles_nbpoints, 5),
c     $  wmles_nbpoints), 1)
c      endif

c     totalarea = glsum(totalarea, 1)
c      nwallnodes = iglsum(nwallnodes, 1)
      
      if (nid .eq. 0) then
c        write(*,*) "        [WMLES] Average predicted u_tau = ",
c     $             sqrt(total)
c        write(*,*) "        [WMLES] Average sampled velocity = ",
c     $             totalvxh/nwallnodes, totalvyh/nwallnodes,
c     $             totalvzh/nwallnodes, totalh/nwallnodes
c        if (ifheat) then
c          write(*,*) "        [WMLES] Average predicted q = ",
c     $               totalq/totalarea
c          write(*,*) "        [WMLES] Average sampled t = ",
c     $               totalt/totalarea
c          write(*,*) "        [WMLES] Surface temprature = ",
c     $               wmles_surface_temp, totalts/totalarea
c        endif
      end if
      end subroutine
!=======================================================================
      subroutine wmles_set_h_from_indices()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'WMLES'

      ! Loop ranges for traversing face nodes
      integer frangex1, frangex2, frangey1, frangey2, frangez1, frangez2

      ! counters
      integer ifacex, ifacey, ifacez, iface, ielem

      ! input h rounded to integer
      integer hind
      
      ! the sampling distace for the current node
      real h

      ! Linear index
      integer i_linear

      i_linear = 0

      do ielem=1, lelv
        do iface=1, 6 
          if (boundaryID(iface, ielem) .eq. wmles_wallbid) then

            ! Grab index limits for traversing the face
            call facind(frangex1, frangex2, frangey1,
     $                  frangey2, frangez1, frangez2,
     $                  lx1, ly1, lz1, iface)
            
            do ifacez=frangez1, frangez2
              do ifacey=frangey1, frangey2
                do ifacex=frangex1, frangex2
                  i_linear = i_linear + 1
                  
                  ! round the input h to integer
                  hind = int(wmles_sampling_h(i_linear))
                  
                  if (hind .lt. 1) then
                    write(*,*) "[WMLES] h < 1!"
                  else if (hind .gt. lx1-1) then
                    write(*,*) "[WMLES] h > lx1-1 "
                  endif

                  call wmles_distance_from_index(h, hind,
     $                                           ifacex, ifacey, ifacez,
     $                                           iface, ielem)
                  wmles_sampling_h(i_linear) = h 

                end do
              end do
            end do

          endif

        enddo
      enddo

      end subroutine
!=======================================================================
!> @brief Compute sampling point distance based on index
!! @param[out]   distance        the value of h
!! @param[in]    ifacex          the x index of the face node
!! @param[in]    ifacey          the y index of the face node
!! @param[in]    ifacez          the z index of the face node
!! @param[in]    iface           the index of the face in an element
!! @param[in]    ielem           the index of the element
!! @param[in]    wmles_sampling_idx     the wall-normal index of the sampling point 
      subroutine wmles_distance_from_index(h, sampling_idx,
     $                                     ifacex, ifacey, ifacez,
     $                                     iface, ielem)
      implicit none
     
      include 'SIZE'
      include 'GEOM'
      include 'WMLES'

      real xw, yw, zw, xh, yh, zh, h
      integer ifacex, ifacey, ifacez
      integer iface, ielem
      
      ! input h rounded to integer
      integer sampling_idx

      ! Timer function and current time place holder
      real dnekclock, ltim
      
      ! Start the timer 
      ltim = dnekclock()
      
      ! The wall node
      xw = xm1(ifacex, ifacey, ifacez, ielem)
      yw = ym1(ifacex, ifacey, ifacez, ielem)
      zw = zm1(ifacex, ifacey, ifacez, ielem)

      ! We figure out the index the sampling point
      ! location based on the plane the face is in
      if (iface .eq. 1) then
        ! Face corresponds to x-z plane at y = -1 
        xh = xm1(ifacex, sampling_idx + 1, ifacez, ielem)
        yh = ym1(ifacex, sampling_idx + 1, ifacez, ielem)
        zh = zm1(ifacex, sampling_idx + 1, ifacez, ielem)
      else if (iface .eq. 2) then
        ! Face corresponds to y-z plane at x = 1
        xh = xm1(lx1 - sampling_idx, ifacey, ifacez, ielem)
        yh = ym1(lx1 - sampling_idx, ifacey, ifacez, ielem)
        zh = zm1(lx1 - sampling_idx, ifacey, ifacez, ielem)
      else if (iface .eq. 3) then
        ! Face corresponds to x-z plane at y = 1
        xh = xm1(ifacex, ly1 - sampling_idx, ifacez, ielem)
        yh = ym1(ifacex, ly1 - sampling_idx, ifacez, ielem)
        zh = zm1(ifacex, ly1 - sampling_idx, ifacez, ielem)
      else if (iface .eq. 4) then
        ! Face corresponds to y-z plane at x = -1
        xh = xm1(sampling_idx + 1, ifacey, ifacez, ielem)
        yh = ym1(sampling_idx + 1, ifacey, ifacez, ielem)
        zh = zm1(sampling_idx + 1, ifacey, ifacez, ielem)
      else if (iface .eq. 5) then
        ! Face corresponds to x-y plane at z = -1
        xh = xm1(ifacex, ifacey, sampling_idx + 1, ielem)
        yh = ym1(ifacex, ifacey, sampling_idx + 1, ielem)
        zh = zm1(ifacex, ifacey, sampling_idx + 1, ielem)
      else if (iface .eq. 6) then
        ! Face corresponds to x-y plane at z = 1
        xh = xm1(ifacex, ifacey, lz1 - sampling_idx, ielem)
        yh = ym1(ifacex, ifacey, lz1 - sampling_idx, ielem)
        zh = zm1(ifacex, ifacey, lz1 - sampling_idx, ielem)
      end if

      h = sqrt((xh - xw)**2 + (yh - yw)**2 + (zh - zw)**2)

      ! Stop the timer and add to total
      ltim = dnekclock() - ltim
      call mntr_tmr_add(wmles_tmr_sampling_id, 1, ltim)
      end subroutine
!=======================================================================
      subroutine wmles_sample_distance_based()
      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'WMLES'

      ! The components of the inward normal to a wall-face + counter
      integer i, j

      ! store scalar product of velocity with normal
      real vdotn

      ! working arrays for interpolation as per
      ! inter_nfld documentation
      real    rwk(NMAX_BOUNDARY_POINTS, ldim+1)
      integer iwk(NMAX_BOUNDARY_POINTS, 3)

      ! sampled velocity and temperature fields
      real sampled_fields(lx1, ly1, lz1, lelt, 5)

      ! temp array to hold the sampled fields
      real temp_solh(NMAX_BOUNDARY_POINTS, 5)

      ! handle variable necessary for the interpolation
      integer intp_h

      ! number of fields to sample, depends on whether we have t
      integer n_fields

      ! weight for time-averaging
      real eps

      ! Timer function and current time place holder
      real dnekclock, ltim

      ! Start the timer 
      ltim = dnekclock()

      if (ISTEP .eq. 0) then
        eps = 1.0
      else
        eps = 1/wmles_navrg
      endif

      if (ifheat) then
          ! temperature and surface temperature
          n_fields = 5
      else 
          n_fields = 3
      endif

      call copy(sampled_fields(:, :, :, :, 1), vx, lx1*ly1*lz1*nelt)
      call copy(sampled_fields(:, :, :, :, 2), vy, lx1*ly1*lz1*nelt)
      call copy(sampled_fields(:, :, :, :, 3), vz, lx1*ly1*lz1*nelt)

      if (ifheat) then
        call copy(sampled_fields(:, :, :, :, 4),
     $    t(:, :, :, :, 1), lx1*ly1*lz1*nelt)
      endif

        call interp_nfld(
     $    temp_solh(:wmles_nbpoints, 1:n_fields), !will store the sampled points
     $    sampled_fields(:, :, :, :, 1:n_fields),     !will sample these
     $    n_fields,                                !number of fields
     $    wmles_sampling_points(1:wmles_nbpoints, 1), ! x
     $    wmles_sampling_points(1:wmles_nbpoints, 2), ! y
     $    wmles_sampling_points(1:wmles_nbpoints, 3), ! z
     $    wmles_nbpoints,                       ! number of points
     $    iwk, rwk, NMAX_BOUNDARY_POINTS,          ! work array stuff
     $    wmles_iffind,                            ! wether to find points
     $    wmles_interpolation_handle)              ! handle

        ! Turn off the flag, actually needs to run once, but whatever
        wmles_iffind = .false.

        ! project sampled velocity on the face plane
        ! and the time-average
        do i=1, wmles_nbpoints
          vdotn = 0
          do j=1, 3
            vdotn = vdotn + temp_solh(i, j)*wmles_normals(i, j)
          enddo

          do j=1, 3
            temp_solh(i, j) = temp_solh(i, j) - vdotn*wmles_normals(i,j)
            wmles_solh(i, j) = (1-eps)*wmles_solh(i, j) +
     $        eps*temp_solh(i, j)
          enddo
          
          if (ifheat) then
            ! time-average the temperature
            wmles_solh(i, 4) = (1-eps)*wmles_solh(i, 4) +
     $        eps*temp_solh(i, 4)
            
            ! sample wall temperature (should we average it?)
            wmles_solh(i, 5) =
     $        t(wmles_indices(i, 1),
     $          wmles_indices(i, 2),
     $          wmles_indices(i, 3),
     $          wmles_indices(i, 4), 1)
          endif  
        enddo
      
      ! Stop the timer and add to total
      ltim = dnekclock() - ltim
      call mntr_tmr_add(wmles_tmr_sampling_id, 1, ltim)
      end subroutine
!=======================================================================
!> @brief Set the artificial viscosity att the wall
      subroutine wmles_set_viscosity()
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'GEOM'
      include 'WMLES'

      ! face normal components
      real normalx, normaly, normalz

      ! the rate of strain tensor
      real sij (lx1,ly1,lz1,6,lelv)

      ! sij components at a node
      real sxx, sxy, sxz, syx, syy, syz, szx, szy, szz
      
      ! S_ij n_j dot product components at a wall node
      real snx, sny, snz, magsij
      
      ! normal wall stress
      real snormal
      
      ! the magnitude of the wall shear stress predicted by the model
      real magtau

      ! stuff necessary to compute sij
      integer lr
      parameter (lr=lx1*ly1*lz1)

      common /scruz/         ur(lr),us(lr),ut(lr)
     $                     , vr(lr),vs(lr),vt(lr)
     $                     , wr(lr),ws(lr),wt(lr)
      real ur, us, ut, vr, vs, vt, wr, ws, wt

      ! counters for boundary traversing
      integer i_linear, ielem, iface, inorm
      integer frangex1, frangex2, frangey1, frangey2, frangez1, frangez2
      integer ifacex, ifacey, ifacez
      
      ! compute the rate of strain. 
      call comp_sij(sij,6,vx,vy,vz,ur,us,ut,vr,vs,vt,wr,ws,wt)

      i_linear = 0
      ! compute the articifial viscosity at each wall face
      do ielem = 1, nelt
        do iface=1, 6
          if (boundaryID(iface, ielem) .eq. wmles_wallbid) then
            call facind(frangex1, frangex2, frangey1,
     $                  frangey2, frangez1, frangez2,
     $                  lx1, ly1, lz1, iface)
            inorm = 0
            do ifacez=frangez1, frangez2
              do ifacey=frangey1, frangey2
                do ifacex=frangex1, frangex2
                  i_linear = i_linear + 1
                  inorm = inorm + 1
                  
                  ! grab the normal
                  normalx = unx(inorm, 1, iface, ielem)
                  normaly = uny(inorm, 1, iface, ielem)
                  normalz = unz(inorm, 1, iface, ielem)

                  ! grab the components of sij
                  ! (some redundancy here due to symmetry, but ok)
                  sxx = sij(ifacex, ifacey, ifacez, 1, ielem)
                  sxy = sij(ifacex, ifacey, ifacez, 4, ielem)
                  sxz = sij(ifacex, ifacey, ifacez, 6, ielem)
                  syx = sij(ifacex, ifacey, ifacez, 4, ielem)
                  syy = sij(ifacex, ifacey, ifacez, 2, ielem)
                  syz = sij(ifacex, ifacey, ifacez, 5, ielem)
                  szx = sij(ifacex, ifacey, ifacez, 6, ielem)
                  szy = sij(ifacex, ifacey, ifacez, 5, ielem)
                  szz = sij(ifacex, ifacey, ifacez, 3, ielem)
                  
                  ! compute the dot product S_ijn_j at the wall
                  snx = sxx*normalx + sxy*normaly + sxz*normalz
                  sny = syx*normalx + syy*normaly + syz*normalz
                  snz = szx*normalx + szy*normaly + szz*normalz

                  ! we want only the shear stress
                  ! first compute the nomrmal stress
                  snormal = (snx*normalx + sny*normaly +snz*normalz)

                  ! and now subtract it
                  snx = snx - normalx*snormal
                  sny = sny - normaly*snormal
                  snz = snz - normalz*snormal
                  
                  ! reconstruct the magnitude of tau from the components
                  magtau =
     $              sqrt
     $              (
     $                  wmles_tau(i_linear, 1)**2
     $                + wmles_tau(i_linear, 2)**2
     $                + wmles_tau(i_linear, 3)**2
     $              )

                  magsij = sqrt(snx**2 + sny**2 + snz**2)
                  vdiff(ifacex, ifacey, ifacez, ielem, 1) = 
     $              magtau/(magsij + 1e-10)
                end do
              end do
            end do
            
          endif
        
        end do
      end do

      end subroutine

