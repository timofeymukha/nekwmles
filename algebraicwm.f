!> @brief Algebraic wall-model based on a law of the wall
      subroutine algebraic_wm()
      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'SOLN'
      include 'GEOM'
      include 'WMLES'
      include 'FRAMELP'

      ! Some counters to traverse the velocity mesh
      integer ielem, iface, igllx, iglly, igllz

      ! Some counters to traversing face nodes
      integer ifacex, ifacey, ifacez

      ! Loop ranges for traversing face nodes
      integer frangex1, frangex2, frangey1, frangey2, frangez1, frangez2

      ! A string for a boundary condition
      character cb*3 
      
      ! Number of wall-nodes
      integer nwallnodes

      ! Values of gll node coordinates
      real xgll, ygll, zgll

      ! The components of the inward normal to a wall-face + counter
      real normalx, normaly, normalz
      integer inorm

      ! The coordinates, velocity, and temprature at the sampling point
      real xh, yh, zh, vxh, vyh, vzh, vhndot, magvh, th, ts
      
      ! For averaging vxh for debug purposes
      real totalvxh, totalvyh, totalvzh, totalh, totalarea, magutau
      real totalq, totalt, totalts

      ! The distance to the sampling point
      real h

      ! Time-averaging parameter
      real eps

      real nu, utau, totalutau, newton
      
      ! A guess for utau sent to the iterative solver
      real guess
      
      ! Timer function and current time place holder
      real dnekclock, ltim
      
      ! Function for reduction across processors
      integer iglsum
      real glsum
      
      ! Start the timer 
      ltim = dnekclock()

      ! Grab various simulation parameters
      nu = param(2)

      nwallnodes = 0
      totalutau = 0
      totalt = 0
      totalts = 0
      totalq = 0
      totalvxh = 0
      totalvyh = 0
      totalvzh = 0
      totalh = 0
      totalarea = 0

      if (ISTEP .eq. 0) then
        eps = 1.0
      else
        eps = 1/wmles_navrg
      endif

      call sample_distance_based
c      call sample_distance_based
c      call sample_distance_based
      
      do ielem=1, lelv
        do iface=1, 6 

          ! Get the velocity bc type for the face
          cb = cbc(iface, ielem, 1)

          if (boundaryID(iface, ielem) .eq. wallbid) then

            ! Grab index limits for traversing the face
            call facind(frangex1, frangex2, frangey1,
     $                  frangey2, frangez1, frangez2,
     $                  lx1, ly1, lz1, iface)
            
            inorm = 0

            do ifacez=frangez1, frangez2
              do ifacey=frangey1, frangey2
                do ifacex=frangex1, frangex2

                  inorm = inorm + 1
                  xgll = xm1(ifacex, ifacey, ifacez, ielem)
                  ygll = ym1(ifacex, ifacey, ifacez, ielem)
                  zgll = zm1(ifacex, ifacey, ifacez, ielem)

                  ! inward face normal
                  normalx =  unx(inorm, 1, iface, ielem)
                  normaly = -uny(inorm, 1, iface, ielem)
                  normalz =  unz(inorm, 1, iface, ielem)

                  totalarea = totalarea + area(inorm, 1, iface, ielem)

                  call sample_index_based(xh, yh, zh, 
     $                                    vxh, vyh, vzh, th, ts,
     $                                    ifacex, ifacey, ifacez,
     $                                    iface, ielem)

                  ! Vector between face node and sampling point
                  xh = xh - xgll
                  yh = yh - ygll
                  zh = zh - zgll
                  
                  ! project onto the local face normal
                  h = abs(xh*normalx + yh*normaly + zh*normalz)
                  
                  ! project sampled velocity on the face plane
                  vhndot =
     $              (vxh*normalx + vyh*normaly + vzh*normalz)
     
                  vxh = vxh - normalx*vhndot
                  vyh = vyh - normaly*vhndot
                  vzh = vzh - normalz*vhndot


                  ! time-averaging the input 
                  vh(1, ifacex, ifacey, ifacez, ielem) =
     $               (1-eps)*vh(1, ifacex, ifacey, ifacez, ielem) +
     $               eps*vxh
                  vh(2, ifacex, ifacey, ifacez, ielem) =
     $               (1-eps)*vh(2, ifacex, ifacey, ifacez, ielem) +
     $               eps*vyh
                  vh(3, ifacex, ifacey, ifacez, ielem) =
     $               (1-eps)*vh(3, ifacex, ifacey, ifacez, ielem) +
     $               eps*vzh
                  temph(ifacex, ifacey, ifacez, ielem) =
     $               (1-eps)*temph(ifacex, ifacey, ifacez, ielem) +
     $               eps*th
                  temps(ifacex, ifacey, ifacez, ielem) =
     $               (1-eps)*temps(ifacex, ifacey, ifacez, ielem) +
     $               eps*ts

                  ! for brevity reassign back
                  vxh = vh(1, ifacex, ifacey, ifacez, ielem)
                  vyh = vh(2, ifacex, ifacey, ifacez, ielem)
                  vzh = vh(3, ifacex, ifacey, ifacez, ielem)
                  th = temph(ifacex, ifacey, ifacez, ielem)
                  ts = temps(ifacex, ifacey, ifacez, ielem)

                  ! Sum up vh for debug output
                  totalvxh = totalvxh + vxh
                  totalvyh = totalvyh + vyh
                  totalvzh = totalvzh + vzh
                  totalh = totalh + h

                  call set_momentum_flux(h, ifacex, ifacey,
     $                                     ifacez, ielem)
                  if (ifheat) then
                    call set_heat_flux(h, ifacex, ifacey, ifacez, ielem)
                    totalq = 
     $                totalq + heat_flux(ifacex, ifacey, ifacez,ielem)*
     $                         area(inorm, 1, iface, ielem) 
                    totalt = totalt + th*area(inorm, 1, iface, ielem)
                    totalts = totalts + ts*area(inorm, 1, iface, ielem)
                  end if

                  magutau = tau(1, ifacex, ifacey, ifacez, ielem)**2 +
     $                      tau(2, ifacex, ifacey, ifacez, ielem)**2 +
     $                      tau(3, ifacex, ifacey, ifacez, ielem)**2
                  magutau = sqrt(sqrt(magutau))
                  
                  ! actually sum tau here, sqrt after averaging
                  totalutau = totalutau + magutau**2*area(inorm, 1,
     $             iface, ielem)

                  nwallnodes = nwallnodes + 1

                end do
              end do
            end do

          endif

        enddo
      enddo

      totalutau = glsum(totalutau, 1)
      totalq = glsum(totalq, 1)
      totalt = glsum(totalt, 1)
      totalts = glsum(totalts, 1)
      totalvxh = glsum(totalvxh, 1)
      totalvyh = glsum(totalvyh, 1)
      totalvzh = glsum(totalvzh, 1)
      totalh = glsum(totalh, 1)
      totalarea = glsum(totalarea, 1)
      nwallnodes = iglsum(nwallnodes, 1)
      
      if (nid .eq. 0) then
        write(*,*) "        [WMLES] Average predicted u_tau = ",
     $             sqrt(totalutau/totalarea)
        write(*,*) "        [WMLES] Average sampled velocity = ",
     $             totalvxh/nwallnodes, totalvyh/nwallnodes,
     $             totalvzh/nwallnodes, totalh/nwallnodes
        if (ifheat) then
          write(*,*) "        [WMLES] Average predicted q = ",
     $               totalq/totalarea
          write(*,*) "        [WMLES] Average sampled t = ",
     $               totalt/totalarea
          write(*,*) "        [WMLES] Surface temprature = ",
     $               wmles_surface_temp, totalts/totalarea
        endif
      end if
      
      ! If we use viscosity to impose the shear stress, then
      ! we expect a no-slip velocity bc and we call the set_viscosity
      ! function. Otherwise, we expect a 'sh ' bc and nothing needs
      ! be done here, isntead traction is to be assigned in userbc.
      if (ifviscosity) then
        call set_viscosity()
      end if

      ! Stop the timer and add to total
      ltim = dnekclock() - ltim
      call mntr_tmr_add(wmles_tmr_tot_id, 1, ltim)
      end subroutine
!=======================================================================
!> @brief Compute sampling point coordinates and sample velocity
!! @param[out]   xh              x coordinate of the sampling point
!! @param[out]   yh              y coordinate of the sampling point
!! @param[out]   zh              z coordinate of the sampling point
!! @param[out]   vxh             x component of the sampled velocity
!! @param[out]   vyh             y component of the sampled velocity
!! @param[out]   vzh             z component of the sampled velocity
!! @param[out]   th              the sampled temperature
!! @param[out]   ts              the sampled surface temperature
!! @param[in]    ifacex          the x index of the face node
!! @param[in]    ifacey          the y index of the face node
!! @param[in]    ifacez          the z index of the face node
!! @param[in]    iface           the index of the face in an element
!! @param[in]    ielem           the index of the element
!! @param[in]    samplingidx     the wall-normal index of the sampling point 
      subroutine sample_index_based(xh, yh, zh,
     $                              vxh, vyh, vzh, th, ts,
     $                              ifacex, ifacey, ifacez,
     $                              iface, ielem)
      implicit none
     
      include 'SIZE'
      include 'GEOM'
      include 'SOLN'
      include 'WMLES'

      real xh, yh, zh
      real vxh, vyh, vzh, th, ts
      integer ifacex, ifacey, ifacez
      integer iface, ielem

      ! Timer function and current time place holder
      real dnekclock, ltim
      
      ! Start the timer 
      ltim = dnekclock()


      ! We figure out the index the sampling point
      ! location based on the plane the face is in
      if (iface .eq. 1) then
        ! Face corresponds to x-z plane at y = -1 
        xh = xm1(ifacex, samplingidx + 1, ifacez, ielem)
        yh = ym1(ifacex, samplingidx + 1, ifacez, ielem)
        zh = zm1(ifacex, samplingidx + 1, ifacez, ielem)
        
        vxh = vx(ifacex, samplingidx + 1, ifacez, ielem)
        vyh = vy(ifacex, samplingidx + 1, ifacez, ielem)
        vzh = vz(ifacex, samplingidx + 1, ifacez, ielem)
        th = t(ifacex, samplingidx + 1, ifacez, ielem, 1)
        ts = t(ifacex, 1, ifacez, ielem, 1)

      else if (iface .eq. 2) then
        ! Face corresponds to y-z plane at x = 1
        xh = xm1(lx1 - samplingidx, ifacey, ifacez, ielem)
        yh = ym1(lx1 - samplingidx, ifacey, ifacez, ielem)
        zh = zm1(lx1 - samplingidx, ifacey, ifacez, ielem)

        vxh = vx(lx1 - samplingidx, ifacey, ifacez, ielem)
        vyh = vy(lx1 - samplingidx, ifacey, ifacez, ielem)
        vzh = vz(lx1 - samplingidx, ifacey, ifacez, ielem)
        th = t(lx1 - samplingidx, ifacey, ifacez, ielem, 1)
        ts = t(lx1, ifacey, ifacez, ielem, 1)
      else if (iface .eq. 3) then
        ! Face corresponds to x-z plane at y = 1
        xh = xm1(ifacex, ly1 - samplingidx, ifacez, ielem)
        yh = ym1(ifacex, ly1 - samplingidx, ifacez, ielem)
        zh = zm1(ifacex, ly1 - samplingidx, ifacez, ielem)

        vxh = vx(ifacex, ly1 - samplingidx, ifacez, ielem)
        vyh = vy(ifacex, ly1 - samplingidx, ifacez, ielem)
        vzh = vz(ifacex, ly1 - samplingidx, ifacez, ielem)
        th = t(ifacex, ly1 - samplingidx, ifacez, ielem, 1)
        ts = t(ifacex, ly1, ifacez, ielem, 1)
      else if (iface .eq. 4) then
        ! Face corresponds to y-z plane at x = -1
        xh = xm1(samplingidx + 1, ifacey, ifacez, ielem)
        yh = ym1(samplingidx + 1, ifacey, ifacez, ielem)
        zh = zm1(samplingidx + 1, ifacey, ifacez, ielem)

        vxh = vx(samplingidx + 1, ifacey, ifacez, ielem)
        vyh = vy(samplingidx + 1, ifacey, ifacez, ielem)
        vzh = vz(samplingidx + 1, ifacey, ifacez, ielem)
        th = t(samplingidx + 1, ifacey, ifacez, ielem, 1)
        ts = t(1, ifacey, ifacez, ielem, 1)
      else if (iface .eq. 5) then
        ! Face corresponds to x-y plane at z = -1
        xh = xm1(ifacex, ifacey, samplingidx + 1, ielem)
        yh = ym1(ifacex, ifacey, samplingidx + 1, ielem)
        zh = zm1(ifacex, ifacey, samplingidx + 1, ielem)

        vxh = vx(ifacex, ifacey, samplingidx + 1, ielem)
        vyh = vy(ifacex, ifacey, samplingidx + 1, ielem)
        vzh = vz(ifacex, ifacey, samplingidx + 1, ielem)
        th = t(ifacex, ifacey, samplingidx + 1, ielem, 1)
        ts = t(ifacex, ifacey, 1, ielem, 1)
      else if (iface .eq. 6) then
        ! Face corresponds to x-y plane at z = 1
        xh = xm1(ifacex, ifacey, lz1 - samplingidx, ielem)
        yh = ym1(ifacex, ifacey, lz1 - samplingidx, ielem)
        zh = zm1(ifacex, ifacey, lz1 - samplingidx, ielem)

        vxh = vx(ifacex, ifacey, lz1 - samplingidx, ielem)
        vyh = vy(ifacex, ifacey, lz1 - samplingidx, ielem)
        vzh = vz(ifacex, ifacey, lz1 - samplingidx, ielem)
        th = t(ifacex, ifacey, lz1 - samplingidx, ielem, 1)
        ts = t(ifacex, ifacey, lz1, ielem, 1)
      end if

      ! Stop the timer and add to total
      ltim = dnekclock() - ltim
      call mntr_tmr_add(wmles_tmr_sampling_id, 1, ltim)
      end subroutine
!=======================================================================
      subroutine find_sampling_points()
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'WMLES'
      ! Values of gll node coordinates
      real xgll, ygll, zgll

      ! The components of the inward normal to a wall-face + counter
      real normalx, normaly, normalz
      integer inorm

      ! Loop ranges for traversing face nodes
      integer frangex1, frangex2, frangey1, frangey2, frangez1, frangez2

      real xh, yh, zh
      real vxh, vyh, vzh, th, ts
      integer ifacex, ifacey, ifacez
      integer iface, ielem

      ! Linear index
      integer n_sampling

      ! Timer function and current time place holder
      real dnekclock, ltim

      ! Start the timer 
      ltim = dnekclock()

      n_sampling = 0

      do ielem=1, lelv
        do iface=1, 6 

          if (boundaryID(iface, ielem) .eq. wallbid) then

            ! Grab index limits for traversing the face
            call facind(frangex1, frangex2, frangey1,
     $                  frangey2, frangez1, frangez2,
     $                  lx1, ly1, lz1, iface)
            
            ! a linear index for the loop below
            inorm = 0

            do ifacez=frangez1, frangez2
              do ifacey=frangey1, frangey2
                do ifacex=frangex1, frangex2
                  n_sampling = n_sampling + 1
c
                  inorm = inorm + 1
                  xgll = xm1(ifacex, ifacey, ifacez, ielem)
                  ygll = ym1(ifacex, ifacey, ifacez, ielem)
                  zgll = zm1(ifacex, ifacey, ifacez, ielem)

                  ! inward face normal
                  normalx =  unx(inorm, 1, iface, ielem)
                  normaly = -uny(inorm, 1, iface, ielem)
                  normalz =  unz(inorm, 1, iface, ielem)

                  sampling_points(n_sampling, 1) = 
     $               xgll + normalx*sampling_h(n_sampling)
                  sampling_points(n_sampling, 2) = 
     $               ygll + normaly*sampling_h(n_sampling)
                  sampling_points(n_sampling, 3) = 
     $               zgll + normalz*sampling_h(n_sampling)

                end do
              end do
            end do

          endif

        enddo
      enddo


      n_boundary_points = n_sampling

c      write(*,*) sampling_points(:n_boundary_points, 2)
      ! Stop the timer and add to total
      ltim = dnekclock() - ltim
      call mntr_tmr_add(wmles_tmr_sampling_id, 1, ltim)
      end subroutine
!=======================================================================
      subroutine sample_distance_based()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'WMLES'

      ! Values of gll node coordinates
      real xgll, ygll, zgll

      ! The components of the inward normal to a wall-face + counter
      real normalx, normaly, normalz
      integer inorm

      real xh, yh, zh
      real vxh, vyh, vzh, th, ts
      integer ifacex, ifacey, ifacez
      integer iface, ielem

      ! working arrays for interpolation as per
      ! inter_nfld documentation
      real    rwk(NMAX_BOUNDARY_POINTS, ldim+1)
      integer iwk(NMAX_BOUNDARY_POINTS, 3)

      ! sampled velocity and temperature fields
      real sampled_fields(lx1, ly1, lz1, lelt, 3)

      ! handle variable necessary for the interpolation
      integer intp_h

      ! Timer function and current time place holder
      real dnekclock, ltim

      ! Start the timer 
      ltim = dnekclock()

c      write(*,*) sampling_points(:n_boundary_points, 2)

      call copy(sampled_fields(:, :, :, :, 1), vx, lx1*ly1*lz1*nelt)
      call copy(sampled_fields(:, :, :, :, 2), vy, lx1*ly1*lz1*nelt)
      call copy(sampled_fields(:, :, :, :, 3), vz, lx1*ly1*lz1*nelt)
c      call copy(sampled_fields(:, :, :, :, 4), t, lx1*ly1*lz1*nelt)

c      write(*,*) nid, "REAL NUMBER OF POINTS", n_boundary_points

c      call outpost(sampled_fields(:,:,:,:,1),
c     $             sampled_fields(:,:,:,:,2),
c     $             sampled_fields(:,:,:,:,3),
c     $             xm2,
c     $             sampled_fields(:,:,:,:,1),
c     $             'kek') 


c          write(*,*) "Running interp_nfld"
        call interp_nfld(
     $    solh(:n_boundary_points, 1:3), !will store the sampled points
     $    sampled_fields,         !will sample these
     $    3,                      !number of fields
     $    sampling_points(1:n_boundary_points, 1), ! x
     $    sampling_points(1:n_boundary_points, 2), ! y
     $    sampling_points(1:n_boundary_points, 3), ! z
     $    n_boundary_points,                       ! number of points
     $    iwk, rwk, NMAX_BOUNDARY_POINTS,      ! work array stuff
     $    wmles_iffind,                            ! wether to find points
     $    wmles_interpolation_handle)                           ! handle

        ! Turn off flag, actually needs to run once, but whatever
        wmles_iffind = .false.

c      if (n_sampling > 0) then
c         write(*,*) "nid", nid, "points", xsampling_points(1:3),
c     $    ysampling_points(1:3),
c     $    "sampled", solh(1:200, 1), solh(1:3, 2)
c       endif

c      write(*,*) "Ran interp_nfld"



      ! Stop the timer and add to total
      ltim = dnekclock() - ltim
      call mntr_tmr_add(wmles_tmr_sampling_id, 1, ltim)
      end subroutine
!=======================================================================
!> @brief Set the artificial viscosity att the wall
      subroutine set_viscosity()
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
      integer ielem, iface, inorm
      integer frangex1, frangex2, frangey1, frangey2, frangez1, frangez2
      integer ifacex, ifacey, ifacez
      
      ! compute the rate of strain. 
      call comp_sij(sij,6,vx,vy,vz,ur,us,ut,vr,vs,vt,wr,ws,wt)

      ! compute the articifial viscosity at each wall face
      do ielem = 1, nelt
        do iface=1, 6
          if (boundaryID(iface, ielem) .eq. wallbid) then
            call facind(frangex1, frangex2, frangey1,
     $                  frangey2, frangez1, frangez2,
     $                  lx1, ly1, lz1, iface)
            inorm = 0
            do ifacez=frangez1, frangez2
              do ifacey=frangey1, frangey2
                do ifacex=frangex1, frangex2
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
     $                  tau(1,ifacex, ifacey, ifacez, ielem)**2
     $                + tau(2,ifacex, ifacey, ifacez, ielem)**2
     $                + tau(3,ifacex, ifacey, ifacez, ielem)**2
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

