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

      ! The components of the inward normal to a wall-face
      real normalx, normaly, normalz

      ! The coordinates and velocity at the sampling point
      real xh, yh, zh, vxh, vyh, vzh, vhndot, magvh
      
      ! For averaging vxh for debug purposes
      real totalvxh, totalvyh, totalvzh

      ! The distance to the sampling point
      real h

      real nu, utau, totalutau, newton
      
      ! A guess for utau sent to the iterative solver
      real guess
      
      ! Laws of the wall
      external spalding_value, spalding_derivative
      
      ! Grab various simulation parameters
      nu = param(2)

      nwallnodes = 0
      totalutau = 0
      totalvxh = 0
      totalvyh = 0
      totalvzh = 0
      
      do ielem=1, lelv
        do iface=1, 6 

          ! Get the velocity bc type for the face
          cb = cbc(iface, ielem, 1)

          if ((cb .eq. 'sh ')) then

            ! Grab index limits for traversing the face
            call facind(frangex1, frangex2, frangey1,
     $                  frangey2, frangez1, frangez2,
     $                  lx1, ly1, lz1, iface)
            do ifacez=frangez1, frangez2
              do ifacey=frangey1, frangey2
                do ifacex=frangex1, frangex2
c                  write(*,*) ielem, iface, ifacez, ifacey, ifacex
                  xgll = xm1(ifacex, ifacey, ifacez, ielem)
                  ygll = ym1(ifacex, ifacey, ifacez, ielem)
                  zgll = zm1(ifacex, ifacey, ifacez, ielem)

                  ! inward face normal
                  normalx = unx(ifacex, ifacey, iface, ielem)
                  normaly = -uny(ifacex, ifacey, iface, ielem)
                  normalz = unz(ifacex, ifacey, iface, ielem)
                  
                  call sample_index_based(xh, yh, zh, 
     $                                    vxh, vyh, vzh,
     $                                    ifacex, ifacey, ifacez,
     $                                    iface, ielem, samplingidx)


                  ! Vector between face node and sampling point
                  xh = xh - xgll
                  yh = yh - ygll
                  zh = zh - zgll
                  
                  ! project onto the local face normal
                  h = abs(xh*normalx + yh*normaly + zh*normalz)
                  
                  ! project sampled velocity on the face plane
                  vhndot = (vxh*normalx + vyh*normaly + vzh*normalz)
                  vxh = vxh - normalx*vhndot
                  vyh = vyh - normaly*vhndot
                  vzh = vzh - normalz*vhndot
                  
                  ! Sum up vh for debug output
                  totalvxh = totalvxh + vxh
                  totalvyh = totalvyh + vyh
                  totalvzh = totalvzh + vzh

                  ! Magnitude of the sampled velocity
                  magvh = sqrt(vxh**2 + vyh**2 + vzh**2)
                 
                  ! take the stress mag at the previous step as a guess
                  ! inital value at simulation start is the GUESS param
                  guess = tau(1, ifacex, ifacey, ifacez, ielem)**2 +
     $                    tau(2, ifacex, ifacey, ifacez, ielem)**2 +
     $                    tau(3, ifacex, ifacey, ifacez, ielem)**2

                  ! Double sqrt to get utau.
                  guess = sqrt(sqrt(guess))
                  
                  utau = newton(spalding_value, spalding_derivative,
     $                          magvh, h, guess,
     $                          1e-4, 50)

                  totalutau = totalutau + utau
                  nwallnodes = nwallnodes + 1
                  
                  ! Assign proportional to the velocity magnitudes at
                  ! the sampling point
                  tau(1, ifacex, ifacey, ifacez, ielem) = 
     $              -utau**2*vxh/magvh
                  tau(2, ifacex, ifacey, ifacez, ielem) = 
     $              -utau**2*vyh/magvh
                  tau(3, ifacex, ifacey, ifacez, ielem) = 
     $              -utau**2*vzh/magvh

                end do
              end do
            end do

          endif

        enddo
      enddo
      
      if (nid .eq. 0) then
      write(*,*) "        [WMLES] Average predicted Re_tau = ",
     $           totalutau/nwallnodes/nu 
      write(*,*) "        [WMLES] Average sampled velocity = ",
     $           totalvxh/nwallnodes, totalvyh/nwallnodes,
     $           totalvzh/nwallnodes, h 
      end if
      end subroutine
!=======================================================================
!> @brief Compute sampling point coordinates and sample velocity
!! @param[out]   xh              x coordinate of the sampling point
!! @param[out]   yh              y coordinate of the sampling point
!! @param[out]   zh              z coordinate of the sampling point
!! @param[out]   vxh             x component of the sampled velocity
!! @param[out]   vyh             y component of the sampled velocity
!! @param[out]   vzh             z component of the sampled velocity
!! @param[in]    ifacex          the x index of the face node
!! @param[in]    ifacey          the y index of the face node
!! @param[in]    ifacez          the z index of the face node
!! @param[in]    iface           the index of the face in an element
!! @param[in]    ielem           the index of the element
!! @param[in]    samplingidx     the wall-normal index of the sampling point 
      subroutine sample_index_based(xh, yh, zh,
     $                              vxh, vyh, vzh,
     $                              ifacex, ifacey, ifacez,
     $                              iface, ielem, samplingidx)
      implicit none
     
      include 'SIZE'
      include 'GEOM'
      include 'SOLN'

      real xh, yh, zh
      real vxh, vyh, vzh
      integer ifacex, ifacey, ifacez
      integer iface, ielem, samplingidx

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
      else if (iface .eq. 2) then
        ! Face corresponds to y-z plane at x = 1
        xh = xm1(lx1 - samplingidx, ifacey, ifacez, ielem)
        yh = ym1(lx1 - samplingidx, ifacey, ifacez, ielem)
        zh = zm1(lx1 - samplingidx, ifacey, ifacez, ielem)

        vxh = vx(lx1 - samplingidx, ifacey, ifacez, ielem)
        vyh = vy(lx1 - samplingidx, ifacey, ifacez, ielem)
        vzh = vz(lx1 - samplingidx, ifacey, ifacez, ielem)
      else if (iface .eq. 3) then
        ! Face corresponds to x-z plane at y = 1
        xh = xm1(ifacex, ly1 - samplingidx, ifacez, ielem)
        yh = ym1(ifacex, ly1 - samplingidx, ifacez, ielem)
        zh = zm1(ifacex, ly1 - samplingidx, ifacez, ielem)

        vxh = vx(ifacex, ly1 - samplingidx, ifacez, ielem)
        vyh = vy(ifacex, ly1 - samplingidx, ifacez, ielem)
        vzh = vz(ifacex, ly1 - samplingidx, ifacez, ielem)
      else if (iface .eq. 4) then
        ! Face corresponds to y-z plane at x = -1
        xh = xm1(samplingidx + 1, ifacey, ifacez, ielem)
        yh = ym1(samplingidx + 1, ifacey, ifacez, ielem)
        zh = zm1(samplingidx + 1, ifacey, ifacez, ielem)

        vxh = vx(samplingidx + 1, ifacey, ifacez, ielem)
        vyh = vy(samplingidx + 1, ifacey, ifacez, ielem)
        vzh = vz(samplingidx + 1, ifacey, ifacez, ielem)
      else if (iface .eq. 5) then
        ! Face corresponds to x-y plane at z = -1
        xh = xm1(ifacex, ifacey, samplingidx + 1, ielem)
        yh = ym1(ifacex, ifacey, samplingidx + 1, ielem)
        zh = zm1(ifacex, ifacey, samplingidx + 1, ielem)

        vxh = vx(ifacex, ifacey, samplingidx + 1, ielem)
        vyh = vy(ifacex, ifacey, samplingidx + 1, ielem)
        vzh = vz(ifacex, ifacey, samplingidx + 1, ielem)
      else if (iface .eq. 6) then
        ! Face corresponds to x-y plane at z = 1
        xh = xm1(ifacex, ifacey, lz1 - samplingidx, ielem)
        yh = ym1(ifacex, ifacey, lz1 - samplingidx, ielem)
        zh = zm1(ifacex, ifacey, lz1 - samplingidx, ielem)

        vxh = vx(ifacex, ifacey, lz1 - samplingidx, ielem)
        vyh = vy(ifacex, ifacey, lz1 - samplingidx, ielem)
        vzh = vz(ifacex, ifacey, lz1 - samplingidx, ielem)
      end if
      end subroutine