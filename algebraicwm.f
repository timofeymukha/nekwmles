
      subroutine algebraic_wm()
      implicit none

      include 'SIZE'
      include 'TSTEP'           ! ISTEP, lastep, time
      include 'INPUT'
      include 'SOLN'
      include 'GEOM'
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
      real xh, yh, zh, vxh, vyh, vzh, magvh

      ! The distance to the sampling point
      real h

      real nu, utau, totalutau, newton

      real guess
      
      ! The array holding the wall shear stress magnitude
      real tau(lx1, ly1, lz1, lelv)

      common /wmles/ tau

      ! Grab various simulation parameters
      nu = param(2)

      ! Fill with a default value for debugging purposes
      !tau = -1

      nwallnodes = 0
      totalutau = 0

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

                  normalx = unx(ifacex, ifacey, iface, ielem)
                  normaly = -uny(ifacex, ifacey, iface, ielem)
                  normalz = unz(ifacex, ifacey, iface, ielem)
                  
c                  write(*,*) ielem, unx(ifacex, ifacey, iface, ielem),
c     $                       uny(ifacex, ifacey, iface, ielem),
c     $                       unz(ifacex, ifacey, iface, ielem),
c     $                       unr(ifacex, iface, ielem),
c     $                       uns(ifacex, iface, ielem),
c     $                       unt(ifacex, iface, ielem)


c                  write(*,*) ielem, iface, xgll, ygll, zgll, ynormal

                  ! Check the direction of the normal to figure out
                  ! how to change the y index to select the sampling
                  ! point.
                  !if (ynormal .gt. 0) then
                    !We are at the bottom wall
                  !  yhind = ly1-1
                  !else
                  !  yhind = 2
                  !endif
                  

                  ! We figure out the index the sampling point
                  ! location based on the plane the face is in
                  if (iface .eq. 1) then
                    ! Face corresponds to x-z plane at y = -1 
                    xh = xm1(ifacex, ly1-1, ifacez, ielem)
                    yh = ym1(ifacex, ly1-1, ifacez, ielem)
                    zh = zm1(ifacex, ly1-1, ifacez, ielem)
                    
                    vxh = vx(ifacex, ly1-1, ifacez, ielem)
                    vyh = vy(ifacex, ly1-1, ifacez, ielem)
                    vzh = vz(ifacex, ly1-1, ifacez, ielem)
                  else if (iface .eq. 2) then
                    ! Face corresponds to y-z plane at x = 1
                    xh = xm1(2, ifacey, ifacez, ielem)
                    yh = ym1(2, ifacey, ifacez, ielem)
                    zh = zm1(2, ifacey, ifacez, ielem)

                    vxh = vx(2, ifacey, ifacez, ielem)
                    vyh = vy(2, ifacey, ifacez, ielem)
                    vzh = vz(2, ifacey, ifacez, ielem)
                  else if (iface .eq. 3) then
                    ! Face corresponds to x-z plane at y = 1
                    xh = xm1(ifacex, 2, ifacez, ielem)
                    yh = ym1(ifacex, 2, ifacez, ielem)
                    zh = zm1(ifacex, 2, ifacez, ielem)

                    vxh = vx(ifacex, 2, ifacez, ielem)
                    vyh = vy(ifacex, 2, ifacez, ielem)
                    vzh = vz(ifacex, 2, ifacez, ielem)
                  else if (iface .eq. 4) then
                    ! Face corresponds to y-z plane at x = -1
                    xh = xm1(lx1-1, ifacey, ifacez, ielem)
                    yh = ym1(lx1-1, ifacey, ifacez, ielem)
                    zh = zm1(lx1-1, ifacey, ifacez, ielem)

                    vxh = vx(lx1-1, ifacey, ifacez, ielem)
                    vyh = vy(lx1-1, ifacey, ifacez, ielem)
                    vzh = vz(lx1-1, ifacey, ifacez, ielem)
                  else if (iface .eq. 5) then
                    ! Face corresponds to x-y plane at z = -1
                    xh = xm1(ifacex, ifacey, lz1-1, ielem)
                    yh = ym1(ifacex, ifacey, lz1-1, ielem)
                    zh = zm1(ifacex, ifacey, lz1-1, ielem)

                    vxh = vx(ifacex, ifacey, lz1-1, ielem)
                    vyh = vy(ifacex, ifacey, lz1-1, ielem)
                    vzh = vz(ifacex, ifacey, lz1-1, ielem)
                  else if (iface .eq. 6) then
                    ! Face corresponds to x-y plane at z = 1
                    xh = xm1(ifacex, ifacey, 2, ielem)
                    yh = ym1(ifacex, ifacey, 2, ielem)
                    zh = zm1(ifacex, ifacey, 2, ielem)

                    vxh = vx(ifacex, ifacey, 2, ielem)
                    vyh = vy(ifacex, ifacey, 2, ielem)
                    vzh = vz(ifacex, ifacey, 2, ielem)
                  end if

                  
c                  write(*,'(I2, I2, 12(F8.5, XX))') ielem, iface,
c     $                       rxm1(ifacex, yhind, ifacez, ielem),
c     $                       sxm1(ifacex, yhind, ifacez, ielem),
c     $                       txm1(ifacex, yhind, ifacez, ielem),
c     $                       rym1(ifacex, yhind, ifacez, ielem),
c     $                       sym1(ifacex, yhind, ifacez, ielem),
c     $                       tym1(ifacex, yhind, ifacez, ielem),
c     $                       rzm1(ifacex, yhind, ifacez, ielem),
c     $                       szm1(ifacex, yhind, ifacez, ielem),
c     $                       tzm1(ifacex, yhind, ifacez, ielem),
c     $                       xm1(ifacex, yhind, ifacez, ielem),
c     $                       ym1(ifacex, yhind, ifacez, ielem),
c     $                       zm1(ifacex, yhind, ifacez, ielem)

                  ! Vector between face node and sampling point
                  xh = xh - xgll
                  yh = yh - ygll
                  zh = zh - zgll
                  
                  ! project onto the local face normal
                  h = abs(xh*normalx + yh*normaly + zh*normalz)
                  
                  magvh = sqrt(vxh**2 + vzh**2)

                  ! take the stress at the previous step as a guess
                  ! inital value is the GUESS parameter
                  guess = sqrt(tau(ifacex, ifacey, ifacez, ielem))
                  utau = newton(magvh, h, nu, guess,
     $                          1e-3, 50, .false.)
                  totalutau = totalutau + utau
                  nwallnodes = nwallnodes + 1
                  tau(ifacex, ifacey, ifacez, ielem) = utau**2

c                 write(*,*) xgll, xh, ygll, yh, zgll, zh, h
c                 write(*,*) vxh, h, tauw
c                 write(*,*) ielem, iface, vxh, vyh, vzh, h

                end do
              end do
            end do

          endif

        enddo
      enddo
      
      write(*,*) "Average predicted Re_tau = ", totalutau/nwallnodes/nu 
      end subroutine