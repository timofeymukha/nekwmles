!> @file wmles.f
!! @ingroup wmles
!! @brief Wall-modelling for LES
!! @details 
!! @author Timofey Mukha
!! @date 2020
!=======================================================================
!> @brief Register the wmles module
!! @ingroup wmles
!! @note This routine should be called in frame_usr_register
      subroutine wmles_register()
      implicit none
 
      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'WMLES'


      ! local variables
      integer lpmid, il
      real ltim
      character*2 str
      ! functions
      real dnekclock
!----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()
 
      ! check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid, wmles_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(wmles_name)//'] already registered')
         return
      endif
      
      ! find parent module
      call mntr_mod_is_name_reg(lpmid,'FRAME')
      if (lpmid.le.0) then
         lpmid = 1
         call mntr_abort(lpmid,
     $        'parent module ['//'FRAME'//'] not registered')
      endif
      
      ! register module
      call mntr_mod_reg(wmles_id,lpmid,wmles_name,
     $      'Wall-Modelled LES')
 
      ! register timers
      call mntr_tmr_is_name_reg(lpmid,'FRM_TOT')

      ! total time
      call mntr_tmr_reg(wmles_tmr_tot_id,lpmid,wmles_id,
     $     'WMLES_TOT','Wall modelling total time',.false.)
      lpmid = wmles_tmr_tot_id

      call mntr_tmr_reg(wmles_tmr_sampling_id,lpmid,wmles_id,
     $     'WMLES_SMP','Sampling the data to the wall model',.true.)
 
 
      ! register and set active section
      call rprm_sec_reg(wmles_sec_id, wmles_id,'_'//adjustl(wmles_name),
     $     'Runtime parameter section for the WMLES module')
      call rprm_sec_set_act(.true.,wmles_sec_id)
 
      ! register parameters
      call rprm_rp_reg(wmles_logkappa_id, wmles_sec_id, 'LOGKAPPA',
     $     'Von Karman coefficient', rpar_real, 0, 0.387, .false.,' ')

      call rprm_rp_reg(wmles_logb_id, wmles_sec_id, 'LOGB',
     $     'The intercept of the log law', rpar_real, 0, 4.21, .false.,
     $     ' ')
      call rprm_rp_reg(wmles_z0_id, wmles_sec_id, 'Z0',
     $     'The roughness height', rpar_real, 0, 0.1, .false.,
     $     ' ')
      call rprm_rp_reg(wmles_z1_id, wmles_sec_id, 'Z1',
     $     'The roughness height', rpar_real, 0, 0.1, .false.,
     $     ' ')

      call rprm_rp_reg(wmles_guess_id, wmles_sec_id, 'GUESS',
     $     'An initial guess for tau_w', rpar_real, 0, 0.002, .false.,
     $     ' ')

      call rprm_rp_reg(wmles_bid_id, wmles_sec_id, 'WALLBID',
     $     'Boundary ID of the wall faces', rpar_int, 1, 0.0, .false.,
     $     ' ')
 
      call rprm_rp_reg(wmles_ifviscosity_id,wmles_sec_id,'IFVISCOSITY',
     $   'Whether tau_w is imposed via wall viscosity',
     $    rpar_log, 0, 0.0, .true.,' ')

      call rprm_rp_reg(wmles_navrg_id,wmles_sec_id,'NAVRG',
     $   'Number of time steps defining the input averaging time scale',
     $    rpar_real, 0, 1.0, .false.,' ')
      call rprm_rp_reg(wmles_surface_temp_id,wmles_sec_id,'TSURFACE',
     $   'Surface temperature',
     $    rpar_real, 0, 273.0, .false.,' ')
      call rprm_rp_reg(wmles_h_is_index_id,wmles_sec_id,'IFHISINDEX',
     $   'Whether h holds consecutive index or distance',
     $    rpar_log, 0, 0.0, .true.,' ')
      call rprm_rp_reg(wmles_theta0_id,wmles_sec_id,'THETA0',
     $   'Reference temperature',
     $    rpar_real, 0, 0.0, .false.,' ')

       ! set initialisation flag
       wmles_ifinit=.false.
 
       ! timing
       ltim = dnekclock() - ltim
       call mntr_tmr_add(wmles_tmr_tot_id,1,ltim)
 
       return
       end subroutine
!=======================================================================
!> @brief Initilise the WMLES module
!! @ingroup wmles
!! @note This routine should be called in frame_usr_init
      subroutine wmles_init()
      implicit none
      
      include 'SIZE'
      include 'FRAMELP'
      include 'TSTEP'
      include 'WMLES'
      
      ! local dummy variables
      integer itmp, il
      real rtmp, ltim
      logical ltmp
      character*20 ctmp
      
      integer i, j, k, ielem
       
      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! check if the module was already initialised
      if (wmles_ifinit) then
        call mntr_warn(wmles_id,
     $      'module ['//trim(wmles_name)//'] already initiaised.')
        return
      endif
      
      ! timing
      ltim = dnekclock()

      ! get the guess for tau and assign
      call rprm_rp_get(itmp, rtmp, ltmp, ctmp, wmles_guess_id,rpar_real)
      
      ! Check that the guess for tau_w magnitude is positive
      if (rtmp .le. 0) then
        if (nid .eq. 0) then 
          write(*,*)
     $      "[WMLES] FATAL ERROR: The GUESS parameter must be positive"
  
        end if
        call exitt
      end if
      
      call cfill(wmles_tau(:, 1), rtmp, NMAX_BOUNDARY_POINTS)
      call cfill(wmles_tau(:, 2), 0.0, NMAX_BOUNDARY_POINTS)
      call cfill(wmles_tau(:, 3), 0.0, NMAX_BOUNDARY_POINTS)
      
      ! get and assign the id of the wall boundary
      call rprm_rp_get(itmp, rtmp, ltmp, ctmp, wmles_bid_id,
     $                 rpar_int)
      wmles_wallbid = itmp
      
      ! get and assign the wall viscosity flag
      call rprm_rp_get(itmp, rtmp, ltmp, ctmp, wmles_ifviscosity_id,
     $                 rpar_log)
      
      wmles_ifviscosity = ltmp

      ! get and assign the h is index flag
      call rprm_rp_get(itmp, rtmp, ltmp, ctmp, wmles_h_is_index_id,
     $                 rpar_log)
      
      wmles_ifhisindex = ltmp

      ! get and assign the number of time-steps for input averaging
      call rprm_rp_get(itmp, rtmp, ltmp, ctmp, wmles_navrg_id,
     $                 rpar_real)
      
      wmles_navrg = rtmp

      ! get and assign surface temperature
      call rprm_rp_get(itmp, rtmp, ltmp, ctmp, wmles_surface_temp_id,
     $                 rpar_real)
      
      wmles_surface_temp = rtmp

      ! get and assign reference temperature
      call rprm_rp_get(itmp, rtmp, ltmp, ctmp, wmles_theta0_id,
     $                 rpar_real)
      
      wmles_theta0 = rtmp

      ! routine that sets h, should be provided by the user
      if (nid .eq. 0) write(*,*) "[WMLES] Setting sampling height"
      call user_set_sampling_height
      
      ! if h is input as an index, compute h as distance
      if (wmles_ifhisindex) call wmles_set_h_from_indices

      if (nid .eq. 0) write(*,*)
     $  "[WMLES] Searching for sampling points"

      call find_sampling_points 
      ! Initialize the sampling
      call interp_setup(wmles_interpolation_handle, 0.0, 0, nelt)
      ! Set flag for finding interpolation points to true
      ! Used by the built-ininterpolation routine
      wmles_iffind = .true.

      ! everything is initialised
      wmles_ifinit = .true.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(wmles_tmr_tot_id, 1, ltim)
      
      return
      end subroutine
!=======================================================================
!> @brief Check if module was initialised
!! @ingroup wmles
!! @return wmles_is_initialised
      logical function wmles_is_initialised()
      implicit none

      include 'SIZE'
      include 'WMLES'
!-----------------------------------------------------------------------
      wmles_is_initialised = wmles_ifinit

      return
      end function
!=======================================================================
!> @brief Outpost h and sampled velocity
!! @ingroup wmles
      subroutine wmles_outpost()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'WMLES'
!-----------------------------------------------------------------------
      real h(lx1, ly1, lz1, nelt)
      real spx(lx1, ly1, lz1, nelt)
      real spy(lx1, ly1, lz1, nelt)
      real spz(lx1, ly1, lz1, nelt)
      real svx(lx1, ly1, lz1, nelt)
      real svy(lx1, ly1, lz1, nelt)
      real svz(lx1, ly1, lz1, nelt)
      real svt(lx1, ly1, lz1, nelt)

      integer ntot
      ! Loop ranges for traversing face nodes
      integer frangex1, frangex2, frangey1, frangey2, frangez1, frangez2

      integer ifacex, ifacey, ifacez
      integer iface, ielem

      ! Linear index
      integer i_linear


      ntot = lx1*ly1*lz1*nelt

      ! Zero everywhere but the walls
      call rzero(h, ntot)
      call rzero(spx, ntot)
      call rzero(spy, ntot)
      call rzero(spz, ntot)
      call rzero(svx, ntot)
      call rzero(svy, ntot)
      call rzero(svz, ntot)
      call rzero(svt, ntot)

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
c
                  h(ifacex, ifacey, ifacez, ielem) = 
     $               wmles_sampling_h(i_linear)
                  svx(ifacex, ifacey, ifacez, ielem) = 
     $               wmles_solh(i_linear, 1)
                  svy(ifacex, ifacey, ifacez, ielem) = 
     $               wmles_solh(i_linear, 2)
                  svz(ifacex, ifacey, ifacez, ielem) = 
     $               wmles_solh(i_linear, 3)
                  spx(ifacex, ifacey, ifacez, ielem) = 
     $               wmles_sampling_points(i_linear, 1)
                  spy(ifacex, ifacey, ifacez, ielem) = 
     $               wmles_sampling_points(i_linear, 2)
                  spz(ifacex, ifacey, ifacez, ielem) = 
     $               wmles_sampling_points(i_linear, 3)
                  if (ifheat) then
                    svt(ifacex, ifacey, ifacez, ielem) = 
     $                wmles_solh(i_linear, 4)
                  endif

                end do
              end do
            end do

          endif

        enddo
      enddo
      call outpost(h, h, h, xm2, h, 'wmh') 
      call outpost(svx, svy, svz, xm2, svt, 'wmv') 
      call outpost(spx, spy, spz, xm2, xm1, 'wmp') 
      
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

      integer iglmax, cpu_max_nbndry

      ! Start the timer 
      ltim = dnekclock()

      n_sampling = 0

      call izero(wmles_inv_indices, lx1*ly1*lz1*nelt)

      do ielem=1, lelv
        do iface=1, 6 

          if (boundaryID(iface, ielem) .eq. wmles_wallbid) then

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

                  if (n_sampling .gt. NMAX_BOUNDARY_POINTS) then
                      write(*,*) "ERROR: Increase N_BOUNDARY_PONTS!"
                  endif

                  inorm = inorm + 1
                  xgll = xm1(ifacex, ifacey, ifacez, ielem)
                  ygll = ym1(ifacex, ifacey, ifacez, ielem)
                  zgll = zm1(ifacex, ifacey, ifacez, ielem)

                  ! inward face normal
                  normalx =  unx(inorm, 1, iface, ielem)
                  normaly = -uny(inorm, 1, iface, ielem)
                  normalz =  unz(inorm, 1, iface, ielem)

                  ! store the normals
                  wmles_normals(n_sampling, 1) = normalx
                  wmles_normals(n_sampling, 2) = normaly
                  wmles_normals(n_sampling, 3) = normalz

                  ! store the area
                  wmles_areas(n_sampling) = area(inorm, 1, iface, ielem) 

                  ! store the wall indices
                  wmles_indices(n_sampling, 1) = ifacex
                  wmles_indices(n_sampling, 2) = ifacey
                  wmles_indices(n_sampling, 3) = ifacez
                  wmles_indices(n_sampling, 4) = ielem
                  
                  ! store the inverse index
                  wmles_inv_indices(ifacex, ifacey, ifacez, ielem) = 
     $              n_sampling


                  wmles_sampling_points(n_sampling, 1) = 
     $               xgll + normalx*wmles_sampling_h(n_sampling)
                  wmles_sampling_points(n_sampling, 2) = 
     $               ygll + normaly*wmles_sampling_h(n_sampling)
                  wmles_sampling_points(n_sampling, 3) = 
     $               zgll + normalz*wmles_sampling_h(n_sampling)

                end do
              end do
            end do

          endif

        enddo
      enddo


      wmles_nbpoints = n_sampling
      cpu_max_nbndry = iglmax(wmles_nbpoints, 1)

      if (nid .eq. 0) then
        write(*,*) "[WMLES] NMAX_BOUNDARY_POINTS must be >= ", 
     $    cpu_max_nbndry
        write(*,*) "[WMLES] NMAX_BOUNDARY_POINTS is", 
     $    NMAX_BOUNDARY_POINTS 
      endif

c      write(*,*) wmles_sampling_points(:wmles_nbpoints, 2)
      ! Stop the timer and add to total
      ltim = dnekclock() - ltim
      call mntr_tmr_add(wmles_tmr_sampling_id, 1, ltim)
      end subroutine

