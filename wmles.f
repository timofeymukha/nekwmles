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

      call rprm_rp_reg(wmles_guess_id, wmles_sec_id, 'GUESS',
     $     'An initial guess for tau_w', rpar_real, 0, 0.002, .false.,
     $     ' ')

      call rprm_rp_reg(wmles_bid_id, wmles_sec_id, 'WALLBID',
     $     'Boundary ID of the wall faces', rpar_int, 1, 0.0, .false.,
     $     ' ')
 
      call rprm_rp_reg(wmles_samplingidx_id,wmles_sec_id,'SAMPLINGIDX',
     $   'Wall-normal sampling point index',rpar_int,2,0.0,.false.,' ')

      call rprm_rp_reg(wmles_ifviscosity_id,wmles_sec_id,'IFVISCOSITY',
     $   'Whether tau_w is imposed via wall viscosity',
     $    rpar_log, 0, 0.0, .true.,' ')

      call rprm_rp_reg(wmles_navrg_id,wmles_sec_id,'NAVRG',
     $   'Number of time steps defining the input averaging time scale',
     $    rpar_real, 0, 1.0, .false.,' ')
       
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
!      ltim = dnekclock()

      vh = 0
      
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
      
      do ielem = 1, lelv
        do k = 1, lz1
          do j = 1, ly1
            do i = 1, lx1
              ! Assign the whole magnitude to the x component.
              ! Does not matter since the magnitude of tau will be
              ! used for the guess.
              tau(1, i, j, k, ielem) = rtmp
              tau(2, i, j, k, ielem) = 0
              tau(3, i, j, k, ielem) = 0
            end do
          end do
        end do
      end do
      
      ! get and assign the wall-normal index of the sampling point
      call rprm_rp_get(itmp, rtmp, ltmp, ctmp, wmles_samplingidx_id,
     $                 rpar_int)
      samplingidx = itmp

      ! get and assign the id of the wall boundary
      call rprm_rp_get(itmp, rtmp, ltmp, ctmp, wmles_bid_id,
     $                 rpar_int)
      wallbid = itmp
      
      ! Check that the sampling index is in [1; lx1-1]
      if (itmp .lt. 1) then
        if (nid .eq. 0) then 
          write(*, *)
     $      "[WMLES] FATAL ERROR: The SAMPLINGIDX must be >= 1"
        end if
        call exitt
      elseif (itmp .ge. lx1) then
        if (nid .eq. 0) then 
          write(*, *)
     $     "[WMLES] FATAL ERROR: The SAMPLINGIDX must be < lx1"
        end if
        call exitt
      end if

      ! get and assign the wall viscosity flag
      call rprm_rp_get(itmp, rtmp, ltmp, ctmp, wmles_ifviscosity_id,
     $                 rpar_log)
      
      ifviscosity = ltmp

      ! get and assign the number of time-steps for input averaging
      call rprm_rp_get(itmp, rtmp, ltmp, ctmp, wmles_navrg_id,
     $                 rpar_real)
      
      wmles_navrg = rtmp
    
      ! everything is initialised
      wmles_ifinit=.true.
    
      ! timing
!      ltim = dnekclock() - ltim
!      call mntr_tmr_add(wmles_tmr_ini_id,1,ltim)
      
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
