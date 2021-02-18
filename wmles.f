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
!      call mntr_tmr_is_name_reg(lpmid,'FRM_TOT')
      ! total time
!      call mntr_tmr_reg(wmles_tmr_tot_id,lpmid,wmles_id,
!     $     'STAT_TOT','Statistics total time',.false.)
!      lpmid = wmles_tmr_tot_id
!      ! initialisation
!      call mntr_tmr_reg(wmles_tmr_ini_id,lpmid,wmles_id,
!     $     'STAT_INI','Statistics initialisation time',.true.)
!      ! averagign
!      call mntr_tmr_reg(wmles_tmr_avg_id,lpmid,wmles_id,
!     $     'STAT_AVG','Statistics averaging time',.true.)
 
!      if (wmles_rdim.eq.1) then
!      ! communication
!         call mntr_tmr_reg(wmles_tmr_cmm_id,lpmid,wmles_id,
!     $        'STAT_CMM','Statistics communication time',.true.)
!      endif
!      ! IO
!      call mntr_tmr_reg(wmles_tmr_io_id,lpmid,wmles_id,
!     $     'STAT_IO','Statistics IO time',.true.)
 
      ! register and set active section
      call rprm_sec_reg(wmles_sec_id, wmles_id,'_'//adjustl(wmles_name),
     $     'Runtime parameter section for the WMLES module')
      call rprm_sec_set_act(.true.,wmles_sec_id)
 
      ! register parameters
      call rprm_rp_reg(wmles_logkappa_id, wmles_sec_id, 'LOGKAPPA',
     $     'von Karman coefficient', rpar_real, 0, 0.41, .false.,' ')

      call rprm_rp_reg(wmles_logb_id, wmles_sec_id, 'LOGB',
     $     'the intercept of the log law', rpar_real, 0, 0.41, .false.,
     $     ' ')

      call rprm_rp_reg(wmles_guess_id, wmles_sec_id, 'GUESS',
     $     'An initial guess for tau_w', rpar_real, 0, 5e-2, .false.,
     $     ' ')
 
!       call rprm_rp_reg(wmles_IOstep_id,wmles_sec_id,'IOSTEP',
!      $     'Frequency of filed saving',rpar_int,100,0.0,.false.,' ')
       
       ! set initialisation flag
       wmles_ifinit=.false.
 
       ! timing
!       ltim = dnekclock() - ltim
!       call mntr_tmr_add(wmles_tmr_tot_id,1,ltim)
 
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
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,wmles_guess_id,rpar_real)
      tau = rtmp
      
    
      ! everything is initialised
      wmles_ifinit=.true.
    
      ! timing
   !   ltim = dnekclock() - ltim
   !   call mntr_tmr_add(wmles_tmr_ini_id,1,ltim)
      
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