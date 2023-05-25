!> @file init.f
!! @ingroup wmles
!! @brief Wall-model LES (WMLES) module
!! @details Register and initialize parameters which control physics of
!!          the wall model and the sub-grid scale model
!! @{
!=======================================================================
!> Register WMLES module
!! @note This routine should be called in frame_usr_register
      subroutine wmles_register()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'WMLES'

      ! local variables
      integer lpmid
!-----------------------------------------------------------------------
      ! check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid, wmles_sec_name)
      if (lpmid.gt.0) then
         call mntr_warn(
     $      lpmid,
     $     'module ['//trim(wmles_sec_name)//'] already registered')
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
      call mntr_mod_reg(wmles_id, lpmid, wmles_sec_name,
     $      'Wall model LES parameters')

      ! register and set active section
      call rprm_sec_reg(
     $    wmles_sec_id,wmles_id,'_'//adjustl(wmles_sec_name),
     $    'Runtime parameter section for WMLES module')
      call rprm_sec_set_act(.true.,wmles_sec_id)

      ! register parameters
      ! subroutine rprm_rp_reg(rpid, mid, pname, pdscr, ptype, ipval, rpval, lpval, cpval)
      call rprm_rp_reg(
     $    wmles_par_id(1), wmles_sec_id,
     $    'BCTEMPFILT', 'Temporal filtering for BC', rpar_log,
     $    1, 0.0, .false., ' '
     $)
      call rprm_rp_reg(
     $    wmles_par_id(2), wmles_sec_id,
     $    'BCZINDEX', 'Z index for BC', rpar_int,
     $    1, 0.0, .false., ' '
     $)
      call rprm_rp_reg(
     $    wmles_par_id(3), wmles_sec_id,
     $    'SGSDELTAMAX', 'Delta max for SGS length scale', rpar_log,
     $    1, 0.0, .false., ' '
     $)
      call rprm_rp_reg(
     $    wmles_par_id(4), wmles_sec_id,
     $    'SGSNPOW', 'Power of SGS wall damping function', rpar_real,
     $    1, 0.5, .false., ' '
     $)
      call rprm_rp_reg(
     $    wmles_par_id(5), wmles_sec_id,
     $    'BCZ0', 'Aerodynamic roughness parameter for BC', rpar_real,
     $    1, 0.1, .false., ' '
     $)
      call rprm_rp_reg(
     $    wmles_par_id(6), wmles_sec_id,
     $    'SGSC0', 'Asymptotic constant in a SGS model', rpar_real,
     $    1, 0.19, .false., ' '
     $)
      call rprm_rp_reg(
     $    wmles_par_id(7), wmles_sec_id,
     $    'SGSBC', 'Set SGS boundary condition', rpar_log,
     $    1, 0.0, .false., ' '
     $)

      ! set initialisation flag
      wmles_ifinit=.false.

      end subroutine
!=======================================================================
!> Initiliase WMLES Module
!! @note This routine should be called in frame_usr_init
      subroutine wmles_init()
      implicit none

      include 'SIZE'
      include 'FRAMELP'
      include 'WMLES'

      ! local variables for reading parameters with appropriate types
      integer itmp
      real rtmp
      logical ltmp
      character*20 ctmp
!-----------------------------------------------------------------------
      ! check if the module was already initialised
      if (wmles_ifinit) then
         call mntr_warn(wmles_id,
     $        'module ['//trim(wmles_sec_name)//'] already initiaised.')
         return
      endif

      ! get runtime parameters
      call rprm_rp_get(itmp, rtmp, ltmp, ctmp,
     $    wmles_par_id(1), rpar_log)
      wmles_bc_temp_filt = ltmp

      call rprm_rp_get(itmp, rtmp, ltmp, ctmp,
     $    wmles_par_id(2), rpar_int)
      wmles_bc_z_index = itmp

      call rprm_rp_get(itmp, rtmp, ltmp, ctmp,
     $    wmles_par_id(3), rpar_log)
      wmles_sgs_delta_max = ltmp

      call rprm_rp_get(itmp, rtmp, ltmp, ctmp,
     $    wmles_par_id(4), rpar_real)
      wmles_sgs_npow = rtmp

      call rprm_rp_get(itmp, rtmp, ltmp, ctmp,
     $    wmles_par_id(5), rpar_real)
      wmles_bc_z0 = rtmp

      call rprm_rp_get(itmp, rtmp, ltmp, ctmp,
     $    wmles_par_id(6), rpar_real)
      wmles_sgs_c0 = rtmp

      call rprm_rp_get(itmp, rtmp, ltmp, ctmp,
     $    wmles_par_id(7), rpar_log)
      wmles_sgs_bc = ltmp

      ! everything is initialised
      wmles_ifinit=.true.
      end subroutine
!=======================================================================
!> Check if module was initialised
!! @return wmles_is_initialised
      logical function wmles_is_initialised()
      implicit none

      include 'SIZE'
      include 'WMLES'
      wmles_is_initialised = wmles_ifinit

      return
      end function
!=======================================================================
!! @}
