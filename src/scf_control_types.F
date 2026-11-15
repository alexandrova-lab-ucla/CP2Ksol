!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief parameters that control an scf iteration
!> \note
!>       not in cp_control_types, to separate operator related parameters from
!>       method related parameters (as suggested by Matthias)
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE scf_control_types

   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE input_constants,                 ONLY: &
        atomic_guess, diag_ot, direct_p_mix, general_roks, high_spin_roks, no_guess, no_mix, &
        ot_algo_taylor_or_diag, outer_scf_basis_center_opt, outer_scf_cdft_constraint, &
        outer_scf_ddapc_constraint, outer_scf_none, outer_scf_optimizer_bisect, &
        outer_scf_optimizer_broyden, outer_scf_optimizer_diis, outer_scf_optimizer_newton, &
        outer_scf_optimizer_newton_ls, outer_scf_optimizer_none, outer_scf_optimizer_sd, &
        outer_scf_optimizer_secant, outer_scf_s2_constraint, smear_energy_window, &
        smear_fermi_dirac, smear_gaussian, smear_list, smear_mp, smear_mv
   USE gce_dyn_types,                   ONLY: gce_dyn_max_chain,&
                                              gce_dyn_type
   USE input_cp2k_scf,                  ONLY: create_scf_section
   USE input_enumeration_types,         ONLY: enum_i2c,&
                                              enumeration_type
   USE input_keyword_types,             ONLY: keyword_get,&
                                              keyword_type
   USE input_section_types,             ONLY: section_get_keyword,&
                                              section_release,&
                                              section_type,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: evolt,&
                                              femtoseconds
   USE outer_scf_control_types,         ONLY: outer_scf_control_type,&
                                              outer_scf_read_parameters
   USE qs_cdft_opt_types,               ONLY: cdft_opt_type_release
   USE qs_ot_types,                     ONLY: ot_readwrite_input,&
                                              qs_ot_settings_init,&
                                              qs_ot_settings_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scf_control_types'
   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.

   ! Public data types

   PUBLIC :: scf_control_type, &
             smear_type, &
             gce_type

   ! Public subroutines

   PUBLIC :: scf_c_create, &
             scf_c_read_parameters, &
             scf_c_release, &
             scf_c_write_parameters

   ! GCE restart state: typed RESTART_* keywords of &GCE, written into the input tree after every ionic step
   PUBLIC :: gce_state_write, &
             gce_state_string, &
             gce_md_report

! **************************************************************************************************
!> \brief contains the parameters needed by a scf run
!> \param density_guess how to choose the initial density
!>        (CORE,RANDOM,RESTART,ATOMIC,FROZEN)
!> \param eps_eigval wanted error on the eigenvalues
!> \param eps_scf whanted error on the whole scf
!> \param level_shift amount of level shift
!> \param p_mix how to mix the new and old densities in non diss iterations
!> \param eps_lumos error on the lumos calculated at the end of the scf
!> \param max_iter_lumus maxumum number of iterations used to calculate
!>        the lumos at the end of the scf
!> \param max_scf max scf iterations
!> \param added_mos additional number of MOs that might be used in the SCF
!> \param step_size the optimizer step size
!> \param cdft_opt_control settings for optimizers that work only together with CDFT constraints
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   TYPE smear_type
      LOGICAL                               :: do_smear = .FALSE.
      LOGICAL                               :: common_mu = .FALSE.
      INTEGER                               :: method = -1
      REAL(KIND=dp)                         :: electronic_temperature = -1.0_dp, &
                                               fixed_mag_mom = -1.0_dp, &
                                               eps_fermi_dirac = -1.0_dp, &
                                               window_size = -1.0_dp, &
                                               smearing_width = -1.0_dp
      REAL(KIND=dp), DIMENSION(:), POINTER  :: list => NULL()
   END TYPE smear_type

   TYPE diagonalization_type
      INTEGER                               :: method = -1
      REAL(KIND=dp)                         :: eps_jacobi = -1.0_dp
      REAL(KIND=dp)                         :: jacobi_threshold = -1.0_dp
      INTEGER                               :: max_iter = -1, nkrylov = -1, nblock_krylov = -1
      ! Maximum Overlap Method
      LOGICAL                               :: mom = .FALSE., mom_didguess = .FALSE.
      INTEGER                               :: mom_proj_formula = -1
      ! indices of de-occupied and newly occupied alpha / beta molecular orbitals
      INTEGER, DIMENSION(:), POINTER        :: mom_deoccA => NULL(), mom_deoccB => NULL(), &
                                               mom_occA => NULL(), mom_occB => NULL()
      ! determines on SCF which iteration MOM will be switched on;
      ! since MOs from the previous iteration should be available, it might be at least
      !  1 when wave-function has been read from restart file, or
      !  2 when the atomic guess method has been used
      INTEGER                               :: mom_start = -1
      INTEGER                               :: mom_type = -1
      REAL(KIND=dp)                         :: eps_iter = -1.0_dp
      REAL(KIND=dp)                         :: eps_adapt = -1.0_dp
      TYPE(qs_ot_settings_type)             :: ot_settings = qs_ot_settings_type()
   END TYPE diagonalization_type

   TYPE scf_control_type
      TYPE(outer_scf_control_type)          :: outer_scf = outer_scf_control_type()
      TYPE(smear_type), POINTER             :: smear => NULL()
      TYPE(diagonalization_type)            :: diagonalization = diagonalization_type()
      TYPE(gce_type), POINTER               :: gce => NULL()
      INTEGER                               :: density_guess = -1, mixing_method = -1
      REAL(KIND=dp)                         :: eps_eigval = -1.0_dp, eps_scf = -1.0_dp, eps_scf_hist = -1.0_dp, &
                                               level_shift = -1.0_dp, &
                                               eps_lumos = -1.0_dp, eps_diis = -1.0_dp
      INTEGER                               :: max_iter_lumos = -1, max_diis = -1, nmixing = -1
      INTEGER                               :: max_scf = -1, max_scf_hist = -1, &
                                               maxl = -1, nkind = -1
      LOGICAL                               :: do_diag_sub = .FALSE., &
                                               use_cholesky = .FALSE., use_ot = .FALSE., &
                                               use_diag = .FALSE., do_outer_scf_reortho = .FALSE., &
                                               ignore_convergence_failure = .FALSE.
      LOGICAL                               :: force_scf_calculation = .FALSE.
      LOGICAL                               :: non_selfconsistent = .FALSE.
      INTEGER, DIMENSION(2)                 :: added_mos = -1
      INTEGER                               :: roks_scheme = -1
      REAL(KIND=dp)                         :: roks_f = -1.0_dp
      REAL(KIND=dp), DIMENSION(0:2, 0:2, 1:2) :: roks_parameter = -1.0_dp
   END TYPE scf_control_type

   TYPE gce_type
      LOGICAL                               :: do_gce = .FALSE.
      REAL(KIND=dp)                         :: target_workfunction = 0.16_dp, &
                                               ref_esp = 0.0_dp, &
                                               mixing_coef = 0.3_dp, &
                                               prev_workfunction = -1001.0_dp, &
                                               nelec = 0.0_dp, &
                                               charge = 0.0_dp
      ! constant potential: capacitance-driven electron-number update
      INTEGER                               :: scheme = 1            ! 1 mu mixing, 2 capacitance
      LOGICAL                               :: update_every_step = .TRUE., &
                                               check_convergence = .FALSE.
      REAL(KIND=dp)                         :: capacitance = 0.0_dp, &   ! e/Hartree (input e/V)
                                               cap_min = 0.0_dp, &
                                               cap_max = 0.0_dp, &
                                               ne_step_max = 0.5_dp, &
                                               fermi_tol = 0.0_dp, &     ! Hartree
                                               fermi_update_tol = 0.0_dp, &
                                               n_target = -1.0_dp, &     ! current fractional electron number (<0: unset)
                                               n_input = -1.0_dp, &      ! electron number of the input charge state
                                               n_excess_init = 0.0_dp, & ! fixed fractional excess (EXCESS_ELECTRONS)
                                               first_update_size = 0.0_dp, &  ! FIRST_UPDATE_SIZE (0 = off)
                                               target_mu = -1000.0_dp, &  ! TARGET_CHEMICAL_POTENTIAL [Hartree]
                                               n_last = 0.0_dp, &
                                               wf_last = 0.0_dp, &
                                               target_potential = -1000.0_dp, &
                                               phi_she = 0.0_dp, &
                                               scf_delta = 1.0_dp, &        ! current SCF density change (set in qs_scf)
                                               update_delta = 1.0E-2_dp, &  ! N_e updated only while scf_delta < this
                                               early_update_delta = 0.0_dp, &   ! >0: early N_e update allowed
                                               keep_mixer_below = 0.0_dp, &     ! |dN_e| below this keeps the mixer history
                                               n_prev_update = 0.0_dp, &     ! N_e before the last update
                                               early_wf_stable = 0.0_dp, &   ! >0: needs |dWF| per step < this
                                               final_tol_frac = 1.0_dp, &    ! converged-SCF update test uses FERMI_TOL*this
                                               carry_gain = 0.0_dp, &        ! per ionic step N_e feed-forward gain
                                               cap_refresh_min = 0.0_dp, &   ! secant refresh needs |dWF| above this
                                               wf_prev_step = -1001.0_dp     ! WF of the previous inner step
      LOGICAL                               :: freeze_cap = .FALSE.   ! no secant refresh in this update
      LOGICAL                               :: have_wf_on_step = .FALSE.   ! wf_on_step was measured in this ionic step
      LOGICAL                               :: pin_mixer_g0 = .FALSE.   ! G=0 of the mixer input follows N_e
      LOGICAL                               :: cap_measured = .FALSE.   ! capacitance came from a secant, not CAPACITANCE_INIT
      TYPE(gce_dyn_type)                    :: dyn                      ! gce_dyn: electron-number dynamics (UPDATE_SCHEME DYNAMICS)
      LOGICAL                               :: early_secant = .FALSE.   ! early updates record / refresh the secant pair
      ! RESPONSE_GUESS: density response d rho/d N_e from two states of one force evaluation (rank-local G-space arrays)
      LOGICAL                               :: response_guess = .FALSE., &
                                               resp_have_rho = .FALSE., &   ! resp_rho, resp_n, resp_step valid
                                               resp_have_f = .FALSE.        ! resp_f valid
      INTEGER                               :: resp_step = -1               ! ion step of resp_rho
      REAL(KIND=dp)                         :: resp_n = 0.0_dp              ! N_e of resp_rho
      COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE :: resp_rho, resp_f       ! density at resp_n; d rho/d N_e
      LOGICAL                               :: constant_potential_requested = .FALSE.   ! steering requested
      LOGICAL                               :: secant_ready = .FALSE., &   ! history (n_last, wf_last) valid
                                               solvent_only_for_charge = .FALSE.   ! SOLVENT_ONLY_FOR_CHARGE
      INTEGER                               :: n_outer = 0, &               ! outer (cycle-mode) updates done
                                               max_outer = 20, &
                                               min_steps = 12, &            ! steps after an N_e update
                                               update_interval = 1, &              ! N_e updated every n-th ionic step
                                               convergence_mode = 3, &              ! CONVERGENCE_MODE 1/2/3
                                               max_outer_now = 20, &        ! outer-update cap of this ionic step
                                               ion_step = 0, &              ! completed SCF runs (ionic steps)
                                               steps_since_update = 1000000, &
                                               iter_base = 0, &             ! SCF steps of the current loop at the last N_e update
                                               secant_step = -1, &          ! ion step of the (n_last, wf_last) pair (-1: none)
                                               early_update_min_steps = 3, &   ! steps at new N_e first
                                               early_update_max = 2, &         ! early updates allowed per SCF
                                               n_early = 0, &                  ! early updates done
                                               ncall_restore = 0               ! Broyden call counter to restore
      ! restart state, per-step Omega, unreachable target
      REAL(KIND=dp)                         :: n_step = -1.0_dp, &          ! N_e of the last ionic step (before CARRY_GAIN)
                                               e_step = 0.0_dp, &           ! energy%total of that step [Ha]
                                               wf_step = -1001.0_dp, &      ! converged W_f of that step [Ha] (the last U)
                                               wf_on_step = 0.0_dp, &       ! solvent-on W_f before the switch-off [Ha]
                                               wf_lo = 0.0_dp, &            ! lowest / highest W_f - target [Ha] before an
                                               wf_hi = 0.0_dp, &            ! update in this ionic step, and their N_e
                                               n_at_lo = 0.0_dp, &
                                               n_at_hi = 0.0_dp
      LOGICAL                               :: enforce_target = .FALSE., &  ! ENFORCE_TARGET
                                               state_pending = .FALSE., &   ! restart state loaded, N_input not yet checked
                                               state_report = .FALSE., &    ! restart state loaded, not yet printed
                                               replay_pending = .FALSE.     ! first force evaluation repeats the last one
   END TYPE gce_type

CONTAINS

! **************************************************************************************************
!> \brief allocates and initializes an scf control object with the default values
!> \param scf_control the object to initialize
!> \par History
!>      09.2002 created [fawzi]
!>      - Default ROKS parameters added (05.04.06,MK)
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE scf_c_create(scf_control)

      TYPE(scf_control_type), INTENT(INOUT)              :: scf_control

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'scf_c_create'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Load the default values

      IF (scf_control%non_selfconsistent) THEN
         scf_control%density_guess = no_guess
      ELSE
         scf_control%density_guess = atomic_guess
      END IF
      scf_control%eps_eigval = 1.0E-5_dp
      scf_control%eps_scf = 1.0E-5_dp
      scf_control%eps_scf_hist = 0.0_dp
      scf_control%eps_lumos = 1.0E-5_dp
      scf_control%max_iter_lumos = 2999
      scf_control%eps_diis = 0.1_dp
      scf_control%level_shift = 0.0_dp
      scf_control%max_diis = 4
      scf_control%max_scf = 50
      scf_control%nmixing = 2
      scf_control%use_cholesky = .TRUE.
      scf_control%use_diag = .TRUE.
      scf_control%do_diag_sub = .FALSE.
      scf_control%use_ot = .FALSE.
      scf_control%ignore_convergence_failure = .FALSE.
      scf_control%force_scf_calculation = .FALSE.
      scf_control%do_outer_scf_reortho = .TRUE.
      scf_control%max_diis = 4
      scf_control%eps_diis = 0.1_dp
      scf_control%added_mos(:) = 0
      scf_control%max_scf_hist = 0

      !Mixing
      IF (scf_control%non_selfconsistent) THEN
         scf_control%mixing_method = no_mix
      ELSE
         scf_control%mixing_method = direct_p_mix
      END IF

      ! Diagonalization
      scf_control%diagonalization%method = 0
      scf_control%diagonalization%eps_jacobi = 0.0_dp
      scf_control%diagonalization%jacobi_threshold = 1.0E-7_dp
      scf_control%diagonalization%max_iter = 0
      scf_control%diagonalization%eps_iter = 0.0_dp
      scf_control%diagonalization%eps_adapt = 0.0_dp
      scf_control%diagonalization%nkrylov = 0
      scf_control%diagonalization%nblock_krylov = 0
      CALL qs_ot_settings_init(scf_control%diagonalization%ot_settings)

      scf_control%diagonalization%mom = .FALSE.
      scf_control%diagonalization%mom_didguess = .FALSE.
      scf_control%diagonalization%mom_proj_formula = 0
      NULLIFY (scf_control%diagonalization%mom_deoccA)
      NULLIFY (scf_control%diagonalization%mom_deoccB)
      NULLIFY (scf_control%diagonalization%mom_occA)
      NULLIFY (scf_control%diagonalization%mom_occB)
      scf_control%diagonalization%mom_start = 0

      ! ROKS

      scf_control%roks_scheme = high_spin_roks
      scf_control%roks_f = 0.5_dp

      ! Initialize the diagonal blocks with the default ROKS parameters
      ! 0 = v)irtual, 1 = o)pen shell, 2 = c)losed shell

      scf_control%roks_parameter(0, 0, 1) = 1.5_dp ! avv
      scf_control%roks_parameter(0, 0, 2) = -0.5_dp ! bvv
      scf_control%roks_parameter(1, 1, 1) = 0.5_dp ! aoo
      scf_control%roks_parameter(1, 1, 2) = 0.5_dp ! boo
      scf_control%roks_parameter(2, 2, 1) = -0.5_dp ! acc
      scf_control%roks_parameter(2, 2, 2) = 1.5_dp ! bcc

      ! Initialize off-diagonal blocks (fixed)

      scf_control%roks_parameter(0, 1, 1) = 1.0_dp ! avo
      scf_control%roks_parameter(0, 1, 2) = 0.0_dp ! bvo
      scf_control%roks_parameter(0, 2, 1) = 0.5_dp ! avc
      scf_control%roks_parameter(0, 2, 2) = 0.5_dp ! bvc
      scf_control%roks_parameter(1, 2, 1) = 0.0_dp ! aoc
      scf_control%roks_parameter(1, 2, 2) = 1.0_dp ! boc

      ! Symmetry enforces

      scf_control%roks_parameter(1, 0, 1) = scf_control%roks_parameter(0, 1, 1) ! aov
      scf_control%roks_parameter(1, 0, 2) = scf_control%roks_parameter(0, 1, 2) ! bov
      scf_control%roks_parameter(2, 0, 1) = scf_control%roks_parameter(0, 2, 1) ! acv
      scf_control%roks_parameter(2, 0, 2) = scf_control%roks_parameter(0, 2, 2) ! bcv
      scf_control%roks_parameter(2, 1, 1) = scf_control%roks_parameter(1, 2, 1) ! aco
      scf_control%roks_parameter(2, 1, 2) = scf_control%roks_parameter(1, 2, 2) ! bco

      ! Outer SCF default settings

      scf_control%outer_scf%have_scf = .FALSE.
      scf_control%outer_scf%max_scf = 0
      scf_control%outer_scf%eps_scf = 0.0_dp
      scf_control%outer_scf%step_size = 0.0_dp
      scf_control%outer_scf%type = -1
      scf_control%outer_scf%optimizer = -1
      scf_control%outer_scf%diis_buffer_length = -1
      NULLIFY (scf_control%outer_scf%cdft_opt_control)

      ! Smearing of the MO occupations

      NULLIFY (scf_control%smear)

      CALL timestop(handle)

   END SUBROUTINE scf_c_create

! **************************************************************************************************
!> \brief releases the given scf_control (see cp2k/doc/ReferenceCounting.html)
!> \param scf_control the object to free
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      at the moment does nothing
! **************************************************************************************************
   SUBROUTINE scf_c_release(scf_control)

      TYPE(scf_control_type), INTENT(INOUT)              :: scf_control

      IF (ASSOCIATED(scf_control%smear%list)) THEN
         DEALLOCATE (scf_control%smear%list)
      END IF
      DEALLOCATE (scf_control%smear)

      IF (ASSOCIATED(scf_control%outer_scf%cdft_opt_control)) &
         CALL cdft_opt_type_release(scf_control%outer_scf%cdft_opt_control)

      IF (ASSOCIATED(scf_control%gce)) THEN
         DEALLOCATE (scf_control%gce)
      END IF

      ! Maximum overlap method orbital indices lists
      ! mom_deoccA, mom_deoccB, mom_occA, mom_occB
      ! points to memory allocated by input file parser,
      ! so they do not have to be deallocated

   END SUBROUTINE scf_c_release

! **************************************************************************************************
!> \brief reads the parameters of the scf section into the given scf_control
!> \param scf_control the object that wil contain the values read
!> \param inp_section ...
!> \par History
!>      05.2001 created [Matthias]
!>      09.2002 creaded separated scf_control type [fawzi]
!> \author Matthias Krack
! **************************************************************************************************
   SUBROUTINE scf_c_read_parameters(scf_control, inp_section)

      TYPE(scf_control_type), INTENT(INOUT)              :: scf_control
      TYPE(section_vals_type), POINTER                   :: inp_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'scf_c_read_parameters'

      INTEGER                                            :: cholesky_flag, handle, ialgo
      INTEGER, DIMENSION(:), POINTER                     :: added_mos
      LOGICAL                                            :: do_mixing, explicit
      REAL(KIND=dp), DIMENSION(:), POINTER               :: roks_parameter
      TYPE(section_vals_type), POINTER                   :: gce_section, mixing_section, &
                                                            outer_scf_section, scf_section, &
                                                            smear_section

      CALL timeset(routineN, handle)

      scf_section => section_vals_get_subs_vals(inp_section, "SCF")
      CALL section_vals_val_get(scf_section, "DIAGONALIZATION%_SECTION_PARAMETERS_", &
                                l_val=scf_control%use_diag)
      IF (scf_control%use_diag) THEN
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%DIAG_SUB_SCF%_SECTION_PARAMETERS_", &
                                   l_val=scf_control%do_diag_sub)
      END IF
      CALL section_vals_val_get(scf_section, "OT%_SECTION_PARAMETERS_", l_val=scf_control%use_ot)
      IF (scf_control%use_diag .AND. scf_control%use_ot) THEN
         ! don't allow both options to be true
         CPABORT("Don't activate OT and Diagonaliztion together")
      ELSEIF (.NOT. (scf_control%use_diag .OR. scf_control%use_ot)) THEN
         ! set default to diagonalization
         scf_control%use_diag = .TRUE.
      END IF
      CALL section_vals_val_get(scf_section, "OT%ALGORITHM", i_val=ialgo)
      scf_control%do_outer_scf_reortho = ialgo == ot_algo_taylor_or_diag
      CALL section_vals_val_get(scf_section, "SCF_GUESS", i_val=scf_control%density_guess)
      CALL section_vals_val_get(scf_section, "eps_eigval", r_val=scf_control%eps_eigval)
      CALL section_vals_val_get(scf_section, "cholesky", i_val=cholesky_flag)
      IF (scf_control%use_ot) THEN
         ! eps_diis default is 0 for OT
         scf_control%eps_diis = 0.0_dp
         CALL section_vals_val_get(scf_section, "EPS_DIIS", explicit=explicit)
         IF (explicit) THEN
            CALL section_vals_val_get(scf_section, "EPS_DIIS", r_val=scf_control%eps_diis)
         END IF
      ELSE
         CALL section_vals_val_get(scf_section, "EPS_DIIS", r_val=scf_control%eps_diis)
      END IF
      IF (cholesky_flag > 0) THEN
         scf_control%use_cholesky = .TRUE.
      END IF
      CALL section_vals_val_get(scf_section, "IGNORE_CONVERGENCE_FAILURE", l_val=scf_control%ignore_convergence_failure)
      CALL section_vals_val_get(scf_section, "FORCE_SCF_CALCULATION", l_val=scf_control%force_scf_calculation)
      CALL section_vals_val_get(scf_section, "eps_scf", r_val=scf_control%eps_scf)
      CALL section_vals_val_get(scf_section, "level_shift", r_val=scf_control%level_shift)
      CALL section_vals_val_get(scf_section, "max_diis", i_val=scf_control%max_diis)
      CALL section_vals_val_get(scf_section, "max_scf", i_val=scf_control%max_scf)

      ! Diagonaliztion section
      IF (scf_control%use_diag) THEN
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%ALGORITHM", &
                                   i_val=scf_control%diagonalization%method)
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%EPS_JACOBI", &
                                   r_val=scf_control%diagonalization%eps_jacobi)
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%JACOBI_THRESHOLD", &
                                   r_val=scf_control%diagonalization%jacobi_threshold)
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%MAX_ITER", &
                                   i_val=scf_control%diagonalization%max_iter)
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%EPS_ITER", &
                                   r_val=scf_control%diagonalization%eps_iter)
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%EPS_ADAPT", &
                                   r_val=scf_control%diagonalization%eps_adapt)
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%KRYLOV%NKRYLOV", &
                                   i_val=scf_control%diagonalization%nkrylov)
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%KRYLOV%NBLOCK", &
                                   i_val=scf_control%diagonalization%nblock_krylov)
         IF (scf_control%diagonalization%method == diag_ot) THEN
            ! read OT section
            CALL ot_diag_read_input(scf_control%diagonalization%ot_settings, scf_section)
         END IF
         ! read maximum overlap method's parameters
         CALL section_vals_val_get(scf_section, "MOM%_SECTION_PARAMETERS_", &
                                   l_val=scf_control%diagonalization%MOM)
         IF (scf_control%diagonalization%mom) THEN
            CALL section_vals_val_get(scf_section, "MOM%MOM_TYPE", &
                                      i_val=scf_control%diagonalization%mom_type)

            CALL section_vals_val_get(scf_section, "MOM%START_ITER", &
                                      i_val=scf_control%diagonalization%mom_start)

            CALL section_vals_val_get(scf_section, "MOM%DEOCC_ALPHA", &
                                      i_vals=scf_control%diagonalization%mom_deoccA)

            CALL section_vals_val_get(scf_section, "MOM%DEOCC_BETA", &
                                      i_vals=scf_control%diagonalization%mom_deoccB)

            CALL section_vals_val_get(scf_section, "MOM%OCC_ALPHA", &
                                      i_vals=scf_control%diagonalization%mom_occA)

            CALL section_vals_val_get(scf_section, "MOM%OCC_BETA", &
                                      i_vals=scf_control%diagonalization%mom_occB)

            CALL section_vals_val_get(scf_section, "MOM%PROJ_FORMULA", &
                                      i_val=scf_control%diagonalization%mom_proj_formula)
         END IF
      END IF

      ! Read ROKS parameters
      CALL section_vals_val_get(scf_section, "ROKS_SCHEME", i_val=scf_control%roks_scheme)

      SELECT CASE (scf_control%roks_scheme)
      CASE (general_roks)
         ! Read parameters for the general ROKS scheme
         CALL section_vals_val_get(scf_section, "ROKS_F", r_val=scf_control%roks_f)
      CASE (high_spin_roks)
         ! Read high-spin ROKS parameters for the diagonal block
         ! 0 = v)irtual, 1 = o)pen shell, 2 = c)losed shell
         NULLIFY (roks_parameter)
         CALL section_vals_val_get(scf_section, "ROKS_PARAMETERS", r_vals=roks_parameter)
         IF (ASSOCIATED(roks_parameter)) THEN
            scf_control%roks_parameter(2, 2, 1) = roks_parameter(1) ! acc
            scf_control%roks_parameter(2, 2, 2) = roks_parameter(2) ! bcc
            scf_control%roks_parameter(1, 1, 1) = roks_parameter(3) ! aoo
            scf_control%roks_parameter(1, 1, 2) = roks_parameter(4) ! boo
            scf_control%roks_parameter(0, 0, 1) = roks_parameter(5) ! avv
            scf_control%roks_parameter(0, 0, 2) = roks_parameter(6) ! bvv
         END IF
      END SELECT

      ! should be moved to printkey
      CALL section_vals_val_get(scf_section, "eps_lumo", r_val=scf_control%eps_lumos)
      CALL section_vals_val_get(scf_section, "max_iter_lumo", i_val=scf_control%max_iter_lumos)

      ! Extra MOs, e.g. for smearing
      CALL section_vals_val_get(scf_section, "added_mos", i_vals=added_mos)
      CPASSERT(ASSOCIATED(added_mos))
      IF (SIZE(added_mos) > 0) THEN
         scf_control%added_mos(1) = added_mos(1)
         IF (SIZE(added_mos) > 1) THEN
            scf_control%added_mos(2) = added_mos(2)
         END IF
      END IF

      CALL section_vals_val_get(scf_section, "max_scf_history", i_val=scf_control%max_scf_hist)
      CALL section_vals_val_get(scf_section, "eps_scf_history", r_val=scf_control%eps_scf_hist)

      IF (scf_control%level_shift /= 0.0_dp) scf_control%use_cholesky = .FALSE.

      ! Outer SCF subsection
      outer_scf_section => section_vals_get_subs_vals(scf_section, "OUTER_SCF")
      CALL outer_scf_read_parameters(scf_control%outer_scf, outer_scf_section)

      smear_section => section_vals_get_subs_vals(scf_section, "SMEAR")
      CALL init_smear(scf_control%smear)
      CALL read_smear_section(scf_control%smear, smear_section)

      do_mixing = .FALSE.
      mixing_section => section_vals_get_subs_vals(scf_section, "MIXING")
      CALL section_vals_val_get(mixing_section, "_SECTION_PARAMETERS_", &
                                l_val=do_mixing)
      IF (do_mixing) THEN
         CALL section_vals_val_get(mixing_section, "METHOD", &
                                   i_val=scf_control%mixing_method)
         CALL section_vals_val_get(mixing_section, "NMIXING", i_val=scf_control%nmixing)
      END IF ! do mixing

      gce_section => section_vals_get_subs_vals(scf_section, "GCE")
      CALL init_gce(scf_control%gce)
      CALL read_gce_section(scf_control%gce, gce_section)

      CALL timestop(handle)

   END SUBROUTINE scf_c_read_parameters

! **************************************************************************************************
!> \brief ...
!> \param smear ...
! **************************************************************************************************
   SUBROUTINE init_smear(smear)
      TYPE(smear_type), POINTER                          :: smear

      CPASSERT(.NOT. ASSOCIATED(smear))
      ALLOCATE (smear)
      smear%do_smear = .FALSE.
      smear%method = smear_energy_window
      smear%electronic_temperature = 0.0_dp
      smear%eps_fermi_dirac = 1.0E-5_dp
      smear%fixed_mag_mom = -100.0_dp
      smear%window_size = 0.0_dp
      smear%smearing_width = 0.0_dp
      NULLIFY (smear%list)
   END SUBROUTINE init_smear

! **************************************************************************************************
!> \brief ...
!> \param smear ...
!> \param smear_section ...
! **************************************************************************************************
   SUBROUTINE read_smear_section(smear, smear_section)
      TYPE(smear_type), POINTER                          :: smear
      TYPE(section_vals_type), POINTER                   :: smear_section

      REAL(KIND=dp), DIMENSION(:), POINTER               :: r_vals

      NULLIFY (r_vals)

      CALL section_vals_val_get(smear_section, "_SECTION_PARAMETERS_", &
                                l_val=smear%do_smear)
      IF (smear%do_smear) THEN
         CALL section_vals_val_get(smear_section, "METHOD", &
                                   i_val=smear%method)
         CALL section_vals_val_get(smear_section, "ELECTRONIC_TEMPERATURE", &
                                   r_val=smear%electronic_temperature)
         CALL section_vals_val_get(smear_section, "EPS_FERMI_DIRAC", &
                                   r_val=smear%eps_fermi_dirac)
         CALL section_vals_val_get(smear_section, "WINDOW_SIZE", &
                                   r_val=smear%window_size)
         IF (smear%method == smear_gaussian .OR. &
             smear%method == smear_mp .OR. &
             smear%method == smear_mv) THEN
            CALL section_vals_val_get(smear_section, "SIGMA", &
                                      r_val=smear%smearing_width)
         END IF
         IF (smear%method == smear_list) THEN
            CALL section_vals_val_get(smear_section, "LIST", &
                                      r_vals=r_vals)
            CPASSERT(ASSOCIATED(r_vals))
            ALLOCATE (smear%list(SIZE(r_vals)))
            smear%list = r_vals
         END IF
         CALL section_vals_val_get(smear_section, "FIXED_MAGNETIC_MOMENT", &
                                   r_val=smear%fixed_mag_mom)
      END IF ! do smear
   END SUBROUTINE read_smear_section

! **************************************************************************************************
!> \brief ...
!> \param gce ...
! **************************************************************************************************
   SUBROUTINE init_gce(gce)
      TYPE(gce_type), POINTER                            :: gce

      CPASSERT(.NOT. ASSOCIATED(gce))
      ALLOCATE (gce)
      gce%do_gce = .FALSE.
      gce%target_workfunction = 0.16_dp
      gce%prev_workfunction = -1001.0_dp
      gce%ref_esp = 0.0_dp
      gce%mixing_coef = 0.3_dp
      gce%scheme = 1
      gce%update_every_step = .TRUE.
      gce%check_convergence = .FALSE.
      gce%capacitance = 1.0_dp*evolt
      gce%cap_min = 0.05_dp*evolt
      gce%cap_max = 2.0_dp*evolt
      gce%ne_step_max = 0.5_dp
      gce%scf_delta = 1.0_dp
      gce%update_delta = 1.0E-2_dp
      gce%early_update_delta = 0.0_dp
      gce%early_update_min_steps = 3
      gce%early_update_max = 2
      gce%early_wf_stable = 0.0_dp
      gce%pin_mixer_g0 = .FALSE.
      gce%response_guess = .FALSE.
      gce%early_secant = .FALSE.
      gce%final_tol_frac = 1.0_dp
      gce%carry_gain = 0.0_dp
      gce%cap_refresh_min = 0.0_dp
      gce%ncall_restore = 0
      gce%keep_mixer_below = 0.0_dp
      gce%n_outer = 0
      gce%max_outer = 20
      gce%min_steps = 12
      gce%update_interval = 1
      gce%convergence_mode = 3
      gce%max_outer_now = 20
      gce%first_update_size = 0.0_dp
      gce%target_mu = -1000.0_dp
      gce%solvent_only_for_charge = .FALSE.
      gce%ion_step = 0
      gce%steps_since_update = 1000000
      gce%iter_base = 0
      gce%secant_step = -1
      gce%secant_ready = .FALSE.
      gce%fermi_tol = 0.01_dp/evolt
      gce%fermi_update_tol = 0.02_dp/evolt
      gce%n_target = -1.0_dp
      gce%n_input = -1.0_dp
      gce%n_excess_init = 0.0_dp
      gce%n_last = 0.0_dp
      gce%wf_last = 0.0_dp
      gce%target_potential = -1000.0_dp
      gce%phi_she = 4.600_dp/evolt
      gce%constant_potential_requested = .FALSE.
      gce%n_step = -1.0_dp
      gce%e_step = 0.0_dp
      gce%wf_step = -1001.0_dp
      gce%wf_on_step = 0.0_dp
      gce%have_wf_on_step = .FALSE.
      gce%wf_lo = 0.0_dp
      gce%wf_hi = 0.0_dp
      gce%n_at_lo = 0.0_dp
      gce%n_at_hi = 0.0_dp
      gce%enforce_target = .FALSE.
      gce%state_pending = .FALSE.
      gce%state_report = .FALSE.
      gce%replay_pending = .FALSE.
   END SUBROUTINE init_gce

! **************************************************************************************************
!> \brief ...
!> \param gce ...
!> \param gce_section ...
! **************************************************************************************************
   SUBROUTINE read_gce_section(gce, gce_section)
      TYPE(gce_type), POINTER                            :: gce
      TYPE(section_vals_type), POINTER                   :: gce_section

      INTEGER                                            :: ival
      LOGICAL                                            :: wf_explicit
      REAL(KIND=dp)                                      :: rval

      CALL section_vals_val_get(gce_section, "_SECTION_PARAMETERS_", &
                                l_val=gce%do_gce)
      IF (gce%do_gce) THEN
         CALL section_vals_val_get(gce_section, "TARGET_WORKFUNCTION", &
                                   r_val=gce%target_workfunction)
         CALL section_vals_val_get(gce_section, "MIXING_COEF", &
                                   r_val=gce%mixing_coef)
         CALL section_vals_val_get(gce_section, "UPDATE_SCHEME", i_val=gce%scheme)
         CALL section_vals_val_get(gce_section, "UPDATE_EVERY", i_val=ival)
         gce%update_every_step = (ival == 1)
         CALL section_vals_val_get(gce_section, "CONVERGENCE_CHECK", l_val=gce%check_convergence)
         CALL section_vals_val_get(gce_section, "CAPACITANCE_INIT", r_val=rval)
         gce%capacitance = rval*evolt   ! e/V -> e/Hartree
         CALL section_vals_val_get(gce_section, "CAPACITANCE_MIN", r_val=rval)
         gce%cap_min = rval*evolt
         CALL section_vals_val_get(gce_section, "CAPACITANCE_MAX", r_val=rval)
         gce%cap_max = rval*evolt
         CALL section_vals_val_get(gce_section, "NE_STEP_MAX", r_val=gce%ne_step_max)
         CALL section_vals_val_get(gce_section, "UPDATE_SCF_DELTA", r_val=gce%update_delta)
         CALL section_vals_val_get(gce_section, "EARLY_UPDATE_DELTA", r_val=gce%early_update_delta)
         CALL section_vals_val_get(gce_section, "EARLY_UPDATE_MIN_STEPS", i_val=gce%early_update_min_steps)
         CALL section_vals_val_get(gce_section, "KEEP_MIXER_BELOW", r_val=gce%keep_mixer_below)
         CALL section_vals_val_get(gce_section, "PIN_MIXER_G0", l_val=gce%pin_mixer_g0)
         CALL section_vals_val_get(gce_section, "RESPONSE_GUESS", l_val=gce%response_guess)
         CALL section_vals_val_get(gce_section, "EARLY_UPDATE_SECANT", l_val=gce%early_secant)
         CALL section_vals_val_get(gce_section, "FINAL_TOL_FRAC", r_val=gce%final_tol_frac)
         CALL section_vals_val_get(gce_section, "CARRY_GAIN", r_val=gce%carry_gain)
         CALL section_vals_val_get(gce_section, "CAP_REFRESH_MIN", r_val=gce%cap_refresh_min)
         CALL section_vals_val_get(gce_section, "ENFORCE_TARGET", l_val=gce%enforce_target)
         IF (gce%final_tol_frac <= 0.0_dp .OR. gce%final_tol_frac > 1.0_dp) &
            CPABORT("GCE: FINAL_TOL_FRAC must be in (0,1]")
         CALL section_vals_val_get(gce_section, "EARLY_UPDATE_MAX", i_val=gce%early_update_max)
         CALL section_vals_val_get(gce_section, "EARLY_UPDATE_WF_STABLE", r_val=gce%early_wf_stable)
         CALL section_vals_val_get(gce_section, "MAX_OUTER_UPDATES", i_val=gce%max_outer)
         CALL section_vals_val_get(gce_section, "EXCESS_ELECTRONS", r_val=gce%n_excess_init)
         CALL section_vals_val_get(gce_section, "MIN_STEPS_AFTER_UPDATE", i_val=gce%min_steps)
         CALL section_vals_val_get(gce_section, "UPDATE_INTERVAL", i_val=gce%update_interval)
         IF (gce%update_interval < 1) CPABORT("GCE: UPDATE_INTERVAL must be >= 1")
         CALL section_vals_val_get(gce_section, "CONVERGENCE_MODE", i_val=gce%convergence_mode)
         IF (gce%convergence_mode < 1 .OR. gce%convergence_mode > 3) &
            CPABORT("GCE: CONVERGENCE_MODE must be 1, 2 or 3")
         gce%max_outer_now = gce%max_outer
         CALL section_vals_val_get(gce_section, "FIRST_UPDATE_SIZE", r_val=gce%first_update_size)
         CALL section_vals_val_get(gce_section, "SOLVENT_ONLY_FOR_CHARGE", l_val=gce%solvent_only_for_charge)
         CALL section_vals_val_get(gce_section, "TARGET_CHEMICAL_POTENTIAL", r_val=gce%target_mu)
         CALL section_vals_val_get(gce_section, "FERMI_TOL", r_val=gce%fermi_tol)
         CALL section_vals_val_get(gce_section, "FERMI_UPDATE_TOL", r_val=gce%fermi_update_tol)
         CALL section_vals_val_get(gce_section, "PHI_SHE", r_val=gce%phi_she)
         CALL section_vals_val_get(gce_section, "TARGET_POTENTIAL", r_val=gce%target_potential)
         IF (gce%target_potential > -999.0_dp/evolt) THEN
            ! TARGET_CHEMICAL_POTENTIAL = -(phi_SHE + U_SHE): the target work function is phi_SHE + U
            gce%target_workfunction = gce%phi_she + gce%target_potential
         END IF
         ! TARGET_CHEMICAL_POTENTIAL (absolute mu_e) overrides both
         IF (gce%target_mu > -999.0_dp/evolt) gce%target_workfunction = -gce%target_mu
         ! the k-point SCF path has no grand-canonical occupation; qs_scf aborts there if any steering is requested.
         ! Requested = outer updates allowed, a target potential / work function / chemical potential given, or a
         ! fractional excess charge. A neutral block (MAX_OUTER_UPDATES 0, EXCESS_ELECTRONS 0) is not a request.
         CALL section_vals_val_get(gce_section, "TARGET_WORKFUNCTION", explicit=wf_explicit)
         gce%constant_potential_requested = (gce%max_outer > 0) .OR. wf_explicit .OR. &
                                            (gce%target_potential > -999.0_dp/evolt) .OR. &
                                            (gce%target_mu > -999.0_dp/evolt) .OR. (gce%n_excess_init /= 0.0_dp)
         ! CAPACITANCE and DYNAMICS steer N_e to a work function: without a target the default (0.16 Ha) would be used silently.
         ! A fixed-N_e block (MAX_OUTER_UPDATES 0 with UPDATE_EVERY SCF_CYCLE, no early update, no carry) needs none.
         IF (gce%scheme == 3 .OR. (gce%scheme == 2 .AND. (gce%update_every_step .OR. gce%max_outer > 0 .OR. &
                                                           gce%early_update_delta > 0.0_dp .OR. gce%carry_gain > 0.0_dp))) THEN
            IF (.NOT. (wf_explicit .OR. gce%target_potential > -999.0_dp/evolt .OR. gce%target_mu > -999.0_dp/evolt)) &
               CALL cp_abort(__LOCATION__, &
                             "GCE: UPDATE_SCHEME CAPACITANCE and DYNAMICS need TARGET_POTENTIAL, TARGET_WORKFUNCTION or "// &
                             "TARGET_CHEMICAL_POTENTIAL (a fixed-N_e block with MAX_OUTER_UPDATES 0 needs none)")
         END IF
         IF (gce%cap_min <= 0.0_dp .OR. gce%cap_max < gce%cap_min .OR. gce%capacitance <= 0.0_dp) THEN
            CPABORT("GCE: CAPACITANCE_INIT/MIN/MAX must be positive with MIN <= MAX")
         END IF
         ! gce_dyn: UPDATE_SCHEME DYNAMICS = the capacitance-scheme occupations at the electron number of the
         ! dynamics, one SCF per ionic step, no N_e update inside the SCF after the first force evaluation.
         ! It restores its own state (read_gce_dyn_section); the capacitance scheme restores gce_state_read's.
         IF (gce%scheme == 3) THEN
            CALL read_gce_dyn_section(gce, gce_section)
         ELSE IF (gce%scheme == 2) THEN
            CALL gce_state_read(gce, gce_section)
         END IF
      END IF ! do gce
   END SUBROUTINE read_gce_section

! **************************************************************************************************
!> \brief the restart state as one printable line
!> \param gce ...
!> \return ...
! **************************************************************************************************
   FUNCTION gce_state_string(gce) RESULT(str)
      TYPE(gce_type), POINTER                            :: gce
      CHARACTER(LEN=500)                                 :: str

      WRITE (str, "(A,ES25.16E3,A,ES25.16E3,A,ES25.16E3,A,I0,A,ES25.16E3,A,ES25.16E3,A,L1,A,I0)") &
         "N_E ", gce%n_target, " N_INPUT ", gce%n_input, " CAPACITANCE ", gce%capacitance, " ION_STEP ", gce%ion_step, &
         " N_LAST ", gce%n_last, " WF_LAST ", gce%wf_last, " SECANT_READY ", gce%secant_ready, &
         " STEPS_SINCE_UPDATE ", gce%steps_since_update
   END FUNCTION gce_state_string

! **************************************************************************************************
!> \brief state of the electron-number control -> RESTART_* keywords of &GCE in the input tree; the restart file that
!>        the MD / geometry optimisation driver writes after this ionic step carries them. To persist another quantity
!>        add its RESTART_<NAME> keyword in input_cp2k_scf.F and one line here and in gce_state_read.
!>        RESTART_VALID, _N_E, _N_INPUT, _CAPACITANCE and _ION_STEP are the keywords that the electron-number
!>        dynamics (UPDATE_SCHEME DYNAMICS) also uses.
!> \param gce ...
!> \param gce_section the SCF%GCE section of the running input
! **************************************************************************************************
   SUBROUTINE gce_state_write(gce, gce_section)
      TYPE(gce_type), POINTER                            :: gce
      TYPE(section_vals_type), POINTER                   :: gce_section

      IF (.NOT. ASSOCIATED(gce)) RETURN
      IF (.NOT. gce%do_gce) RETURN
      IF (gce%scheme /= 2 .OR. gce%dyn%active .OR. gce%n_step < 0.0_dp) RETURN   ! dynamics: gce_dyn_write_state
      ! N_E is the electron number the last evaluation ran at (before CARRY_GAIN): the evaluation that repeats it in
      ! the restarted run is the same one, and its CARRY_GAIN gives the same N_e for the next step
      CALL section_vals_val_set(gce_section, "RESTART_VALID", l_val=.TRUE.)
      CALL section_vals_val_set(gce_section, "RESTART_N_E", r_val=gce%n_step)
      CALL section_vals_val_set(gce_section, "RESTART_N_INPUT", r_val=gce%n_input)
      CALL section_vals_val_set(gce_section, "RESTART_CAPACITANCE", r_val=gce%capacitance)
      CALL section_vals_val_set(gce_section, "RESTART_ION_STEP", i_val=gce%ion_step)
      CALL section_vals_val_set(gce_section, "RESTART_N_LAST", r_val=gce%n_last)
      CALL section_vals_val_set(gce_section, "RESTART_WF_LAST", r_val=gce%wf_last)
      CALL section_vals_val_set(gce_section, "RESTART_SECANT_READY", l_val=gce%secant_ready)
      CALL section_vals_val_set(gce_section, "RESTART_LAST_WF", r_val=gce%wf_step)
      CALL section_vals_val_set(gce_section, "RESTART_STEPS_SINCE_UPDATE", i_val=gce%steps_since_update)
   END SUBROUTINE gce_state_write

! **************************************************************************************************
!> \brief MD hook (md_energy): the conserved quantity in grand-potential form. The MD driver integrates with the
!>        energy E = energy%total, so Cons = K + E + (thermostat energies) is conserved only at fixed N_e; with a
!>        varying N_e the conserved quantity is Cons_Omega = Cons + WF_target n_excess (mu_e = -WF_target),
!>        n_excess = N_e - N_input of the force evaluation of this step. dE is the MD potential energy minus
!>        the energy%total printed in the GCE| omega line of the same step (0 when they are the same quantity).
!>        Not printed with UPDATE_SCHEME DYNAMICS: there the conserved quantity also holds the kinetic energy of N_e
!>        and the thermostat (Hextra of the GCE| dyn line).
!> \param gce ...
!> \param ekin kinetic energy of the MD step [Ha]
!> \param epot potential energy the MD driver uses [Ha]
!> \param econs Cons Qty of the MD step [Ha]
!> \param itimes MD step number
! **************************************************************************************************
   SUBROUTINE gce_md_report(gce, ekin, epot, econs, itimes)
      TYPE(gce_type), POINTER                            :: gce
      REAL(KIND=dp), INTENT(IN)                          :: ekin, epot, econs
      INTEGER, INTENT(IN)                                :: itimes

      INTEGER                                            :: unit_nr
      REAL(KIND=dp)                                      :: n_exc

      IF (.NOT. ASSOCIATED(gce)) RETURN
      IF (.NOT. gce%do_gce) RETURN
      IF (gce%scheme /= 2 .OR. gce%dyn%active .OR. gce%n_step < 0.0_dp .OR. gce%n_input < 0.0_dp) RETURN
      unit_nr = cp_logger_get_default_io_unit()
      IF (unit_nr <= 0) RETURN
      n_exc = gce%n_step - gce%n_input
      WRITE (UNIT=unit_nr, FMT="(T3,A,I8,A,F22.12,A,F22.12,A,ES11.3,A,F10.6,A,F22.12,A,F22.12,A,F22.12,A)") &
         "GCE| conserved  MD step ", itimes, "  Kin ", ekin, " Ha  E ", epot, " Ha  dE ", epot - gce%e_step, &
         " Ha  n_exc ", n_exc, " e  Omega ", epot + gce%target_workfunction*n_exc, " Ha  Cons ", econs, &
         " Ha  Cons_Omega ", econs + gce%target_workfunction*n_exc, " Ha"
   END SUBROUTINE gce_md_report

! **************************************************************************************************
!> \brief reads the state from the RESTART_* keywords (a restart file used as input, or EXT_RESTART with RESTART_GCE).
!>        The state is the one after the last completed ionic step; the first force evaluation of the restarted run
!>        repeats that step (MD, GEO_OPT) at the same N_e. replay_pending makes qs_scf take the step counter back by one
!>        for that evaluation (same update-step phase as in the first run). The state itself is restored unchanged.
!> \param gce ...
!> \param gce_section ...
! **************************************************************************************************
   SUBROUTINE gce_state_read(gce, gce_section)
      TYPE(gce_type), POINTER                            :: gce
      TYPE(section_vals_type), POINTER                   :: gce_section

      LOGICAL                                            :: valid

      CALL section_vals_val_get(gce_section, "RESTART_VALID", l_val=valid)
      IF (.NOT. valid) RETURN
      CALL section_vals_val_get(gce_section, "RESTART_N_E", r_val=gce%n_target)
      CALL section_vals_val_get(gce_section, "RESTART_N_INPUT", r_val=gce%n_input)
      CALL section_vals_val_get(gce_section, "RESTART_CAPACITANCE", r_val=gce%capacitance)
      CALL section_vals_val_get(gce_section, "RESTART_ION_STEP", i_val=gce%ion_step)
      CALL section_vals_val_get(gce_section, "RESTART_N_LAST", r_val=gce%n_last)
      CALL section_vals_val_get(gce_section, "RESTART_WF_LAST", r_val=gce%wf_last)
      CALL section_vals_val_get(gce_section, "RESTART_SECANT_READY", l_val=gce%secant_ready)
      CALL section_vals_val_get(gce_section, "RESTART_LAST_WF", r_val=gce%wf_step)
      CALL section_vals_val_get(gce_section, "RESTART_STEPS_SINCE_UPDATE", i_val=gce%steps_since_update)
      gce%replay_pending = .TRUE.
      gce%state_pending = .TRUE.
      gce%state_report = .TRUE.
   END SUBROUTINE gce_state_read

! **************************************************************************************************
!> \brief gce_dyn: keywords of UPDATE_SCHEME DYNAMICS and the state written by the previous run.
!>        The dynamics is the capacitance scheme (scheme 2) in cycle mode with N_e supplied by the integrator.
!> \param gce ...
!> \param gce_section ...
! **************************************************************************************************
   SUBROUTINE read_gce_dyn_section(gce, gce_section)
      TYPE(gce_type), POINTER                            :: gce
      TYPE(section_vals_type), POINTER                   :: gce_section

      INTEGER                                            :: m
      LOGICAL                                            :: explicit
      REAL(KIND=dp)                                      :: rval
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rvals

      NULLIFY (rvals)
      gce%scheme = 2
      gce%dyn%active = .TRUE.
      gce%constant_potential_requested = .TRUE.
      IF (gce%update_every_step) &
         CPABORT("GCE UPDATE_SCHEME DYNAMICS needs UPDATE_EVERY SCF_CYCLE")
      IF (gce%carry_gain > 0.0_dp) &
         CPABORT("GCE UPDATE_SCHEME DYNAMICS does not combine with CARRY_GAIN")
      CALL section_vals_val_get(gce_section, "POTENTIOSTAT_PERIOD", r_val=gce%dyn%period)
      IF (gce%dyn%period <= 0.0_dp) &
         CALL cp_abort(__LOCATION__, &
               "GCE UPDATE_SCHEME DYNAMICS needs POTENTIOSTAT_PERIOD > 0: M_n = 1/(C omega^2) is set from the period")
      CALL section_vals_val_get(gce_section, "POTENTIOSTAT_THERMOSTAT_PERIOD", r_val=gce%dyn%thermostat_period)
      CALL section_vals_val_get(gce_section, "POTENTIOSTAT_CHAIN_LENGTH", i_val=gce%dyn%chain_length)
      IF (gce%dyn%chain_length < 0 .OR. gce%dyn%chain_length > gce_dyn_max_chain) &
         CPABORT("GCE: POTENTIOSTAT_CHAIN_LENGTH must be in 0..16")
      CALL section_vals_val_get(gce_section, "POTENTIOSTAT_TEMPERATURE", explicit=explicit)
      IF (explicit) CALL section_vals_val_get(gce_section, "POTENTIOSTAT_TEMPERATURE", r_val=gce%dyn%kt_input)
      CALL section_vals_val_get(gce_section, "POTENTIOSTAT_CAPACITANCE", r_val=rval)
      gce%dyn%capacitance_input = rval*evolt   ! e/V -> e/Hartree
      CALL section_vals_val_get(gce_section, "POTENTIOSTAT_INITIAL_RATE", r_val=rval)
      gce%dyn%initial_rate = rval*femtoseconds   ! e/fs -> e/a.u.
      ! the integrator moves N_e between the SCF runs: the G = 0 charge of the mixer input must follow it at once
      ! (otherwise the SCF can converge falsely), and a small |dN_e| keeps the Broyden history
      CALL section_vals_val_get(gce_section, "PIN_MIXER_G0", explicit=explicit)
      IF (.NOT. explicit) gce%pin_mixer_g0 = .TRUE.
      CALL section_vals_val_get(gce_section, "KEEP_MIXER_BELOW", explicit=explicit)
      IF (.NOT. explicit) gce%keep_mixer_below = 0.05_dp
      ! state of the previous run (written by qs_scf_gce_dyn after every ionic step)
      CALL section_vals_val_get(gce_section, "RESTART_VALID", l_val=gce%dyn%restart_valid)
      IF (gce%dyn%restart_valid) THEN
         m = gce%dyn%chain_length
         CALL section_vals_val_get(gce_section, "RESTART_N_E", r_val=gce%dyn%n_e)
         CALL section_vals_val_get(gce_section, "RESTART_P_N", r_val=gce%dyn%p_n)
         CALL section_vals_val_get(gce_section, "RESTART_DU_LAST", r_val=gce%dyn%du_last)
         CALL section_vals_val_get(gce_section, "RESTART_CAPACITANCE", r_val=gce%dyn%capacitance)
         CALL section_vals_val_get(gce_section, "RESTART_ION_STEP", i_val=gce%dyn%ion_step_saved)
         CALL section_vals_val_get(gce_section, "RESTART_NEXT_LEN", i_val=gce%dyn%next_len)
         CALL section_vals_val_get(gce_section, "RESTART_OPEN_PENDING", l_val=gce%dyn%open_pending)
         IF (m > 0) THEN
            CALL section_vals_val_get(gce_section, "RESTART_XI", r_vals=rvals)
            IF (SIZE(rvals) /= m) &
               CALL cp_abort(__LOCATION__, &
                     "GCE restart: RESTART_XI length differs from POTENTIOSTAT_CHAIN_LENGTH")
            gce%dyn%xi(1:m) = rvals(1:m)
            CALL section_vals_val_get(gce_section, "RESTART_P_XI", r_vals=rvals)
            IF (SIZE(rvals) /= m) &
               CALL cp_abort(__LOCATION__, &
                     "GCE restart: RESTART_P_XI length differs from POTENTIOSTAT_CHAIN_LENGTH")
            gce%dyn%p_xi(1:m) = rvals(1:m)
         END IF
         CALL section_vals_val_get(gce_section, "RESTART_N_INPUT", r_val=gce%dyn%n_input)
         gce%n_target = gce%dyn%n_e
         gce%n_input = gce%dyn%n_input
         gce%capacitance = gce%dyn%capacitance
         gce%state_pending = .TRUE.   ! the first occupation call checks RESTART_N_INPUT against this input
      END IF
   END SUBROUTINE read_gce_dyn_section

! **************************************************************************************************
!> \brief writes out the scf parameters
!> \param scf_control the object you want to print
!> \param dft_section ...
!> \par History
!>      05.2001 created [Matthias]
!>      09.2002 created separated scf_control type [fawzi]
!> \author Matthias Krack
! **************************************************************************************************
   SUBROUTINE scf_c_write_parameters(scf_control, dft_section)

      TYPE(scf_control_type), INTENT(IN)                 :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'scf_c_write_parameters'

      INTEGER                                            :: handle, output_unit, roks_scheme
      LOGICAL                                            :: roks
      REAL(KIND=dp)                                      :: elec_temp
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(enumeration_type), POINTER                    :: enum
      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: section
      TYPE(section_vals_type), POINTER                   :: scf_section

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()

      NULLIFY (scf_section)
      NULLIFY (section)

      scf_section => section_vals_get_subs_vals(dft_section, "SCF")
      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".scfLog")

      IF (output_unit > 0) THEN

         IF (scf_control%max_scf > 0) THEN

            CALL create_scf_section(section)

            keyword => section_get_keyword(section, "SCF_GUESS")
            CALL keyword_get(keyword, enum=enum)

            IF (.NOT. scf_control%non_selfconsistent .OR. scf_control%force_scf_calculation) THEN
               WRITE (UNIT=output_unit, &
                      FMT="(/,/,T2,A,T25,A,T51,A30,/,T25,56('-'),3(/,T25,A,T76,I5),/, "// &
                      "T25,56('-'),4(/,T25,A,T72,ES9.2),/,T25,56('-'), "// &
                      "1(/,T25,A,T71,F10.6))") &
                  "SCF PARAMETERS", &
                  "Density guess:     ", ADJUSTR(TRIM(enum_i2c(enum, scf_control%density_guess))), &
                  "max_scf:           ", scf_control%max_scf, &
                  "max_scf_history:   ", scf_control%max_scf_hist, &
                  "max_diis:          ", scf_control%max_diis, &
                  "eps_scf:           ", scf_control%eps_scf, &
                  "eps_scf_history:   ", scf_control%eps_scf_hist, &
                  "eps_diis:          ", scf_control%eps_diis, &
                  "eps_eigval:        ", scf_control%eps_eigval, &
                  "level_shift [a.u.]:", scf_control%level_shift
            END IF

            IF (SUM(ABS(scf_control%added_mos)) > 0) THEN
               WRITE (UNIT=output_unit, FMT="(T25,A,T71,2I5)") &
                  "added MOs          ", scf_control%added_mos
            END IF

            IF (scf_control%diagonalization%mom) THEN
               ! TODO extend the output with further parameters
               WRITE (UNIT=output_unit, FMT="(T25,A)") "MOM enabled"
            END IF

            IF (scf_control%mixing_method > 0 .AND. .NOT. scf_control%use_ot .AND. &
                .NOT. scf_control%non_selfconsistent) THEN
               keyword => section_get_keyword(section, "MIXING%METHOD")
               CALL keyword_get(keyword, enum=enum)
               WRITE (UNIT=output_unit, FMT="(T25,A,/,T25,A,T51,A30)") &
                  REPEAT("-", 56), &
                  "Mixing method:      ", ADJUSTR(TRIM(enum_i2c(enum, scf_control%mixing_method)))
               IF (scf_control%mixing_method > 1) THEN
                  WRITE (UNIT=output_unit, FMT="(T47,A34)") "charge density mixing in g-space"
               END IF
            END IF
            IF (scf_control%smear%do_smear) THEN
               keyword => section_get_keyword(section, "SMEAR%METHOD")
               CALL keyword_get(keyword, enum=enum)
               WRITE (UNIT=output_unit, FMT="(T25,A,/,T25,A,T51,A30)") &
                  REPEAT("-", 56), &
                  "Smear method:      ", ADJUSTR(TRIM(enum_i2c(enum, scf_control%smear%method)))
               SELECT CASE (scf_control%smear%method)
               CASE (smear_fermi_dirac)
                  elec_temp = cp_unit_from_cp2k(scf_control%smear%electronic_temperature, &
                                                "K")
                  WRITE (UNIT=output_unit, FMT="(T25,A,T61,F20.1)") &
                     "Electronic temperature [K]:", elec_temp
                  WRITE (UNIT=output_unit, FMT="(T25,A,T71,ES10.2)") &
                     "Electronic temperature [a.u.]:", scf_control%smear%electronic_temperature, &
                     "Accuracy threshold:", scf_control%smear%eps_fermi_dirac
                  IF (scf_control%smear%fixed_mag_mom > 0.0_dp) WRITE (UNIT=output_unit, FMT="(T25,A,T61,F20.1)") &
                     "Fixed magnetic moment set to:", scf_control%smear%fixed_mag_mom
               CASE (smear_gaussian, smear_mp, smear_mv)
                  WRITE (UNIT=output_unit, FMT="(T25,A,T71,ES10.2)") &
                     "Smearing width (sigma) [a.u.]:", scf_control%smear%smearing_width, &
                     "Accuracy threshold:", scf_control%smear%eps_fermi_dirac
                  IF (scf_control%smear%fixed_mag_mom > 0.0_dp) &
                     WRITE (UNIT=output_unit, FMT="(T25,A,T61,F20.1)") &
                     "Fixed magnetic moment set to:", scf_control%smear%fixed_mag_mom
               CASE (smear_energy_window)
                  WRITE (UNIT=output_unit, FMT="(T25,A,T71,F10.6)") &
                     "Smear window [a.u.]:       ", scf_control%smear%window_size
               END SELECT
            END IF

            IF (scf_control%gce%do_gce) THEN
               WRITE (UNIT=output_unit, FMT="(T25,A,/,T25,A,T61,F20.1)") &
                  REPEAT("-", 56), &
                  "Target workfunction [eV]:", scf_control%gce%target_workfunction
               WRITE (UNIT=output_unit, FMT="(T25,A,T61,F20.1)") "Mixing coefficient:", scf_control%gce%mixing_coef
               WRITE (UNIT=output_unit, FMT="(T25,A)") "Grand canonical SCF is activated."
            END IF

            CALL section_vals_val_get(dft_section, "ROKS", l_val=roks)
            IF (roks .AND. (.NOT. scf_control%use_ot)) THEN
               CALL section_vals_val_get(scf_section, "ROKS_SCHEME", &
                                         i_val=roks_scheme)
               keyword => section_get_keyword(section, "ROKS_SCHEME")
               CALL keyword_get(keyword, enum=enum)
               WRITE (UNIT=output_unit, FMT="(T25,A,/,T25,A,T51,A30)") &
                  REPEAT("-", 56), &
                  "ROKS scheme:", ADJUSTR(TRIM(enum_i2c(enum, roks_scheme)))
               SELECT CASE (roks_scheme)
               CASE (general_roks)
                  WRITE (UNIT=output_unit, FMT="(T25,A,T71,F10.6)") &
                     "ROKS parameter f:", scf_control%roks_f
               CASE (high_spin_roks)
                  WRITE (UNIT=output_unit, &
                         FMT="(T25,A,6(/,T25,A,T71,F10.6))") &
                     "ROKS parameters: a)lpha, b)eta; c)losed, o)pen, v)irtual", &
                     "acc", scf_control%roks_parameter(2, 2, 1), &
                     "bcc", scf_control%roks_parameter(2, 2, 2), &
                     "aoo", scf_control%roks_parameter(1, 1, 1), &
                     "boo", scf_control%roks_parameter(1, 1, 2), &
                     "avv", scf_control%roks_parameter(0, 0, 1), &
                     "bvv", scf_control%roks_parameter(0, 0, 2)
               END SELECT
            END IF
            CALL section_release(section)

            IF (scf_control%outer_scf%have_scf) THEN
               WRITE (output_unit, "(T25,56('-'),/,T25,A)") "Outer loop SCF in use "
               SELECT CASE (scf_control%outer_scf%type)
               CASE (outer_scf_none)
                  WRITE (output_unit, '(T25,A)') "No variables optimised in outer loop"
               CASE (outer_scf_ddapc_constraint)
                  WRITE (output_unit, '(T25,A)') "DDAPC constraint enforced"
               CASE (outer_scf_s2_constraint)
                  WRITE (output_unit, '(T25,A)') "S2 constraint enforced"
               CASE (outer_scf_basis_center_opt)
                  WRITE (output_unit, '(T25,A)') "Floating basis function optimization enforced"
               CASE (outer_scf_cdft_constraint)
                  CPABORT("CDFT constraints must be defined in QS&CDFT")
               CASE DEFAULT
                  CPABORT("Unknown outer SCF type")
               END SELECT
               WRITE (output_unit, '(T25,A,T72,ES9.2)') "eps_scf", scf_control%outer_scf%eps_scf
               WRITE (output_unit, '(T25,A,T72,I9)') "max_scf", scf_control%outer_scf%max_scf
               SELECT CASE (scf_control%outer_scf%optimizer)
               CASE (outer_scf_optimizer_none)
                  WRITE (output_unit, '(T25,A)') "No outer loop optimization"
               CASE (outer_scf_optimizer_sd)
                  WRITE (output_unit, '(T25,A)') "Steepest descent optimization"
               CASE (outer_scf_optimizer_bisect)
                  WRITE (output_unit, '(T25,A)') "Gradient bisection"
                  WRITE (output_unit, '(T25,A,T72,I9)') "bisect_trust_count", scf_control%outer_scf%bisect_trust_count
               CASE (outer_scf_optimizer_diis)
                  WRITE (output_unit, '(T25,A)') "DIIS optimization"
                  WRITE (output_unit, '(T25,A,T72,I9)') "DIIS buffer length", &
                     scf_control%outer_scf%diis_buffer_length
               CASE (outer_scf_optimizer_broyden, outer_scf_optimizer_newton, &
                     outer_scf_optimizer_newton_ls)
                  CPABORT("Selected optimizer only compatible with CDFT")
               CASE (outer_scf_optimizer_secant)
                  WRITE (output_unit, '(T25,A)') "Optimization with the secant method"
               CASE DEFAULT
                  CPABORT("Unknown outer SCF optimizer")
               END SELECT
               WRITE (output_unit, '(T25,A,T72,ES9.2)') "step_size", scf_control%outer_scf%step_size
            ELSE
               WRITE (output_unit, "(T25,56('-'),/,T25,A)") "No outer SCF"
            END IF

         END IF ! max_scf > 0

      END IF ! output_unit > 0

      CALL cp_print_key_finished_output(output_unit, logger, scf_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

      CALL timestop(handle)

   END SUBROUTINE scf_c_write_parameters

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param settings ...
!> \param scf_section ...
! **************************************************************************************************
   SUBROUTINE ot_diag_read_input(settings, scf_section)
      TYPE(qs_ot_settings_type)                          :: settings
      TYPE(section_vals_type), POINTER                   :: scf_section

      CHARACTER(len=*), PARAMETER :: routineN = 'ot_diag_read_input'

      INTEGER                                            :: handle, output_unit
      LOGICAL                                            :: explicit
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: ot_section

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".log")

      ! decide default settings
      CALL qs_ot_settings_init(settings)

      ! use ot input new style
      ot_section => section_vals_get_subs_vals(scf_section, "DIAGONALIZATION%OT")
      CALL section_vals_get(ot_section, explicit=explicit)

      CALL ot_readwrite_input(settings, ot_section, output_unit)

      CALL cp_print_key_finished_output(output_unit, logger, scf_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

      CALL timestop(handle)

   END SUBROUTINE ot_diag_read_input

! **************************************************************************************************

END MODULE scf_control_types
