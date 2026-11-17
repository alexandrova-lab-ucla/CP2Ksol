!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-3.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief SCF-driver side of the electron-number dynamics (&SCF &GCE UPDATE_SCHEME DYNAMICS).
!>
!>        One MD step is velocity Verlet in (R, n, P_n) with a symmetric split of the potentiostat:
!>
!>          start of the step   [opening half]  NHC(h); P_n += h (U_last - U_T)       h = next_len dt / 2
!>                              [drift]         n += dt P_n / M_n                      every MD step
!>          SCF at (R_t+dt, n)  -> forces on the ions and the electrode potential U
!>          end of the step     [closing half]  P_n += h (U - U_T); NHC(h)            h = cycle length dt / 2
!>                              state written to the input section (restart)
!>
!>        The kick uses the U of the SCF that was run at the same (R, n) as the force on the ions, so the
!>        scheme is the time-reversible velocity Verlet of the extended system. Without
!>        SOLVENT_ONLY_FOR_CHARGE the potential is measured every step (cycle length 1). With it the solvent
!>        (and so U) exists only on the update steps; n then drifts every MD step and P_n receives the
!>        impulse of the cycle at the update steps (impulse multiple-time-step; the closing half of one
!>        cycle and the opening half of the next act on the same U).
!>        The very first force evaluation of a run is the initial one of the MD driver: it only converges
!>        n to the target (cold start, secant capacitance) and starts the state; a continuation run
!>        repeats the evaluation of the last step and leaves the restored state untouched.
!> \author  ai-cp2k-lpb (2026-11-08)
! **************************************************************************************************
MODULE qs_scf_gce_dyn

   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE gce_dyn_types,                   ONLY: gce_dyn_drift,&
                                              gce_dyn_kick,&
                                              gce_dyn_kinetic_n,&
                                              gce_dyn_next_len,&
                                              gce_dyn_nhc_propagate,&
                                              gce_dyn_set_masses,&
                                              gce_dyn_thermostat_energy,&
                                              gce_dyn_type
   USE global_types,                    ONLY: global_environment_type
   USE input_constants,                 ONLY: mol_dyn_run
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: evolt
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE scf_control_types,               ONLY: gce_type,&
                                              scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_gce_dyn'

   PUBLIC :: gce_dyn_read_md_input, gce_dyn_begin_step, gce_dyn_end_step, gce_dyn_note_solvent_off

CONTAINS

! **************************************************************************************************
!> \brief Hook of qs_init: the MD time step and temperature come from MOTION&MD, so the dynamics needs
!>        neither a copy of the time step in &GCE nor a change of the MD drivers.
!> \param scf_control ...
!> \param root_section ...
!> \param globenv ...
! **************************************************************************************************
   SUBROUTINE gce_dyn_read_md_input(scf_control, root_section, globenv)
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), OPTIONAL, POINTER         :: root_section
      TYPE(global_environment_type), OPTIONAL, POINTER   :: globenv

      TYPE(section_vals_type), POINTER                   :: md_section

      IF (.NOT. ASSOCIATED(scf_control)) RETURN
      IF (.NOT. scf_control%gce%do_gce) RETURN
      IF (.NOT. scf_control%gce%dyn%active) RETURN
      IF (.NOT. (PRESENT(root_section) .AND. PRESENT(globenv))) RETURN
      IF (globenv%run_type_id /= mol_dyn_run) &
         CALL cp_abort(__LOCATION__, &
               "GCE UPDATE_SCHEME DYNAMICS needs RUN_TYPE MD (the electron number moves with the MD time)")
      md_section => section_vals_get_subs_vals(root_section, "MOTION%MD")
      CALL section_vals_val_get(md_section, "TIMESTEP", r_val=scf_control%gce%dyn%dt)
      CALL section_vals_val_get(md_section, "TEMPERATURE", r_val=scf_control%gce%dyn%md_kt)
      scf_control%gce%dyn%have_md = .TRUE.
   END SUBROUTINE gce_dyn_read_md_input

! **************************************************************************************************
!> \brief Before the SCF of a force evaluation: opening half step and drift of n, n -> the SCF.
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE gce_dyn_begin_step(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      REAL(KIND=dp)                                      :: dt_obs, h, kt
      TYPE(gce_type), POINTER                            :: gce

      gce => qs_env%scf_control%gce
      IF (.NOT. gce%dyn%active) RETURN
      IF (.NOT. gce%dyn%have_md) &
         CALL cp_abort(__LOCATION__, &
               "GCE DYNAMICS: MOTION&MD TIMESTEP was not read (qs_init got no root section)")
      gce%dyn%have_wf_on = .FALSE.
      gce%dyn%is_repeat = .FALSE.

      IF (.NOT. gce%dyn%initialised) THEN
         ! first force evaluation of this run: the MD driver's evaluation at the starting configuration
         gce%dyn%initialised = .TRUE.
         gce%dyn%is_initial_call = .TRUE.
         gce%dyn%last_sim_step = qs_env%sim_step
         gce%dyn%sim_time_last = qs_env%sim_time
         IF (gce%dyn%restart_valid) THEN
            gce%dyn%restored = .TRUE.
            kt = potentiostat_kt(gce%dyn)
            CALL gce_dyn_set_masses(gce%dyn, gce%dyn%capacitance, kt)
            ! the evaluation repeats the last completed step: same solvent phase, ion_step ends where it ended
            gce%ion_step = gce%dyn%ion_step_saved - 1
            gce%dyn%max_outer_cap = 0
         ELSE
            gce%dyn%max_outer_cap = gce%max_outer
         END IF
         RETURN
      END IF

      gce%dyn%is_initial_call = .FALSE.
      gce%dyn%max_outer_cap = 0
      IF (qs_env%sim_step == gce%dyn%last_sim_step) THEN
         ! a second evaluation inside the same MD step: no integration
         gce%dyn%is_repeat = .TRUE.
         RETURN
      END IF
      ! the time step is read off two evaluations of this run that both belong to MD steps; the counters of the initial
      ! evaluation can still be those of a fresh MD driver (zero) when the run starts from an MD restart, so the
      ! first step only sets the baseline
      IF (gce%dyn%dt_baseline) THEN
         dt_obs = (qs_env%sim_time - gce%dyn%sim_time_last)/REAL(qs_env%sim_step - gce%dyn%last_sim_step, KIND=dp)
         IF (ABS(dt_obs - gce%dyn%dt) > 1.0E-6_dp*gce%dyn%dt) &
            CALL cp_abort(__LOCATION__, &
                  "GCE DYNAMICS: the MD time step is not MOTION&MD TIMESTEP (variable time step is not supported)")
      END IF
      gce%dyn%last_sim_step = qs_env%sim_step
      gce%dyn%sim_time_last = qs_env%sim_time
      gce%dyn%dt_baseline = .TRUE.

      IF (gce%dyn%open_pending) THEN
         h = 0.5_dp*REAL(gce%dyn%next_len, KIND=dp)*gce%dyn%dt
         CALL gce_dyn_nhc_propagate(gce%dyn, h)
         CALL gce_dyn_kick(gce%dyn, h, gce%dyn%du_last)
         gce%dyn%open_pending = .FALSE.
      END IF
      CALL gce_dyn_drift(gce%dyn, gce%dyn%dt)
      gce%n_target = gce%dyn%n_e
   END SUBROUTINE gce_dyn_begin_step

! **************************************************************************************************
!> \brief SOLVENT_ONLY_FOR_CHARGE: the solvent is about to be switched off at the converged potential
!>        measurement of an update step; keep that work function and energy for the potentiostat.
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE gce_dyn_note_solvent_off(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      TYPE(qs_energy_type), POINTER                      :: energy

      IF (.NOT. qs_env%scf_control%gce%dyn%active) RETURN
      CALL get_qs_env(qs_env, energy=energy)
      qs_env%scf_control%gce%dyn%wf_on = qs_env%scf_control%gce%prev_workfunction
      qs_env%scf_control%gce%dyn%e_on = energy%total
      qs_env%scf_control%gce%dyn%have_wf_on = .TRUE.
   END SUBROUTINE gce_dyn_note_solvent_off

! **************************************************************************************************
!> \brief After the SCF of a force evaluation (before gce%ion_step is incremented): closing half step,
!>        printout, state to the input section.
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE gce_dyn_end_step(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: output_unit
      LOGICAL                                            :: flash, meas
      REAL(KIND=dp)                                      :: c_use, du, e_extra, e_kin, e_th, h, kt, &
                                                            omega, wf
      TYPE(gce_type), POINTER                            :: gce
      TYPE(qs_energy_type), POINTER                      :: energy

      gce => qs_env%scf_control%gce
      IF (.NOT. gce%dyn%active) RETURN
      IF (gce%dyn%is_repeat) RETURN
      CALL get_qs_env(qs_env, energy=energy)
      output_unit = cp_logger_get_default_io_unit()

      flash = gce%solvent_only_for_charge
      meas = (.NOT. flash) .OR. gce%ion_step == 0 .OR. MOD(gce%ion_step + 1, gce%update_interval) == 0
      ! work function of the potential measurement: the solvent-on one when the solvent was switched off afterwards
      IF (flash .AND. meas .AND. gce%dyn%have_wf_on) THEN
         wf = gce%dyn%wf_on
      ELSE
         wf = gce%ref_esp - energy%efermi
      END IF
      du = wf - gce%target_workfunction

      IF (gce%dyn%is_initial_call .AND. gce%dyn%restored) THEN
         ! continuation: the evaluation repeats the last completed step, the restored state stays as written
         IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(T3,A,I7,A,F14.8,A,F12.6,A)") &
            "GCE| dyn restored  step ", gce%ion_step, "  N_e ", gce%dyn%n_e, " e  U-U_T of the repeated SCF ", du*evolt, " V"
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(T3,A,ES24.16,A,ES24.16,A,ES24.16,A,ES24.16,A,I0,A,I0,A,L1)") &
               "GCE| dyn restored state  N_e ", gce%dyn%n_e, "  N_input ", gce%dyn%n_input, "  P_n ", gce%dyn%p_n, &
               "  dU_last ", gce%dyn%du_last, "  ion_step ", gce%dyn%ion_step_saved, "  next_len ", gce%dyn%next_len, &
               "  open_pending ", gce%dyn%open_pending
            WRITE (UNIT=output_unit, FMT="(T3,A,16ES24.16)") "GCE| dyn restored xi    ", gce%dyn%xi(1:gce%dyn%chain_length)
            WRITE (UNIT=output_unit, FMT="(T3,A,16ES24.16)") "GCE| dyn restored p_xi  ", gce%dyn%p_xi(1:gce%dyn%chain_length)
         END IF
         RETURN
      END IF

      IF (gce%dyn%is_initial_call) THEN
         ! fresh start: n is wherever the first (cold-start) convergence left it
         gce%dyn%n_e = gce%n_target
         gce%dyn%n_input = gce%n_input
         IF (gce%dyn%capacitance_input > 0.0_dp) THEN
            c_use = gce%dyn%capacitance_input
         ELSE IF (gce%cap_measured) THEN
            c_use = gce%capacitance
         ELSE
            CALL cp_abort(__LOCATION__, &
                  "GCE DYNAMICS: no capacitance. Set POTENTIOSTAT_CAPACITANCE or start off target with MAX_OUTER_UPDATES > 1")
         END IF
         kt = potentiostat_kt(gce%dyn)
         CALL gce_dyn_set_masses(gce%dyn, c_use, kt)
         ! the kick is applied once per cycle of CYCLE MD steps: stable only for a period well above that
         IF (gce%dyn%period < 10.0_dp*REAL(MERGE(gce%update_interval, 1, flash), KIND=dp)*gce%dyn%dt) &
            CPWARN("GCE DYNAMICS: POTENTIOSTAT_PERIOD below 10 update cycles")
         gce%dyn%p_n = gce%dyn%mass_n*gce%dyn%initial_rate
         gce%dyn%xi = 0.0_dp
         gce%dyn%p_xi = 0.0_dp
         gce%dyn%du_last = du
         gce%dyn%next_len = gce_dyn_next_len(gce%ion_step, gce%update_interval, flash)
         gce%dyn%open_pending = .TRUE.
      ELSE IF (meas) THEN
         h = 0.5_dp*REAL(gce%dyn%next_len, KIND=dp)*gce%dyn%dt
         CALL gce_dyn_kick(gce%dyn, h, du)
         CALL gce_dyn_nhc_propagate(gce%dyn, h)
         gce%dyn%du_last = du
         gce%dyn%next_len = gce_dyn_next_len(gce%ion_step, gce%update_interval, flash)
         gce%dyn%open_pending = .TRUE.
      END IF
      gce%dyn%n_steps = gce%dyn%n_steps + 1
      gce%dyn%ion_step_saved = gce%ion_step + 1

      ! observability: Omega = E + WF_T (n - n_input); Hextra = Omega - E + K_n + E_th, so that
      ! H_ext = (MD Cons Qty) + Hextra when the energy that drives the ions is the one that defines U
      e_kin = gce_dyn_kinetic_n(gce%dyn)
      e_th = gce_dyn_thermostat_energy(gce%dyn)
      e_extra = gce%target_workfunction*(gce%dyn%n_e - gce%dyn%n_input)
      omega = energy%total + e_extra
      IF (output_unit > 0) WRITE (UNIT=output_unit, &
                                  FMT="(T3,A,I7,A,I7,A,F14.8,A,F12.6,A,L1,A,ES16.8,A,F20.10,A,F18.10,A,ES14.6,A,ES14.6)") &
         "GCE| dyn  step ", gce%ion_step, "  md_step ", qs_env%sim_step, "  N_e ", gce%dyn%n_e, " e  dU ", du*evolt, &
         " V  meas ", meas, "  P_n ", gce%dyn%p_n, "  Omega ", omega, " Ha  Hextra ", e_extra + e_kin + e_th, &
         " Ha  K_n ", e_kin, "  E_th ", e_th

      CALL gce_dyn_write_state(qs_env)
   END SUBROUTINE gce_dyn_end_step

! **************************************************************************************************
!> \brief kT of the potentiostat thermostat: POTENTIOSTAT_TEMPERATURE or MOTION&MD TEMPERATURE.
! **************************************************************************************************
   FUNCTION potentiostat_kt(dyn) RESULT(kt)
      TYPE(gce_dyn_type), INTENT(IN)                     :: dyn
      REAL(KIND=dp)                                      :: kt

      IF (dyn%kt_input > 0.0_dp) THEN
         kt = dyn%kt_input
      ELSE
         kt = dyn%md_kt
      END IF
      IF (dyn%chain_length > 0 .AND. kt <= 0.0_dp) &
         CALL cp_abort(__LOCATION__, &
               "GCE DYNAMICS: POTENTIOSTAT_CHAIN_LENGTH > 0 needs POTENTIOSTAT_TEMPERATURE or MOTION&MD TEMPERATURE > 0")
   END FUNCTION potentiostat_kt

! **************************************************************************************************
!> \brief state of the dynamics -> RESTART_* keywords of &GCE; the restart file the MD driver writes after
!>        this step carries them.
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE gce_dyn_write_state(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: m
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rv
      TYPE(gce_type), POINTER                            :: gce
      TYPE(section_vals_type), POINTER                   :: gce_section

      gce => qs_env%scf_control%gce
      gce_section => section_vals_get_subs_vals(qs_env%input, "DFT%SCF%GCE")
      CALL section_vals_val_set(gce_section, "RESTART_VALID", l_val=.TRUE.)
      CALL section_vals_val_set(gce_section, "RESTART_OPEN_PENDING", l_val=gce%dyn%open_pending)
      CALL section_vals_val_set(gce_section, "RESTART_N_E", r_val=gce%dyn%n_e)
      CALL section_vals_val_set(gce_section, "RESTART_N_INPUT", r_val=gce%dyn%n_input)
      CALL section_vals_val_set(gce_section, "RESTART_P_N", r_val=gce%dyn%p_n)
      CALL section_vals_val_set(gce_section, "RESTART_DU_LAST", r_val=gce%dyn%du_last)
      CALL section_vals_val_set(gce_section, "RESTART_CAPACITANCE", r_val=gce%dyn%capacitance)
      CALL section_vals_val_set(gce_section, "RESTART_ION_STEP", i_val=gce%dyn%ion_step_saved)
      CALL section_vals_val_set(gce_section, "RESTART_NEXT_LEN", i_val=gce%dyn%next_len)
      m = gce%dyn%chain_length
      IF (m > 0) THEN
         ALLOCATE (rv(m))
         rv(1:m) = gce%dyn%xi(1:m)
         CALL section_vals_val_set(gce_section, "RESTART_XI", r_vals_ptr=rv)
         ALLOCATE (rv(m))
         rv(1:m) = gce%dyn%p_xi(1:m)
         CALL section_vals_val_set(gce_section, "RESTART_P_XI", r_vals_ptr=rv)
      END IF
   END SUBROUTINE gce_dyn_write_state

END MODULE qs_scf_gce_dyn
