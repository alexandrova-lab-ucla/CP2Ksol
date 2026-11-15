!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the Quickstep SCF run.
!> \par History
!>      - Joost VandeVondele (02.2002)
!>           added code for: incremental (pab and gvg) update
!>                            initialisation (init_cube, l_info)
!>      - Joost VandeVondele (02.2002)
!>           called the poisson code of the classical part
!>           this takes into account the spherical cutoff and allows for
!>           isolated systems
!>      - Joost VandeVondele (02.2002)
!>           added multiple grid feature
!>           changed to spherical cutoff consistently (?)
!>           therefore removed the gradient correct functionals
!>      - updated with the new QS data structures (10.04.02,MK)
!>      - copy_matrix replaced by transfer_matrix (11.04.02,MK)
!>      - nrebuild_rho and nrebuild_gvg unified (12.04.02,MK)
!>      - set_mo_occupation for smearing of the MO occupation numbers
!>        (17.04.02,MK)
!>      - MO level shifting added (22.04.02,MK)
!>      - Usage of TYPE mo_set_p_type
!>      - Joost VandeVondele (05.2002)
!>            added cholesky based diagonalisation
!>      - 05.2002 added pao method [fawzi]
!>      - parallel FFT (JGH 22.05.2002)
!>      - 06.2002 moved KS matrix construction to qs_build_KS_matrix.F [fawzi]
!>      - started to include more LSD (01.2003,Joost VandeVondele)
!>      - 02.2003 scf_env [fawzi]
!>      - got rid of nrebuild (01.2004, Joost VandeVondele)
!>      - 10.2004 removed pao [fawzi]
!>      - 03.2006 large cleaning action [Joost VandeVondele]
!>      - High-spin ROKS added (05.04.06,MK)
!>      - Mandes (10.2013)
!>        intermediate energy communication with external communicator added
!>      - kpoints (08.2014, JGH)
!>      - unified k-point and gamma-point code (2014.11) [Ole Schuett]
!>      - added extra SCF loop for CDFT constraints (12.2015) [Nico Holmberg]
!> \author Matthias Krack (30.04.2001)
! **************************************************************************************************
MODULE qs_scf
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_get_info,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_add_default_logger,&
                                              cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_release,&
                                              cp_logger_type,&
                                              cp_rm_default_logger,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_p_file,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr,&
                                              cp_rm_iter_level
   USE cp_result_methods,               ONLY: get_results,&
                                              test_for_result
   USE cp_result_types,                 ONLY: cp_result_type
   USE ec_env_types,                    ONLY: energy_correction_type
   USE input_constants,                 ONLY: &
        broyden_type_1, broyden_type_1_explicit, broyden_type_1_explicit_ls, broyden_type_1_ls, &
        broyden_type_2, broyden_type_2_explicit, broyden_type_2_explicit_ls, broyden_type_2_ls, &
        cdft2ot, history_guess, ot2cdft, ot_precond_full_all, ot_precond_full_single, &
        ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse, &
        outer_scf_becke_constraint, outer_scf_hirshfeld_constraint, outer_scf_optimizer_broyden, &
        outer_scf_optimizer_newton_ls, tblite_scc_mixer_tblite
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE kpoint_io,                       ONLY: write_kpoints_restart
   USE kpoint_types,                    ONLY: kpoint_type
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathlib,                         ONLY: invert_matrix
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE preconditioner,                  ONLY: prepare_preconditioner,&
                                              restart_preconditioner
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE qs_block_davidson_types,         ONLY: block_davidson_deallocate
   USE qs_cdft_scf_utils,               ONLY: build_diagonal_jacobian,&
                                              create_tmp_logger,&
                                              initialize_inverse_jacobian,&
                                              prepare_jacobian_stencil,&
                                              print_inverse_jacobian,&
                                              restart_inverse_jacobian
   USE qs_cdft_types,                   ONLY: cdft_control_type
   USE qs_charge_mixing,                ONLY: charge_mixing_scc_error
   USE qs_charges_types,                ONLY: qs_charges_type
   USE qs_density_matrices,             ONLY: calculate_density_matrix
   USE qs_density_mixing_types,         ONLY: broyden_mixing_nr,&
                                              gspace_mixing_nr,&
                                              modified_broyden_mixing_nr
   USE qs_diis,                         ONLY: qs_diis_b_clear,&
                                              qs_diis_b_clear_kp,&
                                              qs_diis_b_create,&
                                              qs_diis_b_create_kp
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_methods,                   ONLY: evaluate_core_matrix_traces,&
                                              qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_did_change,&
                                              qs_ks_env_type
   USE qs_mo_io,                        ONLY: write_mo_set_to_restart
   USE qs_mo_methods,                   ONLY: make_basis_simple,&
                                              make_basis_sm
   USE qs_mo_occupation,                ONLY: gce_early_update_due,&
                                                gce_keep_mixer,&
                                                gce_track_update,&
                                                gce_update_nelectron,&
                                                set_mo_occupation
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              duplicate_mo_set,&
                                              get_mo_set,&
                                              mo_set_type,&
                                              reassign_allocated_mos
   USE qs_ot,                           ONLY: qs_ot_new_preconditioner
   USE qs_ot_scf,                       ONLY: ot_scf_init,&
                                              ot_scf_read_input
   USE qs_sccs_nlpb,                    ONLY: nlpb_new_scf_run,&
                                              nlpb_reset_scf_memory
   USE qs_outer_scf,                    ONLY: outer_loop_gradient,&
                                              outer_loop_optimize,&
                                              outer_loop_purge_history,&
                                              outer_loop_switch,&
                                              outer_loop_update_qs_env
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE pw_methods,                      ONLY: pw_integrate_function,&
                                              pw_transfer
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_gce_dyn,                  ONLY: gce_dyn_begin_step,&
                                              gce_dyn_end_step,&
                                              gce_dyn_note_solvent_off
   USE qs_scf_initialization,           ONLY: qs_scf_env_initialize
   USE qs_scf_loop_utils,               ONLY: qs_scf_check_inner_exit,&
                                              qs_scf_check_outer_exit,&
                                              qs_scf_density_mixing,&
                                              qs_scf_inner_finalize,&
                                              qs_scf_new_mos,&
                                              qs_scf_new_mos_kp,&
                                              qs_scf_rho_update,&
                                              qs_scf_set_loop_flags
   USE qs_scf_output,                   ONLY: gce_excess_electrons,&
                                              qs_scf_cdft_info,&
                                              qs_scf_cdft_initial_info,&
                                              qs_scf_gce_info,&
                                              qs_scf_loop_info,&
                                              qs_scf_loop_print,&
                                              qs_scf_outer_loop_info,&
                                              qs_scf_write_mos
   USE qs_scf_post_scf,                 ONLY: qs_scf_compute_properties
   USE qs_scf_types,                    ONLY: &
        block_davidson_diag_method_nr, block_krylov_diag_method_nr, filter_matrix_diag_method_nr, &
        general_diag_method_nr, ot_diag_method_nr, ot_method_nr, qs_scf_env_type, &
        smeagol_method_nr, special_diag_method_nr
   USE qs_wf_history_methods,           ONLY: wfi_purge_history,&
                                              wfi_update
   USE scf_control_types,               ONLY: gce_state_string,&
                                              gce_state_write,&
                                              gce_type,&
                                              scf_control_type
   USE smeagol_interface,               ONLY: run_smeagol_bulktrans,&
                                              run_smeagol_emtrans
   USE tblite_interface,                ONLY: tb_get_energy,&
                                              tb_native_scc_mixer_active,&
                                              tb_scf_mixer_error,&
                                              tb_update_charges
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf'
   LOGICAL, PRIVATE                     :: reuse_precond = .FALSE.
   LOGICAL, PRIVATE                     :: used_history = .FALSE.

   PUBLIC :: scf, scf_env_cleanup, scf_env_do_scf, cdft_scf, init_scf_loop

CONTAINS

! **************************************************************************************************
!> \brief perform an scf procedure in the given qs_env
!> \param qs_env the qs_environment where to perform the scf procedure
!> \param has_converged ...
!> \param total_scf_steps ...
!> \par History
!>      02.2003 introduced scf_env, moved real work to scf_env_do_scf [fawzi]
!> \author fawzi
!> \note
! **************************************************************************************************
   SUBROUTINE scf(qs_env, has_converged, total_scf_steps)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(OUT), OPTIONAL                     :: has_converged
      INTEGER, INTENT(OUT), OPTIONAL                     :: total_scf_steps

      INTEGER                                            :: ihistory, max_scf_tmp, tsteps
      LOGICAL                                            :: converged, outer_scf_loop, should_stop
      LOGICAL, SAVE                                      :: first_step_flag = .TRUE.
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient_history, variable_history
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, input, scf_section

      NULLIFY (scf_env)
      logger => cp_get_default_logger()
      CPASSERT(ASSOCIATED(qs_env))
      IF (PRESENT(has_converged)) THEN
         has_converged = .FALSE.
      END IF
      IF (PRESENT(total_scf_steps)) THEN
         total_scf_steps = 0
      END IF
      CALL get_qs_env(qs_env, scf_env=scf_env, input=input, &
                      dft_control=dft_control, scf_control=scf_control)
      IF (scf_control%max_scf > 0) THEN

         dft_section => section_vals_get_subs_vals(input, "DFT")
         scf_section => section_vals_get_subs_vals(dft_section, "SCF")

         IF (.NOT. ASSOCIATED(scf_env)) THEN
            CALL qs_scf_env_initialize(qs_env, scf_env)
            ! Moved here from qs_scf_env_initialize to be able to have more scf_env
            CALL set_qs_env(qs_env, scf_env=scf_env)
         ELSE
            CALL qs_scf_env_initialize(qs_env, scf_env)
         END IF

         IF ((scf_control%density_guess == history_guess) .AND. (first_step_flag)) THEN
            max_scf_tmp = scf_control%max_scf
            scf_control%max_scf = 1
            outer_scf_loop = scf_control%outer_scf%have_scf
            scf_control%outer_scf%have_scf = .FALSE.
         END IF

         ! CONVERGENCE_MODE / SOLVENT_ONLY_FOR_CHARGE, per ionic step: outer-update cap of this step; with
         ! SOLVENT_ONLY_FOR_CHARGE the solvent is on only on update steps (force evaluations 1, n, 2n, ... n = UPDATE_INTERVAL)
         IF (scf_control%gce%do_gce) THEN
            ! gce_dyn (UPDATE_SCHEME DYNAMICS): opening half step and drift of the electron number, N_e -> this SCF;
            ! a continuation run also restores ion_step here, before the solvent phase below is decided
            CALL gce_dyn_begin_step(qs_env)
            scf_control%gce%have_wf_on_step = .FALSE.
            ! the first evaluation of a restarted run repeats the last one of the previous run: same update-step phase
            ! (UPDATE_SCHEME CAPACITANCE; the dynamics restores ion_step itself and never sets replay_pending)
            IF (scf_control%gce%replay_pending .AND. scf_control%gce%ion_step > 0) &
               scf_control%gce%ion_step = scf_control%gce%ion_step - 1
            ! ENFORCE_TARGET: in a geometry optimisation every force evaluation is a full update step, so the last one
            ! (where the optimiser stops) is within FERMI_TOL of the target
            IF (scf_control%gce%enforce_target .AND. gce_in_geo_opt()) THEN
               scf_control%gce%update_interval = 1
               scf_control%gce%convergence_mode = 3
            END IF
            scf_control%gce%max_outer_now = scf_control%gce%max_outer
            IF (scf_control%gce%convergence_mode == 2 .OR. &
                (scf_control%gce%convergence_mode == 1 .AND. scf_control%gce%ion_step > 0)) &
               scf_control%gce%max_outer_now = MIN(1, scf_control%gce%max_outer)
            IF (scf_control%gce%dyn%active) scf_control%gce%max_outer_now = scf_control%gce%dyn%max_outer_cap
            IF (scf_control%gce%solvent_only_for_charge .AND. dft_control%do_sccs) THEN
               dft_control%sccs_control%force_off = .NOT. (scf_control%gce%ion_step == 0 .OR. &
                                                           MOD(scf_control%gce%ion_step + 1, scf_control%gce%update_interval) == 0)
            END IF
         END IF
         IF (.NOT. dft_control%qs_control%cdft) THEN
            CALL scf_env_do_scf(scf_env=scf_env, scf_control=scf_control, qs_env=qs_env, &
                                converged=converged, should_stop=should_stop, total_scf_steps=tsteps)
         ELSE
            ! Third SCF loop needed for CDFT with OT to properly restart OT inner loop
            CALL cdft_scf(qs_env=qs_env, should_stop=should_stop)
         END IF

         ! If SCF has not converged, then we should not start MP2
         IF (ASSOCIATED(qs_env%mp2_env)) qs_env%mp2_env%hf_fail = .NOT. converged

         ! cycle-mode updates happen inside scf_env_do_scf; reset the outer counter for the next call
         IF (qs_env%scf_control%gce%do_gce) THEN
            CALL gce_dyn_end_step(qs_env)   ! gce_dyn: closing half step, state to the restart section (ion_step not yet counted)
            qs_env%scf_control%gce%n_outer = 0
            qs_env%scf_control%gce%n_early = 0
            qs_env%scf_control%gce%ion_step = qs_env%scf_control%gce%ion_step + 1   ! counts force evaluations (UPDATE_INTERVAL)
            ! one machine-readable line per ionic step, whether or not N_e was updated, so a
            ! constant-potential trajectory has a potential trace. W_f is the converged value of this ionic step.
            BLOCK
               INTEGER                                            :: gce_unit
               REAL(KIND=dp)                                      :: carry_dn, n_exc, u_now, u_tgt, wf_now
               TYPE(qs_energy_type), POINTER                      :: gce_energy
               TYPE(section_vals_type), POINTER                   :: gce_section
               gce_unit = cp_logger_get_default_io_unit()
               ! W_f from the converged Fermi level (energy%efermi is set on the Gamma and the k-point path;
               ! prev_workfunction is only set by the Gamma occupation routines and is meaningless with k-points)
               CALL get_qs_env(qs_env, energy=gce_energy)
               wf_now = qs_env%scf_control%gce%ref_esp - gce_energy%efermi
               ! SOLVENT_ONLY_FOR_CHARGE: on an update step the converged state is the solvent-off one; the potential the
               ! step was steered to is the solvent-on work function measured before the switch
               IF (qs_env%scf_control%gce%have_wf_on_step) wf_now = qs_env%scf_control%gce%wf_on_step
               IF (gce_unit > 0) THEN
                  u_now = evolt*(wf_now - qs_env%scf_control%gce%phi_she)
                  u_tgt = evolt*(qs_env%scf_control%gce%target_workfunction - qs_env%scf_control%gce%phi_she)
                  WRITE (UNIT=gce_unit, FMT="(T3,A,I7,A,F14.6,A,F10.4,A,F9.4,A,F9.4,A,F9.4,A)") &
                     "GCE| trace  step ", qs_env%scf_control%gce%ion_step, "  N_e ", &
                     qs_env%scf_control%gce%n_target, " e  W_f ", &
                     evolt*wf_now, " eV  U ", u_now, &
                     " V vs SHE  U_target ", u_tgt, "  dU ", u_now - u_tgt, " V"
               END IF
               ! what the MD driver sees for this step: its potential energy is energy%total, and the grand potential
               ! is Omega = E + WF_target n_excess (mu_e = -WF_target); N_e is read before CARRY_GAIN changes it
               qs_env%scf_control%gce%n_step = qs_env%scf_control%gce%n_target
               qs_env%scf_control%gce%e_step = gce_energy%total
               qs_env%scf_control%gce%wf_step = wf_now
               IF (gce_unit > 0) THEN
                  n_exc = gce_excess_electrons(qs_env)
                  WRITE (UNIT=gce_unit, FMT="(T3,A,I7,A,F22.12,A,F12.6,A,F22.12,A,F10.4,A)") &
                     "GCE| omega  step ", qs_env%scf_control%gce%ion_step, "  E ", gce_energy%total, &
                     " Ha  n_exc ", n_exc, " e  Omega ", &
                     gce_energy%total + qs_env%scf_control%gce%target_workfunction*n_exc, &
                     " Ha  WF_target ", evolt*qs_env%scf_control%gce%target_workfunction, " eV"
               END IF
               ! CARRY_GAIN: feed-forward of N_e to the next ionic step
               IF (qs_env%scf_control%gce%carry_gain > 0.0_dp .AND. qs_env%scf_control%gce%scheme == 2 .AND. &
                   (.NOT. qs_env%scf_control%gce%update_every_step) .AND. converged .AND. &
                   (qs_env%scf_control%gce%have_wf_on_step .OR. (.NOT. qs_env%scf_control%gce%solvent_only_for_charge))) THEN
                  carry_dn = qs_env%scf_control%gce%carry_gain*qs_env%scf_control%gce%capacitance* &
                             (wf_now - qs_env%scf_control%gce%target_workfunction)
                  carry_dn = MAX(-qs_env%scf_control%gce%ne_step_max, MIN(qs_env%scf_control%gce%ne_step_max, carry_dn))
                  qs_env%scf_control%gce%n_target = qs_env%scf_control%gce%n_target + carry_dn
                  qs_env%scf_control%gce%steps_since_update = 0
                  IF (gce_unit > 0) WRITE (UNIT=gce_unit, FMT="(T3,A,F12.6,A,F14.6,A)") &
                     "GCE| carry  dN_e ", carry_dn, " e  N_e for the next step ", qs_env%scf_control%gce%n_target, " e"
               END IF
               qs_env%scf_control%gce%replay_pending = .FALSE.
               ! the state after this ionic step goes into the RESTART_* keywords the restart file carries
               gce_section => section_vals_get_subs_vals(scf_section, "GCE")
               CALL gce_state_write(qs_env%scf_control%gce, gce_section)
            END BLOCK
         END IF

         ! Add the converged outer_scf SCF gradient(s)/variable(s) to history
         IF (scf_control%outer_scf%have_scf) THEN
            ihistory = scf_env%outer_scf%iter_count
            CALL get_qs_env(qs_env, gradient_history=gradient_history, &
                            variable_history=variable_history)
            ! We only store the latest two values
            gradient_history(:, 1) = gradient_history(:, 2)
            gradient_history(:, 2) = scf_env%outer_scf%gradient(:, ihistory)
            variable_history(:, 1) = variable_history(:, 2)
            variable_history(:, 2) = scf_env%outer_scf%variables(:, ihistory)
            ! Reset flag
            IF (used_history) used_history = .FALSE.
            ! Update a counter and check if the Jacobian should be deallocated
            IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian)) THEN
               scf_control%outer_scf%cdft_opt_control%ijacobian(2) = scf_control%outer_scf%cdft_opt_control%ijacobian(2) + 1
               IF (scf_control%outer_scf%cdft_opt_control%ijacobian(2) >= &
                   scf_control%outer_scf%cdft_opt_control%jacobian_freq(2) .AND. &
                   scf_control%outer_scf%cdft_opt_control%jacobian_freq(2) > 0) &
                  scf_env%outer_scf%deallocate_jacobian = .TRUE.
            END IF
         END IF
         !   *** add the converged wavefunction to the wavefunction history
         IF ((ASSOCIATED(qs_env%wf_history)) .AND. &
             ((scf_control%density_guess /= history_guess) .OR. &
              (.NOT. first_step_flag))) THEN
            IF (.NOT. dft_control%qs_control%cdft) THEN
               CALL wfi_update(qs_env%wf_history, qs_env=qs_env, dt=1.0_dp)
            ELSE
               IF (dft_control%qs_control%cdft_control%should_purge) THEN
                  CALL wfi_purge_history(qs_env)
                  CALL outer_loop_purge_history(qs_env)
                  dft_control%qs_control%cdft_control%should_purge = .FALSE.
               ELSE
                  CALL wfi_update(qs_env%wf_history, qs_env=qs_env, dt=1.0_dp)
               END IF
            END IF
         ELSE IF ((scf_control%density_guess == history_guess) .AND. &
                  (first_step_flag)) THEN
            scf_control%max_scf = max_scf_tmp
            scf_control%outer_scf%have_scf = outer_scf_loop
            first_step_flag = .FALSE.
         END IF

         ! *** compute properties that depend on the converged wavefunction
         IF (.NOT. (should_stop)) CALL qs_scf_compute_properties(qs_env)

         ! *** SMEAGOL interface ***
         IF (.NOT. (should_stop)) THEN
            ! compute properties that depend on the converged wavefunction ..
            CALL run_smeagol_emtrans(qs_env, last=.TRUE., iter=0)
            ! .. or save matrices related to bulk leads
            CALL run_smeagol_bulktrans(qs_env)
         END IF

         ! *** cleanup
         CALL scf_env_cleanup(scf_env)
         IF (dft_control%qs_control%cdft) &
            CALL cdft_control_cleanup(dft_control%qs_control%cdft_control)

         IF (PRESENT(has_converged)) THEN
            has_converged = converged
         END IF
         IF (PRESENT(total_scf_steps)) THEN
            total_scf_steps = tsteps
         END IF

      END IF

   END SUBROUTINE scf

! **************************************************************************************************
!> \brief perform an scf loop
!> \param scf_env the scf_env where to perform the scf procedure
!> \param scf_control ...
!> \param qs_env the qs_env, the scf_env lives in
!> \param converged will be true / false if converged is reached
!> \param should_stop ...
!> \param total_scf_steps ...
!> \par History
!>      long history, see cvs and qs_scf module history
!>      02.2003 introduced scf_env [fawzi]
!>      09.2005 Frozen density approximation [TdK]
!>      06.2007 Check for SCF iteration count early [jgh]
!>      10.2019 switch_surf_dip [SGh]
!> \author Matthias Krack
!> \note
! **************************************************************************************************
   SUBROUTINE scf_env_do_scf(scf_env, scf_control, qs_env, converged, should_stop, total_scf_steps)

      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(OUT)                               :: converged, should_stop
      INTEGER, INTENT(OUT)                               :: total_scf_steps

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'scf_env_do_scf'

      CHARACTER(LEN=default_string_length)               :: description, name
      INTEGER                                            :: ext_master_id, handle, handle2, i_tmp, &
                                                            ic, ispin, iter_count, output_unit, &
                                                            scf_energy_message_tag, total_steps
      LOGICAL :: density_full_step, diis_step, do_kpoints, energy_only, exit_inner_loop, &
         exit_outer_loop, inner_loop_converged, internal_tblite_density_full_step, &
         internal_tblite_mixer, just_energy, outer_loop_converged, tblite_native_mixer
      LOGICAL                                            :: gce_dbg, gce_resp_done
      REAL(KIND=dp)                                      :: gce_n_before
      LOGICAL, SAVE                                      :: gce_spin_warned = .FALSE.
      CHARACTER(LEN=8)                                   :: gce_dbg_env
      REAL(KIND=dp)                                      :: t1, t2
      REAL(KIND=dp), DIMENSION(3)                        :: res_val_3
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos, mos_last_converged
      TYPE(mp_comm_type)                                 :: external_comm
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_charges_type), POINTER                     :: qs_charges
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: dft_section, input, scf_section

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, rho, energy, &
               logger, qs_charges, ks_env, mos, atomic_kind_set, qs_kind_set, &
               particle_set, dft_section, input, &
               scf_section, para_env, results, kpoints, pw_env, rho_ao_kp, mos_last_converged)

      CPASSERT(ASSOCIATED(scf_env))
      CPASSERT(ASSOCIATED(qs_env))

      logger => cp_get_default_logger()
      t1 = m_walltime()

      CALL get_qs_env(qs_env=qs_env, &
                      energy=energy, &
                      particle_set=particle_set, &
                      qs_charges=qs_charges, &
                      ks_env=ks_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      rho=rho, &
                      mos=mos, &
                      input=input, &
                      dft_control=dft_control, &
                      do_kpoints=do_kpoints, &
                      kpoints=kpoints, &
                      results=results, &
                      pw_env=pw_env, &
                      para_env=para_env)
      tblite_native_mixer = dft_control%qs_control%xtb_control%do_tblite .AND. &
                            scf_env%method /= ot_method_nr .AND. &
                            tb_native_scc_mixer_active(dft_control)
      internal_tblite_mixer = (dft_control%qs_control%dftb .AND. &
                               dft_control%qs_control%dftb_control%tblite_scc_mixer == tblite_scc_mixer_tblite) .OR. &
                              (dft_control%qs_control%xtb .AND. &
                               .NOT. dft_control%qs_control%xtb_control%do_tblite .AND. &
                               dft_control%qs_control%xtb_control%tblite_scc_mixer == tblite_scc_mixer_tblite)
      internal_tblite_density_full_step = dft_control%qs_control%xtb .AND. &
                                          .NOT. dft_control%qs_control%xtb_control%do_tblite .AND. &
                                          dft_control%qs_control%xtb_control%tblite_scc_mixer == tblite_scc_mixer_tblite

      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)

      dft_section => section_vals_get_subs_vals(input, "DFT")
      scf_section => section_vals_get_subs_vals(dft_section, "SCF")

      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".scfLog")

      gce_dbg_env = ' '
      CALL GET_ENVIRONMENT_VARIABLE("GCE_DBG", gce_dbg_env)
      gce_dbg = scf_control%gce%do_gce .AND. LEN_TRIM(gce_dbg_env) > 0
      scf_control%gce%ncall_restore = 0
      ! UKS: the electron number is shared by both spins and the Fermi level is common (set_mo_occupation_3). With
      ! FIXED_MAGNETIC_MOMENT every spin gets its own occupation without the GCE state and N_e would never be steered.
      IF (scf_control%gce%do_gce .AND. scf_control%gce%constant_potential_requested .AND. dft_control%nspins == 2 .AND. &
          ASSOCIATED(scf_control%smear)) THEN
         IF (scf_control%smear%do_smear .AND. (scf_control%smear%fixed_mag_mom >= 0.0_dp .OR. &
                                               scf_control%smear%fixed_mag_mom == -1.0_dp)) &
            CPABORT("GCE with UKS needs one Fermi level for both spins: do not set FIXED_MAGNETIC_MOMENT")
      END IF
      ! one Fermi level for both spins: the multiplicity only sets the initial occupation; the SCF can end in another spin state
      IF (.NOT. gce_spin_warned) THEN
         IF (scf_control%gce%do_gce .AND. scf_control%gce%scheme == 2 .AND. dft_control%nspins == 2 .AND. &
             dft_control%multiplicity > 1) THEN
            CALL cp_warn(__LOCATION__, &
                         "GCE: with UKS, CAPACITANCE and DYNAMICS share one Fermi level between the spins, "// &
                         "so MULTIPLICITY only sets the initial occupation")
            gce_spin_warned = .TRUE.
         END IF
      END IF
      ! PIN_MIXER_G0 sets G = 0 of the mixer input to N_e/V, which is the whole density only for one spin channel
      IF (scf_control%gce%do_gce .AND. scf_control%gce%pin_mixer_g0 .AND. dft_control%nspins /= 1) THEN
         CPWARN("GCE: PIN_MIXER_G0 is implemented for spin-restricted runs only and is switched off for UKS")
         scf_control%gce%pin_mixer_g0 = .FALSE.
      END IF
      ! With the G = 0 charge of the mixer input not pinned (PIN_MIXER_G0 off), a new N_e reaches the density
      ! only after the Broyden buffer has wrapped (the work function stays frozen for NBROYDEN + 1 steps). A
      ! convergence check inside that window sees a frozen work function, so the minimum number of steps after
      ! an update must cover it.
      IF (scf_control%gce%do_gce .AND. scf_control%gce%scheme == 2 .AND. (.NOT. scf_control%gce%update_every_step) .AND. &
          (.NOT. scf_control%gce%pin_mixer_g0) .AND. ASSOCIATED(scf_env%mixing_store)) THEN
         IF (scf_env%mixing_method == broyden_mixing_nr .OR. scf_env%mixing_method == modified_broyden_mixing_nr) THEN
            IF (scf_control%gce%min_steps < scf_env%mixing_store%nbuffer + 2) THEN
               CPWARN("GCE: MIN_STEPS_AFTER_UPDATE raised to NBROYDEN + 2 (Broyden buffer wrap)")
               scf_control%gce%min_steps = scf_env%mixing_store%nbuffer + 2
            END IF
         END IF
      END IF
      ! The k-point path supports the CAPACITANCE scheme (fixed fractional charge and constant potential;
      ! kpoint_set_mo_occupation). The mixing scheme (UPDATE_SCHEME /= CAPACITANCE) has no k-point
      ! implementation, so a requested steering there would be silently ignored: abort. Neutral blocks pass.
      IF (do_kpoints .AND. scf_control%gce%do_gce .AND. scf_control%gce%constant_potential_requested .AND. &
          scf_control%gce%scheme /= 2) THEN
         CPABORT("GCE with k-points needs UPDATE_SCHEME CAPACITANCE")
      END IF

      ! a new SCF run starts a new tolerance history of the NLPB solver (running minimum of the entry residual)
      IF (dft_control%do_sccs) CALL nlpb_reset_scf_memory()
      IF (dft_control%do_sccs) CALL nlpb_new_scf_run()   ! the potential history of the solver restarts with every SCF run

      IF (scf_control%gce%do_gce .AND. output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,78('-'))")
         WRITE (UNIT=output_unit, FMT="(T31,A)") &
            "GRAND-CANONICAL SCF"
         WRITE (UNIT=output_unit, FMT="(T20,A,F12.6,A)") &
            "Target work function (TWF):", &
            evolt*scf_control%gce%target_workfunction, " eV"
         WRITE (UNIT=output_unit, FMT="(T2,78('-'))")
         IF (scf_control%gce%state_report) THEN
            WRITE (UNIT=output_unit, FMT="(T3,A,A)") "GCE| state read from the restart: ", TRIM(gce_state_string(scf_control%gce))
            scf_control%gce%state_report = .FALSE.
         END IF
      END IF

      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
         "SCF WAVEFUNCTION OPTIMIZATION"

      ! when switch_surf_dip is switched on, indicate storing mos from the last converged step
      IF (dft_control%switch_surf_dip) THEN
         CALL get_qs_env(qs_env, mos_last_converged=mos_last_converged)
         DO ispin = 1, dft_control%nspins
            CALL reassign_allocated_mos(mos(ispin), mos_last_converged(ispin))
         END DO
         IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
            "COPIED mos_last_converged ---> mos"
      END IF

      IF ((output_unit > 0) .AND. (.NOT. scf_control%use_ot)) THEN
         WRITE (UNIT=output_unit, &
                FMT="(/,T3,A,T12,A,T31,A,T39,A,T59,A,T75,A,/,T3,A)") &
            "Step", "Update method", "Time", "Convergence", "Total energy", "Change", &
            REPEAT("-", 78)
      END IF
      CALL cp_add_iter_level(logger%iter_info, "QS_SCF")

      ! check for external communicator and if the intermediate energy should be sent
      res_val_3(:) = -1.0_dp
      description = "[EXT_SCF_ENER_COMM]"
      IF (test_for_result(results, description=description)) THEN
         CALL get_results(results, description=description, &
                          values=res_val_3, n_entries=i_tmp)
         CPASSERT(i_tmp == 3)
         IF (ALL(res_val_3(:) <= 0.0)) &
            CALL cp_abort(__LOCATION__, &
                          " Trying to access result ("//TRIM(description)// &
                          ") which is not correctly stored.")
         CALL external_comm%set_handle(NINT(res_val_3(1)))
      END IF
      ext_master_id = NINT(res_val_3(2))
      scf_energy_message_tag = NINT(res_val_3(3))

      ! *** outer loop of the scf, can treat other variables,
      ! *** such as lagrangian multipliers
      scf_env%outer_scf%iter_count = 0
      iter_count = 0
      total_steps = 0
      energy%tot_old = 0.0_dp

      scf_outer_loop: DO

         CALL init_scf_loop(scf_env=scf_env, qs_env=qs_env, &
                            scf_section=scf_section)

         CALL qs_scf_set_loop_flags(scf_env, diis_step, &
                                    energy_only, just_energy, exit_inner_loop)
         scf_control%gce%iter_base = 0   ! MAX_SCF counts from here, and again from every N_e update

         ! decide whether to switch off dipole correction for convergence purposes
         dft_control%surf_dip_correct_switch = dft_control%correct_surf_dip
         IF ((dft_control%correct_surf_dip) .AND. (scf_control%outer_scf%have_scf) .AND. &
             (scf_env%outer_scf%iter_count > FLOOR(scf_control%outer_scf%max_scf/2.0_dp))) THEN
            IF (dft_control%switch_surf_dip) THEN
               dft_control%surf_dip_correct_switch = .FALSE.
               IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
                  "SURFACE DIPOLE CORRECTION switched off"
            END IF
         END IF

         scf_loop: DO

            CALL timeset(routineN//"_inner_loop", handle2)

            IF (.NOT. just_energy) scf_env%iter_count = scf_env%iter_count + 1
            iter_count = iter_count + 1
            ! let the GCE electron-number update see how converged the density is
            IF (scf_control%gce%do_gce) THEN
               scf_control%gce%scf_delta = scf_env%iter_delta
               scf_control%gce%steps_since_update = scf_control%gce%steps_since_update + 1
            END IF
            CALL cp_iterate(logger%iter_info, last=.FALSE., iter_nr=iter_count)

            IF (output_unit > 0) CALL m_flush(output_unit)

            total_steps = total_steps + 1
            just_energy = energy_only

            CALL qs_ks_update_qs_env(qs_env, just_energy=just_energy, &
                                     calculate_forces=.FALSE.)
            IF (gce_dbg) CALL gce_dbg_g0(rho, para_env, output_unit, "after_ks_build", scf_env)

            ! print 'heavy weight' or relatively expensive quantities
            CALL qs_scf_loop_print(qs_env, scf_env, para_env)

            IF (do_kpoints) THEN
               ! kpoints
               IF (dft_control%hairy_probes .EQV. .TRUE.) THEN
                  scf_control%smear%do_smear = .FALSE.
                  CALL qs_scf_new_mos_kp(qs_env, scf_env, scf_control, diis_step, dft_control%probe)
               ELSE
                  CALL qs_scf_new_mos_kp(qs_env, scf_env, scf_control, diis_step)
               END IF
            ELSE
               ! Gamma points only
               IF (dft_control%hairy_probes .EQV. .TRUE.) THEN
                  scf_control%smear%do_smear = .FALSE.
                  CALL qs_scf_new_mos(qs_env, scf_env, scf_control, scf_section, diis_step, energy_only, &
                                      dft_control%probe)
               ELSE
                  CALL qs_scf_new_mos(qs_env, scf_env, scf_control, scf_section, diis_step, energy_only)
               END IF
            END IF

            ! the boundary step after an N_e update was a plain (Kerker) step; the Broyden history is valid again
            IF (scf_control%gce%ncall_restore > 0 .AND. ASSOCIATED(scf_env%mixing_store)) THEN
               scf_env%mixing_store%ncall = scf_control%gce%ncall_restore
               scf_control%gce%ncall_restore = 0
            END IF
            IF (gce_dbg) CALL gce_dbg_g0(rho, para_env, output_unit, "after_new_mos", scf_env)
            ! Print requested MO information (can be computationally expensive with OT)
            CALL qs_scf_write_mos(qs_env, scf_env, final_mos=.FALSE.)

            IF (dft_control%qs_control%xtb_control%do_tblite) THEN
               IF (scf_env%method == ot_method_nr) THEN
                  CALL tb_update_charges(qs_env, dft_control, qs_env%tb_tblite, .FALSE., .TRUE.)
                  CALL evaluate_core_matrix_traces(qs_env)
               ELSE
                  CPASSERT(scf_env%mixing_method > 0)
                  CALL tb_update_charges(qs_env, dft_control, qs_env%tb_tblite, .FALSE., .FALSE.)
                  CALL evaluate_core_matrix_traces(qs_env, rho_ao_ext=scf_env%p_mix_new)
               END IF
               CALL tb_get_energy(qs_env, qs_env%tb_tblite, energy)
            END IF

            density_full_step = diis_step .OR. tblite_native_mixer .OR. internal_tblite_density_full_step
            CALL qs_scf_density_mixing(scf_env, rho, para_env, density_full_step)
            ! GCE CAPACITANCE: the total charge is set by the occupations, not by the mixer:
            ! pin the mixed density's G = 0 coefficient to N_e(target)/V (spin-restricted only)
            IF (scf_control%gce%do_gce .AND. scf_control%gce%scheme == 2 .AND. scf_env%mixing_method > 0 .AND. &
                scf_control%gce%n_target > 0.0_dp .AND. dft_control%nspins == 1) THEN
               CALL gce_pin_g0(rho, scf_env, scf_control)
               IF (gce_dbg) CALL gce_dbg_g0(rho, para_env, output_unit, "after_pin", scf_env)
            END IF
            IF (dft_control%qs_control%xtb_control%do_tblite .AND. &
                .NOT. (dft_control%qs_control%do_ls_scf .OR. scf_control%use_ot)) THEN
               scf_env%iter_delta = MAX(scf_env%iter_delta, &
                                        tb_scf_mixer_error(dft_control, qs_env%tb_tblite, &
                                                           scf_control%eps_scf))
            END IF
            IF (dft_control%qs_control%dftb .OR. &
                (dft_control%qs_control%xtb .AND. .NOT. dft_control%qs_control%xtb_control%do_tblite)) THEN
               scf_env%iter_delta = MAX(scf_env%iter_delta, &
                                        charge_mixing_scc_error(scf_env%mixing_store, scf_control%eps_scf))
            END IF
            IF (tblite_native_mixer) THEN
               scf_env%iter_param = dft_control%qs_control%xtb_control%tblite_mixer_damping
               scf_env%iter_method = "TBLite/Diag"
            ELSEIF (internal_tblite_mixer) THEN
               scf_env%iter_method = "TBLite/Diag"
               IF (dft_control%qs_control%dftb) THEN
                  scf_env%iter_param = dft_control%qs_control%dftb_control%tblite_mixer_damping
               ELSE
                  scf_env%iter_param = dft_control%qs_control%xtb_control%tblite_mixer_damping
               END IF
            END IF

            t2 = m_walltime()

            CALL qs_scf_loop_info(scf_env, output_unit, just_energy, t1, t2, energy)

            IF (scf_control%gce%do_gce) THEN
               CALL qs_scf_gce_info(output_unit, qs_env, just_energy)
            END IF

            IF (.NOT. just_energy) energy%tot_old = energy%total

            ! check for external communicator and if the intermediate energy should be sent
            IF (scf_energy_message_tag > 0) THEN
               CALL external_comm%send(energy%total, ext_master_id, scf_energy_message_tag)
            END IF

            CALL qs_scf_check_inner_exit(qs_env, scf_env, scf_control, should_stop, just_energy, &
                                         exit_inner_loop, inner_loop_converged, output_unit)

            ! GCE CAPACITANCE in cycle mode (UPDATE_EVERY SCF_CYCLE): after a converged
            ! inner SCF at fixed N_e, update N_e once from the converged work function and go on
            IF (exit_inner_loop .AND. inner_loop_converged .AND. (.NOT. just_energy) .AND. &
                scf_control%gce%do_gce .AND. scf_control%gce%scheme == 2 .AND. &
                (.NOT. scf_control%gce%update_every_step)) THEN
               ! UPDATE_INTERVAL n: update steps are force evaluations 1, n, 2n, ... = ion_step == 0 .OR. MOD(ion_step + 1, n) == 0
               IF (.NOT. gce_solvent_is_off(dft_control)) CALL gce_check_target_reached(scf_control%gce, output_unit)
               IF (ABS(scf_control%gce%prev_workfunction - scf_control%gce%target_workfunction) > &
                   scf_control%gce%fermi_tol*scf_control%gce%final_tol_frac .AND. &
                   scf_control%gce%n_outer < scf_control%gce%max_outer_now .AND. &
                   (scf_control%gce%ion_step == 0 .OR. MOD(scf_control%gce%ion_step + 1, scf_control%gce%update_interval) == 0) &
                   .AND. (.NOT. gce_solvent_is_off(dft_control))) THEN
                  scf_control%gce%n_outer = scf_control%gce%n_outer + 1
                  CALL gce_track_update(scf_control%gce)
                  gce_n_before = scf_control%gce%n_target
                  CALL gce_update_nelectron(scf_control%gce, scf_control%gce%prev_workfunction)
                  ! the mixed density must carry the new N_e before the work function is read again:
                  ! restart the mixer history and demand MIN_STEPS_AFTER_UPDATE inner steps
                  scf_control%gce%steps_since_update = 0
                  scf_control%gce%iter_base = scf_env%iter_count
                  CALL gce_mixer_after_update(scf_env, scf_control)
                  ! the new N_e changes the density abruptly, so the NLPB entry residuals of the next steps are
                  ! not bounded by the running minimum of the previous N_e
                  IF (dft_control%do_sccs) CALL nlpb_reset_scf_memory()
                  gce_resp_done = .FALSE.
                  IF (scf_control%gce%response_guess) &
                     CALL gce_response_guess(rho, scf_env, scf_control, gce_n_before, gce_resp_done, output_unit)
                  IF (scf_control%gce%pin_mixer_g0 .AND. .NOT. gce_resp_done) CALL gce_pin_g0(rho, scf_env, scf_control)
                  IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(T3,A,I3,A,F14.6,A,F10.4,A,F10.4,A,F10.4,A)") &
                     "GCE| outer update ", scf_control%gce%n_outer, ": N_e -> ", scf_control%gce%n_target, &
                     " e  C ", scf_control%gce%capacitance/evolt, " e/V  WF-TWF was ", &
                     (scf_control%gce%prev_workfunction - scf_control%gce%target_workfunction)*evolt, &
                     " eV  U was ", (scf_control%gce%prev_workfunction - scf_control%gce%phi_she)*evolt, " V vs SHE"
                  exit_inner_loop = .FALSE.
                  inner_loop_converged = .FALSE.
               ELSE IF (scf_control%gce%solvent_only_for_charge .AND. dft_control%do_sccs .AND. &
                        (.NOT. gce_solvent_is_off(dft_control))) THEN
                  ! SOLVENT_ONLY_FOR_CHARGE: N_e is set for this ionic step -> switch the solvent off and
                  ! re-converge the SCF at fixed N_e without it (mixer restart as after an N_e update)
                  scf_control%gce%wf_on_step = scf_control%gce%prev_workfunction   ! the potential this step is steered to
                  scf_control%gce%have_wf_on_step = .TRUE.
                  CALL gce_dyn_note_solvent_off(qs_env)   ! gce_dyn: U of the solvent-on measurement
                  dft_control%sccs_control%force_off = .TRUE.
                  scf_control%gce%steps_since_update = 0
                  scf_control%gce%iter_base = scf_env%iter_count
                  IF (ASSOCIATED(scf_env%mixing_store)) scf_env%mixing_store%ncall = 0
                  IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(T3,A,F14.6,A)") &
                     "GCE| SOLVENT_ONLY_FOR_CHARGE: N_e set to ", scf_control%gce%n_target, " e, solvent off, re-converging SCF"
                  exit_inner_loop = .FALSE.
                  inner_loop_converged = .FALSE.
               END IF
            END IF

            ! early N_e update (cycle mode) before the inner SCF has converged; off unless EARLY_UPDATE_DELTA > 0
            IF ((.NOT. exit_inner_loop) .AND. (.NOT. just_energy) .AND. scf_control%gce%do_gce .AND. &
                scf_control%gce%scheme == 2 .AND. (.NOT. scf_control%gce%update_every_step) .AND. &
                scf_control%gce%early_update_delta > 0.0_dp) THEN
               IF (gce_early_update_due(scf_control%gce, scf_env%iter_delta) .AND. &
                   (.NOT. gce_solvent_is_off(dft_control))) THEN
                  scf_control%gce%n_early = scf_control%gce%n_early + 1
                  ! EARLY_UPDATE_SECANT: the first early update only records its (N_e, WF) pair, later early updates
                  ! refresh the secant; otherwise the capacitance stays frozen in every early update
                  scf_control%gce%freeze_cap = (scf_control%gce%n_early == 1) .OR. (.NOT. scf_control%gce%early_secant)
                  gce_n_before = scf_control%gce%n_target
                  CALL gce_update_nelectron(scf_control%gce, scf_control%gce%prev_workfunction)
                  scf_control%gce%freeze_cap = .FALSE.
                  scf_control%gce%steps_since_update = 0
                  scf_control%gce%iter_base = scf_env%iter_count
                  CALL gce_mixer_after_update(scf_env, scf_control)
                  IF (dft_control%do_sccs) CALL nlpb_reset_scf_memory()
                  gce_resp_done = .FALSE.
                  IF (scf_control%gce%response_guess) &
                     CALL gce_response_guess(rho, scf_env, scf_control, gce_n_before, gce_resp_done, output_unit)
                  IF (scf_control%gce%pin_mixer_g0 .AND. .NOT. gce_resp_done) CALL gce_pin_g0(rho, scf_env, scf_control)
                  IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(T3,A,I3,A,F14.6,A,F10.4,A,F10.4,A,F10.4,A)") &
                     "GCE| early update ", scf_control%gce%n_early, ": N_e -> ", scf_control%gce%n_target, &
                     " e  C ", scf_control%gce%capacitance/evolt, " e/V  WF-TWF was ", &
                     (scf_control%gce%prev_workfunction - scf_control%gce%target_workfunction)*evolt, &
                     " eV  U was ", (scf_control%gce%prev_workfunction - scf_control%gce%phi_she)*evolt, " V vs SHE"
               END IF
            END IF

            ! In case we decide to exit we perform few more check to see if this one
            ! is really the last SCF step
            IF (exit_inner_loop) THEN

               CALL qs_scf_inner_finalize(scf_env, qs_env, density_full_step, output_unit)

               CALL qs_scf_check_outer_exit(qs_env, scf_env, scf_control, should_stop, &
                                            outer_loop_converged, exit_outer_loop)

               ! Let's tag the last SCF cycle so we can print informations only of the last step
               IF (exit_outer_loop) CALL cp_iterate(logger%iter_info, last=.TRUE., iter_nr=iter_count)

            END IF

            IF (do_kpoints) THEN
               CALL write_kpoints_restart(rho_ao_kp, kpoints, scf_env, dft_section, particle_set, qs_kind_set)
            ELSE
               ! Write wavefunction restart file
               IF (scf_env%method == ot_method_nr) THEN
                  ! With OT: provide the Kohn-Sham matrix for the calculation of the MO eigenvalues
                  CALL get_ks_env(ks_env=ks_env, matrix_ks=matrix_ks)
                  CALL write_mo_set_to_restart(mos, particle_set, dft_section=dft_section, qs_kind_set=qs_kind_set, &
                                               matrix_ks=matrix_ks)
               ELSE
                  CALL write_mo_set_to_restart(mos, particle_set, dft_section=dft_section, qs_kind_set=qs_kind_set)
               END IF

            END IF

            ! Exit if we have finished with the SCF inner loop
            IF (exit_inner_loop) THEN
               CALL timestop(handle2)
               EXIT scf_loop
            END IF

            IF (.NOT. BTEST(cp_print_key_should_output(logger%iter_info, &
                                                       scf_section, "PRINT%ITERATION_INFO/TIME_CUMUL"), cp_p_file)) &
               t1 = m_walltime()

            ! mixing methods have the new density matrix in p_mix_new
            IF (scf_env%mixing_method > 0) THEN
               DO ic = 1, SIZE(rho_ao_kp, 2)
                  DO ispin = 1, dft_control%nspins
                     CALL dbcsr_get_info(rho_ao_kp(ispin, ic)%matrix, name=name) ! keep the name
                     CALL dbcsr_copy(rho_ao_kp(ispin, ic)%matrix, scf_env%p_mix_new(ispin, ic)%matrix, name=name)
                  END DO
               END DO
            END IF

            CALL qs_scf_rho_update(rho, qs_env, scf_env, ks_env, &
                                   mix_rho=scf_env%mixing_method >= gspace_mixing_nr)

            CALL timestop(handle2)

         END DO scf_loop

         IF (.NOT. scf_control%outer_scf%have_scf) EXIT scf_outer_loop

         ! In case we use the OUTER SCF loop let's print some info..
         CALL qs_scf_outer_loop_info(output_unit, scf_control, scf_env, &
                                     energy, total_steps, should_stop, outer_loop_converged)

         ! Save MOs to converged MOs if outer_loop_converged and surf_dip_correct_switch is true
         IF (exit_outer_loop) THEN
            IF ((dft_control%switch_surf_dip) .AND. (outer_loop_converged) .AND. &
                (dft_control%surf_dip_correct_switch)) THEN
               DO ispin = 1, dft_control%nspins
                  CALL reassign_allocated_mos(mos_last_converged(ispin), mos(ispin))
               END DO
               IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
                  "COPIED mos ---> mos_last_converged"
            END IF
         END IF

         IF (exit_outer_loop) EXIT scf_outer_loop

         !
         CALL outer_loop_optimize(scf_env, scf_control)
         CALL outer_loop_update_qs_env(qs_env, scf_env)
         CALL qs_ks_did_change(ks_env, potential_changed=.TRUE.)

      END DO scf_outer_loop

      converged = inner_loop_converged .AND. outer_loop_converged
      total_scf_steps = total_steps

      IF (dft_control%qs_control%cdft) &
         dft_control%qs_control%cdft_control%total_steps = &
         dft_control%qs_control%cdft_control%total_steps + total_steps

      IF (.NOT. converged) THEN
         IF (scf_control%ignore_convergence_failure .OR. should_stop) THEN
            CALL cp_warn(__LOCATION__, "SCF run NOT converged")
         ELSE
            CALL cp_abort(__LOCATION__, &
                          "SCF run NOT converged. To continue the calculation "// &
                          "regardless, please set the keyword IGNORE_CONVERGENCE_FAILURE.")
         END IF
      END IF

      ! Skip Harris functional calculation if ground-state is NOT converged
      IF (qs_env%energy_correction) THEN
         CALL get_qs_env(qs_env, ec_env=ec_env)
         ec_env%do_skip = .FALSE.
         IF (ec_env%skip_ec .AND. .NOT. converged) ec_env%do_skip = .TRUE.
      END IF

      ! if needed copy mo_coeff dbcsr->fm for later use in post_scf!fm->dbcsr
      DO ispin = 1, SIZE(mos) !fm -> dbcsr
         IF (mos(ispin)%use_mo_coeff_b) THEN !fm->dbcsr
            IF (.NOT. ASSOCIATED(mos(ispin)%mo_coeff_b)) & !fm->dbcsr
               CPABORT("mo_coeff_b is not allocated") !fm->dbcsr
            CALL copy_dbcsr_to_fm(mos(ispin)%mo_coeff_b, & !fm->dbcsr
                                  mos(ispin)%mo_coeff) !fm -> dbcsr
         END IF !fm->dbcsr
      END DO !fm -> dbcsr

      CALL cp_rm_iter_level(logger%iter_info, level_name="QS_SCF")
      CALL timestop(handle)

   END SUBROUTINE scf_env_do_scf

! **************************************************************************************************
!> \brief inits those objects needed if you want to restart the scf with, say
!>        only a new initial guess, or different density functional or ...
!>        this will happen just before the scf loop starts
!> \param scf_env ...
!> \param qs_env ...
!> \param scf_section ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE init_scf_loop(scf_env, qs_env, scf_section)

      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: scf_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'init_scf_loop'

      INTEGER                                            :: handle, ispin, nmo, number_of_OT_envs
      LOGICAL                                            :: do_kpoints, do_rotation, &
                                                            has_unit_metric, is_full_all
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dbcsr_type), POINTER                          :: orthogonality_metric
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (scf_control, matrix_s, matrix_ks, dft_control, mos, mo_coeff, kpoints)

      CPASSERT(ASSOCIATED(scf_env))
      CPASSERT(ASSOCIATED(qs_env))

      CALL get_qs_env(qs_env=qs_env, &
                      scf_control=scf_control, &
                      dft_control=dft_control, &
                      do_kpoints=do_kpoints, &
                      kpoints=kpoints, &
                      mos=mos)

      ! if using mo_coeff_b then copy to fm
      DO ispin = 1, SIZE(mos) !fm->dbcsr
         IF (mos(1)%use_mo_coeff_b) THEN !fm->dbcsr
            CALL copy_dbcsr_to_fm(mos(ispin)%mo_coeff_b, mos(ispin)%mo_coeff) !fm->dbcsr
         END IF !fm->dbcsr
      END DO !fm->dbcsr

      ! this just guarantees that all mo_occupations match the eigenvalues, if smear
      DO ispin = 1, dft_control%nspins
         ! do not reset mo_occupations if the maximum overlap method is in use
         IF (.NOT. scf_control%diagonalization%mom) THEN
            !if the hair probes section is present, this sends hairy_probes to set_mo_occupation subroutine
            !and switches off the standard smearing
            IF (dft_control%hairy_probes .EQV. .TRUE.) THEN
               IF (scf_env%outer_scf%iter_count > 0) THEN
                  scf_control%smear%do_smear = .FALSE.
                  CALL set_mo_occupation(mo_set=mos(ispin), &
                                         smear=scf_control%smear, &
                                         probe=dft_control%probe)
               END IF
            ELSE
               ! GCE scheme 2 with two spins: one Fermi level and one N_e for both spins (set_mo_occupation_3, called
               ! with both spin sets in the diagonalisation); a call with one spin set would take that spin's electron
               ! count as N_e. Scheme 1 (mu mixing) uses the per-spin call.
               IF (.NOT. scf_control%gce%do_gce .OR. (dft_control%nspins /= 1 .AND. scf_control%gce%scheme == 2)) THEN
                  CALL set_mo_occupation(mo_set=mos(ispin), &
                                         smear=scf_control%smear)
               ELSE
                  CALL set_mo_occupation(mo_set=mos(ispin), &
                                         smear=scf_control%smear, &
                                         gce=scf_control%gce)
               END IF
            END IF
         END IF
      END DO

      SELECT CASE (scf_env%method)
      CASE DEFAULT

         CPABORT("Unknown SCF method <"//TRIM(cp_to_string(scf_env%method))//"> found. Check the code!")

      CASE (filter_matrix_diag_method_nr)

         IF (.NOT. scf_env%skip_diis) THEN
            IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
               ALLOCATE (scf_env%scf_diis_buffer)
               CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=scf_control%max_diis)
            END IF
            CALL qs_diis_b_clear(scf_env%scf_diis_buffer)
         END IF

      CASE (general_diag_method_nr, special_diag_method_nr, block_krylov_diag_method_nr, smeagol_method_nr)
         IF (.NOT. scf_env%skip_diis) THEN
            IF (do_kpoints) THEN
               IF (.NOT. ASSOCIATED(kpoints%scf_diis_buffer)) THEN
                  ALLOCATE (kpoints%scf_diis_buffer)
                  CALL qs_diis_b_create_kp(kpoints%scf_diis_buffer, nbuffer=scf_control%max_diis)
               END IF
               CALL qs_diis_b_clear_kp(kpoints%scf_diis_buffer)
            ELSE
               IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
                  ALLOCATE (scf_env%scf_diis_buffer)
                  CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=scf_control%max_diis)
               END IF
               CALL qs_diis_b_clear(scf_env%scf_diis_buffer)
            END IF
         END IF

      CASE (ot_diag_method_nr)
         CALL get_qs_env(qs_env, matrix_ks=matrix_ks, matrix_s=matrix_s)

         IF (.NOT. scf_env%skip_diis) THEN
            IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
               ALLOCATE (scf_env%scf_diis_buffer)
               CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=scf_control%max_diis)
            END IF
            CALL qs_diis_b_clear(scf_env%scf_diis_buffer)
         END IF

         ! disable DFTB and SE for now
         IF (dft_control%qs_control%dftb .OR. &
             dft_control%qs_control%xtb .OR. &
             dft_control%qs_control%semi_empirical) THEN
            CPABORT("DFTB and SE not available with OT/DIAG")
         END IF

         ! if an old preconditioner is still around (i.e. outer SCF is active),
         ! remove it if this could be worthwhile
         CALL restart_preconditioner(qs_env, scf_env%ot_preconditioner, &
                                     scf_control%diagonalization%ot_settings%preconditioner_type, &
                                     dft_control%nspins)

         CALL prepare_preconditioner(qs_env, mos, matrix_ks, matrix_s, scf_env%ot_preconditioner, &
                                     scf_control%diagonalization%ot_settings%preconditioner_type, &
                                     scf_control%diagonalization%ot_settings%precond_solver_type, &
                                     scf_control%diagonalization%ot_settings%energy_gap, dft_control%nspins)

      CASE (block_davidson_diag_method_nr)
         ! Preconditioner initialized within the loop, when required
      CASE (ot_method_nr)
         CALL get_qs_env(qs_env, &
                         has_unit_metric=has_unit_metric, &
                         matrix_s=matrix_s, &
                         matrix_ks=matrix_ks)

         ! reortho the wavefunctions if we are having an outer scf and
         ! this is not the first iteration
         ! this is useful to avoid the build-up of numerical noise
         ! however, we can not play this trick if restricted (don't mix non-equivalent orbs)
         IF (scf_control%do_outer_scf_reortho) THEN
            IF (scf_control%outer_scf%have_scf .AND. .NOT. dft_control%restricted) THEN
               IF (scf_env%outer_scf%iter_count > 0) THEN
                  DO ispin = 1, dft_control%nspins
                     CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff, nmo=nmo)
                     IF (has_unit_metric) THEN
                        CALL make_basis_simple(mo_coeff, nmo)
                     ELSE
                        CALL make_basis_sm(mo_coeff, nmo, matrix_s(1)%matrix)
                     END IF
                  END DO
               END IF
            END IF
         ELSE
            ! dont need any dirty trick for the numerically stable irac algorithm.
         END IF

         IF (.NOT. ASSOCIATED(scf_env%qs_ot_env)) THEN

            ! restricted calculations require just one set of OT orbitals
            number_of_OT_envs = dft_control%nspins
            IF (dft_control%restricted) number_of_OT_envs = 1

            ALLOCATE (scf_env%qs_ot_env(number_of_OT_envs))

            ! XXX Joost XXX should disentangle reading input from this part
            IF (scf_env%outer_scf%iter_count > 0) THEN
               IF (scf_env%iter_delta < scf_control%eps_diis) THEN
                  scf_env%qs_ot_env(1)%settings%ot_state = 1
               END IF
            END IF
            !
            CALL ot_scf_read_input(scf_env%qs_ot_env, scf_section)
            !
            IF (scf_env%outer_scf%iter_count > 0) THEN
               IF (scf_env%qs_ot_env(1)%settings%ot_state == 1) THEN
                  scf_control%max_scf = MAX(scf_env%qs_ot_env(1)%settings%max_scf_diis, &
                                            scf_control%max_scf)
               END IF
            END IF

            ! keep a note that we are restricted
            IF (dft_control%restricted) THEN
               scf_env%qs_ot_env(1)%restricted = .TRUE.
               ! requires rotation
               IF (.NOT. scf_env%qs_ot_env(1)%settings%do_rotation) &
                  CALL cp_abort(__LOCATION__, &
                                "Restricted calculation with OT requires orbital rotation. Please "// &
                                "activate the OT%ROTATION keyword!")
            ELSE
               scf_env%qs_ot_env(:)%restricted = .FALSE.
            END IF

            ! this will rotate the MOs to be eigen states, which is not compatible with rotation
            ! e.g. mo_derivs here do not yet include potentially different occupations numbers
            do_rotation = scf_env%qs_ot_env(1)%settings%do_rotation
            ! only full all needs rotation
            is_full_all = scf_env%qs_ot_env(1)%settings%preconditioner_type == ot_precond_full_all
            IF (do_rotation .AND. is_full_all) &
               CPABORT('PRECONDITIONER FULL_ALL is not compatible with ROTATION.')

            ! might need the KS matrix to init properly
            CALL qs_ks_update_qs_env(qs_env, just_energy=.FALSE., &
                                     calculate_forces=.FALSE.)

            ! if an old preconditioner is still around (i.e. outer SCF is active),
            ! remove it if this could be worthwhile
            IF (.NOT. reuse_precond) &
               CALL restart_preconditioner(qs_env, scf_env%ot_preconditioner, &
                                           scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                                           dft_control%nspins)

            !
            ! preconditioning still needs to be done correctly with has_unit_metric
            ! notice that a big part of the preconditioning (S^-1) is fine anyhow
            !
            IF (has_unit_metric) THEN
               NULLIFY (orthogonality_metric)
            ELSE
               orthogonality_metric => matrix_s(1)%matrix
            END IF

            IF (.NOT. reuse_precond) &
               CALL prepare_preconditioner(qs_env, mos, matrix_ks, matrix_s, scf_env%ot_preconditioner, &
                                           scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                                           scf_env%qs_ot_env(1)%settings%precond_solver_type, &
                                           scf_env%qs_ot_env(1)%settings%energy_gap, dft_control%nspins, &
                                           has_unit_metric=has_unit_metric, &
                                           chol_type=scf_env%qs_ot_env(1)%settings%cholesky_type)
            IF (reuse_precond) reuse_precond = .FALSE.

            CALL ot_scf_init(mo_array=mos, matrix_s=orthogonality_metric, &
                             broyden_adaptive_sigma=qs_env%broyden_adaptive_sigma, &
                             qs_ot_env=scf_env%qs_ot_env, matrix_ks=matrix_ks(1)%matrix)

            SELECT CASE (scf_env%qs_ot_env(1)%settings%preconditioner_type)
            CASE (ot_precond_none)
            CASE (ot_precond_full_all, ot_precond_full_single_inverse)
               DO ispin = 1, SIZE(scf_env%qs_ot_env)
                  CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin), &
                                                scf_env%ot_preconditioner(ispin)%preconditioner)
               END DO
            CASE (ot_precond_s_inverse, ot_precond_full_single)
               DO ispin = 1, SIZE(scf_env%qs_ot_env)
                  CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin), &
                                                scf_env%ot_preconditioner(1)%preconditioner)
               END DO
            CASE DEFAULT
               DO ispin = 1, SIZE(scf_env%qs_ot_env)
                  CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin), &
                                                scf_env%ot_preconditioner(1)%preconditioner)
               END DO
            END SELECT
         END IF

         ! if we have non-uniform occupations we should be using rotation
         do_rotation = scf_env%qs_ot_env(1)%settings%do_rotation
         DO ispin = 1, SIZE(mos)
            IF (.NOT. mos(ispin)%uniform_occupation) THEN
               CPASSERT(do_rotation)
            END IF
         END DO
      END SELECT

      ! another safety check
      IF (dft_control%low_spin_roks) THEN
         CPASSERT(scf_env%method == ot_method_nr)
         do_rotation = scf_env%qs_ot_env(1)%settings%do_rotation
         CPASSERT(do_rotation)
      END IF

      CALL timestop(handle)

   END SUBROUTINE init_scf_loop

! **************************************************************************************************
!> \brief perform cleanup operations (like releasing temporary storage)
!>      at the end of the scf
!> \param scf_env ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE scf_env_cleanup(scf_env)
      TYPE(qs_scf_env_type), INTENT(INOUT)               :: scf_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'scf_env_cleanup'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Release SCF work storage
      CALL cp_fm_release(scf_env%scf_work1)

      IF (ASSOCIATED(scf_env%scf_work1_red)) THEN
         CALL cp_fm_release(scf_env%scf_work1_red)
      END IF
      IF (ASSOCIATED(scf_env%scf_work2)) THEN
         CALL cp_fm_release(scf_env%scf_work2)
         DEALLOCATE (scf_env%scf_work2)
         NULLIFY (scf_env%scf_work2)
      END IF
      IF (ASSOCIATED(scf_env%scf_work2_red)) THEN
         CALL cp_fm_release(scf_env%scf_work2_red)
         DEALLOCATE (scf_env%scf_work2_red)
         NULLIFY (scf_env%scf_work2_red)
      END IF
      IF (ASSOCIATED(scf_env%ortho)) THEN
         CALL cp_fm_release(scf_env%ortho)
         DEALLOCATE (scf_env%ortho)
         NULLIFY (scf_env%ortho)
      END IF
      IF (ASSOCIATED(scf_env%ortho_red)) THEN
         CALL cp_fm_release(scf_env%ortho_red)
         DEALLOCATE (scf_env%ortho_red)
         NULLIFY (scf_env%ortho_red)
      END IF
      IF (ASSOCIATED(scf_env%ortho_m1)) THEN
         CALL cp_fm_release(scf_env%ortho_m1)
         DEALLOCATE (scf_env%ortho_m1)
         NULLIFY (scf_env%ortho_m1)
      END IF
      IF (ASSOCIATED(scf_env%ortho_m1_red)) THEN
         CALL cp_fm_release(scf_env%ortho_m1_red)
         DEALLOCATE (scf_env%ortho_m1_red)
         NULLIFY (scf_env%ortho_m1_red)
      END IF

      IF (ASSOCIATED(scf_env%ortho_dbcsr)) THEN
         CALL dbcsr_deallocate_matrix(scf_env%ortho_dbcsr)
      END IF
      IF (ASSOCIATED(scf_env%buf1_dbcsr)) THEN
         CALL dbcsr_deallocate_matrix(scf_env%buf1_dbcsr)
      END IF
      IF (ASSOCIATED(scf_env%buf2_dbcsr)) THEN
         CALL dbcsr_deallocate_matrix(scf_env%buf2_dbcsr)
      END IF

      IF (ASSOCIATED(scf_env%p_mix_new)) THEN
         CALL dbcsr_deallocate_matrix_set(scf_env%p_mix_new)
      END IF

      IF (ASSOCIATED(scf_env%p_delta)) THEN
         CALL dbcsr_deallocate_matrix_set(scf_env%p_delta)
      END IF

      ! Method dependent cleanup
      SELECT CASE (scf_env%method)
      CASE (ot_method_nr)
         !
      CASE (ot_diag_method_nr)
         !
      CASE (general_diag_method_nr)
         !
      CASE (special_diag_method_nr)
         !
      CASE (block_krylov_diag_method_nr)
      CASE (block_davidson_diag_method_nr)
         CALL block_davidson_deallocate(scf_env%block_davidson_env)
      CASE (filter_matrix_diag_method_nr)
         !
      CASE (smeagol_method_nr)
         !
      CASE DEFAULT
         CPABORT("unknown scf method method:"//cp_to_string(scf_env%method))
      END SELECT

      IF (ASSOCIATED(scf_env%outer_scf%variables)) THEN
         DEALLOCATE (scf_env%outer_scf%variables)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%count)) THEN
         DEALLOCATE (scf_env%outer_scf%count)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%gradient)) THEN
         DEALLOCATE (scf_env%outer_scf%gradient)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%energy)) THEN
         DEALLOCATE (scf_env%outer_scf%energy)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian) .AND. &
          scf_env%outer_scf%deallocate_jacobian) THEN
         DEALLOCATE (scf_env%outer_scf%inv_jacobian)
      END IF

      CALL timestop(handle)

   END SUBROUTINE scf_env_cleanup

! **************************************************************************************************
!> \brief perform a CDFT scf procedure in the given qs_env
!> \param qs_env the qs_environment where to perform the scf procedure
!> \param should_stop flag determining if calculation should stop
!> \par History
!>      12.2015 Created
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE cdft_scf(qs_env, should_stop)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(OUT)                               :: should_stop

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cdft_scf'

      INTEGER                                            :: handle, iatom, ispin, ivar, nmo, nvar, &
                                                            output_unit, tsteps
      LOGICAL                                            :: cdft_loop_converged, converged, &
                                                            exit_cdft_loop, first_iteration, &
                                                            my_uocc, uniform_occupation
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_occupations
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, input, scf_section

      NULLIFY (scf_env, ks_env, energy, rho, matrix_s, rho_ao, cdft_control, logger, &
               dft_control, pw_env, auxbas_pw_pool, energy, ks_env, scf_env, dft_section, &
               input, scf_section, scf_control, mos, mo_occupations)
      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, scf_env=scf_env, energy=energy, &
                      dft_control=dft_control, scf_control=scf_control, &
                      ks_env=ks_env, input=input)

      CALL timeset(routineN//"_loop", handle)
      dft_section => section_vals_get_subs_vals(input, "DFT")
      scf_section => section_vals_get_subs_vals(dft_section, "SCF")
      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".scfLog")
      first_iteration = .TRUE.

      cdft_control => dft_control%qs_control%cdft_control

      scf_env%outer_scf%iter_count = 0
      cdft_control%total_steps = 0

      ! Write some info about the CDFT calculation
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
            "CDFT EXTERNAL SCF WAVEFUNCTION OPTIMIZATION"
         CALL qs_scf_cdft_initial_info(output_unit, cdft_control)
      END IF
      IF (cdft_control%reuse_precond) THEN
         reuse_precond = .FALSE.
         cdft_control%nreused = 0
      END IF
      cdft_outer_loop: DO
         ! Change outer_scf settings to OT settings
         CALL outer_loop_switch(scf_env, scf_control, cdft_control, cdft2ot)
         ! Solve electronic structure with fixed value of constraint
         CALL scf_env_do_scf(scf_env=scf_env, scf_control=scf_control, qs_env=qs_env, &
                             converged=converged, should_stop=should_stop, total_scf_steps=tsteps)
         ! Decide whether to reuse the preconditioner on the next iteration
         IF (cdft_control%reuse_precond) THEN
            ! For convergence in exactly one step, the preconditioner is always reused (assuming max_reuse > 0)
            ! usually this means that the electronic structure has already converged to the correct state
            ! but the constraint optimizer keeps jumping over the optimal solution
            IF (scf_env%outer_scf%iter_count == 1 .AND. scf_env%iter_count == 1 &
                .AND. cdft_control%total_steps /= 1) &
               cdft_control%nreused = cdft_control%nreused - 1
            ! SCF converged in less than precond_freq steps
            IF (scf_env%outer_scf%iter_count == 1 .AND. scf_env%iter_count <= cdft_control%precond_freq .AND. &
                cdft_control%total_steps /= 1 .AND. cdft_control%nreused < cdft_control%max_reuse) THEN
               reuse_precond = .TRUE.
               cdft_control%nreused = cdft_control%nreused + 1
            ELSE
               reuse_precond = .FALSE.
               cdft_control%nreused = 0
            END IF
         END IF
         ! Update history purging counters
         IF (first_iteration .AND. cdft_control%purge_history) THEN
            cdft_control%istep = cdft_control%istep + 1
            IF (scf_env%outer_scf%iter_count > 1) THEN
               cdft_control%nbad_conv = cdft_control%nbad_conv + 1
               IF (cdft_control%nbad_conv >= cdft_control%purge_freq .AND. &
                   cdft_control%istep >= cdft_control%purge_offset) THEN
                  cdft_control%nbad_conv = 0
                  cdft_control%istep = 0
                  cdft_control%should_purge = .TRUE.
               END IF
            END IF
         END IF
         first_iteration = .FALSE.
         ! Change outer_scf settings to CDFT settings
         CALL outer_loop_switch(scf_env, scf_control, cdft_control, ot2cdft)
         CALL qs_scf_check_outer_exit(qs_env, scf_env, scf_control, should_stop, &
                                      cdft_loop_converged, exit_cdft_loop)
         CALL qs_scf_cdft_info(output_unit, scf_control, scf_env, cdft_control, &
                               energy, cdft_control%total_steps, &
                               should_stop, cdft_loop_converged, cdft_loop=.TRUE.)
         IF (exit_cdft_loop) EXIT cdft_outer_loop
         ! Check if the inverse Jacobian needs to be calculated
         CALL qs_calculate_inverse_jacobian(qs_env)
         ! Check if a line search should be performed to find an optimal step size for the optimizer
         CALL qs_cdft_line_search(qs_env)
         ! Optimize constraint
         CALL outer_loop_optimize(scf_env, scf_control)
         CALL outer_loop_update_qs_env(qs_env, scf_env)
         CALL qs_ks_did_change(ks_env, potential_changed=.TRUE.)
      END DO cdft_outer_loop

      cdft_control%ienergy = cdft_control%ienergy + 1

      ! Store needed arrays for ET coupling calculation
      IF (cdft_control%do_et) THEN
         CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, mos=mos)
         nvar = SIZE(cdft_control%target)
         ! Matrix representation of weight function
         ALLOCATE (cdft_control%wmat(nvar))
         DO ivar = 1, nvar
            CALL dbcsr_init_p(cdft_control%wmat(ivar)%matrix)
            CALL dbcsr_copy(cdft_control%wmat(ivar)%matrix, matrix_s(1)%matrix, &
                            name="ET_RESTRAINT_MATRIX")
            CALL dbcsr_set(cdft_control%wmat(ivar)%matrix, 0.0_dp)
            CALL integrate_v_rspace(cdft_control%group(ivar)%weight, &
                                    hmat=cdft_control%wmat(ivar), qs_env=qs_env, &
                                    calculate_forces=.FALSE., &
                                    gapw=dft_control%qs_control%gapw)
         END DO
         ! Overlap matrix
         CALL dbcsr_init_p(cdft_control%matrix_s%matrix)
         CALL dbcsr_copy(cdft_control%matrix_s%matrix, matrix_s(1)%matrix, &
                         name="OVERLAP")
         ! Molecular orbital coefficients
         NULLIFY (cdft_control%mo_coeff)
         ALLOCATE (cdft_control%mo_coeff(dft_control%nspins))
         DO ispin = 1, dft_control%nspins
            CALL cp_fm_create(matrix=cdft_control%mo_coeff(ispin), &
                              matrix_struct=qs_env%mos(ispin)%mo_coeff%matrix_struct, &
                              name="MO_COEFF_A"//TRIM(ADJUSTL(cp_to_string(ispin)))//"MATRIX")
            CALL cp_fm_to_fm(qs_env%mos(ispin)%mo_coeff, &
                             cdft_control%mo_coeff(ispin))
         END DO
         ! Density matrix
         IF (cdft_control%calculate_metric) THEN
            CALL get_qs_env(qs_env, rho=rho)
            CALL qs_rho_get(rho, rho_ao=rho_ao)
            ALLOCATE (cdft_control%matrix_p(dft_control%nspins))
            DO ispin = 1, dft_control%nspins
               NULLIFY (cdft_control%matrix_p(ispin)%matrix)
               CALL dbcsr_init_p(cdft_control%matrix_p(ispin)%matrix)
               CALL dbcsr_copy(cdft_control%matrix_p(ispin)%matrix, rho_ao(ispin)%matrix, &
                               name="DENSITY MATRIX")
            END DO
         END IF
         ! Copy occupation numbers if non-uniform occupation
         uniform_occupation = .TRUE.
         DO ispin = 1, dft_control%nspins
            CALL get_mo_set(mo_set=mos(ispin), uniform_occupation=my_uocc)
            uniform_occupation = uniform_occupation .AND. my_uocc
         END DO
         IF (.NOT. uniform_occupation) THEN
            ALLOCATE (cdft_control%occupations(dft_control%nspins))
            DO ispin = 1, dft_control%nspins
               CALL get_mo_set(mo_set=mos(ispin), &
                               nmo=nmo, &
                               occupation_numbers=mo_occupations)
               ALLOCATE (cdft_control%occupations(ispin)%array(nmo))
               cdft_control%occupations(ispin)%array(1:nmo) = mo_occupations(1:nmo)
            END DO
         END IF
      END IF

      ! Deallocate constraint storage if forces are not needed
      ! In case of a simulation with multiple force_evals,
      ! deallocate only if weight function should not be copied to different force_evals
      IF (.NOT. (cdft_control%save_pot .OR. cdft_control%transfer_pot)) THEN
         CALL get_qs_env(qs_env, pw_env=pw_env)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         DO iatom = 1, SIZE(cdft_control%group)
            CALL auxbas_pw_pool%give_back_pw(cdft_control%group(iatom)%weight)
            DEALLOCATE (cdft_control%group(iatom)%weight)
         END DO
         IF (cdft_control%atomic_charges) THEN
            DO iatom = 1, cdft_control%natoms
               CALL auxbas_pw_pool%give_back_pw(cdft_control%charge(iatom))
            END DO
            DEALLOCATE (cdft_control%charge)
         END IF
         IF (cdft_control%type == outer_scf_becke_constraint .AND. &
             cdft_control%becke_control%cavity_confine) THEN
            IF (.NOT. ASSOCIATED(cdft_control%becke_control%cavity_mat)) THEN
               CALL auxbas_pw_pool%give_back_pw(cdft_control%becke_control%cavity)
            ELSE
               DEALLOCATE (cdft_control%becke_control%cavity_mat)
            END IF
         ELSE IF (cdft_control%type == outer_scf_hirshfeld_constraint) THEN
            IF (ASSOCIATED(cdft_control%hirshfeld_control%hirshfeld_env%fnorm)) THEN
               CALL auxbas_pw_pool%give_back_pw(cdft_control%hirshfeld_control%hirshfeld_env%fnorm)
            END IF
         END IF
         IF (ASSOCIATED(cdft_control%charges_fragment)) DEALLOCATE (cdft_control%charges_fragment)
         cdft_control%need_pot = .TRUE.
         cdft_control%external_control = .FALSE.
      END IF

      CALL timestop(handle)

   END SUBROUTINE cdft_scf

! **************************************************************************************************
!> \brief perform cleanup operations for cdft_control
!> \param cdft_control container for the external CDFT SCF loop variables
!> \par History
!>      12.2015 created [Nico Holmberg]
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE cdft_control_cleanup(cdft_control)
      TYPE(cdft_control_type), POINTER                   :: cdft_control

      IF (ASSOCIATED(cdft_control%constraint%variables)) &
         DEALLOCATE (cdft_control%constraint%variables)
      IF (ASSOCIATED(cdft_control%constraint%count)) &
         DEALLOCATE (cdft_control%constraint%count)
      IF (ASSOCIATED(cdft_control%constraint%gradient)) &
         DEALLOCATE (cdft_control%constraint%gradient)
      IF (ASSOCIATED(cdft_control%constraint%energy)) &
         DEALLOCATE (cdft_control%constraint%energy)
      IF (ASSOCIATED(cdft_control%constraint%inv_jacobian) .AND. &
          cdft_control%constraint%deallocate_jacobian) &
         DEALLOCATE (cdft_control%constraint%inv_jacobian)

   END SUBROUTINE cdft_control_cleanup

! **************************************************************************************************
!> \brief Calculates the finite difference inverse Jacobian
!> \param qs_env the qs_environment_type where to compute the Jacobian
!> \par History
!>      01.2017 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE qs_calculate_inverse_jacobian(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_calculate_inverse_jacobian'

      CHARACTER(len=default_path_length)                 :: project_name
      INTEGER                                            :: counter, handle, i, ispin, iter_count, &
                                                            iwork, j, max_scf, nspins, nsteps, &
                                                            nvar, nwork, output_unit, pwork, &
                                                            tsteps, twork
      LOGICAL                                            :: converged, explicit_jacobian, &
                                                            should_build, should_stop, &
                                                            use_md_history
      REAL(KIND=dp)                                      :: inv_error, step_size
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coeff, dh, step_multiplier
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: jacobian
      REAL(KIND=dp), DIMENSION(:), POINTER               :: energy
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient, inv_jacobian
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(cp_logger_type), POINTER                      :: logger, tmp_logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_rmpv
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), ALLOCATABLE, DIMENSION(:)       :: mos_stashed
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy_qs
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      NULLIFY (energy, gradient, p_rmpv, rho_ao_kp, mos, rho, &
               ks_env, scf_env, scf_control, dft_control, cdft_control, &
               inv_jacobian, para_env, tmp_logger, energy_qs)
      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, scf_env=scf_env, ks_env=ks_env, &
                      scf_control=scf_control, mos=mos, rho=rho, &
                      dft_control=dft_control, &
                      para_env=para_env, energy=energy_qs)
      explicit_jacobian = .FALSE.
      should_build = .FALSE.
      use_md_history = .FALSE.
      iter_count = scf_env%outer_scf%iter_count
      ! Quick exit if optimizer does not require Jacobian
      IF (.NOT. ASSOCIATED(scf_control%outer_scf%cdft_opt_control)) RETURN
      ! Check if Jacobian should be calculated and initialize
      CALL timeset(routineN, handle)
      CALL initialize_inverse_jacobian(scf_control, scf_env, explicit_jacobian, should_build, used_history)
      IF (scf_control%outer_scf%cdft_opt_control%jacobian_restart) THEN
         ! Restart from previously calculated inverse Jacobian
         should_build = .FALSE.
         CALL restart_inverse_jacobian(qs_env)
      END IF
      IF (should_build) THEN
         scf_env%outer_scf%deallocate_jacobian = .FALSE.
         ! Actually need to (re)build the Jacobian
         IF (explicit_jacobian) THEN
            ! Build Jacobian with finite differences
            cdft_control => dft_control%qs_control%cdft_control
            IF (.NOT. ASSOCIATED(cdft_control)) &
               CALL cp_abort(__LOCATION__, &
                             "Optimizers that need the explicit Jacobian can"// &
                             " only be used together with a valid CDFT constraint.")
            ! Redirect output from Jacobian calculation to a new file by creating a temporary logger
            project_name = logger%iter_info%project_name
            CALL create_tmp_logger(para_env, project_name, "-JacobianInfo.out", output_unit, tmp_logger)
            ! Save last converged state so we can roll back to it (mo_coeff and some outer_loop variables)
            nspins = dft_control%nspins
            ALLOCATE (mos_stashed(nspins))
            DO ispin = 1, nspins
               CALL duplicate_mo_set(mos_stashed(ispin), mos(ispin))
            END DO
            CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
            p_rmpv => rho_ao_kp(:, 1)
            ! Allocate work
            nvar = SIZE(scf_env%outer_scf%variables, 1)
            max_scf = scf_control%outer_scf%max_scf + 1
            ALLOCATE (gradient(nvar, max_scf))
            gradient = scf_env%outer_scf%gradient
            ALLOCATE (energy(max_scf))
            energy = scf_env%outer_scf%energy
            ALLOCATE (jacobian(nvar, nvar))
            jacobian = 0.0_dp
            nsteps = cdft_control%total_steps
            ! Setup finite difference scheme
            CALL prepare_jacobian_stencil(qs_env, output_unit, nwork, pwork, coeff, step_multiplier, dh)
            twork = pwork - nwork
            DO i = 1, nvar
               jacobian(i, :) = coeff(0)*scf_env%outer_scf%gradient(i, iter_count)
            END DO
            ! Calculate the Jacobian by perturbing each Lagrangian and recalculating the energy self-consistently
            CALL cp_add_default_logger(tmp_logger)
            DO i = 1, nvar
               IF (output_unit > 0) THEN
                  WRITE (output_unit, FMT="(A)") " "
                  WRITE (output_unit, FMT="(A)") " #####################################"
                  WRITE (output_unit, '(A,I3,A,I3,A)') &
                     " ###  Constraint        ", i, " of ", nvar, " ###"
                  WRITE (output_unit, FMT="(A)") " #####################################"
               END IF
               counter = 0
               DO iwork = nwork, pwork
                  IF (iwork == 0) CYCLE
                  counter = counter + 1
                  IF (output_unit > 0) THEN
                     WRITE (output_unit, FMT="(A)") " #####################################"
                     WRITE (output_unit, '(A,I3,A,I3,A)') &
                        " ###  Energy evaluation ", counter, " of ", twork, " ###"
                     WRITE (output_unit, FMT="(A)") " #####################################"
                  END IF
                  IF (SIZE(scf_control%outer_scf%cdft_opt_control%jacobian_step) == 1) THEN
                     step_size = scf_control%outer_scf%cdft_opt_control%jacobian_step(1)
                  ELSE
                     step_size = scf_control%outer_scf%cdft_opt_control%jacobian_step(i)
                  END IF
                  scf_env%outer_scf%variables(:, iter_count + 1) = scf_env%outer_scf%variables(:, iter_count)
                  scf_env%outer_scf%variables(i, iter_count + 1) = scf_env%outer_scf%variables(i, iter_count) + &
                                                                   step_multiplier(iwork)*step_size
                  CALL outer_loop_update_qs_env(qs_env, scf_env)
                  CALL qs_ks_did_change(ks_env, potential_changed=.TRUE.)
                  CALL outer_loop_switch(scf_env, scf_control, cdft_control, cdft2ot)
                  CALL scf_env_do_scf(scf_env=scf_env, scf_control=scf_control, qs_env=qs_env, &
                                      converged=converged, should_stop=should_stop, total_scf_steps=tsteps)
                  CALL outer_loop_switch(scf_env, scf_control, cdft_control, ot2cdft)
                  ! Update (iter_count + 1) element of gradient and print constraint info
                  scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count + 1
                  CALL outer_loop_gradient(qs_env, scf_env)
                  CALL qs_scf_cdft_info(output_unit, scf_control, scf_env, cdft_control, &
                                        energy_qs, cdft_control%total_steps, &
                                        should_stop=.FALSE., outer_loop_converged=.FALSE., cdft_loop=.FALSE.)
                  scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count - 1
                  ! Update Jacobian
                  DO j = 1, nvar
                     jacobian(j, i) = jacobian(j, i) + coeff(iwork)*scf_env%outer_scf%gradient(j, iter_count + 1)
                  END DO
                  ! Reset everything to last converged state
                  scf_env%outer_scf%variables(:, iter_count + 1) = 0.0_dp
                  scf_env%outer_scf%gradient = gradient
                  scf_env%outer_scf%energy = energy
                  cdft_control%total_steps = nsteps
                  DO ispin = 1, nspins
                     CALL deallocate_mo_set(mos(ispin))
                     CALL duplicate_mo_set(mos(ispin), mos_stashed(ispin))
                     CALL calculate_density_matrix(mos(ispin), &
                                                   p_rmpv(ispin)%matrix)
                  END DO
                  CALL qs_rho_update_rho(rho, qs_env=qs_env)
                  CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
               END DO
            END DO
            CALL cp_rm_default_logger()
            CALL cp_logger_release(tmp_logger)
            ! Finalize and invert Jacobian
            DO j = 1, nvar
               DO i = 1, nvar
                  jacobian(i, j) = jacobian(i, j)/dh(j)
               END DO
            END DO
            IF (.NOT. ASSOCIATED(scf_env%outer_scf%inv_jacobian)) &
               ALLOCATE (scf_env%outer_scf%inv_jacobian(nvar, nvar))
            inv_jacobian => scf_env%outer_scf%inv_jacobian
            CALL invert_matrix(jacobian, inv_jacobian, inv_error)
            scf_control%outer_scf%cdft_opt_control%broyden_update = .FALSE.
            ! Release temporary storage
            DO ispin = 1, nspins
               CALL deallocate_mo_set(mos_stashed(ispin))
            END DO
            DEALLOCATE (mos_stashed, jacobian, gradient, energy, coeff, step_multiplier, dh)
            IF (output_unit > 0) THEN
               WRITE (output_unit, FMT="(/,A)") &
                  " ================================== JACOBIAN CALCULATED =================================="
               CALL close_file(unit_number=output_unit)
            END IF
         ELSE
            ! Build a strictly diagonal Jacobian from history and invert it
            CALL build_diagonal_jacobian(qs_env, used_history)
         END IF
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian) .AND. para_env%is_source()) THEN
         ! Write restart file for inverse Jacobian
         CALL print_inverse_jacobian(logger, scf_env%outer_scf%inv_jacobian, iter_count)
      END IF
      ! Update counter
      scf_control%outer_scf%cdft_opt_control%ijacobian(1) = scf_control%outer_scf%cdft_opt_control%ijacobian(1) + 1
      CALL timestop(handle)

   END SUBROUTINE qs_calculate_inverse_jacobian

! **************************************************************************************************
!> \brief Perform backtracking line search to find the optimal step size for the CDFT constraint
!>        optimizer. Assumes that the CDFT gradient function is a smooth function of the constraint
!>        variables.
!> \param qs_env the qs_environment_type where to perform the line search
!> \par History
!>      02.2017 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE qs_cdft_line_search(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_cdft_line_search'

      CHARACTER(len=default_path_length)                 :: project_name
      INTEGER                                            :: handle, i, ispin, iter_count, &
                                                            max_linesearch, max_scf, nspins, &
                                                            nsteps, nvar, output_unit, tsteps
      LOGICAL :: continue_ls, continue_ls_exit, converged, do_linesearch, found_solution, &
         reached_maxls, should_exit, should_stop, sign_changed
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: positive_sign
      REAL(KIND=dp)                                      :: alpha, alpha_ls, factor, norm_ls
      REAL(KIND=dp), DIMENSION(:), POINTER               :: energy
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient, inv_jacobian
      REAL(KIND=dp), EXTERNAL                            :: dnrm2
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(cp_logger_type), POINTER                      :: logger, tmp_logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_rmpv
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy_qs
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (energy, gradient, p_rmpv, rho_ao_kp, mos, rho, &
               ks_env, scf_env, scf_control, dft_control, &
               cdft_control, inv_jacobian, para_env, &
               tmp_logger, energy_qs)
      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, scf_env=scf_env, ks_env=ks_env, &
                      scf_control=scf_control, mos=mos, rho=rho, &
                      dft_control=dft_control, &
                      para_env=para_env, energy=energy_qs)
      do_linesearch = .FALSE.
      SELECT CASE (scf_control%outer_scf%optimizer)
      CASE DEFAULT
         do_linesearch = .FALSE.
      CASE (outer_scf_optimizer_newton_ls)
         do_linesearch = .TRUE.
      CASE (outer_scf_optimizer_broyden)
         SELECT CASE (scf_control%outer_scf%cdft_opt_control%broyden_type)
         CASE (broyden_type_1, broyden_type_2, broyden_type_1_explicit, broyden_type_2_explicit)
            do_linesearch = .FALSE.
         CASE (broyden_type_1_ls, broyden_type_1_explicit_ls, broyden_type_2_ls, broyden_type_2_explicit_ls)
            cdft_control => dft_control%qs_control%cdft_control
            IF (.NOT. ASSOCIATED(cdft_control)) &
               CALL cp_abort(__LOCATION__, &
                             "Optimizers that perform a line search can"// &
                             " only be used together with a valid CDFT constraint")
            IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian)) &
               do_linesearch = .TRUE.
         END SELECT
      END SELECT
      IF (do_linesearch) THEN
         BLOCK
            TYPE(mo_set_type), DIMENSION(:), ALLOCATABLE :: mos_ls, mos_stashed
            cdft_control => dft_control%qs_control%cdft_control
            IF (.NOT. ASSOCIATED(cdft_control)) &
               CALL cp_abort(__LOCATION__, &
                             "Optimizers that perform a line search can"// &
                             " only be used together with a valid CDFT constraint")
            CPASSERT(ASSOCIATED(scf_env%outer_scf%inv_jacobian))
            CPASSERT(ASSOCIATED(scf_control%outer_scf%cdft_opt_control))
            alpha = scf_control%outer_scf%cdft_opt_control%newton_step_save
            iter_count = scf_env%outer_scf%iter_count
            ! Redirect output from line search procedure to a new file by creating a temporary logger
            project_name = logger%iter_info%project_name
            CALL create_tmp_logger(para_env, project_name, "-LineSearch.out", output_unit, tmp_logger)
            ! Save last converged state so we can roll back to it (mo_coeff and some outer_loop variables)
            nspins = dft_control%nspins
            ALLOCATE (mos_stashed(nspins))
            DO ispin = 1, nspins
               CALL duplicate_mo_set(mos_stashed(ispin), mos(ispin))
            END DO
            CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
            p_rmpv => rho_ao_kp(:, 1)
            nsteps = cdft_control%total_steps
            ! Allocate work
            nvar = SIZE(scf_env%outer_scf%variables, 1)
            max_scf = scf_control%outer_scf%max_scf + 1
            max_linesearch = scf_control%outer_scf%cdft_opt_control%max_ls
            continue_ls = scf_control%outer_scf%cdft_opt_control%continue_ls
            factor = scf_control%outer_scf%cdft_opt_control%factor_ls
            continue_ls_exit = .FALSE.
            found_solution = .FALSE.
            ALLOCATE (gradient(nvar, max_scf))
            gradient = scf_env%outer_scf%gradient
            ALLOCATE (energy(max_scf))
            energy = scf_env%outer_scf%energy
            reached_maxls = .FALSE.
            ! Broyden optimizers: perform update of inv_jacobian if necessary
            IF (scf_control%outer_scf%cdft_opt_control%broyden_update) THEN
               CALL outer_loop_optimize(scf_env, scf_control)
               ! Reset the variables and prevent a reupdate of inv_jacobian
               scf_env%outer_scf%variables(:, iter_count + 1) = 0
               scf_control%outer_scf%cdft_opt_control%broyden_update = .FALSE.
            END IF
            ! Print some info
            IF (output_unit > 0) THEN
               WRITE (output_unit, FMT="(/,A)") &
                  " ================================== LINE SEARCH STARTED  =================================="
               WRITE (output_unit, FMT="(A,I5,A)") &
                  " Evaluating optimal step size for optimizer using a maximum of", max_linesearch, " steps"
               IF (continue_ls) THEN
                  WRITE (output_unit, FMT="(A)") &
                     " Line search continues until best step size is found or max steps are reached"
               END IF
               WRITE (output_unit, '(/,A,F5.3)') &
                  " Initial step size: ", alpha
               WRITE (output_unit, '(/,A,F5.3)') &
                  " Step size update factor: ", factor
               WRITE (output_unit, '(/,A,I10,A,I10)') &
                  " Energy evaluation: ", cdft_control%ienergy, ", CDFT SCF iteration: ", iter_count
            END IF
            ! Perform backtracking line search
            CALL cp_add_default_logger(tmp_logger)
            DO i = 1, max_linesearch
               IF (output_unit > 0) THEN
                  WRITE (output_unit, FMT="(A)") " "
                  WRITE (output_unit, FMT="(A)") " #####################################"
                  WRITE (output_unit, '(A,I10,A)') &
                     " ###  Line search step: ", i, " ###"
                  WRITE (output_unit, FMT="(A)") " #####################################"
               END IF
               inv_jacobian => scf_env%outer_scf%inv_jacobian
               ! Newton update of CDFT variables with a step size of alpha
               scf_env%outer_scf%variables(:, iter_count + 1) = scf_env%outer_scf%variables(:, iter_count) - alpha* &
                                                                MATMUL(inv_jacobian, scf_env%outer_scf%gradient(:, iter_count))
               ! With updated CDFT variables, perform SCF
               CALL outer_loop_update_qs_env(qs_env, scf_env)
               CALL qs_ks_did_change(ks_env, potential_changed=.TRUE.)
               CALL outer_loop_switch(scf_env, scf_control, cdft_control, cdft2ot)
               CALL scf_env_do_scf(scf_env=scf_env, scf_control=scf_control, qs_env=qs_env, &
                                   converged=converged, should_stop=should_stop, total_scf_steps=tsteps)
               CALL outer_loop_switch(scf_env, scf_control, cdft_control, ot2cdft)
               ! Update (iter_count + 1) element of gradient and print constraint info
               scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count + 1
               CALL outer_loop_gradient(qs_env, scf_env)
               CALL qs_scf_cdft_info(output_unit, scf_control, scf_env, cdft_control, &
                                     energy_qs, cdft_control%total_steps, &
                                     should_stop=.FALSE., outer_loop_converged=.FALSE., cdft_loop=.FALSE.)
               scf_env%outer_scf%iter_count = scf_env%outer_scf%iter_count - 1
               ! Store sign of initial gradient for each variable for continue_ls
               IF (continue_ls .AND. .NOT. ALLOCATED(positive_sign)) THEN
                  ALLOCATE (positive_sign(nvar))
                  DO ispin = 1, nvar
                     positive_sign(ispin) = scf_env%outer_scf%gradient(ispin, iter_count + 1) >= 0.0_dp
                  END DO
               END IF
               ! Check if the L2 norm of the gradient decreased
               inv_jacobian => scf_env%outer_scf%inv_jacobian
               IF (dnrm2(nvar, scf_env%outer_scf%gradient(:, iter_count + 1), 1) < &
                   dnrm2(nvar, scf_env%outer_scf%gradient(:, iter_count), 1)) THEN
                  ! Optimal step size found
                  IF (.NOT. continue_ls) THEN
                     should_exit = .TRUE.
                  ELSE
                     ! But line search continues for at least one more iteration in an attempt to find a better solution
                     ! if max number of steps is not exceeded
                     IF (found_solution) THEN
                        ! Check if the norm also decreased w.r.t. to previously found solution
                        IF (dnrm2(nvar, scf_env%outer_scf%gradient(:, iter_count + 1), 1) > norm_ls) THEN
                           ! Norm increased => accept previous solution and exit
                           continue_ls_exit = .TRUE.
                        END IF
                     END IF
                     ! Store current state and the value of alpha
                     IF (.NOT. continue_ls_exit) THEN
                        should_exit = .FALSE.
                        alpha_ls = alpha
                        found_solution = .TRUE.
                        norm_ls = dnrm2(nvar, scf_env%outer_scf%gradient(:, iter_count + 1), 1)
                        ! Check if the sign of the gradient has changed for all variables (w.r.t initial gradient)
                        ! In this case we should exit because further line search steps will just increase the norm
                        sign_changed = .TRUE.
                        DO ispin = 1, nvar
                           sign_changed = sign_changed .AND. (positive_sign(ispin) .NEQV. &
                                                              scf_env%outer_scf%gradient(ispin, iter_count + 1) >= 0.0_dp)
                        END DO
                        IF (.NOT. ALLOCATED(mos_ls)) THEN
                           ALLOCATE (mos_ls(nspins))
                        ELSE
                           DO ispin = 1, nspins
                              CALL deallocate_mo_set(mos_ls(ispin))
                           END DO
                        END IF
                        DO ispin = 1, nspins
                           CALL duplicate_mo_set(mos_ls(ispin), mos(ispin))
                        END DO
                        alpha = alpha*factor
                        ! Exit on last iteration
                        IF (i == max_linesearch) continue_ls_exit = .TRUE.
                        ! Exit if constraint target is satisfied to requested tolerance
                        IF (SQRT(MAXVAL(scf_env%outer_scf%gradient(:, scf_env%outer_scf%iter_count + 1)**2)) < &
                            scf_control%outer_scf%eps_scf) &
                           continue_ls_exit = .TRUE.
                        ! Exit if line search jumped over the optimal step length
                        IF (sign_changed) continue_ls_exit = .TRUE.
                     END IF
                  END IF
               ELSE
                  ! Gradient increased => alpha is too large (if the gradient function is smooth)
                  should_exit = .FALSE.
                  ! Update alpha using Armijo's scheme
                  alpha = alpha*factor
               END IF
               IF (continue_ls_exit) THEN
                  ! Continuation of line search did not yield a better alpha, use previously located solution and exit
                  alpha = alpha_ls
                  DO ispin = 1, nspins
                     CALL deallocate_mo_set(mos(ispin))
                     CALL duplicate_mo_set(mos(ispin), mos_ls(ispin))
                     CALL calculate_density_matrix(mos(ispin), &
                                                   p_rmpv(ispin)%matrix)
                     CALL deallocate_mo_set(mos_ls(ispin))
                  END DO
                  CALL qs_rho_update_rho(rho, qs_env=qs_env)
                  CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
                  DEALLOCATE (mos_ls)
                  should_exit = .TRUE.
               END IF
               ! Reached max steps and SCF converged: continue with last iterated step size
               IF (.NOT. should_exit .AND. &
                   (i == max_linesearch .AND. converged .AND. .NOT. found_solution)) THEN
                  should_exit = .TRUE.
                  reached_maxls = .TRUE.
                  alpha = alpha*(1.0_dp/factor)
               END IF
               ! Reset outer SCF environment to last converged state
               scf_env%outer_scf%variables(:, iter_count + 1) = 0.0_dp
               scf_env%outer_scf%gradient = gradient
               scf_env%outer_scf%energy = energy
               ! Exit line search if a suitable step size was found
               IF (should_exit) EXIT
               ! Reset the electronic structure
               cdft_control%total_steps = nsteps
               DO ispin = 1, nspins
                  CALL deallocate_mo_set(mos(ispin))
                  CALL duplicate_mo_set(mos(ispin), mos_stashed(ispin))
                  CALL calculate_density_matrix(mos(ispin), &
                                                p_rmpv(ispin)%matrix)
               END DO
               CALL qs_rho_update_rho(rho, qs_env=qs_env)
               CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
            END DO
            scf_control%outer_scf%cdft_opt_control%newton_step = alpha
            IF (.NOT. should_exit) THEN
               CALL cp_warn(__LOCATION__, &
                            "Line search did not converge. CDFT SCF proceeds with fixed step size.")
               scf_control%outer_scf%cdft_opt_control%newton_step = scf_control%outer_scf%cdft_opt_control%newton_step_save
            END IF
            IF (reached_maxls) &
               CALL cp_warn(__LOCATION__, &
                            "Line search did not converge. CDFT SCF proceeds with lasted iterated step size.")
            CALL cp_rm_default_logger()
            CALL cp_logger_release(tmp_logger)
            ! Release temporary storage
            DO ispin = 1, nspins
               CALL deallocate_mo_set(mos_stashed(ispin))
            END DO
            DEALLOCATE (mos_stashed, gradient, energy)
            IF (ALLOCATED(positive_sign)) DEALLOCATE (positive_sign)
            IF (output_unit > 0) THEN
               WRITE (output_unit, FMT="(/,A)") &
                  " ================================== LINE SEARCH COMPLETE =================================="
               CALL close_file(unit_number=output_unit)
            END IF
         END BLOCK
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_cdft_line_search


! **************************************************************************************************
!> \brief debug trace (env GCE_DBG): electron count carried by G = 0 of rho_g and of the mixer input
! **************************************************************************************************
   SUBROUTINE gce_dbg_g0(rho, para_env, output_unit, tag, scf_env)
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: output_unit
      CHARACTER(LEN=*), INTENT(IN)                       :: tag
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      REAL(KIND=dp)                                      :: q(2)
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: g

      NULLIFY (g)
      CALL qs_rho_get(rho, rho_g=g)
      q = 0.0_dp
      IF (g(1)%pw_grid%have_g0) THEN
         q(1) = REAL(g(1)%array(1), KIND=dp)*g(1)%pw_grid%vol
         IF (ASSOCIATED(scf_env%mixing_store)) THEN
            IF (ASSOCIATED(scf_env%mixing_store%rhoin)) THEN
               IF (ASSOCIATED(scf_env%mixing_store%rhoin(1)%cc)) &
                  q(2) = REAL(scf_env%mixing_store%rhoin(1)%cc(1), KIND=dp)*g(1)%pw_grid%vol
            END IF
         END IF
      END IF
      CALL para_env%sum(q)
      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(T3,A,A,A,F16.8,A,F16.8)") &
         "GCE| dbg ", tag, " Ne(rho_g,G0) ", q(1), "  Ne(rhoin,G0) ", q(2)
   END SUBROUTINE gce_dbg_g0


! **************************************************************************************************
!> \brief GCE CAPACITANCE: the total charge is set by the occupations, not by the mixer: set the G = 0
!>        coefficient of the (mixed) density to N_e(target)/V (spin-restricted only). With PIN_MIXER_G0 the
!>        G = 0 coefficient of the mixer input is set too, so the Broyden/Pulay residual has no charge part.
! **************************************************************************************************
   SUBROUTINE gce_pin_g0(rho, scf_env, scf_control)
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      REAL(KIND=dp), DIMENSION(:), POINTER               :: gce_tot_rho_r
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: gce_rho_g
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: gce_rho_r

      NULLIFY (gce_rho_g, gce_rho_r, gce_tot_rho_r)
      CALL qs_rho_get(rho, rho_g=gce_rho_g, rho_r=gce_rho_r, tot_rho_r=gce_tot_rho_r)
      IF (SIZE(gce_rho_g) /= 1) RETURN
      IF (gce_rho_g(1)%pw_grid%have_g0) &
         gce_rho_g(1)%array(1) = CMPLX(scf_control%gce%n_target/gce_rho_g(1)%pw_grid%vol, 0.0_dp, KIND=dp)
      CALL pw_transfer(gce_rho_g(1), gce_rho_r(1))
      IF (ASSOCIATED(gce_tot_rho_r)) gce_tot_rho_r(1) = pw_integrate_function(gce_rho_r(1), isign=-1)
      IF (scf_control%gce%pin_mixer_g0 .AND. ASSOCIATED(scf_env%mixing_store)) THEN
         IF (ASSOCIATED(scf_env%mixing_store%rhoin)) THEN
            IF (ASSOCIATED(scf_env%mixing_store%rhoin(1)%cc) .AND. gce_rho_g(1)%pw_grid%have_g0) &
               scf_env%mixing_store%rhoin(1)%cc(1) = CMPLX(scf_control%gce%n_target/gce_rho_g(1)%pw_grid%vol, 0.0_dp, KIND=dp)
         END IF
      END IF
   END SUBROUTINE gce_pin_g0


! **************************************************************************************************
!> \brief RESPONSE_GUESS: start density of the inner SCF after an N_e update, rho + dN_e (d rho/d N_e).
!>        Called after the update with rho still the density at n_old. Two states of the same force evaluation
!>        (same geometry) give the difference quotient f = (rho - rho_stored)/(n_old - n_stored); its G = 0
!>        coefficient is 1/V because the G = 0 coefficient of every state is N_e/V. f is kept for later
!>        updates and force evaluations. The density, its real-space copy, the electron count and the mixer
!>        input are set to rho + (n_target - n_old) f, so the mixer residual has no part from the change of N_e.
!>        done = F (no f yet, more than one spin, or a mixer without an input copy of the same size): the
!>        caller falls back to the G = 0 shift.
! **************************************************************************************************
   SUBROUTINE gce_response_guess(rho, scf_env, scf_control, n_old, done, output_unit)
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      REAL(KIND=dp), INTENT(IN)                          :: n_old
      LOGICAL, INTENT(OUT)                               :: done
      INTEGER, INTENT(IN)                                :: output_unit

      REAL(KIND=dp), PARAMETER                           :: dn_min = 1.0E-3_dp   ! smallest |dN_e| of a pair [e]

      INTEGER                                            :: ng
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(gce_type), POINTER                            :: gce
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r

      done = .FALSE.
      NULLIFY (rho_g, rho_r, tot_rho_r)
      CALL qs_rho_get(rho, rho_g=rho_g, rho_r=rho_r, tot_rho_r=tot_rho_r)
      IF (SIZE(rho_g) /= 1) RETURN
      ng = SIZE(rho_g(1)%array)
      gce => scf_control%gce
      IF (ALLOCATED(gce%resp_rho)) THEN
         IF (SIZE(gce%resp_rho) /= ng) THEN   ! the grid changed
            DEALLOCATE (gce%resp_rho)
            IF (ALLOCATED(gce%resp_f)) DEALLOCATE (gce%resp_f)
            gce%resp_have_rho = .FALSE.
            gce%resp_have_f = .FALSE.
         END IF
      END IF
      IF (.NOT. ALLOCATED(gce%resp_rho)) ALLOCATE (gce%resp_rho(ng))
      IF (gce%resp_have_rho .AND. gce%resp_step == gce%ion_step .AND. ABS(n_old - gce%resp_n) >= dn_min) THEN
         IF (.NOT. ALLOCATED(gce%resp_f)) ALLOCATE (gce%resp_f(ng))
         gce%resp_f(:) = (rho_g(1)%array(:) - gce%resp_rho(:))/(n_old - gce%resp_n)
         gce%resp_have_f = .TRUE.
      END IF
      gce%resp_rho(:) = rho_g(1)%array(:)
      gce%resp_n = n_old
      gce%resp_step = gce%ion_step
      gce%resp_have_rho = .TRUE.
      IF (.NOT. gce%resp_have_f) RETURN
      IF (.NOT. ASSOCIATED(scf_env%mixing_store)) RETURN
      IF (.NOT. ASSOCIATED(scf_env%mixing_store%rhoin)) RETURN
      IF (.NOT. ASSOCIATED(scf_env%mixing_store%rhoin(1)%cc)) RETURN
      IF (SIZE(scf_env%mixing_store%rhoin(1)%cc) /= ng) RETURN
      rho_g(1)%array(:) = rho_g(1)%array(:) + (gce%n_target - n_old)*gce%resp_f(:)
      CALL pw_transfer(rho_g(1), rho_r(1))
      IF (ASSOCIATED(tot_rho_r)) tot_rho_r(1) = pw_integrate_function(rho_r(1), isign=-1)
      scf_env%mixing_store%rhoin(1)%cc(:) = rho_g(1)%array(:)
      done = .TRUE.
      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(T3,A,F12.6,A)") &
         "GCE| response guess: density moved along d rho/d N_e by ", gce%n_target - n_old, " e"
   END SUBROUTINE gce_response_guess

! **************************************************************************************************
!> \brief What the density mixer does after an N_e update. Default: restart the history (ncall = 0).
!>        KEEP_MIXER_BELOW with PIN_MIXER_G0 and |dN_e| small: the history is kept and only the first step
!>        after the update is a plain Kerker step, so no (dRho, dF) pair spans the change of the fixed-point map.
! **************************************************************************************************
   SUBROUTINE gce_mixer_after_update(scf_env, scf_control)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      IF (.NOT. ASSOCIATED(scf_env%mixing_store)) RETURN
      IF (.NOT. gce_keep_mixer(scf_control%gce)) THEN
         scf_env%mixing_store%ncall = 0
      ELSE IF (scf_control%gce%pin_mixer_g0 .AND. scf_env%mixing_store%ncall > 1) THEN
         scf_control%gce%ncall_restore = scf_env%mixing_store%ncall
         scf_env%mixing_store%ncall = 0
      END IF
   END SUBROUTINE gce_mixer_after_update

! **************************************************************************************************
!> \brief is the force evaluation part of a geometry optimisation? (the GEO_OPT iteration level is set)
! **************************************************************************************************
!> \brief true when the GCE has switched the solvent off for this SCF phase (SOLVENT_ONLY_FOR_CHARGE); without &SCCS
!>        there is no sccs_control to ask
!> \param dft_control ...
!> \return ...
! **************************************************************************************************
   FUNCTION gce_solvent_is_off(dft_control) RESULT(off)
      TYPE(dft_control_type), POINTER                    :: dft_control
      LOGICAL                                            :: off

      off = .FALSE.
      IF (dft_control%do_sccs) off = dft_control%sccs_control%force_off

   END FUNCTION gce_solvent_is_off

! **************************************************************************************************
   FUNCTION gce_in_geo_opt() RESULT(in_geo_opt)
      LOGICAL                                            :: in_geo_opt

      INTEGER                                            :: i
      TYPE(cp_logger_type), POINTER                      :: logger

      in_geo_opt = .FALSE.
      logger => cp_get_default_logger()
      IF (.NOT. ASSOCIATED(logger%iter_info)) RETURN
      DO i = 1, logger%iter_info%n_rlevel
         IF (TRIM(logger%iter_info%level_name(i)) == "GEO_OPT") in_geo_opt = .TRUE.
      END DO
   END FUNCTION gce_in_geo_opt

! **************************************************************************************************
!> \brief ENFORCE_TARGET: an update step that has used all MAX_OUTER_UPDATES updates and is still off the target by
!>        more than FERMI_TOL ends the run with the reason. Called on a converged inner SCF with the solvent on.
!>        Steps with a reduced budget (CONVERGENCE_MODE 1 / 2 after the first step) and steps without an update are
!>        not meant to reach the target and are not checked.
!> \param gce ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE gce_check_target_reached(gce, output_unit)
      TYPE(gce_type), POINTER                            :: gce
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(LEN=600)                                 :: msg
      REAL(KIND=dp)                                      :: dwf, hi, lo, n_hi, n_lo

      IF (.NOT. gce%enforce_target) RETURN
      IF (gce%max_outer <= 0 .OR. gce%max_outer_now < gce%max_outer) RETURN
      IF (.NOT. (gce%ion_step == 0 .OR. MOD(gce%ion_step + 1, gce%update_interval) == 0)) RETURN
      IF (gce%n_outer < gce%max_outer_now) RETURN
      dwf = gce%prev_workfunction - gce%target_workfunction
      IF (ABS(dwf) <= gce%fermi_tol) RETURN
      ! extremes of W_f - target over the updates of this step, including the value at the present N_e
      lo = gce%wf_lo; hi = gce%wf_hi; n_lo = gce%n_at_lo; n_hi = gce%n_at_hi
      IF (dwf < lo) THEN
         lo = dwf; n_lo = gce%n_target
      END IF
      IF (dwf > hi) THEN
         hi = dwf; n_hi = gce%n_target
      END IF
      IF (lo < 0.0_dp .AND. hi > 0.0_dp) THEN
         WRITE (msg, "(A,I0,A,F9.4,A,F9.4,A,F9.4,A,F9.4,A)") &
            "GCE ENFORCE_TARGET: the target is not reached after ", gce%n_outer, " updates. W_f - target was ", &
            lo*evolt, " eV at N_e ", n_lo, " e and ", hi*evolt, " eV at N_e ", n_hi, " e: it crosses the target without "// &
            "hitting it, so the target lies in a gap or in a region with too few states for this smearing "// &
            "(raise ELECTRONIC_TEMPERATURE or ADDED_MOS, or choose another potential)."
      ELSE
         WRITE (msg, "(A,I0,A,F9.4,A,F9.4,A,F9.4,A,F9.4,A,F8.4,A)") &
            "GCE ENFORCE_TARGET: the target is not reached after ", gce%n_outer, " updates. W_f - target stays on one side, ", &
            dwf*evolt, " eV at N_e ", gce%n_target, " e (range ", lo*evolt, " to ", hi*evolt, " eV, capacitance ", &
            gce%capacitance/evolt, " e/V): N_e moves too slowly (NE_STEP_MAX, MAX_OUTER_UPDATES) or the target is beyond "// &
            "what N_e can reach."
      END IF
      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(T3,A)") TRIM(msg)
      CPABORT(TRIM(msg))
   END SUBROUTINE gce_check_target_reached

END MODULE qs_scf
