!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Self-consistent continuum solvation (SCCS) model implementation
!> \author  Matthias Krack (MK)
!> \version 1.0
!> \par Literature:
!>      - J.-L. Fattebert and F. Gygi,
!>        Density functional theory for efficient ab initio molecular dynamics
!>        simulations in solution, J. Comput. Chem. 23, 662-666 (2002)
!>      - O. Andreussi, I. Dabo, and N. Marzari,
!>        Revised self-consistent continuum solvation in electronic-structure
!>        calculations, J. Chem. Phys. 136, 064102-20 (2012)
!> \par History:
!>      - Creation (10.10.2013,MK)
!>      - Derivatives using finite differences (26.11.2013,MK)
!>      - Cube file dump of the dielectric function (19.12.2013,MK)
!>      - Cube file dump of the polarisation potential (20.12.2013,MK)
!>      - Calculation of volume and surface of the cavity (21.12.2013,MK)
!>      - Functional derivative of the cavitation energy (28.12.2013,MK)
!>      - Update printout (11.11.2022,MK)
! **************************************************************************************************

MODULE qs_sccs

   USE cp_control_types,                ONLY: dft_control_type,&
                                              sccs_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE machine,                         ONLY: m_walltime
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr,&
                                              low_print_level
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE cp_realspace_grid_init,          ONLY: init_input_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE input_constants,                 ONLY: sccs_andreussi,&
                                              sccs_derivative_cd3,&
                                              sccs_derivative_cd5,&
                                              sccs_derivative_cd7,&
                                              sccs_derivative_fft,&
                                              sccs_fattebert_gygi,&
                                              sccs_saa_andreussi,&
                                              sccs_solver_newton_cg,&
                                              sccs_isodensity,&
                                              sccs_solvent_excluded
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_get_lval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              int_8
   USE mathconstants,                   ONLY: fourpi,&
                                              pi,&
                                              twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_list_types,             ONLY: particle_list_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: HALFSPACE
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_derive,&
                                              pw_integral_ab,&
                                              pw_integrate_function,&
                                              pw_laplace,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_func_u_convolution,&
                                              pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_analytic,&
                                              pw_poisson_mt,&
                                              pw_poisson_periodic,&
                                              pw_poisson_type
   USE qs_sccs_pw,                      ONLY: sol_poisson,&
                                              sol_prolong,&
                                              sol_restrict
   USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_rho_core
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_sccs_cavity,                  ONLY: cavity_build,&
                                              cavity_params_type,&
                                              cavity_potential_sum,&
                                              cavity_release,&
                                              cavity_state_type
   USE qs_sccs_nldiel,                  ONLY: nldiel_setup
   USE qs_sccs_nlpb,                    ONLY: nlpb_ion_density,&
                                              nlpb_last_psi_q_min,&
                                              nlpb_params_type,&
                                              nlpb_solve
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE realspace_grid_types,            ONLY: realspace_grid_desc_type,&
                                              realspace_grid_input_type,&
                                              realspace_grid_type,&
                                              rs_grid_create,&
                                              rs_grid_create_descriptor,&
                                              rs_grid_release,&
                                              rs_grid_release_descriptor
   USE rs_methods,                      ONLY: derive_fdm_cd3,&
                                              derive_fdm_cd5,&
                                              derive_fdm_cd7
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_sccs'

   PUBLIC :: print_sccs_results, sccs, sccs_pseudo_core_kinds

   REAL(KIND=dp), SAVE :: t_sccs_cum = 0.0_dp   ! cumulative solvation wall time [s]
   INTEGER, SAVE :: n_nlpb_calls = 0             ! SOLVER_STRIDE: activated NLPB calls so far
   INTEGER, SAVE :: n_nlpb_fail = 0              ! consecutive SCCS calls whose NLPB solve did not converge

   ! Masks that do not depend on the density. The REGION window depends on the grid and its parameters; the
   ! EXCLUSION masks on the grid, the excluded atoms' positions, radii and width. Each is rebuilt whenever its key (kept
   ! as an exact copy, compared with ==) differs from the previous call, so ionic steps, cell changes and input changes
   ! invalidate it. Plain arrays: the key and the data never reference a pw pool.
   TYPE mask_cache_type
      LOGICAL :: valid = .FALSE.
      INTEGER :: id_nr = -1, nex = 0, iax = 0
      INTEGER, DIMENSION(2, 3) :: bounds = 0, bounds_local = 0
      REAL(KIND=dp), DIMENSION(3, 3) :: dh = 0.0_dp
      REAL(KIND=dp), DIMENSION(4) :: par = 0.0_dp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: pos
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: rd, ri
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :) :: m1, m2
   END TYPE mask_cache_type
   TYPE(mask_cache_type), SAVE :: region_cache, excl_cache, flex_cache

CONTAINS

! **************************************************************************************************
!> \brief   Self-consistent continuum solvation (SCCS) model implementation
!> \param qs_env ...
!> \param rho_tot_gspace ...
!> \param v_hartree_gspace ...
!> \param v_sccs ...
!> \param h_stress ...
!> \par History:
!>      - Creation (10.10.2013,MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! **************************************************************************************************

   SUBROUTINE sccs(qs_env, rho_tot_gspace_io, v_hartree_gspace_out, v_sccs_out, h_stress, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rho_tot_gspace_io, v_hartree_gspace_out
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: v_sccs_out
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT), &
         OPTIONAL                                        :: h_stress
      LOGICAL, INTENT(IN), OPTIONAL                      :: calculate_forces   ! forces are requested in this call

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'sccs'
      REAL(KIND=dp), PARAMETER                           :: epstol = 1.0E-8_dp

      CHARACTER(LEN=4*default_string_length)             :: message, my_pos_cube
      CHARACTER(LEN=default_path_length)                 :: cube_path, filename, mpi_filename, &
                                                            print_path
      INTEGER                                            :: cube_unit, handle, i, ig, ispin, iter, j, k, &
                                                            l, ncg_nl, niter_nl, nspin, output_unit, print_level, iw
      INTEGER(KIND=int_8)                                :: ngpts
      INTEGER, DIMENSION(3)                              :: lb, ub
      LOGICAL                                            :: use_psi, append_cube, calculate_stress_tensor, &
                                                            conv_nl, do_electrolyte, do_kpoints, &
                                                            do_newton, do_nlpb, do_ppcav, mpi_io, &
                                                            should_output
      REAL(KIND=dp) :: alpha_zeta, cavity_surface, cavity_volume, cell_volume, delta_eta, &
         delta_zeta, dphi2, dvol, e_ion, e_tot, epsilon_solvent, f, f0, f_ion, kappa0sq, &
         a_nl, d_ion, e_diel_nl, e_ion_nl, e_ion_phi, q_ion_nl, rms_nl, theta_b, &
         polarisation_charge, q_ion, r, R_solv, rho_delta, rho_delta_avg, rho_delta_max, &
         rho_iter_new, tot_rho_elec, tot_rho_solute, zeta
      REAL(KIND=dp), DIMENSION(3)                        :: abc, lxyz, uxyz
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool, fine_pool
      ! solver-grid locals; the dummies above stay on the fine (auxbas) grid
      TYPE(pw_c1d_gs_type)                               :: rho_tot_gspace, v_hartree_gspace, phi_vac_f
      TYPE(pw_r3d_rs_type)                               :: v_sccs
      TYPE(pw_r3d_rs_type), ALLOCATABLE, DIMENSION(:), TARGET :: rho_r_c
      TYPE(pw_r3d_rs_type), TARGET                       :: phi_c
      TYPE(pw_r3d_rs_type)                               :: rmask   ! &SCCS&REGION mask m(z)
      REAL(KIND=dp)                                      :: zpos, zlen, wsq2, rg_in, rg_tot, rg_in2, rg_tot2, mval, zimg
      INTEGER                                            :: iimg
      INTEGER                                            :: iax, idx
      INTEGER, DIMENSION(3)                              :: n3
      LOGICAL                                            :: region_on, region_apply, region_hit, n_e_moves
      ! diagnostics (polarisation energy terms) and failure handling
      REAL(KIND=dp)                                      :: eps_min_diag, rho_phi_diag, e_grad_diag, e_grad_raw, phi_avg_diag
      REAL(KIND=dp)                                      :: q_net_diag
      LOGICAL                                            :: pol_neutral
      LOGICAL, SAVE                                      :: warned_charged_pol = .FALSE.
      CHARACTER(LEN=24)                                  :: region_tag
      ! &SCCS&EXCLUSION, atom-following exclusion mask
      LOGICAL                                            :: excl_on, calc_forces
      INTEGER                                            :: nex, iex
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ex_atom, ex_ikind, ex_iat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ex_rd, ex_ri
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: fex
      TYPE(pw_r3d_rs_type)                               :: excl_md, excl_mi, excl_sd0, excl_si0
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force_ex
      ! &SCCS&REGION&FLEXIBLE, region mask from the coarse-grained density of the listed atoms
      LOGICAL                                            :: flex_on, flex_frc, flex_cav, region_planar
      INTEGER                                            :: nfx
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: fx_atom, fx_ikind, fx_iat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: fx_w, fx_w2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ffx
      TYPE(pw_r3d_rs_type)                               :: flex_m, flex_dm, flex_sd0, flex_si0, flex_b
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_pw_r_sccs_fine
      REAL(KIND=dp)                                      :: e_h_c, e_h_fine
      INTEGER                                            :: level
      LOGICAL, SAVE                                      :: warned_se_diel = .FALSE.
      TYPE(nlpb_params_type)                             :: nlpb
      TYPE(cavity_params_type)                           :: cpar
      TYPE(cavity_state_type)                            :: cav
      TYPE(pw_r3d_rs_type)                               :: acav_fld, afld, s_diel_nl, vchain
      TYPE(pw_r3d_rs_type), DIMENSION(3)                 :: dsc
      TYPE(pw_c1d_gs_type)                               :: greenfn_scr
      TYPE(pw_r3d_rs_type)                               :: adiel_unit, aion_unit, geff_unit, d_s_rhoel, deps_elec, &
                                                            dkap2e, &
                                                            dtf, eps_elec, eps_nl, ff, kap2, &
                                                            p_e_interface, s, s_ion, t, &
                                                            tem_convolution, tem_func, u
      TYPE(pw_r3d_rs_type), DIMENSION(3)                 :: dln_eps_elec, dphi_tot
      ! fields kept for the strain derivative of the local-cavity functional
      TYPE(pw_r3d_rs_type)                               :: st_eps, st_neps
      REAL(KIND=dp), DIMENSION(3, 3)                     :: sig_cav, sig_el
      REAL(KIND=dp)                                      :: st_e2, st_f
      INTEGER                                            :: st_a, st_b
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_pw_r
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_pw_r_sccs
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: pscore_alpha, pscore_charge
      REAL(KIND=dp)                                      :: tot_rho_pscore
      REAL(KIND=dp)                                      :: t_sccs, t_start
      ! timers; surface/density-gradient skips for SOLVENT_EXCLUDED; polarisation-charge transform on demand
      INTEGER                                            :: h_cav, h_kspot, h_nlpb, h_surf
      LOGICAL                                            :: need_drho, need_rho_tot, se_skip_surf
      TYPE(pw_r3d_rs_type)                               :: aion_f
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)
      t_start = m_walltime()

      NULLIFY (auxbas_pw_pool)
      NULLIFY (cp_subsys)
      NULLIFY (dft_control)
      NULLIFY (energy)
      NULLIFY (input)
      NULLIFY (logger)
      NULLIFY (para_env)
      NULLIFY (particles)
      NULLIFY (poisson_env)
      NULLIFY (pw_env)
      NULLIFY (pw_pools)
      NULLIFY (rho)
      NULLIFY (sccs_control)
      NULLIFY (scf_env)

      ! Load data from Quickstep environment
      CALL get_qs_env(qs_env=qs_env, &
                      cp_subsys=cp_subsys, &
                      do_kpoints=do_kpoints, &
                      dft_control=dft_control, &
                      energy=energy, &
                      input=input, &
                      para_env=para_env, &
                      pw_env=pw_env, &
                      rho=rho, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      scf_env=scf_env)
      CALL cp_subsys_get(cp_subsys, particles=particles)

      sccs_control => dft_control%sccs_control

      CPASSERT(ASSOCIATED(qs_env))

      excl_on = .FALSE.
      nex = 0
      calc_forces = .FALSE.
      IF (PRESENT(calculate_forces)) calc_forces = calculate_forces

      IF (do_kpoints) THEN
         CPWARN("SCCS with k-points has not yet been fully validated")
      END IF

      ! Implicit electrolyte (linearised Poisson-Boltzmann ionic screening)
      do_electrolyte = sccs_control%do_electrolyte
      kappa0sq = 0.0_dp
      e_ion = 0.0_dp
      q_ion = 0.0_dp
      eps_min_diag = 1.0_dp
      IF (do_electrolyte) THEN
         kappa0sq = 1.0_dp/sccs_control%debye_length**2
         IF (sccs_control%method_id == sccs_saa_andreussi) &
            CPABORT("SCCS DEBYE_LENGTH is not implemented for METHOD SAA_ANDREUSSI")
         IF (ABS(sccs_control%epsilon_solvent - 1.0_dp) <= epstol) &
            CPABORT("SCCS DEBYE_LENGTH requires RELATIVE_PERMITTIVITY > 1")
         IF (PRESENT(h_stress)) THEN
            CPABORT("SCCS: no analytic stress with DEBYE_LENGTH, use NUMERICAL")
         END IF
      END IF
      do_ppcav = (sccs_control%method_id == sccs_solvent_excluded)
      ! analytic stress: the strain derivative exists for the LOCAL ISODENSITY cavity with a linear dielectric
      ! (Maxwell stress with eps(n) + cavitation/dispersion/repulsion); everything else aborts.
      IF (PRESENT(h_stress)) THEN
         ! the analytic stress is not released: on the ISODENSITY path it carries an isotropic offset of about E_pol/2
         ! per diagonal element. The checks below keep the intended scope of the implemented strain derivative.
         CPABORT("SCCS: analytic stress not validated, use NUMERICAL (STRESS_TENSOR NUMERICAL)")
         ! every abort names the working alternative (STRESS_TENSOR NUMERICAL)
         IF (do_ppcav .OR. sccs_control%pseudo_core) THEN
            CPABORT("SCCS: no analytic stress for SE/ADD_CORE_DENSITY, use NUMERICAL")
         END IF
         IF (sccs_control%method_id /= sccs_isodensity) THEN
            CPABORT("SCCS: analytic stress only for ISODENSITY, use NUMERICAL")
         END IF
         IF (sccs_control%nldiel .OR. sccs_control%psi_form .OR. sccs_control%region_on .OR. sccs_control%excl_on) THEN
            CPABORT("SCCS: no analytic stress for NLD/PSI/REGION/EXCL, use NUMERICAL")
         END IF
      END IF
      IF (do_ppcav) THEN
         cpar%nc = sccs_control%rho_cavity
         cpar%sigma = sccs_control%switch_width
         cpar%a = sccs_control%cav_a
         cpar%radius(1) = sccs_control%r_cav_solv
         cpar%radius(2) = sccs_control%r_cav_ion
         cpar%radius(3) = sccs_control%r_cav_diel
         cpar%radius(4) = sccs_control%r_cav_cav
      END IF
      ! R_diel = 0 bypasses the dielectric kernel, S_diel = S_solv, i.e. the dielectric starts at the
      ! solvent-accessible surface R_solv (1.4 A) outside the solute surface, not at the solvent-excluded one.
      ! A linear dielectric with R_diel >= 1.0 has no saturation and the SCF may run away.
      IF (do_ppcav .AND. (.NOT. warned_se_diel)) THEN
         IF (sccs_control%r_cav_diel <= 0.0_dp) THEN
            CPWARN("SCCS SE: DIELECTRIC_RADIUS 0 puts eps on the SAS (S_diel = S_solv)")
            warned_se_diel = .TRUE.
         ELSE IF (.NOT. sccs_control%nldiel) THEN
            CPWARN("SCCS SE + linear dielectric: SCF may run away at R_diel >= 1 A")
            warned_se_diel = .TRUE.
         END IF
      END IF
      do_newton = (sccs_control%solver_id == sccs_solver_newton_cg)
      do_nlpb = do_newton .OR. (do_electrolyte .AND. sccs_control%nlpb)
      IF (sccs_control%nlpb .AND. .NOT. do_electrolyte) &
         CPABORT("SCCS NONLINEAR_IONS requires DEBYE_LENGTH > 0")
      IF (do_newton .AND. (ABS(sccs_control%epsilon_solvent - 1.0_dp) <= epstol)) &
         CPABORT("SCCS SOLVER NEWTON_CG needs RELATIVE_PERMITTIVITY > 1")
      IF (do_newton .AND. (sccs_control%method_id == sccs_saa_andreussi)) &
         CPABORT("SCCS SOLVER NEWTON_CG needs METHOD /= SAA_ANDREUSSI")
      IF (do_nlpb) THEN
         ! no electrolyte -> pure dielectric Poisson problem; electrolyte without NONLINEAR_IONS -> linear PB
         nlpb%ions = do_electrolyte
         nlpb%linear_ions = do_electrolyte .AND. (.NOT. sccs_control%nlpb)
         nlpb%kappa0sq = kappa0sq
         nlpb%eps_b = sccs_control%epsilon_solvent
         nlpb%beta = 1.0_dp/sccs_control%temperature   ! TEMPERATURE is stored in internal units (k_B T in Hartree)
         nlpb%zion = sccs_control%zion
         theta_b = 0.0_dp
         IF (nlpb%ions .AND. sccs_control%nlpb) THEN
            ! theta_b = 2 c0/n_max with 2 c0 = eps_b kappa0^2/(4 pi beta z^2) and n_max = 1/d_ion^3,
            ! d_ion = ION_PACKING_DIAMETER or the close-packed 2^(5/6) R_ion (Islam et al. 2023)
            d_ion = sccs_control%ion_packing_diameter
            IF (d_ion <= 0.0_dp) d_ion = 2.0_dp**(5.0_dp/6.0_dp)*sccs_control%ion_radius
            theta_b = nlpb%eps_b*kappa0sq*d_ion**3/(fourpi*nlpb%beta*nlpb%zion**2)
            IF (theta_b >= 1.0_dp) &
               CPABORT("SCCS NLPB: theta_b >= 1 (ion packing diameter too large for c)")
         END IF
         nlpb%theta_b = theta_b
         nlpb%eps_resid = sccs_control%eps_nlpb
         nlpb%max_iter = sccs_control%max_iter_nlpb
         nlpb%nldiel = sccs_control%nldiel
         nlpb%a_smooth = sccs_control%a_smooth
         nlpb%dev_r(:) = sccs_control%dev_r(:)
         nlpb%dev_i(:) = sccs_control%dev_i(:)
         nlpb%dev_l(:) = sccs_control%dev_l(:)
         IF (sccs_control%nldiel) THEN
            nlpb%dl%eps_b = sccs_control%epsilon_solvent
            nlpb%dl%eps_inf = sccs_control%eps_inf
            nlpb%dl%n_mol = sccs_control%n_mol
            nlpb%dl%p_mol = sccs_control%p_mol
            nlpb%dl%beta = nlpb%beta
            CALL nldiel_setup(nlpb%dl)
         END IF
      END IF
      IF (sccs_control%nldiel .AND. .NOT. do_nlpb) &
         CPABORT("SCCS NONLINEAR_DIELECTRIC requires NONLINEAR_IONS and DEBYE_LENGTH > 0")

      sig_cav = 0.0_dp
      IF (PRESENT(h_stress)) THEN
         calculate_stress_tensor = .TRUE.
         h_stress(:, :) = 0.0_dp
      ELSE
         calculate_stress_tensor = .FALSE.
      END IF

      ! Get access to the PW grid pool
      CALL pw_env_get(pw_env, &
                      auxbas_pw_pool=auxbas_pw_pool, &
                      pw_pools=pw_pools, &
                      poisson_env=poisson_env)

      IF (do_electrolyte .OR. do_nlpb) THEN
         IF (poisson_env%parameters%solver /= pw_poisson_periodic) &
            CPABORT("SCCS DEBYE_LENGTH/NEWTON_CG need POISSON_SOLVER PERIODIC")
      END IF

      CALL pw_zero(v_sccs_out)

      ! Calculate no SCCS contribution, if the requested SCF convergence threshold is not reached yet
      ! SOLVENT_ONLY_FOR_CHARGE: force_off = the GCE has switched the solvent off for this SCF phase
      IF ((.NOT. sccs_control%sccs_activated) .OR. sccs_control%force_off) THEN
         IF (sccs_control%eps_scf > 0.0_dp .OR. sccs_control%force_off) THEN
            IF (sccs_control%force_off .OR. (scf_env%iter_delta > sccs_control%eps_scf) .OR. &
                ((qs_env%scf_env%outer_scf%iter_count == 0) .AND. &
                 (qs_env%scf_env%iter_count <= 1))) THEN
               IF (calculate_stress_tensor) THEN
                  ! Request also the calculation of the stress tensor contribution
                  CALL pw_poisson_solve(poisson_env=poisson_env, &
                                        density=rho_tot_gspace_io, &
                                        ehartree=energy%hartree, &
                                        vhartree=v_hartree_gspace_out, &
                                        h_stress=h_stress)
               ELSE
                  CALL pw_poisson_solve(poisson_env=poisson_env, &
                                        density=rho_tot_gspace_io, &
                                        ehartree=energy%hartree, &
                                        vhartree=v_hartree_gspace_out)
               END IF
               energy%sccs_pol = 0.0_dp
               energy%sccs_cav = 0.0_dp
               energy%sccs_dis = 0.0_dp
               energy%sccs_rep = 0.0_dp
               energy%sccs_sol = 0.0_dp
               energy%sccs_hartree = energy%hartree
               CALL timestop(handle)
               RETURN
            END IF
         END IF
         sccs_control%sccs_activated = .TRUE.
      END IF

      ! optional coarser solver grid; every SCCS field below lives on pw_pools(level)
      level = MAX(1, sccs_control%solver_grid_level)
      fine_pool => auxbas_pw_pool
      IF (level > 1) THEN
         IF (.NOT. do_nlpb) CPABORT("SCCS SOLVER_GRID_LEVEL > 1 requires SOLVER NEWTON_CG")
         IF (calculate_stress_tensor) CPABORT("SCCS SOLVER_GRID_LEVEL > 1: stress tensor not available")
         IF (level > SIZE(pw_pools)) CPABORT("SCCS SOLVER_GRID_LEVEL exceeds MGRID NGRIDS")
         auxbas_pw_pool => pw_pools(level)%pool
      END IF
      CALL auxbas_pw_pool%create_pw(rho_tot_gspace)
      CALL pw_zero(rho_tot_gspace)
      CALL pw_copy(rho_tot_gspace_io, rho_tot_gspace)
      CALL auxbas_pw_pool%create_pw(v_hartree_gspace)
      CALL pw_zero(v_hartree_gspace)
      CALL auxbas_pw_pool%create_pw(v_sccs)
      CALL pw_zero(v_sccs)
      use_psi = .FALSE.   ! set after the NLPB solve
      e_h_fine = 0.0_dp
      IF (level > 1) THEN
         ! exact fine-grid vacuum Hartree potential and energy of the solute; the coarse solve supplies only the
         ! polarisation / ionic response and the solvation free energy
         CALL fine_pool%create_pw(phi_vac_f)
         CALL pw_poisson_solve(poisson_env=poisson_env, density=rho_tot_gspace_io, ehartree=e_h_fine, &
                               vhartree=phi_vac_f)
      END IF

      nspin = dft_control%nspins

      ! Manage print output control
      logger => cp_get_default_logger()
      print_level = logger%iter_info%print_level
      print_path = "DFT%PRINT%SCCS"
      should_output = (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                                        TRIM(print_path)), cp_p_file))
      output_unit = cp_print_key_unit_nr(logger, input, TRIM(print_path), &
                                         extension=".sccs", &
                                         ignore_should_output=should_output, &
                                         log_filename=.FALSE.)
      ! Under NEWTON_CG the real-space polarisation charge is only printed or written to a cube
      need_rho_tot = should_output .OR. &
                     BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                                      TRIM(print_path)//"%POLARISATION_CHARGE_DENSITY"), cp_p_file)

      ! Get rho
      CALL qs_rho_get(rho_struct=rho, &
                      rho_r=rho_pw_r, &
                      rho_r_sccs=rho_pw_r_sccs)

      ! Retrieve the last rho_iter from the previous SCCS cycle if available
      CPASSERT(ASSOCIATED(rho_pw_r_sccs))
      rho_pw_r_sccs_fine => rho_pw_r_sccs
      IF (level > 1) THEN
         ! restricted electron density / stored potential on the solver grid
         ALLOCATE (rho_r_c(SIZE(rho_pw_r)))
         DO ispin = 1, SIZE(rho_pw_r)
            CALL auxbas_pw_pool%create_pw(rho_r_c(ispin))
            CALL sol_restrict(fine_pool, auxbas_pw_pool, rho_pw_r(ispin), rho_r_c(ispin))
         END DO
         rho_pw_r => rho_r_c
         CALL auxbas_pw_pool%create_pw(phi_c)
         CALL sol_restrict(fine_pool, auxbas_pw_pool, rho_pw_r_sccs_fine, phi_c)
         rho_pw_r_sccs => phi_c
      END IF

      ! Retrieve the total electronic density in r-space
      BLOCK
         TYPE(pw_r3d_rs_type) :: rho_elec
         CALL auxbas_pw_pool%create_pw(rho_elec)

         ! Retrieve grid parameters
         ngpts = rho_elec%pw_grid%ngpts
         dvol = rho_elec%pw_grid%dvol
         cell_volume = rho_elec%pw_grid%vol
         abc(1:3) = REAL(rho_elec%pw_grid%npts(1:3), KIND=dp)*rho_elec%pw_grid%dr(1:3)
         lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
         ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

         CALL pw_copy(rho_pw_r(1), rho_elec)
         DO ispin = 2, nspin
            CALL pw_axpy(rho_pw_r(ispin), rho_elec)
         END DO
         tot_rho_elec = pw_integrate_function(rho_elec)

         ! The cavity is built from n_valence + a pseudo-core density of Z - Zeff electrons per atom.
         ! Reuse CP2K's own core-charge collocation with explicit amplitudes: ccore = Z - Zeff
         ! (scaled by CORE_DENSITY_SCALE), alpha = alpha_core_charge of the kind.
         IF (sccs_control%pseudo_core .AND. &
             ((sccs_control%method_id == sccs_isodensity) .OR. &
              (sccs_control%method_id == sccs_solvent_excluded))) THEN
            BLOCK
               TYPE(pw_r3d_rs_type) :: rho_pscore
               ALLOCATE (pscore_alpha(SIZE(atomic_kind_set)), pscore_charge(SIZE(atomic_kind_set)))
               CALL sccs_pseudo_core_kinds(qs_kind_set, sccs_control, pscore_alpha, pscore_charge)
               IF (ANY(pscore_charge /= 0.0_dp)) THEN
                  CALL fine_pool%create_pw(rho_pscore)
                  CALL calculate_rho_core(rho_pscore, tot_rho_pscore, qs_env, &
                                          calpha=pscore_alpha, ccore=pscore_charge)
                  ! calculate_rho_core returns a negative (nuclear-convention) density
                  IF (level > 1) THEN
                     BLOCK
                        TYPE(pw_r3d_rs_type) :: rho_pscore_c
                        CALL auxbas_pw_pool%create_pw(rho_pscore_c)
                        CALL sol_restrict(fine_pool, auxbas_pw_pool, rho_pscore, rho_pscore_c)
                        CALL pw_axpy(rho_pscore_c, rho_elec, -1.0_dp)
                        CALL auxbas_pw_pool%give_back_pw(rho_pscore_c)
                     END BLOCK
                  ELSE
                     CALL pw_axpy(rho_pscore, rho_elec, -1.0_dp)
                  END IF
                  CALL fine_pool%give_back_pw(rho_pscore)
                  IF (should_output .AND. (output_unit > 0)) &
                     WRITE (UNIT=output_unit, FMT="(T2,A,T61,F20.10)") &
                     "SCCS| Pseudo-core charge added to the cavity density", -tot_rho_pscore
               END IF
               DEALLOCATE (pscore_alpha, pscore_charge)
            END BLOCK
         END IF

         ! Calculate the dielectric (smoothed) function of rho_elec in r-space
         CALL auxbas_pw_pool%create_pw(eps_elec)
         CALL auxbas_pw_pool%create_pw(deps_elec)
         IF (sccs_control%method_id == sccs_saa_andreussi) THEN
            CALL auxbas_pw_pool%create_pw(s)
            CALL auxbas_pw_pool%create_pw(d_s_rhoel)
            CALL auxbas_pw_pool%create_pw(u)
            CALL auxbas_pw_pool%create_pw(ff)
            CALL auxbas_pw_pool%create_pw(t)
            CALL auxbas_pw_pool%create_pw(dtf)
            f0 = sccs_control%f0
            delta_eta = sccs_control%delta_eta
            alpha_zeta = sccs_control%alpha_zeta
            R_solv = sccs_control%R_solv
            delta_zeta = sccs_control%delta_zeta
         END IF

         ! Relative permittivity or dielectric constant of the solvent (medium)
         epsilon_solvent = sccs_control%epsilon_solvent
         SELECT CASE (sccs_control%method_id)
         CASE (sccs_andreussi)
            CALL andreussi(rho_elec, eps_elec, deps_elec, epsilon_solvent, sccs_control%rho_max, &
                           sccs_control%rho_min)
         CASE (sccs_fattebert_gygi)
            CALL fattebert_gygi(rho_elec, eps_elec, deps_elec, epsilon_solvent, sccs_control%beta, &
                                sccs_control%rho_zero)
         CASE (sccs_isodensity)
            CALL isodensity_switch(rho_elec, eps_elec, deps_elec, epsilon_solvent, sccs_control%rho_cavity, &
                              sccs_control%switch_width)
            IF (calculate_stress_tensor) THEN
               CALL auxbas_pw_pool%create_pw(st_eps)
               CALL auxbas_pw_pool%create_pw(st_neps)
               CALL pw_copy(eps_elec, st_eps)
               st_neps%array(:, :, :) = rho_elec%array(:, :, :)*deps_elec%array(:, :, :)
            END IF
         CASE (sccs_solvent_excluded)
            ! nonlocal cavities: eps = 1 + (eps_b - 1) S_diel; all local density derivatives (deps_elec) vanish,
            ! the Kohn-Sham potential is assembled through the chained cavity derivatives below
            CALL timeset('sccs_cavity', h_cav)
            CALL cavity_build(cav, cpar, rho_elec, auxbas_pw_pool)
            CALL timestop(h_cav)
            eps_elec%array(:, :, :) = 1.0_dp + (epsilon_solvent - 1.0_dp)*cav%s_diel%array(:, :, :)
            CALL pw_zero(deps_elec)
         CASE (sccs_saa_andreussi)
            CALL sa_andreussi(rho_elec, eps_elec, deps_elec, epsilon_solvent, sccs_control%rho_max, &
                              sccs_control%rho_min, s, d_s_rhoel)
         CASE DEFAULT
            CPABORT("Invalid method specified for SCCS model")
         END SELECT

         ! &SCCS&REGION: planar confinement of the implicit solvent for hybrid explicit-water cells.
         ! m(z) = sum_n 1/2[1+erf((z+nL-z_lo)/(w sqrt2))] 1/2[1+erf((z_hi-z-nL)/(w sqrt2))], n = -1, 0, 1, along one cell
         ! axis, z wrapped into [0, L): the image sum keeps m smooth across the cell boundary (a one-sided window
         ! would jump from 1 to 0 at z = L). UPPER_EDGE may exceed L for a region that wraps through the boundary. Report how
         ! much of the ionic and dielectric cavity lies outside the region (the leakage into the explicit water), then
         ! (APPLY T) multiply the cavities by m. The mask does not depend on the density, so the Kohn-Sham potential
         ! only needs the same factor on the chained A fields (SOLVENT_EXCLUDED) or on deps
         ! (local cavities); the cavitation surface S_cav is left unmasked.
         region_on = sccs_control%region_on
         region_apply = region_on .AND. sccs_control%region_apply
         region_planar = region_on .AND. sccs_control%region_has_lo   ! LOWER_EDGE and UPPER_EDGE given (optional with FLEXIBLE)
         flex_on = region_on .AND. sccs_control%flex_on
         flex_cav = flex_on .AND. sccs_control%flex_cavity   ! mask from the density (1 - S_solv), not from atoms
         ! flex_frc: dE/dm is needed, for the force on the listed atoms or for the Kohn-Sham potential (flex_cav)
         flex_frc = flex_on .AND. region_apply .AND. (calc_forces .OR. flex_cav)
         IF (region_on) THEN
            CALL auxbas_pw_pool%create_pw(rmask)
            iax = sccs_control%region_axis
            zlen = REAL(rmask%pw_grid%npts(iax), KIND=dp)*rmask%pw_grid%dr(iax)
            wsq2 = sccs_control%region_width*SQRT(2.0_dp)
            IF (region_planar) THEN
               IF (sccs_control%region_zlo < 0.0_dp .OR. sccs_control%region_zlo >= zlen) &
                  CPABORT("SCCS REGION: LOWER_EDGE must lie in [0, L) of the cell axis")
               IF (sccs_control%region_zhi - sccs_control%region_zlo > zlen - 6.0_dp*sccs_control%region_width) &
                  CPABORT("SCCS REGION: UPPER_EDGE - LOWER_EDGE must be < L - 6 WIDTH (the window has to close)")
               region_hit = mask_key_match(region_cache, rmask, iax, 0, [sccs_control%region_zlo, sccs_control%region_zhi, &
                                                                         sccs_control%region_width, 0.0_dp])
               IF (region_hit) rmask%array(:, :, :) = region_cache%m1(:, :, :)
            ELSE
               region_hit = .TRUE.   ! no planes: the mask is the flexible one alone
               rmask%array(:, :, :) = 1.0_dp
            END IF
            IF (.NOT. region_hit) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,idx,zpos,mval,iimg,zimg) SHARED(iax,lb,rmask,sccs_control,ub,wsq2,zlen)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     SELECT CASE (iax)
                     CASE (1)
                        idx = i
                     CASE (2)
                        idx = j
                     CASE DEFAULT
                        idx = k
                     END SELECT
                     ! grid index -> position: the point with the global lower bound sits at the cell origin (the bounds
                     ! are centred, bounds(1) ~ -n/2)
                     zpos = MODULO(REAL(idx - rmask%pw_grid%bounds(1, iax), KIND=dp)*rmask%pw_grid%dr(iax), zlen)
                     mval = 0.0_dp
                     DO iimg = -1, 1
                        zimg = zpos + REAL(iimg, KIND=dp)*zlen
                        mval = mval + 0.25_dp*(1.0_dp + ERF((zimg - sccs_control%region_zlo)/wsq2))* &
                               (1.0_dp + ERF((sccs_control%region_zhi - zimg)/wsq2))
                     END DO
                     rmask%array(i, j, k) = mval
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
            CALL mask_key_store(region_cache, rmask, iax, 0, [sccs_control%region_zlo, sccs_control%region_zhi, &
                                                              sccs_control%region_width, 0.0_dp], rmask)
            END IF
            ! &REGION&FLEXIBLE: m_flex(r) from the coarse-grained density of the listed atoms multiplies the planar
            ! mask (or stands alone). flex_dm = d(rmask)/d(rho_bar) is kept for the force on the listed atoms.
            IF (flex_on) THEN
               IF (.NOT. do_ppcav) &
                  CPABORT("SCCS REGION FLEXIBLE is implemented only for METHOD SOLVENT_EXCLUDED")
               CALL auxbas_pw_pool%create_pw(flex_m)
               CALL auxbas_pw_pool%create_pw(flex_dm)
               IF (flex_cav) THEN
                  CALL flex_build_cavity(flex_m, flex_dm, auxbas_pw_pool, sccs_control, cav%s_solv, &
                                         cp_logger_get_default_io_unit(logger))
               ELSE
                  IF (ALLOCATED(sccs_control%flex_atoms)) THEN
                     CALL exclusion_resolve(atomic_kind_set, SIZE(particles%els), "REGION FLEXIBLE", &
                                            sccs_control%flex_kind_name, sccs_control%flex_weight, &
                                            sccs_control%flex_weight, nfx, fx_atom, fx_ikind, fx_iat, fx_w, fx_w2, &
                                            atoms=sccs_control%flex_atoms)
                  ELSE
                     CALL exclusion_resolve(atomic_kind_set, SIZE(particles%els), "REGION FLEXIBLE", &
                                            sccs_control%flex_kind_name, sccs_control%flex_weight, &
                                            sccs_control%flex_weight, nfx, fx_atom, fx_ikind, fx_iat, fx_w, fx_w2)
                  END IF
                  IF (nfx < 1) CPABORT("SCCS REGION FLEXIBLE: no atom selected (KIND / ATOM_LIST)")
                  CALL flex_cached(flex_m, flex_dm, auxbas_pw_pool, sccs_control, particles, nfx, fx_atom, fx_w, &
                                   cp_logger_get_default_io_unit(logger))
               END IF
               flex_dm%array(:, :, :) = flex_dm%array(:, :, :)*rmask%array(:, :, :)
               rmask%array(:, :, :) = rmask%array(:, :, :)*flex_m%array(:, :, :)
               CALL auxbas_pw_pool%give_back_pw(flex_m)
            END IF
            ! leakage report: cavity volume outside the region, for the ionic (S_ion) and dielectric (S_diel) cavities
            IF (do_ppcav) THEN
               rg_tot = pw_integrate_function(cav%s_ion); rg_in = pw_integral_ab(cav%s_ion, rmask)
               rg_tot2 = pw_integrate_function(cav%s_diel); rg_in2 = pw_integral_ab(cav%s_diel, rmask)
            ELSE
               ! local cavities: S = (eps - 1)/(eps_b - 1) carries both the dielectric and the ions
               rg_tot = (pw_integrate_function(eps_elec) - cell_volume)/(epsilon_solvent - 1.0_dp)
               rg_in = (pw_integral_ab(eps_elec, rmask) - pw_integrate_function(rmask))/(epsilon_solvent - 1.0_dp)
               rg_tot2 = rg_tot; rg_in2 = rg_in
            END IF
            iw = cp_logger_get_default_io_unit(logger)   ! like "SCCS| t_solv": every call, not gated by the print key
            IF (iw > 0) THEN
               WRITE (UNIT=iw, FMT="(T2,A,ES12.4,A,ES12.4,A,F8.4)") &
                  "SCCS| region: S_ion volume outside [bohr^3]", rg_tot - rg_in, "  total", rg_tot, &
                  "  pre-mask fraction", (rg_tot - rg_in)/MAX(rg_tot, 1.0E-30_dp)
               WRITE (UNIT=iw, FMT="(T2,A,ES12.4,A,ES12.4,A,F8.4)") &
                  "SCCS| region: S_diel volume outside [bohr^3]", rg_tot2 - rg_in2, "  total", rg_tot2, &
                  "  pre-mask fraction", (rg_tot2 - rg_in2)/MAX(rg_tot2, 1.0E-30_dp)
            END IF
            IF (region_apply) THEN
               IF (do_ppcav) THEN
                  IF (flex_frc) THEN   ! cavities before the region mask, for the force of the flexible mask
                     CALL auxbas_pw_pool%create_pw(flex_sd0)
                     CALL auxbas_pw_pool%create_pw(flex_si0)
                     CALL pw_copy(cav%s_diel, flex_sd0)
                     CALL pw_copy(cav%s_ion, flex_si0)
                  END IF
                  cav%s_ion%array(:, :, :) = cav%s_ion%array(:, :, :)*rmask%array(:, :, :)
                  cav%s_diel%array(:, :, :) = cav%s_diel%array(:, :, :)*rmask%array(:, :, :)
                  eps_elec%array(:, :, :) = 1.0_dp + (epsilon_solvent - 1.0_dp)*cav%s_diel%array(:, :, :)
               ELSE
                  eps_elec%array(:, :, :) = 1.0_dp + (eps_elec%array(:, :, :) - 1.0_dp)*rmask%array(:, :, :)
                  deps_elec%array(:, :, :) = deps_elec%array(:, :, :)*rmask%array(:, :, :)
               END IF
            END IF
         END IF

         ! &SCCS&EXCLUSION. M_c(r) = prod_a [1 - E_a(r)] per channel c (dielectric, ions); the cavities are
         ! multiplied by M_c AFTER the REGION mask (so the total mask of a channel is rmask*M_c) and eps is rebuilt.
         ! Unmasked-by-M copies (S_c*rmask) are kept for the force. Nothing here runs when EXCLUSION is off.
         excl_on = sccs_control%excl_on
         IF (excl_on) THEN
            IF (.NOT. do_ppcav) &
               CPABORT("SCCS EXCLUSION is implemented only for METHOD SOLVENT_EXCLUDED")
            IF (ALLOCATED(sccs_control%excl_atoms)) THEN
               CALL exclusion_resolve(atomic_kind_set, SIZE(particles%els), "EXCLUSION", sccs_control%excl_kind_name, &
                                      sccs_control%excl_r_diel, sccs_control%excl_r_ion, nex, ex_atom, ex_ikind, &
                                      ex_iat, ex_rd, ex_ri, atoms=sccs_control%excl_atoms)
            ELSE
               CALL exclusion_resolve(atomic_kind_set, SIZE(particles%els), "EXCLUSION", sccs_control%excl_kind_name, &
                                      sccs_control%excl_r_diel, sccs_control%excl_r_ion, nex, ex_atom, ex_ikind, &
                                      ex_iat, ex_rd, ex_ri)
            END IF
            CALL auxbas_pw_pool%create_pw(excl_md)
            CALL auxbas_pw_pool%create_pw(excl_mi)
            CALL exclusion_cached(excl_md, excl_mi, sccs_control, particles, nex, ex_atom, ex_rd, ex_ri)
            IF (calc_forces) THEN
               CALL auxbas_pw_pool%create_pw(excl_sd0)
               CALL auxbas_pw_pool%create_pw(excl_si0)
               CALL pw_copy(cav%s_diel, excl_sd0)
               CALL pw_copy(cav%s_ion, excl_si0)
               ALLOCATE (fex(3, MAX(nex, 1)))
               fex(:, :) = 0.0_dp
            END IF
            cav%s_ion%array(:, :, :) = cav%s_ion%array(:, :, :)*excl_mi%array(:, :, :)
            cav%s_diel%array(:, :, :) = cav%s_diel%array(:, :, :)*excl_md%array(:, :, :)
            eps_elec%array(:, :, :) = 1.0_dp + (epsilon_solvent - 1.0_dp)*cav%s_diel%array(:, :, :)
            IF (sccs_control%excl_print) THEN
               iw = cp_logger_get_default_io_unit(logger)
               rg_tot = excl_md%pw_grid%vol - pw_integrate_function(excl_md)
               rg_tot2 = excl_mi%pw_grid%vol - pw_integrate_function(excl_mi)
               IF (iw > 0) WRITE (UNIT=iw, FMT="(T2,A,I8,A,ES12.4,A,ES12.4)") &
                  "SCCS| exclusion: excluded atoms", nex, "  int(1-M_diel) dV [bohr^3]", rg_tot, &
                  "  int(1-M_ion) dV [bohr^3]", rg_tot2
            END IF
         END IF

         ! the readings are taken after the exclusion mask too (they run once both masks are applied)
         IF (region_on) THEN
            ! readings AFTER the mask, over the user's window, the whole cell and the PROBE_ATOMS spheres.
            ! Reporting only; nothing below reads these values. (With APPLY F the cavities are unmasked and the readings
            ! show what the explicit region would receive.)
            region_tag = "post-mask"
            IF (.NOT. region_apply) region_tag = "unmasked (APPLY F)"
            IF (.NOT. region_apply .AND. excl_on) region_tag = "post-mask (EXCLUSION)"
            CALL region_readings(auxbas_pw_pool, para_env, sccs_control, particles, iw, region_tag, "eps-1", &
                                 eps_elec, 1.0_dp, 1.0_dp, .TRUE.)
            IF (do_ppcav) THEN
               CALL region_readings(auxbas_pw_pool, para_env, sccs_control, particles, iw, region_tag, "S_ion", &
                                    cav%s_ion, 0.0_dp, 1.0_dp, .TRUE.)
               CALL region_readings(auxbas_pw_pool, para_env, sccs_control, particles, iw, region_tag, "S_diel", &
                                    cav%s_diel, 0.0_dp, 1.0_dp, .TRUE.)
            ELSE IF (ABS(epsilon_solvent - 1.0_dp) > epstol) THEN
               ! local cavity: S_ion = S_diel = (eps - 1)/(eps_b - 1)
               CALL region_readings(auxbas_pw_pool, para_env, sccs_control, particles, iw, region_tag, "S_ion=S_diel", &
                                    eps_elec, 1.0_dp, 1.0_dp/(epsilon_solvent - 1.0_dp), .TRUE.)
            END IF
         END IF

         IF (do_electrolyte) THEN
            ! kap2   = eps_b*zeta(r)/lambda_D^2           (full LPB screening coefficient, zeta = (eps - 1)/(eps_b - 1):
            !                                              zeta = 1 in the bulk electrolyte, 0 inside the solute; NOTE the CP2K
            !                                              cavity indicator theta = (eps_b - eps)/(eps_b - 1) is 1 - zeta)
            ! dkap2e = kap2/eps - kappa0^2                (deviation from the bulk value that is iterated in the SCCS cycle;
            !                                              the bulk part kappa0^2 lives in the screened Green function)
            CALL auxbas_pw_pool%create_pw(kap2)
            CALL auxbas_pw_pool%create_pw(dkap2e)
            f = epsilon_solvent*kappa0sq/(epsilon_solvent - 1.0_dp)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k,zeta) &
!$OMP             SHARED(cav,dkap2e,do_ppcav,eps_elec,epsilon_solvent,f,kap2,kappa0sq,lb,ub)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     zeta = eps_elec%array(i, j, k) - 1.0_dp
                     IF (do_ppcav) zeta = (epsilon_solvent - 1.0_dp)*cav%s_ion%array(i, j, k)
                     kap2%array(i, j, k) = f*zeta
                     dkap2e%array(i, j, k) = f*zeta/eps_elec%array(i, j, k) - kappa0sq
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
            IF (do_nlpb) THEN
               ! keep the real permittivity and the ionic cavity S_ion = (eps-1)/(eps_b-1) for the NLPB solver
               CALL auxbas_pw_pool%create_pw(eps_nl)
               CALL auxbas_pw_pool%create_pw(s_ion)
               CALL pw_copy(eps_elec, eps_nl)
               CALL auxbas_pw_pool%create_pw(s_diel_nl)
               f = 1.0_dp/(epsilon_solvent - 1.0_dp)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(eps_elec,f,lb,s_ion,ub)
               DO k = lb(3), ub(3)
                  DO j = lb(2), ub(2)
                     DO i = lb(1), ub(1)
                        s_ion%array(i, j, k) = f*(eps_elec%array(i, j, k) - 1.0_dp)
                     END DO
                  END DO
               END DO
!$OMP END PARALLEL DO
               CALL pw_copy(s_ion, s_diel_nl)
               IF (do_ppcav) THEN
                  CALL pw_copy(cav%s_ion, s_ion)
                  CALL pw_copy(cav%s_diel, s_diel_nl)
               END IF
            END IF
         END IF

         IF (do_nlpb .AND. (.NOT. do_electrolyte)) THEN
            ! SOLVER NEWTON_CG on the pure dielectric problem: S_ion = 0, S_diel = (eps - 1)/(eps_b - 1)
            CALL auxbas_pw_pool%create_pw(eps_nl)
            CALL auxbas_pw_pool%create_pw(s_ion)
            CALL auxbas_pw_pool%create_pw(s_diel_nl)
            CALL pw_copy(eps_elec, eps_nl)
            CALL pw_zero(s_ion)
            f = 1.0_dp/(epsilon_solvent - 1.0_dp)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(eps_elec,f,lb,s_diel_nl,ub)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     s_diel_nl%array(i, j, k) = f*(eps_elec%array(i, j, k) - 1.0_dp)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
            IF (do_ppcav) CALL pw_copy(cav%s_diel, s_diel_nl)
         END IF

         IF (sccs_control%method_id /= sccs_saa_andreussi) THEN
            ! Optional output of the dielectric function in cube file format
            filename = "DIELECTRIC_FUNCTION"
            cube_path = TRIM(print_path)//"%"//TRIM(filename)
            IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                      cp_p_file)) THEN
               append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
               my_pos_cube = "REWIND"
               IF (append_cube) my_pos_cube = "APPEND"
               mpi_io = .TRUE.
               cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                                extension=".cube", middle_name=TRIM(filename), &
                                                file_position=my_pos_cube, log_filename=.FALSE., &
                                                mpi_io=mpi_io, fout=mpi_filename)
               IF (output_unit > 0) THEN
                  IF (.NOT. mpi_io) THEN
                     INQUIRE (UNIT=cube_unit, NAME=filename)
                  ELSE
                     filename = mpi_filename
                  END IF
                  WRITE (UNIT=output_unit, FMT="(T3,A)") &
                     "SCCS| The dielectric function is written in cube file format to the file:", &
                     "SCCS| "//TRIM(filename)
               END IF
               CALL cp_pw_to_cube(eps_elec, cube_unit, TRIM(filename), particles=particles, &
                                  stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), &
                                  mpi_io=mpi_io)
               CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
            END IF
         END IF

         ! Calculate the (quantum) volume and surface of the solute cavity
         cavity_surface = 0.0_dp
         cavity_volume = 0.0_dp

         IF (ABS(epsilon_solvent - 1.0_dp) > epstol) THEN

            BLOCK
               TYPE(pw_r3d_rs_type) :: theta, norm_drho_elec
               LOGICAL :: need_surf
               TYPE(pw_r3d_rs_type), DIMENSION(3)                        :: drho_elec
               CALL auxbas_pw_pool%create_pw(theta)
               CALL pw_zero(theta)

               ! Calculate the (quantum) volume of the solute cavity
               f = 1.0_dp/(epsilon_solvent - 1.0_dp)
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(epsilon_solvent,eps_elec,f,lb,theta,ub)
               DO k = lb(3), ub(3)
                  DO j = lb(2), ub(2)
                     DO i = lb(1), ub(1)
                        theta%array(i, j, k) = f*(epsilon_solvent - eps_elec%array(i, j, k))
                     END DO
                  END DO
               END DO
!$OMP    END PARALLEL DO
               IF (do_ppcav) theta%array(:, :, :) = 1.0_dp - cav%s_cav%array(:, :, :)
               cavity_volume = pw_integrate_function(theta)

               ! |grad n|, the cavity surface and the cavitation potential are needed only when a
               ! non-electrostatic term is active, the SOLVENT_EXCLUDED chain needs its cavitation field, or SCCS prints
               need_surf = do_ppcav .OR. should_output .OR. &
                           (ABS(sccs_control%alpha_solvent) + ABS(sccs_control%beta_solvent) + &
                            ABS(sccs_control%gamma_solvent) > 0.0_dp)
               ! SOLVENT_EXCLUDED builds its cavitation field from S_cav, never from the density gradient, and the
               ! surface integral only multiplies ALPHA and GAMMA (or is printed). With ALPHA = GAMMA = 0, no SCCS print and no
               ! DENSITY_GRADIENT cube, the field is the constant -BETA that the general path gives (g*div(...) = 0*x = 0).
               filename = "DENSITY_GRADIENT"
               cube_path = TRIM(print_path)//"%"//TRIM(filename)
               need_drho = (.NOT. do_ppcav) .OR. BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                                                      cp_p_file)
               se_skip_surf = do_ppcav .AND. (.NOT. should_output) .AND. (.NOT. need_drho) .AND. &
                              (sccs_control%alpha_solvent == 0.0_dp) .AND. (sccs_control%gamma_solvent == 0.0_dp)
               CALL timeset('sccs_surface', h_surf)
               IF (se_skip_surf) THEN
                  CALL auxbas_pw_pool%create_pw(acav_fld)
                  acav_fld%array(:, :, :) = -sccs_control%beta_solvent
               END IF
               IF (need_surf .AND. (.NOT. se_skip_surf)) THEN
               ! Calculate the derivative of the electronic density in r-space
               ! TODO: Could be retrieved from the QS environment
               IF (need_drho) THEN
               DO i = 1, 3
                  CALL auxbas_pw_pool%create_pw(drho_elec(i))
               END DO
               CALL derive(rho_elec, drho_elec, sccs_derivative_fft, pw_env, input, pool=auxbas_pw_pool)

               CALL auxbas_pw_pool%create_pw(norm_drho_elec)

               ! Calculate the norm of the gradient of the electronic density in r-space
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(drho_elec,lb,norm_drho_elec,ub)
               DO k = lb(3), ub(3)
                  DO j = lb(2), ub(2)
                     DO i = lb(1), ub(1)
                        norm_drho_elec%array(i, j, k) = SQRT(drho_elec(1)%array(i, j, k)* &
                                                             drho_elec(1)%array(i, j, k) + &
                                                             drho_elec(2)%array(i, j, k)* &
                                                             drho_elec(2)%array(i, j, k) + &
                                                             drho_elec(3)%array(i, j, k)* &
                                                             drho_elec(3)%array(i, j, k))
                     END DO
                  END DO
               END DO
!$OMP    END PARALLEL DO
               END IF ! need_drho

               ! Optional output of the norm of the density gradient in cube file format
               IF (need_drho .AND. BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                         cp_p_file)) THEN
                  append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
                  my_pos_cube = "REWIND"
                  IF (append_cube) my_pos_cube = "APPEND"
                  mpi_io = .TRUE.
                  cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                                   extension=".cube", middle_name=TRIM(filename), &
                                                   file_position=my_pos_cube, log_filename=.FALSE., &
                                                   mpi_io=mpi_io, fout=mpi_filename)
                  IF (output_unit > 0) THEN
                     IF (.NOT. mpi_io) THEN
                        INQUIRE (UNIT=cube_unit, NAME=filename)
                     ELSE
                        filename = mpi_filename
                     END IF
                     WRITE (UNIT=output_unit, FMT="(T3,A)") &
                        "SCCS| The norm of the density gradient is written in cube file format to the file:", &
                        "SCCS| "//TRIM(filename)
                  END IF
                  CALL cp_pw_to_cube(norm_drho_elec, cube_unit, TRIM(filename), particles=particles, &
                                     stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), &
                                     mpi_io=mpi_io)
                  CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
               END IF

               ! Calculate the (quantum) surface of the solute cavity
               SELECT CASE (sccs_control%method_id)
               CASE (sccs_andreussi)
                  CALL surface_andreussi(rho_elec, norm_drho_elec, theta, epsilon_solvent, &
                                         sccs_control%rho_max, sccs_control%rho_min, &
                                         sccs_control%delta_rho)
               CASE (sccs_isodensity)
                  CALL isodensity_surface(rho_elec, norm_drho_elec, theta, sccs_control%rho_cavity, &
                                       sccs_control%switch_width)
               CASE (sccs_solvent_excluded)
                  ! |grad S_cav| and the cavitation field A_cav = -(alpha+gamma) div(grad S_cav/|grad S_cav|) - beta
                  BLOCK
                     REAL(KIND=dp), PARAMETER :: gtol = 1.0E-12_dp
                     REAL(KIND=dp) :: g
                     TYPE(pw_c1d_gs_type) :: dg, ag
                     DO i = 1, 3
                        CALL auxbas_pw_pool%create_pw(dsc(i))
                     END DO
                     CALL derive(cav%s_cav, dsc, sccs_derivative_fft, pw_env, input, pool=auxbas_pw_pool)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,g) SHARED(dsc,lb,theta,ub)
                     DO k = lb(3), ub(3)
                        DO j = lb(2), ub(2)
                           DO i = lb(1), ub(1)
                              g = SQRT(dsc(1)%array(i, j, k)**2 + dsc(2)%array(i, j, k)**2 + dsc(3)%array(i, j, k)**2)
                              theta%array(i, j, k) = g
                              IF (g > gtol) THEN
                                 dsc(1)%array(i, j, k) = dsc(1)%array(i, j, k)/g
                                 dsc(2)%array(i, j, k) = dsc(2)%array(i, j, k)/g
                                 dsc(3)%array(i, j, k) = dsc(3)%array(i, j, k)/g
                              ELSE
                                 dsc(1)%array(i, j, k) = 0.0_dp; dsc(2)%array(i, j, k) = 0.0_dp; dsc(3)%array(i, j, k) = 0.0_dp
                              END IF
                           END DO
                        END DO
                     END DO
!$OMP END PARALLEL DO
                     CALL auxbas_pw_pool%create_pw(acav_fld)
                     ! div u in G space: three forward transforms, one inverse
                     CALL auxbas_pw_pool%create_pw(dg)
                     CALL auxbas_pw_pool%create_pw(ag)
                     CALL pw_zero(ag)
                     DO l = 1, 3
                        n3(:) = 0
                        n3(l) = 1
                        CALL pw_transfer(dsc(l), dg)
                        CALL pw_derive(dg, n3)
                        CALL pw_axpy(dg, ag)
                     END DO
                     CALL pw_transfer(ag, acav_fld)
                     CALL auxbas_pw_pool%give_back_pw(ag)
                     CALL auxbas_pw_pool%give_back_pw(dg)
                     g = -(sccs_control%alpha_solvent + sccs_control%gamma_solvent)
                     acav_fld%array(:, :, :) = g*acav_fld%array(:, :, :) - sccs_control%beta_solvent
                     DO i = 1, 3
                        CALL auxbas_pw_pool%give_back_pw(dsc(i))
                     END DO
                  END BLOCK
               CASE (sccs_fattebert_gygi)
                  CALL surface_fattebert_gygi(rho_elec, norm_drho_elec, theta, epsilon_solvent, &
                                              sccs_control%beta, sccs_control%rho_zero, &
                                              sccs_control%delta_rho)
               CASE (sccs_saa_andreussi)
               CASE DEFAULT
                  CPABORT("Invalid method specified for SCCS model")
               END SELECT
               cavity_surface = pw_integrate_function(theta)

               IF (calculate_stress_tensor .AND. sccs_control%method_id == sccs_isodensity) THEN
                  ! strain derivative of  (alpha+gamma) int f(n)|grad n| + beta int theta(n), f = |dS/dn|:
                  !   surface: d/de_ab = -int [ n f'(n) delta_ab |grad n| + f d_a n d_b n/|grad n| ]
                  !   volume : d/de_ab = delta_ab int [ theta - n f ]       (n -> n/J, grid-index measure fixed)
                  BLOCK
                     REAL(KIND=dp), PARAMETER :: rhotol = 1.0E-12_dp, gtol = 1.0E-12_dp
                     REAL(KIND=dp) :: rr, xx, ff, nfp, gg, th, ssurf(3, 3), svol
                     REAL(KIND=dp) :: grd(3)
                     ssurf = 0.0_dp; svol = 0.0_dp
                     DO k = lb(3), ub(3)
                        DO j = lb(2), ub(2)
                           DO i = lb(1), ub(1)
                              rr = rho_elec%array(i, j, k)
                              IF (rr < rhotol) CYCLE
                              xx = LOG(rr) - LOG(sccs_control%rho_cavity)
                              ff = EXP(-0.5_dp*(xx/sccs_control%switch_width)**2)/(SQRT(twopi)*sccs_control%switch_width*rr)
                              nfp = ff*(-xx/sccs_control%switch_width**2 - 1.0_dp)
                              grd(1) = drho_elec(1)%array(i, j, k); grd(2) = drho_elec(2)%array(i, j, k)
                              grd(3) = drho_elec(3)%array(i, j, k)
                              gg = SQRT(grd(1)**2 + grd(2)**2 + grd(3)**2)
                              th = (epsilon_solvent - eps_elec%array(i, j, k))/(epsilon_solvent - 1.0_dp)
                              svol = svol + (th - rr*ff)
                              DO st_a = 1, 3
                                 ssurf(st_a, st_a) = ssurf(st_a, st_a) - nfp*gg
                                 IF (gg > gtol) THEN
                                    DO st_b = 1, 3
                                       ssurf(st_a, st_b) = ssurf(st_a, st_b) - ff*grd(st_a)*grd(st_b)/gg
                                    END DO
                                 END IF
                              END DO
                           END DO
                        END DO
                     END DO
                     ssurf = ssurf*dvol; svol = svol*dvol
                     CALL para_env%sum(ssurf); CALL para_env%sum(svol)
                     sig_cav = (sccs_control%alpha_solvent + sccs_control%gamma_solvent)*ssurf
                     DO st_a = 1, 3
                        sig_cav(st_a, st_a) = sig_cav(st_a, st_a) + sccs_control%beta_solvent*svol
                     END DO
                  END BLOCK
               END IF

               IF (sccs_control%method_id == sccs_isodensity) THEN
                  ! Variational cavitation/dispersion/repulsion potential (Mathew 2014 Eq 13, 2nd term):
                  ! v = dvol*|S'|*(beta - (alpha+gamma)*div u), u = grad(n)/|grad(n)|, |S'| = -deps/(eps_b-1)
                  BLOCK
                     TYPE(pw_r3d_rs_type) :: div_u
                     TYPE(pw_r3d_rs_type), DIMENSION(3) :: du
                     REAL(KIND=dp), PARAMETER :: gtol = 1.0E-12_dp
                     REAL(KIND=dp) :: ab, g, sp
                     ! u_i = drho_i/|drho| in place
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k,g) &
!$OMP                SHARED(drho_elec,lb,norm_drho_elec,ub)
                     DO k = lb(3), ub(3)
                        DO j = lb(2), ub(2)
                           DO i = lb(1), ub(1)
                              g = norm_drho_elec%array(i, j, k)
                              IF (g > gtol) THEN
                                 drho_elec(1)%array(i, j, k) = drho_elec(1)%array(i, j, k)/g
                                 drho_elec(2)%array(i, j, k) = drho_elec(2)%array(i, j, k)/g
                                 drho_elec(3)%array(i, j, k) = drho_elec(3)%array(i, j, k)/g
                              ELSE
                                 drho_elec(1)%array(i, j, k) = 0.0_dp
                                 drho_elec(2)%array(i, j, k) = 0.0_dp
                                 drho_elec(3)%array(i, j, k) = 0.0_dp
                              END IF
                           END DO
                        END DO
                     END DO
!$OMP END PARALLEL DO
                     CALL auxbas_pw_pool%create_pw(div_u)
                     CALL pw_zero(div_u)
                     DO i = 1, 3
                        CALL auxbas_pw_pool%create_pw(du(i))
                     END DO
                     DO l = 1, 3
                        CALL derive(drho_elec(l), du, sccs_derivative_fft, pw_env, input, pool=auxbas_pw_pool)
                        CALL pw_axpy(du(l), div_u)
                     END DO
                     DO i = 1, 3
                        CALL auxbas_pw_pool%give_back_pw(du(i))
                     END DO
                     ab = sccs_control%alpha_solvent + sccs_control%gamma_solvent
                     sp = -dvol/(epsilon_solvent - 1.0_dp)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(ab,deps_elec,div_u,lb,sccs_control,sp,ub,v_sccs)
                     DO k = lb(3), ub(3)
                        DO j = lb(2), ub(2)
                           DO i = lb(1), ub(1)
                              v_sccs%array(i, j, k) = v_sccs%array(i, j, k) + sp*deps_elec%array(i, j, k)* &
                                                      (sccs_control%beta_solvent - ab*div_u%array(i, j, k))
                           END DO
                        END DO
                     END DO
!$OMP END PARALLEL DO
                     CALL auxbas_pw_pool%give_back_pw(div_u)
                  END BLOCK
               END IF

               ! Release storage
               IF (need_drho) THEN
                  CALL auxbas_pw_pool%give_back_pw(norm_drho_elec)
                  DO i = 1, 3
                     CALL auxbas_pw_pool%give_back_pw(drho_elec(i))
                  END DO
               END IF
               END IF ! need_surf
               CALL timestop(h_surf)
               CALL auxbas_pw_pool%give_back_pw(theta)
            END BLOCK

         END IF ! epsilon_solvent > 1

         CALL auxbas_pw_pool%give_back_pw(rho_elec)
      END BLOCK

      BLOCK
         TYPE(pw_r3d_rs_type) :: rho_tot, phi_tot, rho_solute, rho_tot_zero
         ! Retrieve the total charge density (core + elec) of the solute in r-space
         CALL auxbas_pw_pool%create_pw(rho_solute)
         CALL pw_zero(rho_solute)
         CALL pw_transfer(rho_tot_gspace, rho_solute)
         tot_rho_solute = pw_integrate_function(rho_solute)

         ! Check total charge
         IF (ABS(tot_rho_solute) >= 1.0E-6_dp) THEN
            IF ((poisson_env%parameters%solver /= pw_poisson_analytic) .AND. &
                (poisson_env%parameters%solver /= pw_poisson_mt)) THEN
               WRITE (UNIT=message, FMT="(A,SP,F0.6,A)") &
                  "The system (solute) has a non-negligible charge of ", -tot_rho_solute, &
                  ". It is recommended to use non-periodic boundary conditions (PERIODIC none) "// &
                  "combined with an appropriate Poisson solver (POISSON_SOLVER MT or analytic)"
               CPWARN(message)
            END IF
         END IF

         ! Reassign work storage to rho_tot_zero, because rho_elec is no longer needed
         CALL auxbas_pw_pool%create_pw(rho_tot_zero)

         IF (sccs_control%method_id /= sccs_saa_andreussi) THEN

            ! Build the initial (rho_iter = 0) total charge density (solute plus polarisation) in r-space
            ! eps_elec <- ln(eps_elec)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k) &
!$OMP             SHARED(eps_elec,lb,message,output_unit,para_env,ub) &
!$OMP             SHARED(rho_solute,rho_tot_zero)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     IF (eps_elec%array(i, j, k) < 1.0_dp) THEN
                        WRITE (UNIT=message, FMT="(A,ES12.3,A,3(I0,A))") &
                           "SCCS| Invalid dielectric function value ", eps_elec%array(i, j, k), &
                           " encountered at grid point (", i, ",", j, ",", k, ")"
                        CPABORT(message)
                     END IF
                     rho_tot_zero%array(i, j, k) = rho_solute%array(i, j, k)/eps_elec%array(i, j, k)
                     eps_elec%array(i, j, k) = LOG(eps_elec%array(i, j, k))
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

         ELSE

            lxyz(:) = REAL(u%pw_grid%bounds(1, :), dp)
            uxyz(:) = REAL(u%pw_grid%bounds(2, :), dp)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k,r) &
!$OMP             SHARED(lb,message,output_unit,para_env,ub,dtf,uxyz,lxyz) &
!$OMP             SHARED(u, alpha_zeta, R_solv, delta_zeta, eps_elec)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     IF (eps_elec%array(i, j, k) < 1.0_dp) THEN
                        WRITE (UNIT=message, FMT="(A,ES12.3,A,3(I0,A))") &
                           "SCCS| Invalid dielectric function value ", eps_elec%array(i, j, k), &
                           " encountered at grid point (", i, ",", j, ",", k, ")"
                        CPABORT(message)
                     END IF
                     IF (i > -1) r = ((REAL(uxyz(1), dp) - REAL(i, dp) + 0.5_dp)*u%pw_grid%dr(1))**2.0_dp
                     IF (i < 0) r = ((REAL(i, dp) - REAL(lxyz(1), dp) + 0.5_dp)*u%pw_grid%dr(1))**2.0_dp
                     IF (j > -1) r = ((REAL(uxyz(2), dp) - REAL(j, dp) + 0.5_dp)*u%pw_grid%dr(2))**2.0_dp + r
                     IF (j < 0) r = ((REAL(j, dp) - REAL(lxyz(2), dp) + 0.5_dp)*u%pw_grid%dr(2))**2.0_dp + r
                     IF (k > -1) r = ((REAL(uxyz(3), dp) - REAL(k, dp) + 0.5_dp)*u%pw_grid%dr(3))**2.0_dp + r
                     IF (k < 0) r = ((REAL(k, dp) - REAL(lxyz(3), dp) + 0.5_dp)*u%pw_grid%dr(3))**2.0_dp + r
                     r = SQRT(r)
                     u%array(i, j, k) = (1.0_dp/2.0_dp)*erfc((r - alpha_zeta*R_solv)/delta_zeta)
                     dtf%array(i, j, k) = 1.0_dp
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

            CALL pw_func_u_convolution(poisson_env=poisson_env, func=dtf, convolution=ff, u=u)
            CALL pw_scale(u, 1.0_dp/ff%array(0, 0, 0))
            CALL pw_func_u_convolution(poisson_env=poisson_env, func=s, convolution=ff, u=u)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k) &
!$OMP             SHARED(lb,para_env,ub, delta_eta, f0, epsilon_solvent) &
!$OMP             SHARED(ff, t, dtf, eps_elec, rho_tot_zero, rho_solute, s)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     t%array(i, j, k) = 0.5_dp*(1.0_dp + erf((ff%array(i, j, k) - f0)/delta_eta))
                     dtf%array(i, j, k) = (1.0_dp/SQRT(pi))*EXP(-(ff%array(i, j, k) - f0)**2.0/(delta_eta**2.0))
                     eps_elec%array(i, j, k) = EXP(LOG(epsilon_solvent)*(1 - (s%array(i, j, k) &
                                                                              + (1 - s%array(i, j, k))*t%array(i, j, k))))
                     rho_tot_zero%array(i, j, k) = rho_solute%array(i, j, k)/eps_elec%array(i, j, k)
!                     eps_elec%array(i, j, k) = LOG(eps_elec%array(i, j, k))
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

            ! Optional output of the dielectric function in cube file format (solvent aware algorithm case)
            filename = "DIELECTRIC_FUNCTION"
            cube_path = TRIM(print_path)//"%"//TRIM(filename)
            IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                      cp_p_file)) THEN
               append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
               my_pos_cube = "REWIND"
               IF (append_cube) my_pos_cube = "APPEND"
               mpi_io = .TRUE.
               cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                                extension=".cube", middle_name=TRIM(filename), &
                                                file_position=my_pos_cube, log_filename=.FALSE., &
                                                mpi_io=mpi_io, fout=mpi_filename)
               IF (output_unit > 0) THEN
                  IF (.NOT. mpi_io) THEN
                     INQUIRE (UNIT=cube_unit, NAME=filename)
                  ELSE
                     filename = mpi_filename
                  END IF
                  WRITE (UNIT=output_unit, FMT="(T3,A)") &
                     "SCCS| The dielectric function is written in cube file format to the file:", &
                     "SCCS| "//TRIM(filename)
               END IF
               CALL cp_pw_to_cube(eps_elec, cube_unit, TRIM(filename), particles=particles, &
                                  stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), &
                                  mpi_io=mpi_io)
               CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
            END IF

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k) &
!$OMP             SHARED(lb,para_env,ub, delta_eta, f0, epsilon_solvent) &
!$OMP             SHARED(ff, t, dtf, eps_elec, rho_tot_zero, rho_solute, s)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     eps_elec%array(i, j, k) = LOG(eps_elec%array(i, j, k))
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

         END IF

         ! Build the derivative of LOG(eps_elec): only the Jacobi iteration uses it
         IF (.NOT. do_nlpb) THEN
         DO i = 1, 3
            CALL auxbas_pw_pool%create_pw(dln_eps_elec(i))
            CALL pw_zero(dln_eps_elec(i))
         END DO
         CALL derive(eps_elec, dln_eps_elec, sccs_control%derivative_method, pw_env, input, pool=auxbas_pw_pool)
         END IF
         IF (sccs_control%method_id /= sccs_saa_andreussi) THEN
            CALL auxbas_pw_pool%give_back_pw(eps_elec)
         END IF

         ! Print header for the SCCS cycle
         IF (should_output .AND. (output_unit > 0)) THEN
            IF (print_level > low_print_level) THEN
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.12)") &
                  "SCCS| Total electronic charge density ", -tot_rho_elec, &
                  "SCCS| Total charge density (solute)   ", -tot_rho_solute
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.3)") &
                  "SCCS| Volume of the cell           [bohr^3]", cell_volume, &
                  "SCCS|                              [angstrom^3]", &
                  cp_unit_from_cp2k(cell_volume, "angstrom^3")
               IF (ABS(epsilon_solvent - 1.0_dp) > epstol) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.3)") &
                     "SCCS| Volume of the solute cavity  [bohr^3]", cavity_volume, &
                     "SCCS|                              [angstrom^3]", &
                     cp_unit_from_cp2k(cavity_volume, "angstrom^3"), &
                     "SCCS| Surface of the solute cavity [bohr^2]", cavity_surface, &
                     "SCCS|                              [angstrom^2]", &
                     cp_unit_from_cp2k(cavity_surface, "angstrom^2")
               END IF
               WRITE (UNIT=output_unit, FMT="(T3,A)") &
                  "SCCS|", &
                  "SCCS|   Step    Average residual    Maximum residual         E(Hartree) [a.u.]"
            END IF
         END IF

         ! Get storage for the derivative of the total potential (field) in r-space
         DO i = 1, 3
            CALL auxbas_pw_pool%create_pw(dphi_tot(i))
         END DO

         ! Initialise the total charge density in r-space rho_tot with rho_tot_zero + rho_iter_zero
         CALL auxbas_pw_pool%create_pw(rho_tot)
         CALL pw_copy(rho_tot_zero, rho_tot)
         CALL pw_axpy(rho_pw_r_sccs, rho_tot)

         CALL auxbas_pw_pool%create_pw(phi_tot)
         CALL pw_zero(phi_tot)

         IF (do_electrolyte) THEN
            ! Screened (Debye-Hueckel) Green function 4*pi/(G^2 + kappa0^2); finite at G = 0
            CALL auxbas_pw_pool%create_pw(greenfn_scr)
            DO ig = 1, greenfn_scr%pw_grid%ngpts_cut_local
               greenfn_scr%array(ig) = fourpi/(greenfn_scr%pw_grid%gsq(ig) + kappa0sq)
            END DO
         ELSE IF (do_nlpb) THEN
            ! No ions: preconditioner = plain periodic Green function 4*pi/G^2 with G = 0 -> 0
            CALL auxbas_pw_pool%create_pw(greenfn_scr)
            CALL pw_zero(greenfn_scr)
            DO ig = greenfn_scr%pw_grid%first_gne0, greenfn_scr%pw_grid%ngpts_cut_local
               greenfn_scr%array(ig) = fourpi/greenfn_scr%pw_grid%gsq(ig)
            END DO
         END IF

         IF (sccs_control%method_id == sccs_saa_andreussi) THEN
            CALL auxbas_pw_pool%create_pw(p_e_interface)
            CALL pw_zero(p_e_interface)
         END IF

         IF (do_nlpb) THEN
            ! Nonlinear Poisson-Boltzmann: Newton-CG on phi; rho_pw_r_sccs stores phi between SCF steps in this mode
            nlpb%output_unit = -1
            IF (should_output) nlpb%output_unit = output_unit
            nlpb%io_unit = cp_logger_get_default_io_unit(logger)   ! guard / line-search lines always visible
            nlpb%scf_factor = sccs_control%nlpb_scf_factor
            nlpb%e_tol = sccs_control%nlpb_energy_tol
            nlpb%cg_diag = sccs_control%print_pol_diag
            nlpb%selfcheck = sccs_control%nlpb_selfcheck
            nlpb%coarse = (level > 1)
            nlpb%psi_form = sccs_control%psi_form
            nlpb%coarse_direct_sum = sccs_control%coarse_direct_sum
            nlpb%lap_carry = sccs_control%lap_carry
            nlpb%linear_single_cg = sccs_control%linear_single_cg
            nlpb%history_guess = sccs_control%history_guess
            nlpb%linear_cg_floor = sccs_control%linear_cg_floor
            ! COARSE_SPACE_IONS and SCF_TOL_JUMP (Debye ions only) belong to the grand-canonical SCF that can change N_e during
            ! the SCF (MAX_OUTER_UPDATES > 0, UPDATE_EVERY SCF_STEP, scheme 1 or EARLY_UPDATE_DELTA > 0); an ion-free run or a
            ! fixed-N_e run (no &GCE, or MAX_OUTER_UPDATES 0) uses the plain preconditioner and behaves as SCF_TOL_JUMP 0
            n_e_moves = .FALSE.
            IF (ASSOCIATED(qs_env%scf_control)) THEN
               IF (qs_env%scf_control%gce%do_gce) THEN
                  IF (qs_env%scf_control%gce%max_outer > 0 .OR. qs_env%scf_control%gce%update_every_step .OR. &
                      qs_env%scf_control%gce%scheme == 1 .OR. qs_env%scf_control%gce%early_update_delta > 0.0_dp) &
                     n_e_moves = .TRUE.
               END IF
            END IF
            nlpb%coarse_space_ions = sccs_control%coarse_space_ions .AND. n_e_moves
            nlpb%scf_tol_jump = 0.0_dp
            IF (do_electrolyte .AND. n_e_moves) nlpb%scf_tol_jump = sccs_control%scf_tol_jump
            nlpb%cold_start_from_last = sccs_control%cold_start_from_last
            nlpb%newton_cg_forcing = sccs_control%newton_cg_forcing
            IF (sccs_control%psi_form .AND. sccs_control%nldiel) &
               CPABORT("SCCS PSI_FORM needs NONLINEAR_DIELECTRIC F")
            ! SOLVER_STRIDE: solve on calls 1, 1+n, 1+2n, ...; a cold start (phi = 0) always solves
            nlpb%skip = (MOD(n_nlpb_calls, sccs_control%solver_stride) /= 0) .AND. &
                        (pw_integral_ab(rho_pw_r_sccs, rho_pw_r_sccs) > 0.0_dp)
            n_nlpb_calls = n_nlpb_calls + 1
            CALL pw_copy(rho_pw_r_sccs, phi_tot)
            CALL auxbas_pw_pool%create_pw(aion_unit)
            CALL auxbas_pw_pool%create_pw(adiel_unit)
            CALL auxbas_pw_pool%create_pw(geff_unit)
            ! ionic and dielectric cavities: identical for the local cavities, separate for SOLVENT_EXCLUDED
            CALL timeset('sccs_nlpb_call', h_nlpb)
            CALL nlpb_solve(nlpb, phi_tot, rho_solute, eps_nl, s_ion, s_diel_nl, greenfn_scr, auxbas_pw_pool, &
                            para_env, a_nl, q_ion_nl, e_ion_nl, e_ion_phi, e_diel_nl, aion_unit, adiel_unit, geff_unit, &
                            rms_nl, niter_nl, conv_nl, ncg_nl)
            CALL timestop(h_nlpb)
            CALL pw_copy(phi_tot, rho_pw_r_sccs)
            use_psi = sccs_control%psi_form   ! g_eff replaces |grad phi|^2 in the KS terms below
            IF (.NOT. conv_nl) CPWARN("SCCS NLPB: Newton iteration did not converge")
            ! diagnostics: smallest permittivity seen by this solve (collective, cheap)
            eps_min_diag = MINVAL(eps_nl%array)
            CALL para_env%min(eps_min_diag)
            ! hand the convergence flag to the SCF. A frozen potential (SOLVER_STRIDE) tells nothing about
            ! convergence and leaves the count alone; SOLVER_MAX_FAILURES = 0 keeps the warn-only behaviour.
            IF (.NOT. nlpb%skip) THEN
               IF (conv_nl) THEN
                  n_nlpb_fail = 0
               ELSE
                  n_nlpb_fail = n_nlpb_fail + 1
                  IF (nlpb%io_unit > 0) WRITE (UNIT=nlpb%io_unit, FMT="(T3,A,I6,A,I6,A,ES12.4)") &
                     "NLPB| not converged: consecutive failed SCF steps ", n_nlpb_fail, "  limit ", &
                     sccs_control%max_nlpb_failures, "  final RMS ", rms_nl
                  IF (sccs_control%max_nlpb_failures > 0 .AND. n_nlpb_fail >= sccs_control%max_nlpb_failures) THEN
                     WRITE (UNIT=message, FMT="(A,I0,A,ES11.4,A,ES11.4,A,ES11.4)") &
                        "SCCS NLPB did not converge in ", n_nlpb_fail, " consecutive SCF steps (SOLVER_MAX_FAILURES): "// &
                        "final RMS residual ", rms_nl, " e/bohr^3, V_ion ", pw_integrate_function(s_ion), &
                        " bohr^3, min eps ", eps_min_diag
                     CPABORT(TRIM(message))
                  END IF
               END IF
            END IF
            iter = niter_nl
         ELSE
         ! Main SCCS iteration loop
         iter = 0

         iter_loop: DO

            ! Increment iteration counter
            iter = iter + 1

            ! Check if the requested maximum number of SCCS iterations is reached
            IF (iter > sccs_control%max_iter) THEN
               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,A,/,T3,A,I0,A)") &
                     "SCCS| Maximum number of SCCS iterations reached", &
                     "SCCS| Iteration cycle did not converge in ", sccs_control%max_iter, " steps"
               ELSE
                  WRITE (UNIT=message, FMT="(A,I0,A)") &
                     "The SCCS iteration cycle did not converge in ", sccs_control%max_iter, " steps"
                  CPWARN(message)
               END IF
               EXIT iter_loop
            END IF

            ! Calculate derivative of the current total potential in r-space
            IF (do_electrolyte) THEN
               CALL pw_poisson_solve(poisson_env=poisson_env, &
                                     density=rho_tot, &
                                     vhartree=phi_tot, &
                                     dvhartree=dphi_tot, &
                                     greenfn=greenfn_scr)
            ELSE
               CALL pw_poisson_solve(poisson_env=poisson_env, &
                                     density=rho_tot, &
                                     vhartree=phi_tot, &
                                     dvhartree=dphi_tot)
            END IF
            energy%sccs_hartree = 0.5_dp*pw_integral_ab(rho_solute, phi_tot)

            ! Update total charge density (solute plus polarisation) in r-space
            ! based on the iterated polarisation charge density
            f = 1.0_dp/fourpi
            rho_delta_avg = 0.0_dp
            rho_delta_max = 0.0_dp
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k,rho_delta,rho_iter_new) &
!$OMP                SHARED(dln_eps_elec,dphi_tot,f,lb,rho_pw_r_sccs,ub) &
!$OMP                SHARED(rho_tot,rho_tot_zero,sccs_control) &
!$OMP                SHARED(dkap2e,do_electrolyte,phi_tot) &
!$OMP                REDUCTION(+:rho_delta_avg) &
!$OMP                REDUCTION(MAX:rho_delta_max)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     rho_iter_new = (dln_eps_elec(1)%array(i, j, k)*dphi_tot(1)%array(i, j, k) + &
                                     dln_eps_elec(2)%array(i, j, k)*dphi_tot(2)%array(i, j, k) + &
                                     dln_eps_elec(3)%array(i, j, k)*dphi_tot(3)%array(i, j, k))*f
                     IF (do_electrolyte) THEN
                        ! LPB ionic screening: -(kap2/eps - kappa0^2)*phi/(4*pi)
                        rho_iter_new = rho_iter_new - dkap2e%array(i, j, k)*phi_tot%array(i, j, k)*f
                     END IF
                     rho_iter_new = rho_pw_r_sccs%array(i, j, k) + &
                                    sccs_control%mixing*(rho_iter_new - rho_pw_r_sccs%array(i, j, k))
                     rho_delta = ABS(rho_iter_new - rho_pw_r_sccs%array(i, j, k))
                     rho_delta_max = MAX(rho_delta, rho_delta_max)
                     rho_delta_avg = rho_delta_avg + rho_delta
                     rho_tot%array(i, j, k) = rho_tot_zero%array(i, j, k) + rho_iter_new
                     rho_pw_r_sccs%array(i, j, k) = rho_iter_new
                  END DO
               END DO
            END DO
!$OMP    END PARALLEL DO

            CALL para_env%sum(rho_delta_avg)
            rho_delta_avg = rho_delta_avg/REAL(ngpts, KIND=dp)
            CALL para_env%max(rho_delta_max)

            IF (should_output .AND. (output_unit > 0)) THEN
               IF (print_level > low_print_level) THEN
                  IF ((ABS(rho_delta_avg) < 1.0E-8_dp) .OR. &
                      (ABS(rho_delta_avg) >= 1.0E5_dp)) THEN
                     WRITE (UNIT=output_unit, FMT="(T3,A,I6,4X,ES16.4,4X,ES16.4,1X,F25.12)") &
                        "SCCS| ", iter, rho_delta_avg, rho_delta_max, energy%sccs_hartree
                  ELSE
                     WRITE (UNIT=output_unit, FMT="(T3,A,I6,4X,F16.8,4X,F16.8,1X,F25.12)") &
                        "SCCS| ", iter, rho_delta_avg, rho_delta_max, energy%sccs_hartree
                  END IF
               END IF
            END IF

            ! Check if the SCCS iteration cycle is converged to the requested tolerance
            IF (rho_delta_max <= sccs_control%eps_sccs) THEN
               IF (should_output .AND. (output_unit > 0)) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,A,I0,A)") &
                     "SCCS| Iteration cycle converged in ", iter, " steps"
               END IF
               EXIT iter_loop
            END IF

         END DO iter_loop
         END IF ! do_nlpb

         IF (sccs_control%method_id == sccs_saa_andreussi) THEN
            CALL auxbas_pw_pool%create_pw(tem_convolution)
            CALL auxbas_pw_pool%create_pw(tem_func)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k,dphi2) &
!$OMP                SHARED(lb,ub, tem_func, s, dtf) &
!$OMP                SHARED(dphi_tot, p_e_interface, eps_elec, epsilon_solvent)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     p_e_interface%array(i, j, k) = 1.0_dp/(4.0_dp*twopi)*(EXP(eps_elec%array(i, j, k))*LOG(epsilon_solvent))
                     dphi2 = dphi_tot(1)%array(i, j, k)*dphi_tot(1)%array(i, j, k) + &
                             dphi_tot(2)%array(i, j, k)*dphi_tot(2)%array(i, j, k) + &
                             dphi_tot(3)%array(i, j, k)*dphi_tot(3)%array(i, j, k)
                     p_e_interface%array(i, j, k) = p_e_interface%array(i, j, k)*dphi2
                     tem_func%array(i, j, k) = (1.0_dp - s%array(i, j, k))*p_e_interface%array(i, j, k)*dtf%array(i, j, k)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

            CALL pw_func_u_convolution(poisson_env=poisson_env, func=tem_func, convolution=tem_convolution, u=u)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(lb,ub) &
!$OMP                SHARED(p_e_interface, t, tem_convolution, phi_tot, d_s_rhoel)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     p_e_interface%array(i, j, k) = p_e_interface%array(i, j, k)*(1.0_dp - t%array(i, j, k)) &
                                                    + tem_convolution%array(i, j, k)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

         END IF

         ! Release work storage which is no longer needed
         CALL auxbas_pw_pool%give_back_pw(rho_tot_zero)
         IF (.NOT. do_nlpb) THEN
         DO i = 1, 3
            CALL auxbas_pw_pool%give_back_pw(dln_eps_elec(i))
         END DO
         END IF

         ! Optional output of the total charge density in cube file format
         filename = "TOTAL_CHARGE_DENSITY"
         cube_path = TRIM(print_path)//"%"//TRIM(filename)
         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), cp_p_file)) THEN
            append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
            my_pos_cube = "REWIND"
            IF (append_cube) my_pos_cube = "APPEND"
            mpi_io = .TRUE.
            cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                             extension=".cube", middle_name=TRIM(filename), &
                                             file_position=my_pos_cube, log_filename=.FALSE., &
                                             mpi_io=mpi_io, fout=mpi_filename)
            IF (output_unit > 0) THEN
               IF (.NOT. mpi_io) THEN
                  INQUIRE (UNIT=cube_unit, NAME=filename)
               ELSE
                  filename = mpi_filename
               END IF
               WRITE (UNIT=output_unit, FMT="(T3,A)") &
                  "SCCS| The total SCCS charge density is written in cube file format to the file:", &
                  "SCCS| "//TRIM(filename)
            END IF
            IF (do_nlpb) THEN
               ! rho_tot is a work field here. Assemble the total charge in a
               ! temporary: -lap(phi)/4pi (= solute + polarisation + ions, zero mean) plus the uniform part
               ! (Q_solute + Q_ion)/V that the Laplacian drops. Energies and potentials are not touched.
               BLOCK
                  TYPE(pw_r3d_rs_type) :: rho_cube_t, rho_ion_t
                  TYPE(pw_c1d_gs_type) :: gs_cube_t
                  REAL(KIND=dp)        :: q_off_t
                  CALL auxbas_pw_pool%create_pw(gs_cube_t)
                  CALL auxbas_pw_pool%create_pw(rho_cube_t)
                  CALL pw_transfer(phi_tot, gs_cube_t)
                  CALL pw_laplace(gs_cube_t)
                  CALL pw_scale(gs_cube_t, -1.0_dp/fourpi)
                  CALL pw_transfer(gs_cube_t, rho_cube_t)
                  q_off_t = tot_rho_solute
                  IF (do_electrolyte) THEN
                     CALL auxbas_pw_pool%create_pw(rho_ion_t)
                     CALL nlpb_ion_density(nlpb, phi_tot, s_ion, rho_ion_t)
                     q_off_t = q_off_t + pw_integrate_function(rho_ion_t)
                     CALL auxbas_pw_pool%give_back_pw(rho_ion_t)
                  END IF
                  rho_cube_t%array(:, :, :) = rho_cube_t%array(:, :, :) + q_off_t/rho_cube_t%pw_grid%vol
                  CALL cp_pw_to_cube(rho_cube_t, cube_unit, TRIM(filename), particles=particles, &
                                     stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), mpi_io=mpi_io)
                  CALL auxbas_pw_pool%give_back_pw(rho_cube_t)
                  CALL auxbas_pw_pool%give_back_pw(gs_cube_t)
               END BLOCK
            ELSE
               CALL cp_pw_to_cube(rho_tot, cube_unit, TRIM(filename), particles=particles, &
                                  stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), mpi_io=mpi_io)
            END IF
            CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
         END IF

         ! Optional output of the ionic (electrolyte) charge density (print key IONIC_CHARGE_DENSITY)
         IF (do_electrolyte) THEN
            filename = "IONIC_CHARGE_DENSITY"
            cube_path = TRIM(print_path)//"%"//TRIM(filename)
            IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), cp_p_file)) THEN
               BLOCK
                  TYPE(pw_r3d_rs_type) :: rho_ion_cube
                  CALL auxbas_pw_pool%create_pw(rho_ion_cube)
                  IF (do_nlpb) THEN
                     CALL nlpb_ion_density(nlpb, phi_tot, s_ion, rho_ion_cube)
                  ELSE
                     ! linear ions: rho_ion = -eps_b kappa0^2 S phi/(4 pi), kap2 = eps_b kappa0^2 S
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(kap2,lb,phi_tot,rho_ion_cube,ub)
                     DO k = lb(3), ub(3)
                        DO j = lb(2), ub(2)
                           DO i = lb(1), ub(1)
                              rho_ion_cube%array(i, j, k) = -kap2%array(i, j, k)*phi_tot%array(i, j, k)/fourpi
                           END DO
                        END DO
                     END DO
!$OMP END PARALLEL DO
                  END IF
                  append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
                  my_pos_cube = "REWIND"
                  IF (append_cube) my_pos_cube = "APPEND"
                  mpi_io = .TRUE.
                  cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                                   extension=".cube", middle_name=TRIM(filename), &
                                                   file_position=my_pos_cube, log_filename=.FALSE., &
                                                   mpi_io=mpi_io, fout=mpi_filename)
                  IF (output_unit > 0) THEN
                     IF (.NOT. mpi_io) THEN
                        INQUIRE (UNIT=cube_unit, NAME=filename)
                     ELSE
                        filename = mpi_filename
                     END IF
                     WRITE (UNIT=output_unit, FMT="(T3,A)") &
                        "SCCS| The ionic charge density is written in cube file format to the file:", &
                        "SCCS| "//TRIM(filename)
                  END IF
                  CALL cp_pw_to_cube(rho_ion_cube, cube_unit, TRIM(filename), particles=particles, &
                                     stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), mpi_io=mpi_io)
                  CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
                  CALL auxbas_pw_pool%give_back_pw(rho_ion_cube)
               END BLOCK
            END IF
            ! implicit ionic charge inside the region after the mask, |rho_ion| (reporting only). It needs the
            ! converged potential, hence it sits here and not with the cavity readings above.
            IF (region_on) THEN
               BLOCK
                  TYPE(pw_r3d_rs_type) :: rho_ion_abs
                  CALL auxbas_pw_pool%create_pw(rho_ion_abs)
                  IF (do_nlpb) THEN
                     CALL nlpb_ion_density(nlpb, phi_tot, s_ion, rho_ion_abs)
                  ELSE
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(kap2,lb,phi_tot,rho_ion_abs,ub)
                     DO k = lb(3), ub(3)
                        DO j = lb(2), ub(2)
                           DO i = lb(1), ub(1)
                              rho_ion_abs%array(i, j, k) = -kap2%array(i, j, k)*phi_tot%array(i, j, k)/fourpi
                           END DO
                        END DO
                     END DO
!$OMP END PARALLEL DO
                  END IF
                  rho_ion_abs%array(:, :, :) = ABS(rho_ion_abs%array(:, :, :))
                  CALL region_readings(auxbas_pw_pool, para_env, sccs_control, particles, iw, region_tag, &
                                       "|rho_ion| [e]", rho_ion_abs, 0.0_dp, 1.0_dp, .FALSE.)
                  CALL auxbas_pw_pool%give_back_pw(rho_ion_abs)
               END BLOCK
            END IF
         END IF

         ! Calculate the total SCCS Hartree energy, potential, and its
         ! derivatives of the solute and the implicit solvent
         IF (do_nlpb) THEN
            ! phi from the Newton solver; effective vacuum source rho_sol + rho_solv = -nabla^2 phi/(4 pi)
            CALL pw_transfer(phi_tot, v_hartree_gspace)
            CALL pw_copy(v_hartree_gspace, rho_tot_gspace)
            CALL pw_laplace(rho_tot_gspace)
            CALL pw_scale(rho_tot_gspace, -1.0_dp/fourpi)
            ! the real-space copy is only read by the polarisation-charge print / cube below
            IF (need_rho_tot) CALL pw_transfer(rho_tot_gspace, rho_tot)
            ! with PSI_FORM every Kohn-Sham term uses g_eff, grad phi is never read
            IF (sccs_control%psi_form) THEN
               DO i = 1, 3
                  CALL pw_zero(dphi_tot(i))
               END DO
            ELSE
            CALL derive(phi_tot, dphi_tot, sccs_derivative_fft, pw_env, input, pool=auxbas_pw_pool)
            END IF
            ! electrostatic + solvent free energy at the stationary point = A[phi] of the solver
            energy%sccs_hartree = a_nl
         ELSE
         CALL pw_transfer(rho_tot, rho_tot_gspace)
         IF (do_electrolyte) THEN
            CALL pw_poisson_solve(poisson_env=poisson_env, &
                                  density=rho_tot_gspace, &
                                  ehartree=e_tot, &
                                  vhartree=v_hartree_gspace, &
                                  dvhartree=dphi_tot, &
                                  greenfn=greenfn_scr)
         ELSE
            CALL pw_poisson_solve(poisson_env=poisson_env, &
                                  density=rho_tot_gspace, &
                                  ehartree=e_tot, &
                                  vhartree=v_hartree_gspace, &
                                  dvhartree=dphi_tot)
         END IF
         CALL pw_transfer(v_hartree_gspace, phi_tot)
         energy%sccs_hartree = 0.5_dp*pw_integral_ab(rho_solute, phi_tot)
         END IF ! do_nlpb

         IF (calculate_stress_tensor) THEN
            ! strain derivative of the medium free energy
            !   A[phi] = int rho phi - (1/8pi) int eps(n)|grad phi|^2   at fixed phi, rho -> rho/J, n -> n/J, grad -> (1-e)grad:
            !   dA/de_ab = -delta_ab int (eps - n eps') |grad phi|^2/8pi + (1/4pi) int eps d_a phi d_b phi
            ! h_stress = -dE/de (CP2K virial sign, same as the Hartree h_stress); cavitation part from sig_cav.
            sig_el = 0.0_dp
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     st_e2 = dphi_tot(1)%array(i, j, k)**2 + dphi_tot(2)%array(i, j, k)**2 + dphi_tot(3)%array(i, j, k)**2
                     st_f = st_eps%array(i, j, k)
                     DO st_a = 1, 3
                        sig_el(st_a, st_a) = sig_el(st_a, st_a) - (st_f - st_neps%array(i, j, k))*st_e2/(2.0_dp*fourpi)
                        DO st_b = 1, 3
                           sig_el(st_a, st_b) = sig_el(st_a, st_b) + st_f*dphi_tot(st_a)%array(i, j, k)* &
                                                dphi_tot(st_b)%array(i, j, k)/fourpi
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            sig_el = sig_el*dvol
            CALL para_env%sum(sig_el)
            h_stress(:, :) = -(sig_el(:, :) + sig_cav(:, :))
            CALL auxbas_pw_pool%give_back_pw(st_eps)
            CALL auxbas_pw_pool%give_back_pw(st_neps)
         END IF

         IF (do_ppcav) THEN
            ! Kohn-Sham potential through the chained cavity derivatives (Islam et al. 2023 Eq. 62-66); dvol like v_sccs.
            ! the dielectric, ionic and cavitation A fields are summed inside one chain pass (cavity_potential_sum)
            CALL timeset('sccs_kspot', h_kspot)
            CALL auxbas_pw_pool%create_pw(afld)
            CALL auxbas_pw_pool%create_pw(vchain)
            ! dielectric: A = n_mol lambda_diel (nonlinear) or -(eps_b-1)|grad phi|^2/(8 pi) (linear)
            IF (do_nlpb .AND. sccs_control%nldiel) THEN
               CALL pw_copy(adiel_unit, afld)
            ELSE
               f = -(epsilon_solvent - 1.0_dp)/(2.0_dp*fourpi)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(afld,dphi_tot,f,geff_unit,lb,ub,use_psi)
               DO k = lb(3), ub(3)
                  DO j = lb(2), ub(2)
                     DO i = lb(1), ub(1)
                        afld%array(i, j, k) = f*(dphi_tot(1)%array(i, j, k)**2 + dphi_tot(2)%array(i, j, k)**2 + &
                                                 dphi_tot(3)%array(i, j, k)**2)
                        IF (use_psi) afld%array(i, j, k) = f*geff_unit%array(i, j, k)
                     END DO
                  END DO
               END DO
!$OMP END PARALLEL DO
            END IF
            ! explicit forces of the flexible region mask and of the exclusion mask, from A_diel BEFORE any mask is
            ! applied to afld. dE/dm = A_diel M_diel S_diel + A_ion M_ion S_ion (S before the region mask) is gathered
            ! in flex_b.
            IF (flex_frc) THEN
               CALL auxbas_pw_pool%create_pw(flex_b)
               flex_b%array(:, :, :) = afld%array(:, :, :)*flex_sd0%array(:, :, :)
               IF (excl_on) flex_b%array(:, :, :) = flex_b%array(:, :, :)*excl_md%array(:, :, :)
            END IF
            IF (excl_on .AND. calc_forces) &
               CALL exclusion_force(fex, afld, excl_sd0, excl_md, sccs_control, particles, nex, ex_atom, ex_rd, dvol)
            IF (region_apply) afld%array(:, :, :) = afld%array(:, :, :)*rmask%array(:, :, :)   ! dA/dS_eff, S_eff = m S
            IF (excl_on) afld%array(:, :, :) = afld%array(:, :, :)*excl_md%array(:, :, :)      ! total mask rmask*M_diel
            ! ionic: A = a_ion(phi) (nonlinear) or -eps_b kappa0^2 phi^2/(8 pi) (linear)
            IF (do_electrolyte) THEN
               CALL auxbas_pw_pool%create_pw(aion_f)
               IF (do_nlpb) THEN
                  CALL pw_copy(aion_unit, aion_f)
               ELSE
                  f = -epsilon_solvent*kappa0sq/(2.0_dp*fourpi)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(aion_f,f,lb,phi_tot,ub)
                  DO k = lb(3), ub(3)
                     DO j = lb(2), ub(2)
                        DO i = lb(1), ub(1)
                           aion_f%array(i, j, k) = f*phi_tot%array(i, j, k)**2
                        END DO
                     END DO
                  END DO
!$OMP END PARALLEL DO
               END IF
               IF (flex_frc .AND. excl_on) THEN
                  flex_b%array(:, :, :) = flex_b%array(:, :, :) + &
                                          aion_f%array(:, :, :)*flex_si0%array(:, :, :)*excl_mi%array(:, :, :)
               ELSE IF (flex_frc) THEN
                  flex_b%array(:, :, :) = flex_b%array(:, :, :) + aion_f%array(:, :, :)*flex_si0%array(:, :, :)
               END IF
               IF (excl_on .AND. calc_forces) &
                  CALL exclusion_force(fex, aion_f, excl_si0, excl_mi, sccs_control, particles, nex, ex_atom, ex_ri, dvol)
               IF (region_apply) aion_f%array(:, :, :) = aion_f%array(:, :, :)*rmask%array(:, :, :)
               IF (excl_on) aion_f%array(:, :, :) = aion_f%array(:, :, :)*excl_mi%array(:, :, :)   ! total mask rmask*M_ion
            END IF
            IF (flex_frc) THEN
               flex_b%array(:, :, :) = flex_b%array(:, :, :)*flex_dm%array(:, :, :)   ! dE/d(rho_bar)
               IF (flex_cav) THEN
                  ! rho_bar = g * (1 - S_solv): dE/dS_solv = -(g * dE/d(rho_bar)) joins the cavity chain
                  CALL auxbas_pw_pool%create_pw(cav%a_solv)
                  CALL flex_smooth(auxbas_pw_pool, sccs_control%flex_sigma, flex_b, cav%a_solv)
                  cav%a_solv%array(:, :, :) = -cav%a_solv%array(:, :, :)
                  cav%has_a_solv = .TRUE.
               ELSE
                  ALLOCATE (ffx(3, MAX(nfx, 1)))
                  ffx(:, :) = 0.0_dp
                  CALL flex_force(ffx, flex_b, auxbas_pw_pool, sccs_control, particles, nfx, fx_atom, fx_w)
               END IF
               CALL auxbas_pw_pool%give_back_pw(flex_b)
            END IF
            ! cavitation / dispersion / repulsion: A_cav field computed with the surface above; with ALPHA = BETA = GAMMA = 0
            ! it is identically zero and drops out of the sum
            IF (ABS(epsilon_solvent - 1.0_dp) > epstol) THEN
               IF (ABS(sccs_control%alpha_solvent) + ABS(sccs_control%beta_solvent) + &
                   ABS(sccs_control%gamma_solvent) > 0.0_dp) THEN
                  IF (do_electrolyte) THEN
                     CALL cavity_potential_sum(cav, afld, vchain, auxbas_pw_pool, a_ion=aion_f, a_cav=acav_fld)
                  ELSE
                     CALL cavity_potential_sum(cav, afld, vchain, auxbas_pw_pool, a_cav=acav_fld)
                  END IF
               ELSE
                  IF (do_electrolyte) THEN
                     CALL cavity_potential_sum(cav, afld, vchain, auxbas_pw_pool, a_ion=aion_f)
                  ELSE
                     CALL cavity_potential_sum(cav, afld, vchain, auxbas_pw_pool)
                  END IF
               END IF
               CALL auxbas_pw_pool%give_back_pw(acav_fld)
            ELSE
               IF (do_electrolyte) THEN
                  CALL cavity_potential_sum(cav, afld, vchain, auxbas_pw_pool, a_ion=aion_f)
               ELSE
                  CALL cavity_potential_sum(cav, afld, vchain, auxbas_pw_pool)
               END IF
            END IF
            CALL pw_axpy(vchain, v_sccs, dvol)
            IF (cav%has_a_solv) THEN
               CALL auxbas_pw_pool%give_back_pw(cav%a_solv)
               cav%has_a_solv = .FALSE.
            END IF
            IF (do_electrolyte) CALL auxbas_pw_pool%give_back_pw(aion_f)
            CALL timestop(h_kspot)
            CALL auxbas_pw_pool%give_back_pw(vchain)
            CALL auxbas_pw_pool%give_back_pw(afld)
            ! add the local partial sums of the exclusion force to the pseudo-core force bucket force(ikind)%rho_core.
            ! That bucket is accumulated per MPI rank from the rank's own grid slab (integrate_v_core_rspace) and summed
            ! over ranks later (replicate_qs_force in qs_force_types / qs_forces), so the rank-local partial sum is
            ! added here WITHOUT a para_env%sum (a sum here would multiply the force by the number of ranks).
            IF (excl_on .AND. calc_forces) THEN
               NULLIFY (force_ex)
               CALL get_qs_env(qs_env=qs_env, force=force_ex)
               IF (ASSOCIATED(force_ex)) THEN
                  DO iex = 1, nex
                     force_ex(ex_ikind(iex))%rho_core(1:3, ex_iat(iex)) = &
                        force_ex(ex_ikind(iex))%rho_core(1:3, ex_iat(iex)) + fex(1:3, iex)
                  END DO
               END IF
            END IF
            ! the force of the flexible region mask goes into the same bucket (rank-local partial sums over the
            ! rank's own G vectors, summed over ranks later)
            IF (flex_frc .AND. .NOT. flex_cav) THEN
               NULLIFY (force_ex)
               CALL get_qs_env(qs_env=qs_env, force=force_ex)
               IF (ASSOCIATED(force_ex)) THEN
                  DO iex = 1, nfx
                     force_ex(fx_ikind(iex))%rho_core(1:3, fx_iat(iex)) = &
                        force_ex(fx_ikind(iex))%rho_core(1:3, fx_iat(iex)) + ffx(1:3, iex)
                  END DO
               END IF
            END IF
         END IF

         IF (do_electrolyte .AND. .NOT. do_nlpb) THEN
            ! Ionic charge q_ion = -1/(4 pi) int kap2*phi, ionic (linearised) free energy e_ion = -1/(8 pi) int kap2*phi^2,
            ! and the SCCS contribution of the ionic term to the Kohn-Sham potential:
            ! d/dn [-1/(8 pi) eps_b zeta/lambda_D^2 phi^2] = -1/(8 pi) eps_b/(lambda_D^2 (eps_b-1)) (d eps/dn) phi^2
            ! (d zeta/dn = (d eps/dn)/(eps_b - 1))
            q_ion = 0.0_dp
            e_ion = 0.0_dp
            f_ion = -0.5_dp*dvol*epsilon_solvent*kappa0sq/(fourpi*(epsilon_solvent - 1.0_dp))
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k,dphi2) &
!$OMP             SHARED(deps_elec,dvol,f_ion,kap2,lb,phi_tot,ub,v_sccs) &
!$OMP             REDUCTION(+:e_ion,q_ion)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     dphi2 = phi_tot%array(i, j, k)*phi_tot%array(i, j, k)
                     q_ion = q_ion - kap2%array(i, j, k)*phi_tot%array(i, j, k)
                     e_ion = e_ion - kap2%array(i, j, k)*dphi2
                     v_sccs%array(i, j, k) = v_sccs%array(i, j, k) + f_ion*deps_elec%array(i, j, k)*dphi2
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
            q_ion = q_ion*dvol/fourpi
            e_ion = 0.5_dp*e_ion*dvol/fourpi
            CALL para_env%sum(q_ion)
            CALL para_env%sum(e_ion)
            CALL auxbas_pw_pool%give_back_pw(kap2)
            CALL auxbas_pw_pool%give_back_pw(dkap2e)
            CALL auxbas_pw_pool%give_back_pw(greenfn_scr)
         END IF
         IF (do_nlpb) THEN
            ! Kohn-Sham potential of the nonlinear ionic term: v = (dS_ion/dn) a_ion(phi), dS_ion/dn = (d eps/dn)/(eps_b-1)
            q_ion = q_ion_nl
            e_ion = e_ion_nl
            f_ion = dvol/(epsilon_solvent - 1.0_dp)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(aion_unit,deps_elec,f_ion,lb,ub,v_sccs)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     v_sccs%array(i, j, k) = v_sccs%array(i, j, k) + f_ion*deps_elec%array(i, j, k)*aion_unit%array(i, j, k)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
            IF (sccs_control%nldiel) THEN
               ! nonlinear dielectric: v_diel = (dS_diel/dn) n_mol lambda_diel replaces -(1/8pi)(d eps/dn)|grad phi|^2
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(adiel_unit,deps_elec,f_ion,lb,ub,v_sccs)
               DO k = lb(3), ub(3)
                  DO j = lb(2), ub(2)
                     DO i = lb(1), ub(1)
                        v_sccs%array(i, j, k) = v_sccs%array(i, j, k) + f_ion*deps_elec%array(i, j, k)*adiel_unit%array(i, j, k)
                     END DO
                  END DO
               END DO
!$OMP END PARALLEL DO
            END IF
            CALL auxbas_pw_pool%give_back_pw(adiel_unit)
            CALL auxbas_pw_pool%give_back_pw(aion_unit)
            CALL auxbas_pw_pool%give_back_pw(s_diel_nl)
            CALL auxbas_pw_pool%give_back_pw(s_ion)
            CALL auxbas_pw_pool%give_back_pw(eps_nl)
            IF (do_electrolyte) THEN
               CALL auxbas_pw_pool%give_back_pw(kap2)
               CALL auxbas_pw_pool%give_back_pw(dkap2e)
            END IF
            CALL auxbas_pw_pool%give_back_pw(greenfn_scr)
         END IF

         ! Calculate the Hartree energy and potential of the solute only
         BLOCK
            TYPE(pw_r3d_rs_type) :: phi_solute
            CALL auxbas_pw_pool%create_pw(phi_solute)
            CALL pw_zero(phi_solute)
            IF (level > 1) THEN
               ! E_H = E_H^fine + (A_coarse - E_H^coarse); v_H(fine) = phi_vac(fine) + P^T phi_pol(coarse)
               CALL sol_poisson(auxbas_pw_pool, rho_solute, phi_solute, ehartree=e_h_c)
               energy%hartree = e_h_fine
               energy%sccs_hartree = a_nl + (e_h_fine - e_h_c)
            ELSE
               CALL pw_poisson_solve(poisson_env=poisson_env, &
                                     density=rho_solute, &
                                     ehartree=energy%hartree, &
                                     vhartree=phi_solute)
            END IF

            ! Calculate the polarisation potential (store it in phi_tot)
            ! phi_pol = phi_tot - phi_solute
            CALL pw_axpy(phi_solute, phi_tot, alpha=-1.0_dp)
            CALL auxbas_pw_pool%give_back_pw(phi_solute)
         END BLOCK
         IF (level > 1) THEN
            BLOCK
               TYPE(pw_c1d_gs_type) :: g_c
               CALL auxbas_pw_pool%create_pw(g_c)
               CALL pw_transfer(phi_tot, g_c)
               CALL pw_zero(v_hartree_gspace_out)
               CALL pw_copy(g_c, v_hartree_gspace_out)
               CALL pw_axpy(phi_vac_f, v_hartree_gspace_out)
               CALL auxbas_pw_pool%give_back_pw(g_c)
            END BLOCK
         END IF

         ! Optional output of the SCCS polarisation potential in cube file format
         filename = "POLARISATION_POTENTIAL"
         cube_path = TRIM(print_path)//"%"//TRIM(filename)
         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                   cp_p_file)) THEN
            append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
            my_pos_cube = "REWIND"
            IF (append_cube) my_pos_cube = "APPEND"
            mpi_io = .TRUE.
            cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                             extension=".cube", middle_name=TRIM(filename), &
                                             file_position=my_pos_cube, log_filename=.FALSE., &
                                             mpi_io=mpi_io, fout=mpi_filename)
            IF (output_unit > 0) THEN
               IF (.NOT. mpi_io) THEN
                  INQUIRE (UNIT=cube_unit, NAME=filename)
               ELSE
                  filename = mpi_filename
               END IF
               WRITE (UNIT=output_unit, FMT="(T3,A)") &
                  "SCCS| The SCCS polarisation potential is written in cube file format to the file:", &
                  "SCCS| "//TRIM(filename)
            END IF
            CALL cp_pw_to_cube(phi_tot, cube_unit, TRIM(filename), particles=particles, &
                               stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
         END IF

         ! Calculate the polarisation charge (store it in rho_tot)
         ! rho_pol = rho_tot - rho_solute
         polarisation_charge = 0.0_dp
         IF (need_rho_tot .OR. (.NOT. do_nlpb)) THEN
            CALL pw_axpy(rho_solute, rho_tot, alpha=-1.0_dp)
            polarisation_charge = pw_integrate_function(rho_tot)
         END IF

         ! Optional output of the SCCS polarisation charge in cube file format
         filename = "POLARISATION_CHARGE_DENSITY"
         cube_path = TRIM(print_path)//"%"//TRIM(filename)
         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                   cp_p_file)) THEN
            append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
            my_pos_cube = "REWIND"
            IF (append_cube) my_pos_cube = "APPEND"
            mpi_io = .TRUE.
            cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                             extension=".cube", middle_name=TRIM(filename), &
                                             file_position=my_pos_cube, log_filename=.FALSE., &
                                             mpi_io=mpi_io, fout=mpi_filename)
            IF (output_unit > 0) THEN
               IF (.NOT. mpi_io) THEN
                  INQUIRE (UNIT=cube_unit, NAME=filename)
               ELSE
                  filename = mpi_filename
               END IF
               WRITE (UNIT=output_unit, FMT="(T3,A)") &
                  "SCCS| The SCCS polarisation charge density is written in cube file format to the file:", &
                  "SCCS| "//TRIM(filename)
            END IF
            CALL cp_pw_to_cube(rho_tot, cube_unit, TRIM(filename), particles=particles, &
                               stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
         END IF

         ! Calculate SCCS polarisation energy
         energy%sccs_pol = 0.5_dp*pw_integral_ab(rho_solute, phi_tot)
         IF (do_nlpb) energy%sccs_pol = energy%sccs_hartree - energy%hartree
         ! For a net solute charge Q the reference is the periodic vacuum with a uniform neutralising background,
         ! whose energy grows linearly with the cell length. E_pol of a charged cell is therefore reference
         ! dependent and its sign carries no information about the solver; the solver-state warning is kept for
         ! (nearly) neutral solutes only. Energies are not touched.
         q_net_diag = pw_integrate_function(rho_solute)
         pol_neutral = ABS(q_net_diag) < 1.0E-3_dp
         IF (energy%sccs_pol > 1.0E-6_dp .AND. pol_neutral) THEN
            CPWARN("SCCS: positive polarisation energy, spurious solver state")
         END IF
         IF (.NOT. pol_neutral .AND. .NOT. warned_charged_pol) THEN
            CPWARN("SCCS: charged solute, E_pol depends on the cell length")
            warned_charged_pol = .TRUE.
         END IF
         ! terms of the polarisation energy, printed when the warning above fires and, with
         ! PRINT_POL_DIAGNOSTICS, at every call. Separates a lost operator positivity (min eps, min q of the PSI_FORM
         ! operator), a mismatch of the smoothed and the raw gradient energy (e_grad from A vs the raw field energy) and
         ! a background/charge convention of the G = 0 component (<phi>, rho.phi, ionic charge). Energies are not touched.
         IF (sccs_control%print_pol_diag .OR. (energy%sccs_pol > 1.0E-6_dp .AND. pol_neutral)) THEN
            rho_phi_diag = pw_integral_ab(rho_solute, phi_tot)
            phi_avg_diag = pw_integrate_function(phi_tot)/cell_volume
            e_grad_diag = 0.0_dp
            e_grad_raw = 0.0_dp
            IF (do_nlpb) THEN
               ! A = int rho phi - e_grad + e_diel + e_ion (qs_sccs_nlpb residual), so e_grad follows from A
               e_grad_diag = rho_phi_diag - a_nl + e_diel_nl + e_ion_nl
               IF (.NOT. sccs_control%psi_form) THEN
                  DO l = 1, 3
                     e_grad_raw = e_grad_raw + pw_integral_ab(dphi_tot(l), dphi_tot(l))
                  END DO
                  e_grad_raw = e_grad_raw/(2.0_dp*fourpi)
               END IF
            END IF
            iw = cp_logger_get_default_io_unit(logger)
            ! G = 0 part of the solute-polarisation coupling. phi_tot is the polarisation potential here
            IF (iw > 0) WRITE (UNIT=iw, FMT="(T2,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)") &
               "SCCS| pol-diag: Q_solute [e]", q_net_diag, "  Q*<phi> [Ha]", q_net_diag*phi_avg_diag, &
               "  int rho*phi - Q*<phi>", rho_phi_diag - q_net_diag*phi_avg_diag, &
               "  E_pol - Q*<phi>/2", energy%sccs_pol - 0.5_dp*q_net_diag*phi_avg_diag
            IF (iw > 0) THEN
               WRITE (UNIT=iw, FMT="(T2,A,ES14.6,A,ES14.6,A,ES14.6)") &
                  "SCCS| pol-diag: E_pol [Ha]", energy%sccs_pol, "  A[phi]", energy%sccs_hartree, &
                  "  E_H(vacuum)", energy%hartree
               WRITE (UNIT=iw, FMT="(T2,A,ES14.6,A,ES14.6,A,ES14.6)") &
                  "SCCS| pol-diag: int rho_sol*phi [Ha]", rho_phi_diag, "  <phi> [Ha]", phi_avg_diag, &
                  "  ionic charge [e]", q_ion
               IF (do_nlpb) THEN
                  WRITE (UNIT=iw, FMT="(T2,A,ES14.6,A,ES14.6,A,ES14.6)") &
                     "SCCS| pol-diag: e_grad(from A) [Ha]", e_grad_diag, "  e_diel", e_diel_nl, "  e_ion", e_ion_nl
                  IF (sccs_control%psi_form) THEN
                     WRITE (UNIT=iw, FMT="(T2,A,ES14.6,A,ES14.6)") &
                        "SCCS| pol-diag: min eps", eps_min_diag, "  min q(psi-form operator) [bohr^-2]", nlpb_last_psi_q_min()
                  ELSE
                     WRITE (UNIT=iw, FMT="(T2,A,ES14.6,A,ES14.6)") &
                        "SCCS| pol-diag: min eps", eps_min_diag, "  raw field energy (1/8pi) int |grad phi|^2 [Ha]", e_grad_raw
                  END IF
               END IF
            END IF
         END IF
         CALL auxbas_pw_pool%give_back_pw(rho_solute)
         CALL auxbas_pw_pool%give_back_pw(phi_tot)
         CALL auxbas_pw_pool%give_back_pw(rho_tot)
         IF (sccs_control%method_id == sccs_saa_andreussi) THEN
            CALL auxbas_pw_pool%give_back_pw(s)
            CALL auxbas_pw_pool%give_back_pw(u)
            CALL auxbas_pw_pool%give_back_pw(ff)
            CALL auxbas_pw_pool%give_back_pw(t)
            CALL auxbas_pw_pool%give_back_pw(dtf)
            CALL auxbas_pw_pool%give_back_pw(tem_convolution)
            CALL auxbas_pw_pool%give_back_pw(tem_func)
            CALL auxbas_pw_pool%give_back_pw(eps_elec)
         END IF
      END BLOCK

      ! Calculate additional solvation terms
      energy%sccs_cav = sccs_control%gamma_solvent*cavity_surface
      energy%sccs_dis = sccs_control%beta_solvent*cavity_volume
      energy%sccs_rep = sccs_control%alpha_solvent*cavity_surface
      ! Calculate solvation free energy: \delta G^el + (alpha + gamma)*S + beta*V
      energy%sccs_sol = energy%sccs_pol + energy%sccs_rep + energy%sccs_cav + energy%sccs_dis

      IF (should_output .AND. (output_unit > 0)) THEN
         WRITE (UNIT=output_unit, FMT="(T3,A)") &
            "SCCS|"
         WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.12)") &
            "SCCS| Polarisation charge", polarisation_charge
         IF (do_electrolyte) THEN
            WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.12)") &
               "SCCS| Ionic (electrolyte) charge", q_ion, &
               "SCCS| Ionic free energy (linearised) [a.u.]", e_ion
         END IF
         IF (do_nlpb) THEN
            IF (sccs_control%psi_form) WRITE (UNIT=output_unit, FMT="(T3,A)") &
               "SCCS| electrostatics: psi-form discretisation"
            WRITE (UNIT=output_unit, FMT="(T3,A,T56,I25)") &
               "SCCS| NLPB Newton steps", niter_nl
            WRITE (UNIT=output_unit, FMT="(T3,A,T56,ES25.6)") &
               "SCCS| NLPB final RMS residual [e/bohr^3]", rms_nl
            WRITE (UNIT=output_unit, FMT="(T3,A,T56,I25)") &
               "SCCS| NLPB consecutive non-converged SCF steps", n_nlpb_fail
            WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.12)") &
               "SCCS| NLPB free energy A[phi] [a.u.]", a_nl, &
               "SCCS| NLPB int rho_ion*phi [a.u.]", e_ion_phi, &
               "SCCS| NLPB dielectric free energy [a.u.]", e_diel_nl
         END IF
         !MK   "SCCS| Total interaction energy [a.u.]", e_tot
         WRITE (UNIT=output_unit, FMT="(T3,A)") &
            "SCCS|"
         CALL print_sccs_results(energy, sccs_control, output_unit)
      END IF

      ! Calculate SCCS contribution to the Kohn-Sham potential
      IF (sccs_control%method_id /= sccs_saa_andreussi) THEN

         f = -0.5_dp*dvol/fourpi
         ! nonlinear dielectric: the linear -(1/8pi)(d eps/dn)|grad phi|^2 term is replaced by n_mol lambda_diel dS/dn
         IF (do_nlpb .AND. sccs_control%nldiel) f = 0.0_dp
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(dphi2,i,j,k) &
!$OMP             SHARED(f,deps_elec,dphi_tot,geff_unit,use_psi) &
!$OMP             SHARED(lb,ub,v_sccs)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  dphi2 = dphi_tot(1)%array(i, j, k)*dphi_tot(1)%array(i, j, k) + &
                          dphi_tot(2)%array(i, j, k)*dphi_tot(2)%array(i, j, k) + &
                          dphi_tot(3)%array(i, j, k)*dphi_tot(3)%array(i, j, k)
                  IF (use_psi) dphi2 = geff_unit%array(i, j, k)
                  v_sccs%array(i, j, k) = v_sccs%array(i, j, k) + f*deps_elec%array(i, j, k)*dphi2
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         IF (do_nlpb) CALL auxbas_pw_pool%give_back_pw(geff_unit)

      ELSE

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(lb,ub, dvol) &
!$OMP                SHARED(p_E_interface, v_sccs, d_s_rhoel)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  v_sccs%array(i, j, k) = dvol*d_s_rhoel%array(i, j, k)*p_E_interface%array(i, j, k)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO

         CALL auxbas_pw_pool%give_back_pw(d_s_rhoel)
         CALL auxbas_pw_pool%give_back_pw(p_E_interface)

      END IF

      CALL auxbas_pw_pool%give_back_pw(deps_elec)
      DO i = 1, 3
         CALL auxbas_pw_pool%give_back_pw(dphi_tot(i))
      END DO
      IF (do_ppcav) CALL cavity_release(cav, auxbas_pw_pool)
      IF (region_on) CALL auxbas_pw_pool%give_back_pw(rmask)
      IF (flex_on) CALL auxbas_pw_pool%give_back_pw(flex_dm)
      IF (flex_frc .AND. do_ppcav) THEN
         CALL auxbas_pw_pool%give_back_pw(flex_sd0)
         CALL auxbas_pw_pool%give_back_pw(flex_si0)
      END IF
      IF (excl_on) THEN
         CALL auxbas_pw_pool%give_back_pw(excl_md)
         CALL auxbas_pw_pool%give_back_pw(excl_mi)
         IF (calc_forces) THEN
            CALL auxbas_pw_pool%give_back_pw(excl_sd0)
            CALL auxbas_pw_pool%give_back_pw(excl_si0)
         END IF
      END IF

      ! Release the SCCS printout environment
      CALL cp_print_key_finished_output(output_unit, logger, input, TRIM(print_path), &
                                        ignore_should_output=should_output)

      ! results back on the fine grid (zero padding = P^T; the dvol convention is level independent)
      IF (level > 1) THEN
         CALL sol_prolong(auxbas_pw_pool, fine_pool, v_sccs, v_sccs_out)
         ! v_sccs is the partial derivative dA/dn(r_j) (dvol-scaled); on the fine grid dA/dn_f(r_i) = (N_c/N_f) x its
         ! interpolation (chain rule through the projector), i.e. a factor dvol_f/dvol_c; potentials (v_H) need none
         CALL pw_scale(v_sccs_out, v_sccs_out%pw_grid%dvol/v_sccs%pw_grid%dvol)
         IF (do_nlpb) CALL sol_prolong(auxbas_pw_pool, fine_pool, phi_c, rho_pw_r_sccs_fine)
         CALL pw_zero(rho_tot_gspace_io)
         CALL pw_copy(rho_tot_gspace, rho_tot_gspace_io)
         DO ispin = 1, SIZE(rho_r_c)
            CALL auxbas_pw_pool%give_back_pw(rho_r_c(ispin))
         END DO
         DEALLOCATE (rho_r_c)
         CALL auxbas_pw_pool%give_back_pw(phi_c)
         CALL fine_pool%give_back_pw(phi_vac_f)
      ELSE
         CALL pw_copy(v_sccs, v_sccs_out)
         CALL pw_copy(v_hartree_gspace, v_hartree_gspace_out)
         CALL pw_copy(rho_tot_gspace, rho_tot_gspace_io)
      END IF
      CALL auxbas_pw_pool%give_back_pw(v_sccs)
      CALL auxbas_pw_pool%give_back_pw(v_hartree_gspace)
      CALL auxbas_pw_pool%give_back_pw(rho_tot_gspace)

      ! per-SCF-step solvation cost on the main output; the SCF table's Time column right after it is the
      ! whole step, so the difference is the electronic part
      t_sccs = m_walltime() - t_start
      t_sccs_cum = t_sccs_cum + t_sccs
      iw = cp_logger_get_default_io_unit(logger)
      IF (iw > 0) THEN
         IF (do_nlpb) THEN
            WRITE (UNIT=iw, FMT="(T3,A,F9.2,A,F10.1,A,I3,A,I5,A,ES10.2)") &
               "SCCS| t_solv", t_sccs, " s  cum", t_sccs_cum, " s  Newton", niter_nl, "  CG", ncg_nl, "  RMS", rms_nl
         ELSE
            WRITE (UNIT=iw, FMT="(T3,A,F9.2,A,F10.1,A,I5)") &
               "SCCS| t_solv", t_sccs, " s  cum", t_sccs_cum, " s  iter", iter
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE sccs

! **************************************************************************************************
!> \brief      Calculate the smoothed dielectric function of Andreussi et al.
!> \param rho_elec ...
!> \param eps_elec ...
!> \param deps_elec ...
!> \param epsilon_solvent ...
!> \param rho_max ...
!> \param rho_min ...
!> \par History:
!>      - Creation (16.10.2013,MK)
!>      - Finite difference of isosurfaces implemented (21.12.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.1
! **************************************************************************************************
   SUBROUTINE andreussi(rho_elec, eps_elec, deps_elec, epsilon_solvent, rho_max, &
                        rho_min)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, eps_elec, deps_elec
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, rho_max, rho_min

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'andreussi'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: diff, dq, dt, f, ln_rho_max, ln_rho_min, &
                                                            q, rho, t, x, y

      CALL timeset(routineN, handle)

      f = LOG(epsilon_solvent)/twopi
      diff = rho_max - rho_min
      IF (diff < SQRT(rhotol)) CPABORT("SCCS: Difference between rho(min) and rho(max) is too small")
      IF (rho_min >= rhotol) THEN
         ln_rho_max = LOG(rho_max)
         ln_rho_min = LOG(rho_min)
         q = twopi/(ln_rho_max - ln_rho_min)
         dq = -f*q
      END IF

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

      ! Calculate the dielectric function and its derivative
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(dt,i,j,k,rho,t,x,y) &
!$OMP             SHARED(deps_elec,dq,eps_elec,epsilon_solvent,f,lb,ub) &
!$OMP             SHARED(ln_rho_max,rho_elec,q,rho_max,rho_min)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               rho = rho_elec%array(i, j, k)
               IF (rho < rho_min) THEN
                  eps_elec%array(i, j, k) = epsilon_solvent
                  deps_elec%array(i, j, k) = 0.0_dp
               ELSE IF (rho <= rho_max) THEN
                  x = LOG(rho)
                  y = q*(ln_rho_max - x)
                  t = f*(y - SIN(y))
                  eps_elec%array(i, j, k) = EXP(t)
                  dt = dq*(1.0_dp - COS(y))
                  deps_elec%array(i, j, k) = eps_elec%array(i, j, k)*dt/rho
               ELSE
                  eps_elec%array(i, j, k) = 1.0_dp
                  deps_elec%array(i, j, k) = 0.0_dp
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE andreussi

! **************************************************************************************************
!> \brief      erfc isodensity switching function (Mathew et al., J. Chem. Phys. 140, 084106 (2014), Eq. 10-11):
!>             S(n) = 1/2 erfc(ln(n/nc)/(sigma sqrt 2)), eps = 1 + (eps_b - 1) S, deps/dn = (eps_b - 1) dS/dn with
!>             dS/dn = -exp(-x^2/(2 sigma^2))/(sqrt(2 pi) sigma n), x = ln(n/nc)
!> \param rho_elec ...
!> \param eps_elec ...
!> \param deps_elec ...
!> \param epsilon_solvent ...
!> \param nc ...
!> \param sigma ...
! **************************************************************************************************
   SUBROUTINE isodensity_switch(rho_elec, eps_elec, deps_elec, epsilon_solvent, nc, sigma)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, eps_elec, deps_elec
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, nc, sigma

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'isodensity_switch'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: e, fd, fs, ln_nc, rho, s, x

      CALL timeset(routineN, handle)

      e = epsilon_solvent - 1.0_dp
      ln_nc = LOG(nc)
      fs = 1.0_dp/(sigma*SQRT(2.0_dp))
      fd = -e/(SQRT(twopi)*sigma)

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k,rho,s,x) &
!$OMP             SHARED(deps_elec,e,eps_elec,epsilon_solvent,fd,fs,lb,ln_nc,rho_elec,sigma,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               rho = rho_elec%array(i, j, k)
               IF (rho < rhotol) THEN
                  eps_elec%array(i, j, k) = epsilon_solvent
                  deps_elec%array(i, j, k) = 0.0_dp
               ELSE
                  x = LOG(rho) - ln_nc
                  s = 0.5_dp*ERFC(x*fs)
                  eps_elec%array(i, j, k) = 1.0_dp + e*s
                  deps_elec%array(i, j, k) = fd*EXP(-0.5_dp*(x/sigma)**2)/rho
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE isodensity_switch

! **************************************************************************************************
!> \brief      Quantum surface |grad theta| for the isodensity erfc cavity, theta = 1 - S (analytic dS/dn)
!> \param rho_elec ...
!> \param norm_drho_elec ...
!> \param dtheta ...
!> \param nc ...
!> \param sigma ...
! **************************************************************************************************
   SUBROUTINE isodensity_surface(rho_elec, norm_drho_elec, dtheta, nc, sigma)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, norm_drho_elec, dtheta
      REAL(KIND=dp), INTENT(IN)                          :: nc, sigma

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'isodensity_surface'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: f, ln_nc, rho, x

      CALL timeset(routineN, handle)

      ln_nc = LOG(nc)
      f = 1.0_dp/(SQRT(twopi)*sigma)

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k,rho,x) &
!$OMP             SHARED(dtheta,f,lb,ln_nc,norm_drho_elec,rho_elec,sigma,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               rho = rho_elec%array(i, j, k)
               IF (rho < rhotol) THEN
                  dtheta%array(i, j, k) = 0.0_dp
               ELSE
                  x = LOG(rho) - ln_nc
                  ! |d theta/d n| = |dS/dn| = exp(-x^2/(2 sigma^2))/(sqrt(2 pi) sigma n)
                  dtheta%array(i, j, k) = f*EXP(-0.5_dp*(x/sigma)**2)/rho*norm_drho_elec%array(i, j, k)
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE isodensity_surface

! **************************************************************************************************
!> \brief      Calculate the smoothed dielectric function of Fattebert and Gygi
!> \param rho_elec ...
!> \param eps_elec ...
!> \param deps_elec ...
!> \param epsilon_solvent ...
!> \param beta ...
!> \param rho_zero ...
!> \par History:
!>      - Creation (15.10.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE fattebert_gygi(rho_elec, eps_elec, deps_elec, epsilon_solvent, beta, &
                             rho_zero)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, eps_elec, deps_elec
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, beta, rho_zero

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fattebert_gygi'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: df, f, p, q, rho, s, t, twobeta

      CALL timeset(routineN, handle)

      df = (1.0_dp - epsilon_solvent)/rho_zero
      f = 0.5_dp*(epsilon_solvent - 1.0_dp)
      q = 1.0_dp/rho_zero
      twobeta = 2.0_dp*beta

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

      ! Calculate the smoothed dielectric function and its derivative
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k,p,rho,s,t) &
!$OMP             SHARED(df,deps_elec,eps_elec,epsilon_solvent,f,lb,ub) &
!$OMP             SHARED(q,rho_elec,twobeta)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               rho = rho_elec%array(i, j, k)
               IF (rho < rhotol) THEN
                  eps_elec%array(i, j, k) = epsilon_solvent
                  deps_elec%array(i, j, k) = 0.0_dp
               ELSE
                  s = rho*q
                  p = s**twobeta
                  t = 1.0_dp/(1.0_dp + p)
                  eps_elec%array(i, j, k) = 1.0_dp + f*(1.0_dp + (1.0_dp - p)*t)
                  deps_elec%array(i, j, k) = df*twobeta*t*t*p/s
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE fattebert_gygi

! **************************************************************************************************
!> \brief      Calculate the smoothed dielectric function of solvent aware algorithm
!> \param rho_elec ...
!> \param eps_elec ...
!> \param deps_elec ...
!> \param epsilon_solvent ...
!> \param rho_max ...
!> \param rho_min ...
!> \param s ...
!> \param d_s_rhoel ...
!> \author     Ziwei Chai
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE sa_andreussi(rho_elec, eps_elec, deps_elec, epsilon_solvent, rho_max, &
                           rho_min, s, d_s_rhoel)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, eps_elec, deps_elec
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, rho_max, rho_min
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: s, d_s_rhoel

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'sa_andreussi'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: diff, dq, dt, f, ln_rho_max, ln_rho_min, &
                                                            q, rho, t, x, y

      CALL timeset(routineN, handle)

      f = LOG(epsilon_solvent)/twopi
      diff = rho_max - rho_min
      IF (diff < SQRT(rhotol)) CPABORT("SCCS: Difference between rho(min) and rho(max) is too small")
      IF (rho_min >= rhotol) THEN
         ln_rho_max = LOG(rho_max)
         ln_rho_min = LOG(rho_min)
         q = twopi/(ln_rho_max - ln_rho_min)
         dq = -f*q
      END IF

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

      ! Calculate the dielectric function and its derivative
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(dt,i,j,k,rho,t,x,y) &
!$OMP             SHARED(deps_elec,dq,eps_elec,epsilon_solvent,f,lb,ub) &
!$OMP             SHARED(ln_rho_max,rho_elec,q,rho_max,rho_min, s, d_s_rhoel)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               rho = rho_elec%array(i, j, k)
               IF (rho < rho_min) THEN
                  eps_elec%array(i, j, k) = epsilon_solvent
                  s%array(i, j, k) = 0.0_dp
                  deps_elec%array(i, j, k) = 0.0_dp
                  d_s_rhoel%array(i, j, k) = 0.0_dp
               ELSE IF (rho <= rho_max) THEN
                  x = LOG(rho)
                  y = q*(ln_rho_max - x)
                  t = f*(y - SIN(y))
                  s%array(i, j, k) = 1 - t/(f*twopi)
                  eps_elec%array(i, j, k) = EXP(t)
                  dt = dq*(1.0_dp - COS(y))
                  deps_elec%array(i, j, k) = eps_elec%array(i, j, k)*dt/rho
                  d_s_rhoel%array(i, j, k) = deps_elec%array(i, j, k)/(-eps_elec%array(i, j, k)*LOG(epsilon_solvent))
               ELSE
                  eps_elec%array(i, j, k) = 1.0_dp
                  deps_elec%array(i, j, k) = 0.0_dp
                  s%array(i, j, k) = 1.0_dp
                  d_s_rhoel%array(i, j, k) = 0.0_dp
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE sa_andreussi

! **************************************************************************************************
!> \brief      Build the numerical derivative of a function on realspace grid
!> \param f ...
!> \param df ...
!> \param method ...
!> \param pw_env ...
!> \param input ...
!> \par History:
!>      - Creation (15.11.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE derive(f, df, method, pw_env, input, pool)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: f
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(INOUT)  :: df
      INTEGER, INTENT(IN)                                :: method
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(pw_pool_type), POINTER, OPTIONAL              :: pool   ! grid of f

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'derive'

      INTEGER                                            :: border_points, handle, i
      INTEGER, DIMENSION(3)                              :: lb, n, ub
      TYPE(pw_c1d_gs_type), DIMENSION(2)                 :: work_g1d
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(realspace_grid_desc_type), POINTER            :: rs_desc
      TYPE(realspace_grid_input_type)                    :: input_settings
      TYPE(realspace_grid_type), POINTER                 :: rs_grid
      TYPE(section_vals_type), POINTER                   :: rs_grid_section

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(pw_env))

      ! Perform method specific setup
      SELECT CASE (method)
      CASE (sccs_derivative_cd3, sccs_derivative_cd5, sccs_derivative_cd7)
         NULLIFY (rs_desc)
         rs_grid_section => section_vals_get_subs_vals(input, "DFT%MGRID%RS_GRID")
         SELECT CASE (method)
         CASE (sccs_derivative_cd3)
            border_points = 1
         CASE (sccs_derivative_cd5)
            border_points = 2
         CASE (sccs_derivative_cd7)
            border_points = 3
         END SELECT
         CALL init_input_type(input_settings, 2*border_points + 1, rs_grid_section, &
                              1, [-1, -1, -1])
         CALL rs_grid_create_descriptor(rs_desc, f%pw_grid, input_settings, &
                                        border_points=border_points)
         ALLOCATE (rs_grid)
         CALL rs_grid_create(rs_grid, rs_desc)
!MK      CALL rs_grid_print(rs_grid, 6)
      CASE (sccs_derivative_fft)
         lb(1:3) = f%pw_grid%bounds_local(1, 1:3)
         ub(1:3) = f%pw_grid%bounds_local(2, 1:3)
         NULLIFY (auxbas_pw_pool)
         IF (PRESENT(pool)) THEN
            auxbas_pw_pool => pool
         ELSE
            CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         END IF
         ! Get work storage for the 1d grids in g-space (derivative calculation)
         DO i = 1, SIZE(work_g1d)
            CALL auxbas_pw_pool%create_pw(work_g1d(i))
         END DO
      END SELECT

      ! Calculate the derivatives
      SELECT CASE (method)
      CASE (sccs_derivative_cd3)
         CALL derive_fdm_cd3(f, df, rs_grid)
      CASE (sccs_derivative_cd5)
         CALL derive_fdm_cd5(f, df, rs_grid)
      CASE (sccs_derivative_cd7)
         CALL derive_fdm_cd7(f, df, rs_grid)
      CASE (sccs_derivative_fft)
         ! FFT
         CALL pw_transfer(f, work_g1d(1))
         DO i = 1, 3
            n(:) = 0
            n(i) = 1
            CALL pw_copy(work_g1d(1), work_g1d(2))
            CALL pw_derive(work_g1d(2), n(:))
            CALL pw_transfer(work_g1d(2), df(i))
         END DO
      CASE DEFAULT
         CPABORT("Invalid derivative method for SCCS specified")
      END SELECT

      ! Perform method specific cleanup
      SELECT CASE (method)
      CASE (sccs_derivative_cd3, sccs_derivative_cd5, sccs_derivative_cd7)
         CALL rs_grid_release(rs_grid)
         DEALLOCATE (rs_grid)
         CALL rs_grid_release_descriptor(rs_desc)
      CASE (sccs_derivative_fft)
         DO i = 1, SIZE(work_g1d)
            CALL auxbas_pw_pool%give_back_pw(work_g1d(i))
         END DO
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE derive

! **************************************************************************************************
!> \brief      Calculate the finite difference between two isosurfaces of the
!>             electronic density. The smoothed dielectric function of
!>             Andreussi et al. is used as switching function eventually
!>             defining the quantum volume and surface of the cavity.
!> \param rho_elec ...
!> \param norm_drho_elec ...
!> \param dtheta ...
!> \param epsilon_solvent ...
!> \param rho_max ...
!> \param rho_min ...
!> \param delta_rho ...
!> \par History:
!>      - Creation (21.12.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE surface_andreussi(rho_elec, norm_drho_elec, dtheta, &
                                epsilon_solvent, rho_max, rho_min, delta_rho)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, norm_drho_elec, dtheta
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, rho_max, rho_min, &
                                                            delta_rho

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'surface_andreussi'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k, l
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: diff, e, eps_elec, f, ln_rho_max, &
                                                            ln_rho_min, q, rho, t, x, y
      REAL(KIND=dp), DIMENSION(2)                        :: theta

      CALL timeset(routineN, handle)

      e = epsilon_solvent - 1.0_dp
      f = LOG(epsilon_solvent)/twopi
      diff = rho_max - rho_min
      IF (diff < SQRT(rhotol)) CPABORT("SCCS: Difference between rho(min) and rho(max) is too small")
      IF (rho_min >= rhotol) THEN
         ln_rho_max = LOG(rho_max)
         ln_rho_min = LOG(rho_min)
         q = twopi/(ln_rho_max - ln_rho_min)
      END IF

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

      ! Calculate finite difference between two isosurfaces
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(eps_elec,i,j,k,l,rho,t,theta,x,y) &
!$OMP             SHARED(delta_rho,dtheta,e,epsilon_solvent,f,lb) &
!$OMP             SHARED(ln_rho_max,norm_drho_elec,rho_elec,q,rho_max,rho_min,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               DO l = 1, 2
                  rho = rho_elec%array(i, j, k) + (REAL(l, KIND=dp) - 1.5_dp)*delta_rho
                  IF (rho < rho_min) THEN
                     eps_elec = epsilon_solvent
                  ELSE IF (rho <= rho_max) THEN
                     x = LOG(rho)
                     y = q*(ln_rho_max - x)
                     t = f*(y - SIN(y))
                     eps_elec = EXP(t)
                  ELSE
                     eps_elec = 1.0_dp
                  END IF
                  theta(l) = (epsilon_solvent - eps_elec)/e
               END DO
               dtheta%array(i, j, k) = (theta(2) - theta(1))*norm_drho_elec%array(i, j, k)/delta_rho
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE surface_andreussi

! **************************************************************************************************
!> \brief      Calculate the finite difference between two isosurfaces of the
!>             the electronic density. The smoothed dielectric function of
!>             Fattebert and Gygi is used as switching function eventually
!>             defining the quantum volume and surface of the cavity.
!> \param rho_elec ...
!> \param norm_drho_elec ...
!> \param dtheta ...
!> \param epsilon_solvent ...
!> \param beta ...
!> \param rho_zero ...
!> \param delta_rho ...
!> \par History:
!>      - Creation (21.12.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE surface_fattebert_gygi(rho_elec, norm_drho_elec, dtheta, &
                                     epsilon_solvent, beta, rho_zero, delta_rho)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, norm_drho_elec, dtheta
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, beta, rho_zero, &
                                                            delta_rho

      CHARACTER(LEN=*), PARAMETER :: routineN = 'surface_fattebert_gygi'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k, l
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: e, eps_elec, f, p, q, rho, s, t, twobeta
      REAL(KIND=dp), DIMENSION(2)                        :: theta

      CALL timeset(routineN, handle)

      e = epsilon_solvent - 1.0_dp
      f = 0.5_dp*e
      q = 1.0_dp/rho_zero
      twobeta = 2.0_dp*beta

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

      ! Calculate finite difference between two isosurfaces
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(eps_elec,i,j,k,l,p,rho,s,t,theta) &
!$OMP             SHARED(delta_rho,dtheta,e,epsilon_solvent,f,lb) &
!$OMP             SHARED(norm_drho_elec,q,rho_elec,twobeta,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               DO l = 1, 2
                  rho = rho_elec%array(i, j, k) + (REAL(l, KIND=dp) - 1.5_dp)*delta_rho
                  IF (rho < rhotol) THEN
                     eps_elec = epsilon_solvent
                  ELSE
                     s = rho*q
                     p = s**twobeta
                     t = 1.0_dp/(1.0_dp + p)
                     eps_elec = 1.0_dp + f*(1.0_dp + (1.0_dp - p)*t)
                  END IF
                  theta(l) = (epsilon_solvent - eps_elec)/e
               END DO
               dtheta%array(i, j, k) = (theta(2) - theta(1))*norm_drho_elec%array(i, j, k)/delta_rho
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE surface_fattebert_gygi

! **************************************************************************************************
!> \brief      Print SCCS results
!> \param energy ...
!> \param sccs_control ...
!> \param output_unit ...
!> \par History:
!>      - Creation (11.11.2022,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE print_sccs_results(energy, sccs_control, output_unit)

      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      INTEGER, INTENT(IN)                                :: output_unit

      IF (output_unit > 0) THEN
         CPASSERT(ASSOCIATED(energy))
         CPASSERT(ASSOCIATED(sccs_control))
         WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
            "SCCS| Hartree energy of solute and solvent [Hartree]", energy%sccs_hartree, &
            "SCCS| Hartree energy of the solute only    [Hartree]", energy%hartree
         WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14,/,T3,A,T61,F20.3)") &
            "SCCS| Polarisation energy                  [Hartree]", energy%sccs_pol, &
            "SCCS|                                      [kcal/mol]", &
            cp_unit_from_cp2k(energy%sccs_pol, "kcalmol"), &
            "SCCS| Cavitation energy                    [Hartree]", energy%sccs_cav, &
            "SCCS|                                      [kcal/mol]", &
            cp_unit_from_cp2k(energy%sccs_cav, "kcalmol"), &
            "SCCS| Dispersion free energy               [Hartree]", energy%sccs_dis, &
            "SCCS|                                      [kcal/mol]", &
            cp_unit_from_cp2k(energy%sccs_dis, "kcalmol"), &
            "SCCS| Repulsion free energy                [Hartree]", energy%sccs_rep, &
            "SCCS|                                      [kcal/mol]", &
            cp_unit_from_cp2k(energy%sccs_rep, "kcalmol"), &
            "SCCS| Solvation free energy                [Hartree]", energy%sccs_sol, &
            "SCCS|                                      [kcal/mol]", &
            cp_unit_from_cp2k(energy%sccs_sol, "kcalmol")
      END IF

   END SUBROUTINE print_sccs_results

! **************************************************************************************************
!> \brief Per-kind Gaussian exponent and electron count of the pseudo-core charge added to
!>        the cavity density (ADD_CORE_DENSITY). GPW mapping: exponent = alpha_core_charge of the GTH
!>        pseudopotential (scaled by 1/CORE_DENSITY_WIDTH_SCALE^2), electrons = CORE_DENSITY_CHARGE, or
!>        CORE_DENSITY_SCALE*(Z - Zeff) when that is negative; zero below CORE_DENSITY_ZMIN or for kinds
!>        without a core charge. Used by sccs (density) and by the force integration (qs_ks_methods).
! **************************************************************************************************
   SUBROUTINE sccs_pseudo_core_kinds(qs_kind_set, sccs_control, calpha, ccharge)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: calpha, ccharge

      INTEGER                                            :: ikind, zatom_kind
      REAL(KIND=dp)                                      :: zeff_kind

      DO ikind = 1, SIZE(qs_kind_set)
         CALL get_qs_kind(qs_kind_set(ikind), zatom=zatom_kind, zeff=zeff_kind, &
                          alpha_core_charge=calpha(ikind))
         IF (sccs_control%pseudo_core_charge >= 0.0_dp) THEN
            ccharge(ikind) = sccs_control%pseudo_core_charge
         ELSE
            ccharge(ikind) = sccs_control%pseudo_core_scale*(REAL(zatom_kind, KIND=dp) - zeff_kind)
         END IF
         IF (calpha(ikind) <= 0.0_dp) ccharge(ikind) = 0.0_dp
         IF (zatom_kind < sccs_control%pseudo_core_zmin) ccharge(ikind) = 0.0_dp
         IF (sccs_control%pseudo_core_sigma > 0.0_dp) THEN
            calpha(ikind) = 0.5_dp/sccs_control%pseudo_core_sigma**2   ! alpha = 1/(2 sigma^2), sigma in bohr
         ELSE
            calpha(ikind) = calpha(ikind)/sccs_control%pseudo_core_sigma_scale**2
         END IF
         ! CP2K core-charge amplitudes are the normalised Gaussian prefactors (external_potential_types.F:
         ! ccore_charge = zeff*sqrt((alpha/pi)**3)), so the pseudo-core integrates to the charge set above
         ccharge(ikind) = ccharge(ikind)*SQRT((calpha(ikind)/pi)**3)
      END DO
   END SUBROUTINE sccs_pseudo_core_kinds

! **************************************************************************************************
!> \brief Readings of one field over the REGION window, the whole cell and the PROBE_ATOMS spheres.
!>        q(r) = scale*(fld(r) - offset); reported per domain: integral of q, domain volume, mean = integral/volume and
!>        (want_max) max q. The window is the sharp periodic z-interval [LOWER_EDGE, UPPER_EDGE] with the grid coordinate
!>        (idx - bounds(1))*dr, the same convention as the mask; a sphere uses the minimum-image distance in grid
!>        index space. Collective on all ranks; only the printing rank (iw > 0) writes. No effect on the model.
!> \param pool ...
!> \param para_env ...
!> \param sc ...
!> \param particles ...
!> \param iw ...
!> \param tag ...
!> \param label ...
!> \param fld ...
!> \param offset ...
!> \param scale ...
!> \param want_max ...
! **************************************************************************************************
   SUBROUTINE region_readings(pool, para_env, sc, particles, iw, tag, label, fld, offset, scale, want_max)
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(sccs_control_type), POINTER                   :: sc
      TYPE(particle_list_type), POINTER                  :: particles
      INTEGER, INTENT(IN)                                :: iw
      CHARACTER(LEN=*), INTENT(IN)                       :: tag, label
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: fld
      REAL(KIND=dp), INTENT(IN)                          :: offset, scale
      LOGICAL, INTENT(IN)                                :: want_max

      INTEGER                                            :: i, iax, ia, idom, iimg, j, k, ndom, nprobe
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: dn(3), dvec(3), na(3), npt(3), rr, rsq, vmax, vol, &
                                                            xint, zlen, zpos, zimg
      TYPE(pw_r3d_rs_type)                               :: mask
      CHARACTER(LEN=40)                                  :: dname
      INTEGER                                            :: handle

      CALL timeset('sccs_region_read', handle)

      lb(1:3) = fld%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = fld%pw_grid%bounds_local(2, 1:3)
      iax = sc%region_axis
      zlen = REAL(fld%pw_grid%npts(iax), KIND=dp)*fld%pw_grid%dr(iax)
      npt(1:3) = REAL(fld%pw_grid%npts(1:3), KIND=dp)
      nprobe = 0
      IF (ALLOCATED(sc%region_probe_atoms)) nprobe = SIZE(sc%region_probe_atoms)
      rsq = sc%region_probe_radius**2
      ndom = 1 + nprobe
      na(1:3) = 0.0_dp

      ! whole cell
      vol = fld%pw_grid%vol
      xint = scale*(pw_integrate_function(fld) - offset*vol)
      IF (iw > 0) WRITE (UNIT=iw, FMT="(T2,A,A,A,A,ES12.4,A,ES12.4,A,ES12.4)") &
         "SCCS| region ", TRIM(tag), " ", TRIM(label)//" cell: integral", xint, "  volume [bohr^3]", vol, &
         "  mean", xint/vol

      CALL pool%create_pw(mask)
      DO idom = MERGE(1, 2, sc%region_has_lo), ndom   ! the z window exists only with LOWER_EDGE/UPPER_EDGE
         IF (idom == 1) THEN
            dname = "window"
         ELSE
            ia = sc%region_probe_atoms(idom - 1)
            IF (ia > SIZE(particles%els)) &
               CPABORT("SCCS REGION PROBE_ATOMS: atom index exceeds the number of atoms")
            WRITE (UNIT=dname, FMT="(A,I0)") "probe atom ", ia
            ! atom position in grid-index units (r = dh * n)
            na(1:3) = MATMUL(fld%pw_grid%dh_inv(1:3, 1:3), particles%els(ia)%r(1:3))
         END IF
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,dn,dvec,rr,zpos,zimg,iimg) &
!$OMP             SHARED(fld,iax,idom,lb,mask,na,npt,rsq,sc,ub,zlen)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  mask%array(i, j, k) = 0.0_dp
                  IF (idom == 1) THEN
                     SELECT CASE (iax)
                     CASE (1)
                        rr = REAL(i - fld%pw_grid%bounds(1, 1), KIND=dp)
                     CASE (2)
                        rr = REAL(j - fld%pw_grid%bounds(1, 2), KIND=dp)
                     CASE DEFAULT
                        rr = REAL(k - fld%pw_grid%bounds(1, 3), KIND=dp)
                     END SELECT
                     zpos = MODULO(rr*fld%pw_grid%dr(iax), zlen)
                     DO iimg = -1, 1
                        zimg = zpos + REAL(iimg, KIND=dp)*zlen
                        IF (zimg >= sc%region_zlo .AND. zimg <= sc%region_zhi) mask%array(i, j, k) = 1.0_dp
                     END DO
                  ELSE
                     dn(1) = REAL(i - fld%pw_grid%bounds(1, 1), KIND=dp) - na(1)
                     dn(2) = REAL(j - fld%pw_grid%bounds(1, 2), KIND=dp) - na(2)
                     dn(3) = REAL(k - fld%pw_grid%bounds(1, 3), KIND=dp) - na(3)
                     dn(1:3) = dn(1:3) - npt(1:3)*ANINT(dn(1:3)/npt(1:3))
                     dvec(1:3) = MATMUL(fld%pw_grid%dh(1:3, 1:3), dn(1:3))
                     IF (dvec(1)*dvec(1) + dvec(2)*dvec(2) + dvec(3)*dvec(3) <= rsq) mask%array(i, j, k) = 1.0_dp
                  END IF
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         vol = pw_integrate_function(mask)
         xint = scale*(pw_integral_ab(fld, mask) - offset*vol)
         vmax = -HUGE(1.0_dp)
         IF (want_max) THEN
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     IF (mask%array(i, j, k) > 0.5_dp) vmax = MAX(vmax, scale*(fld%array(i, j, k) - offset))
                  END DO
               END DO
            END DO
            CALL para_env%max(vmax)
         END IF
         IF (iw > 0) THEN
            IF (want_max) THEN
               WRITE (UNIT=iw, FMT="(T2,A,A,A,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)") &
                  "SCCS| region ", TRIM(tag), " ", TRIM(label)//" "//TRIM(dname)//": integral", xint, &
                  "  volume [bohr^3]", vol, "  mean", xint/MAX(vol, 1.0E-30_dp), "  max", vmax
            ELSE
               WRITE (UNIT=iw, FMT="(T2,A,A,A,A,ES12.4,A,ES12.4,A,ES12.4)") &
                  "SCCS| region ", TRIM(tag), " ", TRIM(label)//" "//TRIM(dname)//": integral", xint, &
                  "  volume [bohr^3]", vol, "  mean", xint/MAX(vol, 1.0E-30_dp)
            END IF
         END IF
      END DO
      CALL pool%give_back_pw(mask)
      CALL timestop(handle)

   END SUBROUTINE region_readings

! **************************************************************************************************
!> \brief Resolve the EXCLUSION kind names (and ATOM_LIST) to the list of excluded atoms.
!>        Output per excluded atom: global atom index, kind index, index inside the kind, radii [bohr].
! **************************************************************************************************
   SUBROUTINE exclusion_resolve(atomic_kind_set, natom, what, kind_name, val1, val2, nex, ex_atom, ex_ikind, ex_iat, &
                                ex_rd, ex_ri, atoms)
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      INTEGER, INTENT(IN)                                :: natom
      CHARACTER(LEN=*), INTENT(IN)                       :: what   ! section name for the messages
      CHARACTER(LEN=*), DIMENSION(:), INTENT(IN)         :: kind_name
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: val1, val2   ! per KIND line (radii or weights)
      INTEGER, INTENT(OUT)                               :: nex
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)   :: ex_atom, ex_ikind, ex_iat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: ex_rd, ex_ri
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: atoms   ! ATOM_LIST

      CHARACTER(LEN=default_string_length)               :: kname
      INTEGER                                            :: iat, ikind, ikl, natom_of_kind, nkl, pass
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: kind_used
      LOGICAL                                            :: take

      nkl = SIZE(kind_name)
      ALLOCATE (kind_used(nkl))
      kind_used(:) = .FALSE.
      IF (PRESENT(atoms)) THEN
         IF (ANY(atoms > natom)) &
            CPABORT("SCCS "//what//" ATOM_LIST: atom index exceeds the number of atoms")
      END IF
      nex = 0
      DO pass = 1, 2
         IF (pass == 2) THEN
            ALLOCATE (ex_atom(MAX(nex, 1)), ex_ikind(MAX(nex, 1)), ex_iat(MAX(nex, 1)))
            ALLOCATE (ex_rd(MAX(nex, 1)), ex_ri(MAX(nex, 1)))
            nex = 0
         END IF
         DO ikind = 1, SIZE(atomic_kind_set)
            NULLIFY (atom_list)
            CALL get_atomic_kind(atomic_kind_set(ikind), name=kname, natom=natom_of_kind, atom_list=atom_list)
            kname = exclusion_upper(kname)
            DO ikl = 1, nkl
               IF (TRIM(kname) /= TRIM(kind_name(ikl))) CYCLE
               kind_used(ikl) = .TRUE.
               DO iat = 1, natom_of_kind
                  take = .TRUE.
                  IF (PRESENT(atoms)) take = ANY(atoms == atom_list(iat))
                  IF (.NOT. take) CYCLE
                  nex = nex + 1
                  IF (pass == 2) THEN
                     ex_atom(nex) = atom_list(iat)
                     ex_ikind(nex) = ikind
                     ex_iat(nex) = iat
                     ex_rd(nex) = val1(ikl)
                     ex_ri(nex) = val2(ikl)
                  END IF
               END DO
            END DO
         END DO
      END DO
      DO ikl = 1, nkl
         IF (.NOT. kind_used(ikl)) &
            CPABORT("SCCS "//what//": unknown atomic kind name "//TRIM(kind_name(ikl)))
      END DO
      DEALLOCATE (kind_used)

   END SUBROUTINE exclusion_resolve

   FUNCTION exclusion_upper(str) RESULT(res)
      CHARACTER(LEN=*), INTENT(IN)                       :: str
      CHARACTER(LEN=LEN(str))                            :: res
      INTEGER                                            :: i, ic

      res = str
      DO i = 1, LEN(str)
         ic = IACHAR(str(i:i))
         IF (ic >= IACHAR('a') .AND. ic <= IACHAR('z')) res(i:i) = ACHAR(ic - 32)
      END DO
   END FUNCTION exclusion_upper

! **************************************************************************************************
!> \brief Local grid indices (and their minimum-image offsets from the atom, grid units) of one axis that can
!>        lie within the distance rcut of the atom. Same index convention as the REGION / PROBE_ATOMS code:
!>        point = (idx - bounds(1)) * dh, atom in grid units na = dh_inv * r.
! **************************************************************************************************
   SUBROUTINE exclusion_axis(lo, hi, b1, npt, na, nb, idxl, dnl, n)
      INTEGER, INTENT(IN)                                :: lo, hi, b1, npt, nb
      REAL(KIND=dp), INTENT(IN)                          :: na
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: idxl
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: dnl
      INTEGER, INTENT(OUT)                               :: n

      INTEGER                                            :: ii
      REAL(KIND=dp)                                      :: dn, rn

      n = 0
      rn = REAL(npt, KIND=dp)
      DO ii = lo, hi
         dn = REAL(ii - b1, KIND=dp) - na
         dn = dn - rn*ANINT(dn/rn)
         IF (ABS(dn) <= REAL(nb, KIND=dp)) THEN
            n = n + 1
            idxl(n) = ii
            dnl(n) = dn
         END IF
      END DO
   END SUBROUTINE exclusion_axis

! **************************************************************************************************
!> \brief Does the cache hold the data of exactly this grid / integer / real key (and optionally these positions)?
!>        All comparisons are exact. Called with identical arguments on every rank.
! **************************************************************************************************
   FUNCTION mask_key_match(c, ref, iax, nex, par, pos, rd, ri) RESULT(hit)
      TYPE(mask_cache_type), INTENT(IN)                  :: c
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: ref
      INTEGER, INTENT(IN)                                :: iax, nex
      REAL(KIND=dp), DIMENSION(4), INTENT(IN)            :: par
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), OPTIONAL :: pos
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: rd, ri
      LOGICAL                                            :: hit

      hit = c%valid
      IF (.NOT. hit) RETURN
      hit = (c%id_nr == ref%pw_grid%id_nr) .AND. ALL(c%bounds == ref%pw_grid%bounds) .AND. &
            ALL(c%bounds_local == ref%pw_grid%bounds_local) .AND. ALL(c%dh == ref%pw_grid%dh) .AND. &
            (c%iax == iax) .AND. (c%nex == nex) .AND. ALL(c%par == par)
      IF (hit .AND. PRESENT(pos)) THEN
         hit = ALLOCATED(c%pos)
         IF (hit) hit = ALL(SHAPE(c%pos) == SHAPE(pos))
         IF (hit) hit = ALL(c%pos == pos)
         IF (hit) hit = ALL(c%rd == rd) .AND. ALL(c%ri == ri)
      END IF
   END FUNCTION mask_key_match

! **************************************************************************************************
!> \brief Store key and data (m1 always, m2 optional) in the cache
! **************************************************************************************************
   SUBROUTINE mask_key_store(c, ref, iax, nex, par, m1, m2, pos, rd, ri)
      TYPE(mask_cache_type), INTENT(INOUT)               :: c
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: ref, m1
      INTEGER, INTENT(IN)                                :: iax, nex
      REAL(KIND=dp), DIMENSION(4), INTENT(IN)            :: par
      TYPE(pw_r3d_rs_type), INTENT(IN), OPTIONAL         :: m2
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), OPTIONAL :: pos
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: rd, ri

      IF (ALLOCATED(c%m1)) DEALLOCATE (c%m1)
      IF (ALLOCATED(c%m2)) DEALLOCATE (c%m2)
      IF (ALLOCATED(c%pos)) DEALLOCATE (c%pos)
      IF (ALLOCATED(c%rd)) DEALLOCATE (c%rd)
      IF (ALLOCATED(c%ri)) DEALLOCATE (c%ri)
      ALLOCATE (c%m1(LBOUND(m1%array, 1):UBOUND(m1%array, 1), LBOUND(m1%array, 2):UBOUND(m1%array, 2), &
                     LBOUND(m1%array, 3):UBOUND(m1%array, 3)))
      c%m1(:, :, :) = m1%array(:, :, :)
      IF (PRESENT(m2)) THEN
         ALLOCATE (c%m2(LBOUND(m2%array, 1):UBOUND(m2%array, 1), LBOUND(m2%array, 2):UBOUND(m2%array, 2), &
                        LBOUND(m2%array, 3):UBOUND(m2%array, 3)))
         c%m2(:, :, :) = m2%array(:, :, :)
      END IF
      IF (PRESENT(pos)) THEN
         ALLOCATE (c%pos(SIZE(pos, 1), SIZE(pos, 2)), c%rd(SIZE(rd)), c%ri(SIZE(ri)))
         c%pos(:, :) = pos(:, :); c%rd(:) = rd(:); c%ri(:) = ri(:)
      END IF
      c%id_nr = ref%pw_grid%id_nr; c%bounds = ref%pw_grid%bounds; c%bounds_local = ref%pw_grid%bounds_local
      c%dh = ref%pw_grid%dh; c%iax = iax; c%nex = nex; c%par = par
      c%valid = .TRUE.
   END SUBROUTINE mask_key_store

! **************************************************************************************************
!> \brief Exclusion masks from the cache when the excluded atoms have not moved (exact), else rebuilt
! **************************************************************************************************
   SUBROUTINE exclusion_cached(md, mi, sc, particles, nex, ex_atom, ex_rd, ex_ri)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: md, mi
      TYPE(sccs_control_type), POINTER                   :: sc
      TYPE(particle_list_type), POINTER                  :: particles
      INTEGER, INTENT(IN)                                :: nex
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ex_atom
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: ex_rd, ex_ri

      INTEGER                                            :: ia
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: pos
      REAL(KIND=dp), DIMENSION(4)                        :: par

      ALLOCATE (pos(3, MAX(nex, 1)))
      pos = 0.0_dp
      DO ia = 1, nex
         pos(1:3, ia) = particles%els(ex_atom(ia))%r(1:3)
      END DO
      par = [sc%excl_width, 0.0_dp, 0.0_dp, 0.0_dp]
      IF (mask_key_match(excl_cache, md, 0, nex, par, pos, ex_rd, ex_ri)) THEN
         md%array(:, :, :) = excl_cache%m1(:, :, :)
         mi%array(:, :, :) = excl_cache%m2(:, :, :)
      ELSE
         CALL exclusion_build(md, mi, sc, particles, nex, ex_atom, ex_rd, ex_ri)
         CALL mask_key_store(excl_cache, md, 0, nex, par, md, mi, pos, ex_rd, ex_ri)
      END IF
      DEALLOCATE (pos)
   END SUBROUTINE exclusion_cached

! **************************************************************************************************
!> \brief M_diel(r) = prod_a [1 - E_a(r)] with E_a = 1/2 erfc[(d_a - r_a)/(sqrt(2) W)], and the same for the
!>        ionic radii (M_ion). 1 - E_a is evaluated as 1/2 erfc(-x). Only grid points with d_a <= r_a + 6 W are visited
!>        (outside, E_a = 0 to 6E-10); a radius <= 0 switches that channel off for the atom.
! **************************************************************************************************
   SUBROUTINE exclusion_build(md, mi, sc, particles, nex, ex_atom, ex_rd, ex_ri)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: md, mi
      TYPE(sccs_control_type), POINTER                   :: sc
      TYPE(particle_list_type), POINTER                  :: particles
      INTEGER, INTENT(IN)                                :: nex
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ex_atom
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: ex_rd, ex_ri

      INTEGER                                            :: i, ia, ii, j, jj, k, kk, ax, ni, nj, nk
      INTEGER, DIMENSION(3)                              :: lb, nb, ub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: il, jl, kl
      REAL(KIND=dp)                                      :: d, rcut, w, wsq2, x
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: din, djn, dkn
      REAL(KIND=dp), DIMENSION(3)                        :: dv, na
      INTEGER                                            :: handle

      CALL timeset('sccs_excl_build', handle)

      lb(1:3) = md%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = md%pw_grid%bounds_local(2, 1:3)
      md%array(:, :, :) = 1.0_dp
      mi%array(:, :, :) = 1.0_dp
      w = sc%excl_width
      wsq2 = w*SQRT(2.0_dp)
      ALLOCATE (il(ub(1) - lb(1) + 1), din(ub(1) - lb(1) + 1))
      ALLOCATE (jl(ub(2) - lb(2) + 1), djn(ub(2) - lb(2) + 1))
      ALLOCATE (kl(ub(3) - lb(3) + 1), dkn(ub(3) - lb(3) + 1))
      DO ia = 1, nex
         rcut = MAX(ex_rd(ia), ex_ri(ia)) + 6.0_dp*w
         na(1:3) = MATMUL(md%pw_grid%dh_inv(1:3, 1:3), particles%els(ex_atom(ia))%r(1:3))
         DO ax = 1, 3
            nb(ax) = CEILING(rcut*SQRT(SUM(md%pw_grid%dh_inv(ax, 1:3)**2))) + 1
         END DO
         CALL exclusion_axis(lb(1), ub(1), md%pw_grid%bounds(1, 1), md%pw_grid%npts(1), na(1), nb(1), il, din, ni)
         CALL exclusion_axis(lb(2), ub(2), md%pw_grid%bounds(1, 2), md%pw_grid%npts(2), na(2), nb(2), jl, djn, nj)
         CALL exclusion_axis(lb(3), ub(3), md%pw_grid%bounds(1, 3), md%pw_grid%npts(3), na(3), nb(3), kl, dkn, nk)
         DO kk = 1, nk
            k = kl(kk)
            DO jj = 1, nj
               j = jl(jj)
               DO ii = 1, ni
                  i = il(ii)
                  dv(1:3) = MATMUL(md%pw_grid%dh(1:3, 1:3), [din(ii), djn(jj), dkn(kk)])
                  d = SQRT(dv(1)*dv(1) + dv(2)*dv(2) + dv(3)*dv(3))
                  IF (ex_rd(ia) > 0.0_dp .AND. d <= ex_rd(ia) + 6.0_dp*w) THEN
                     x = (d - ex_rd(ia))/wsq2
                     md%array(i, j, k) = md%array(i, j, k)*0.5_dp*ERFC(-x)
                  END IF
                  IF (ex_ri(ia) > 0.0_dp .AND. d <= ex_ri(ia) + 6.0_dp*w) THEN
                     x = (d - ex_ri(ia))/wsq2
                     mi%array(i, j, k) = mi%array(i, j, k)*0.5_dp*ERFC(-x)
                  END IF
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE (il, jl, kl, din, djn, dkn)
      CALL timestop(handle)

   END SUBROUTINE exclusion_build

! **************************************************************************************************
!> \brief Explicit force of one channel's exclusion mask on the excluded atoms (rank-local partial sums).
!>        E = E[S_eff], S_eff = rmask M S, afld = A = dE/dS_eff (per unit volume, integrated with dvol), sr = rmask S
!>        (unmasked by M). dM/dR_a = -prod_{b/=a}(1-E_b) dE_a/dR_a, dE_a/dR_a = exp(-x^2)/(sqrt(2 pi) W) (r-R_a)/d,
!>        x = (d - r_a)/(sqrt(2) W)  (d/dR_a = -(r-R_a)/d d/dd, dE_a/dd = -exp(-x^2)/(sqrt(2 pi) W)).
!>        F_a = -dE/dR_a = - int A sr dM/dR_a dV = + int A sr prod_{b/=a}(1-E_b) dE_a/dR_a dV,
!>        prod_{b/=a} = M/(1-E_a); a point is skipped if 1-E_a < 1E-12 or d < 1E-10 bohr. Shell |d - r_a| <= 6 W.
!>        fex is incremented by the GRADIENT dE/dR_a = -F_a (INTENT(INOUT)): the qs_force_type arrays hold
!>        gradients (qs_force sets f = -total at the end). The caller adds fex to the force array without an MPI sum.
! **************************************************************************************************
   SUBROUTINE exclusion_force(fex, afld, sr, mfld, sc, particles, nex, ex_atom, rad, dvol)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: fex
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: afld, sr, mfld
      TYPE(sccs_control_type), POINTER                   :: sc
      TYPE(particle_list_type), POINTER                  :: particles
      INTEGER, INTENT(IN)                                :: nex
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ex_atom
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rad
      REAL(KIND=dp), INTENT(IN)                          :: dvol

      INTEGER                                            :: i, ia, ii, j, jj, k, kk, ax, ni, nj, nk
      INTEGER, DIMENSION(3)                              :: lb, nb, ub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: il, jl, kl
      REAL(KIND=dp)                                      :: coef, d, omega, rcut, w, wsq2, x
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: din, djn, dkn
      REAL(KIND=dp), DIMENSION(3)                        :: dv, fsum, na
      INTEGER                                            :: handle

      CALL timeset('sccs_excl_force', handle)

      lb(1:3) = mfld%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = mfld%pw_grid%bounds_local(2, 1:3)
      w = sc%excl_width
      wsq2 = w*SQRT(2.0_dp)
      ALLOCATE (il(ub(1) - lb(1) + 1), din(ub(1) - lb(1) + 1))
      ALLOCATE (jl(ub(2) - lb(2) + 1), djn(ub(2) - lb(2) + 1))
      ALLOCATE (kl(ub(3) - lb(3) + 1), dkn(ub(3) - lb(3) + 1))
      DO ia = 1, nex
         IF (rad(ia) <= 0.0_dp) CYCLE
         rcut = rad(ia) + 6.0_dp*w
         na(1:3) = MATMUL(mfld%pw_grid%dh_inv(1:3, 1:3), particles%els(ex_atom(ia))%r(1:3))
         DO ax = 1, 3
            nb(ax) = CEILING(rcut*SQRT(SUM(mfld%pw_grid%dh_inv(ax, 1:3)**2))) + 1
         END DO
         CALL exclusion_axis(lb(1), ub(1), mfld%pw_grid%bounds(1, 1), mfld%pw_grid%npts(1), na(1), nb(1), il, din, ni)
         CALL exclusion_axis(lb(2), ub(2), mfld%pw_grid%bounds(1, 2), mfld%pw_grid%npts(2), na(2), nb(2), jl, djn, nj)
         CALL exclusion_axis(lb(3), ub(3), mfld%pw_grid%bounds(1, 3), mfld%pw_grid%npts(3), na(3), nb(3), kl, dkn, nk)
         fsum(:) = 0.0_dp
         DO kk = 1, nk
            k = kl(kk)
            DO jj = 1, nj
               j = jl(jj)
               DO ii = 1, ni
                  i = il(ii)
                  dv(1:3) = MATMUL(mfld%pw_grid%dh(1:3, 1:3), [din(ii), djn(jj), dkn(kk)])
                  d = SQRT(dv(1)*dv(1) + dv(2)*dv(2) + dv(3)*dv(3))
                  IF (ABS(d - rad(ia)) > 6.0_dp*w .OR. d < 1.0E-10_dp) CYCLE
                  x = (d - rad(ia))/wsq2
                  omega = 0.5_dp*ERFC(-x)   ! 1 - E_a
                  IF (omega < 1.0E-12_dp) CYCLE
                  coef = afld%array(i, j, k)*sr%array(i, j, k)*(mfld%array(i, j, k)/omega)* &
                         EXP(-x*x)/(SQRT(twopi)*w*d)
                  fsum(1:3) = fsum(1:3) + coef*dv(1:3)
               END DO
            END DO
         END DO
         fex(1:3, ia) = fex(1:3, ia) - fsum(1:3)*dvol   ! fsum is -dE/dR_a; the force arrays hold gradients
      END DO
      DEALLOCATE (il, jl, kl, din, djn, dkn)
      CALL timestop(handle)

   END SUBROUTINE exclusion_force

! **************************************************************************************************
!> \brief Flexible region mask and its density derivative from the cache when the listed atoms have not moved
!>        (exact comparison), else rebuilt
! **************************************************************************************************
   SUBROUTINE flex_cached(m, dm, pool, sc, particles, nfx, fx_atom, fx_w, iw)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: m, dm
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(sccs_control_type), POINTER                   :: sc
      TYPE(particle_list_type), POINTER                  :: particles
      INTEGER, INTENT(IN)                                :: nfx
      INTEGER, DIMENSION(:), INTENT(IN)                  :: fx_atom
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: fx_w
      INTEGER, INTENT(IN)                                :: iw

      INTEGER                                            :: ia
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: pos
      REAL(KIND=dp), DIMENSION(4)                        :: par

      ALLOCATE (pos(3, MAX(nfx, 1)))
      pos = 0.0_dp
      DO ia = 1, nfx
         pos(1:3, ia) = particles%els(fx_atom(ia))%r(1:3)
      END DO
      par = [sc%flex_sigma, sc%flex_rho0, sc%flex_offset, sc%region_width]
      IF (mask_key_match(flex_cache, m, 0, nfx, par, pos, fx_w, fx_w)) THEN
         m%array(:, :, :) = flex_cache%m1(:, :, :)
         dm%array(:, :, :) = flex_cache%m2(:, :, :)
      ELSE
         CALL flex_build(m, dm, pool, sc, particles, nfx, fx_atom, fx_w, iw)
         CALL mask_key_store(flex_cache, m, 0, nfx, par, m, dm, pos, fx_w, fx_w)
      END IF
      DEALLOCATE (pos)
   END SUBROUTINE flex_cached

! **************************************************************************************************
!> \brief Flexible region mask. rho_bar(r) = sum_i w_i g_sigma(r - R_i) over all periodic images is built in
!>        G space, rho_bar(G) = exp(-sigma^2 G^2/2) sum_i w_i exp(-i G.R_i)/V (one inverse FFT; G vectors with
!>        sigma^2 G^2/2 >= 40 are dropped, exp(-40) = 4E-18). With t = erfc^-1(2 rho_bar/rho0),
!>        a = offset/(sqrt2 w), b = sigma/w:  m = 1/2 erfc(a - b t),  dm/drho_bar = -(b/rho0) exp(t^2 - (a - b t)^2).
!>        A uniform half space of density rho0 has rho_bar = rho0/2 erfc(z/(sqrt2 sigma)), z the height above its
!>        Gibbs surface, so t = z/(sqrt2 sigma) and m = 1/2 [1 + erf((z - offset)/(sqrt2 w))]: the planar REGION
!>        edge. b > 1 keeps dm bounded. rho_bar >= rho0 gives m = 0, rho_bar <= 0 gives m = 1 (both dm = 0).
! **************************************************************************************************
   SUBROUTINE flex_build(m, dm, pool, sc, particles, nfx, fx_atom, fx_w, iw)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: m, dm
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(sccs_control_type), POINTER                   :: sc
      TYPE(particle_list_type), POINTER                  :: particles
      INTEGER, INTENT(IN)                                :: nfx
      INTEGER, DIMENSION(:), INTENT(IN)                  :: fx_atom
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: fx_w
      INTEGER, INTENT(IN)                                :: iw

      COMPLEX(KIND=dp)                                   :: sf
      INTEGER                                            :: ia, ig
      REAL(KIND=dp)                                      :: arg, s2h, x, xex, xint
      TYPE(pw_c1d_gs_type)                               :: rg
      INTEGER                                            :: handle

      CALL timeset('sccs_flex_build', handle)

      CALL pool%create_pw(rg)
      s2h = 0.5_dp*sc%flex_sigma**2
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(arg,ia,ig,sf,x) SHARED(fx_atom,fx_w,nfx,particles,rg,s2h)
      DO ig = 1, SIZE(rg%array)
         x = s2h*rg%pw_grid%gsq(ig)
         sf = (0.0_dp, 0.0_dp)
         IF (x < 40.0_dp) THEN
            DO ia = 1, nfx
               arg = DOT_PRODUCT(rg%pw_grid%g(1:3, ig), particles%els(fx_atom(ia))%r(1:3))
               sf = sf + fx_w(ia)*CMPLX(COS(arg), -SIN(arg), KIND=dp)
            END DO
            sf = sf*EXP(-x)/rg%pw_grid%vol
         END IF
         rg%array(ig) = sf
      END DO
!$OMP END PARALLEL DO
      CALL pw_transfer(rg, m)   ! m holds rho_bar
      CALL pool%give_back_pw(rg)
      xint = pw_integrate_function(m)
      CALL flex_map(m, dm, sc)
      xex = m%pw_grid%vol - pw_integrate_function(m)   ! parallel sum: on every rank, not only on the output rank
      IF (iw > 0) WRITE (UNIT=iw, FMT="(T2,A,I8,A,F12.6,A,F12.6,A,ES12.4)") &
         "SCCS| flexible region: atoms", nfx, "  int rho_bar dV", xint, "  sum of weights", SUM(fx_w(1:nfx)), &
         "  int(1-m) dV [bohr^3]", xex
      CALL timestop(handle)

   END SUBROUTINE flex_build

! **************************************************************************************************
!> \brief Pointwise map of the coarse-grained field (in m on entry) to the mask m and to dm/d(field), see flex_build
! **************************************************************************************************
   SUBROUTINE flex_map(m, dm, sc)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: m, dm
      TYPE(sccs_control_type), POINTER                   :: sc

      INTEGER                                            :: i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: a, b, rho0, t, u

      lb(1:3) = m%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = m%pw_grid%bounds_local(2, 1:3)
      rho0 = sc%flex_rho0
      a = sc%flex_offset/(SQRT(2.0_dp)*sc%region_width)
      b = sc%flex_sigma/sc%region_width
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,t,u) SHARED(a,b,dm,lb,m,rho0,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               u = 2.0_dp*m%array(i, j, k)/rho0
               IF (u >= 2.0_dp - 1.0E-14_dp) THEN
                  m%array(i, j, k) = 0.0_dp
                  dm%array(i, j, k) = 0.0_dp
               ELSE IF (u <= 1.0E-280_dp) THEN
                  m%array(i, j, k) = 1.0_dp
                  dm%array(i, j, k) = 0.0_dp
               ELSE
                  t = erfc_inv(u)
                  m%array(i, j, k) = 0.5_dp*ERFC(a - b*t)
                  dm%array(i, j, k) = -(b/rho0)*EXP(t*t - (a - b*t)**2)
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
   END SUBROUTINE flex_map

! **************************************************************************************************
!> \brief g = (normalised Gaussian of width sigma) * f, periodic convolution in G space (two FFTs)
! **************************************************************************************************
   SUBROUTINE flex_smooth(pool, sigma, f, g)
      TYPE(pw_pool_type), POINTER                        :: pool
      REAL(KIND=dp), INTENT(IN)                          :: sigma
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: f
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: g

      INTEGER                                            :: ig
      REAL(KIND=dp)                                      :: s2h
      TYPE(pw_c1d_gs_type)                               :: fg

      CALL pool%create_pw(fg)
      CALL pw_transfer(f, fg)
      s2h = 0.5_dp*sigma**2
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ig) SHARED(fg,s2h)
      DO ig = 1, SIZE(fg%array)
         fg%array(ig) = fg%array(ig)*EXP(-MIN(s2h*fg%pw_grid%gsq(ig), 700.0_dp))
      END DO
!$OMP END PARALLEL DO
      CALL pw_transfer(fg, g)
      CALL pool%give_back_pw(fg)
   END SUBROUTINE flex_smooth

! **************************************************************************************************
!> \brief Flexible region mask from the electron density (FROM_DENSITY): the coarse-grained field is
!>        g * (1 - S_solv), S_solv the solvent cavity function (0 where no solvent molecule fits, 1 in the open
!>        solvent). It is 1 deep inside any matter, so the reference density is 1 and the map of flex_build puts
!>        m = 1/2 at OFFSET from the coarse-grained surface S_solv = 1/2. Rebuilt at every call (the density moves).
! **************************************************************************************************
   SUBROUTINE flex_build_cavity(m, dm, pool, sc, s_solv, iw)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: m, dm
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(sccs_control_type), POINTER                   :: sc
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: s_solv
      INTEGER, INTENT(IN)                                :: iw

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: xex

      CALL timeset('sccs_flex_cav', handle)
      dm%array(:, :, :) = 1.0_dp - s_solv%array(:, :, :)
      CALL flex_smooth(pool, sc%flex_sigma, dm, m)
      CALL flex_map(m, dm, sc)
      xex = m%pw_grid%vol - pw_integrate_function(m)   ! parallel sum: on every rank, not only on the output rank
      IF (iw > 0) WRITE (UNIT=iw, FMT="(T2,A,ES12.4)") &
         "SCCS| flexible region (density): int(1-m) dV [bohr^3]", xex
      CALL timestop(handle)

   END SUBROUTINE flex_build_cavity

! **************************************************************************************************
!> \brief Inverse complementary error function, 0 < p < 2: rational start (Abramowitz & Stegun 26.2.23) and three
!>        Halley steps on erfc(x) - p.
! **************************************************************************************************
   ELEMENTAL FUNCTION erfc_inv(p) RESULT(x)
      REAL(KIND=dp), INTENT(IN)                          :: p
      REAL(KIND=dp)                                      :: x

      INTEGER                                            :: it
      REAL(KIND=dp)                                      :: err, pp, t

      pp = MERGE(p, 2.0_dp - p, p < 1.0_dp)
      t = SQRT(-2.0_dp*LOG(0.5_dp*pp))
      x = -0.70711_dp*((2.30753_dp + t*0.27061_dp)/(1.0_dp + t*(0.99229_dp + t*0.04481_dp)) - t)
      DO it = 1, 3
         err = ERFC(x) - pp
         x = x + err/(2.0_dp/SQRT(pi)*EXP(-x*x) - x*err)
      END DO
      IF (p >= 1.0_dp) x = -x
   END FUNCTION erfc_inv

! **************************************************************************************************
!> \brief Force of the flexible region mask on the listed atoms (rank-local partial sums over the rank's G vectors).
!>        bfld = dE/d(rho_bar) per unit volume, B(G) its plane-wave coefficients. With
!>        rho_bar(G) = exp(-sigma^2 G^2/2) sum_i w_i exp(-i G.R_i)/V and E-variation V sum_G conj(B(G)) d rho_bar(G):
!>        dE/dR_i = -w_i sum_G exp(-sigma^2 G^2/2) G Im[B(G) exp(i G.R_i)]  (sum over all G; a half-space
!>        grid holds one of each +-G pair, hence the factor 2). ffx is incremented by the GRADIENT dE/dR_i: the
!>        qs_force_type arrays hold gradients (qs_force sets f = -total at the end). The caller adds ffx to the
!>        force array without an MPI sum.
! **************************************************************************************************
   SUBROUTINE flex_force(ffx, bfld, pool, sc, particles, nfx, fx_atom, fx_w)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: ffx
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: bfld
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(sccs_control_type), POINTER                   :: sc
      TYPE(particle_list_type), POINTER                  :: particles
      INTEGER, INTENT(IN)                                :: nfx
      INTEGER, DIMENSION(:), INTENT(IN)                  :: fx_atom
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: fx_w

      INTEGER                                            :: ia, ig
      REAL(KIND=dp)                                      :: arg, c, fac, s2h, x
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: fl
      TYPE(pw_c1d_gs_type)                               :: bg
      INTEGER                                            :: handle

      CALL timeset('sccs_flex_force', handle)

      CALL pool%create_pw(bg)
      CALL pw_transfer(bfld, bg)
      fac = 1.0_dp
      IF (bg%pw_grid%grid_span == HALFSPACE) fac = 2.0_dp
      s2h = 0.5_dp*sc%flex_sigma**2
      ALLOCATE (fl(3, MAX(nfx, 1)))
      fl(:, :) = 0.0_dp
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(arg,c,ia,ig,x) SHARED(bg,fac,fx_atom,fx_w,nfx,particles,s2h) REDUCTION(+:fl)
      DO ig = 1, SIZE(bg%array)
         x = s2h*bg%pw_grid%gsq(ig)
         IF (x >= 40.0_dp .OR. bg%pw_grid%gsq(ig) < 1.0E-12_dp) CYCLE
         DO ia = 1, nfx
            arg = DOT_PRODUCT(bg%pw_grid%g(1:3, ig), particles%els(fx_atom(ia))%r(1:3))
            c = fac*fx_w(ia)*EXP(-x)*AIMAG(bg%array(ig)*CMPLX(COS(arg), SIN(arg), KIND=dp))
            fl(1:3, ia) = fl(1:3, ia) + c*bg%pw_grid%g(1:3, ig)
         END DO
      END DO
!$OMP END PARALLEL DO
      ffx(1:3, 1:nfx) = ffx(1:3, 1:nfx) - fl(1:3, 1:nfx)   ! gradient = -force
      DEALLOCATE (fl)
      CALL pool%give_back_pw(bg)
      CALL timestop(handle)

   END SUBROUTINE flex_force

END MODULE qs_sccs
