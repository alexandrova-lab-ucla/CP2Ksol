!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Nonlinear Poisson-Boltzmann (NLPB) solver for the SCCS implicit electrolyte with optional
!>          nonlinear dielectric response and Gaussian smoothing of the polarisation
!>          (equations of Islam et al., J. Chem. Phys. 159, 234117 (2023) Eq. 14, 27-38, 49-57, 70-85 and
!>          Borukhov, Andelman, Orland, Phys. Rev. Lett. 79, 435 (1997)).
!>
!>          -(1/4pi) nabla^2 phi = rho_sol + rho_b{phi} + rho_ion{phi}          (Hartree a.u.)
!>            linear dielectric:    rho_b = (1/4pi) nabla.((eps-1) nabla phi)
!>            nonlinear dielectric: rho_b = -w_b * nabla.(P),  P = n_mol S_diel p(E),  E from qs_sccs_nldiel
!>            rho_ion = -(eps_b kappa0^2/4pi) S_ion g_ion(beta z phi) phi,  g_ion = sinh x/(x D), D = 1-theta_b+theta_b cosh x
!>          Newton: J dphi = R,  J = -(1/4pi)[nabla.(eps_r nabla) - eps_b kappa0^2 S_ion a_ion(phi)],
!>            eps_r = eps (linear) or 1 + w_b * chi{phi} . w_b *  with chi = chi_perp I + (chi_par-chi_perp) g g^T,
!>          inner preconditioned CG (preconditioner = screened Poisson 4pi/(G^2+kappa0^2)), backtracking line search
!>          on the concave free energy
!>            A[phi] = int rho_sol phi - (1/8pi) int |grad phi|^2 (- (1/8pi) int (eps-1)|grad phi|^2 if linear diel.)
!>                     + n_mol int S_diel lambda_diel (nonlinear diel.) + int a_ion(phi),
!>            a_ion = -(2c0/beta) ln(1+theta_b(cosh x-1))/theta_b, 2c0 = eps_b kappa0^2/(4 pi beta z^2).
!> \author  ai-cp2k-lpb (2026-09-10)
! **************************************************************************************************
MODULE qs_sccs_nlpb

   USE ISO_C_BINDING,                   ONLY: C_F_POINTER, C_LOC
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi,&
                                              twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_grid_types,                   ONLY: HALFSPACE
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_derive,&
                                              pw_integral_ab,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE qs_sccs_pw,                      ONLY: sol_convolution,&
                                              sol_poisson
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_sccs_nldiel,                  ONLY: nldiel_params_type,&
                                              nldiel_point
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_sccs_nlpb'

   PUBLIC :: nlpb_params_type, nlpb_solve, nlpb_ion_density, nlpb_reset_scf_memory, nlpb_last_psi_q_min
   PUBLIC :: nlpb_new_scf_run

   TYPE nlpb_params_type
      REAL(KIND=dp) :: kappa0sq = 0.0_dp      ! 1/lambda_D^2 [bohr^-2]
      REAL(KIND=dp) :: eps_b = 1.0_dp         ! bulk relative permittivity
      REAL(KIND=dp) :: beta = 0.0_dp          ! 1/(k_B T) [1/Hartree]
      REAL(KIND=dp) :: zion = 1.0_dp          ! |z| of the symmetric z:z electrolyte
      REAL(KIND=dp) :: theta_b = 0.0_dp       ! bulk site occupation of ions (0 = point ions)
      REAL(KIND=dp) :: eps_resid = 1.0E-6_dp  ! absolute RMS residual tolerance [e/bohr^3]
      REAL(KIND=dp) :: scf_factor = 0.0_dp    ! tol = scf_factor*min_SCF(RMS(R0)) on restarts (typically 0.1); 0 = off
      REAL(KIND=dp) :: eps_floor = 1.0E-11_dp ! hard lower bound of the coupled tolerance (FFT-derivative noise)
      REAL(KIND=dp) :: e_tol = 0.0_dp         ! linear-path CG stops when last 20 energy gains of A sum < e_tol [Ha]; 0 = off
      LOGICAL :: cg_diag = .FALSE.            ! print one line per linear-path CG solve
      REAL(KIND=dp), DIMENSION(8) :: dev_r = 0.0_dp   ! &SCCS&DEVEL developer slots (see input description)
      INTEGER, DIMENSION(8) :: dev_i = 0
      LOGICAL, DIMENSION(8) :: dev_l = .FALSE.
      ! options of the PB solver, input keywords of &SCCS (see their descriptions)
      LOGICAL :: coarse_direct_sum = .TRUE.   ! COARSE_DIRECT_SUM: FFT-free transforms of the coarse correction
      LOGICAL :: lap_carry = .TRUE.           ! LAP_CARRY: lap(s phi) carried from the CG into the residuals
      LOGICAL :: linear_single_cg = .TRUE.    ! LINEAR_SINGLE_CG: ion-free linear Newton step solved in one CG
      INTEGER :: history_guess = 2            ! HISTORY_GUESS: potential increments of the initial guess (0 = off)
      REAL(KIND=dp) :: linear_cg_floor = 0.8_dp   ! LINEAR_CG_FLOOR: CG stops at MAX(0.1 rms, floor*tol); 0 = off
      LOGICAL :: coarse_space_ions = .TRUE.   ! COARSE_SPACE_IONS: two-level preconditioner with linearised ions
      REAL(KIND=dp) :: scf_tol_jump = 100.0_dp ! SCF_TOL_JUMP: restart the tolerance minimum on an entry-RMS jump by this factor
      LOGICAL :: cold_start_from_last = .TRUE. ! COLD_START_FROM_LAST: first solve of an SCF run starts from the last potential
      REAL(KIND=dp) :: newton_cg_forcing = 0.0_dp ! NEWTON_CG_FORCING: nonlinear paths, CG stops at MAX(0.1 rms, f*tol); 0 = off
      REAL(KIND=dp) :: a_smooth = 0.0_dp      ! Gaussian smoothing width of the polarisation, 0 = off [bohr]
      LOGICAL :: nldiel = .FALSE.             ! nonlinear dielectric response
      LOGICAL :: ions = .TRUE.                ! include the ionic term
      LOGICAL :: linear_ions = .FALSE.        ! linearised (Debye-Hueckel) ions instead of the Boltzmann ones
      LOGICAL :: selfcheck = .FALSE.          ! finite-difference checks of J and A on the first Newton step (all ranks)
      LOGICAL :: skip = .FALSE.               ! SOLVER_STRIDE: evaluate at the frozen phi, no iteration
      LOGICAL :: coarse = .FALSE.             ! solver grid coarser than auxbas -> Newton loop for linear models
      LOGICAL :: psi_form = .FALSE.           ! psi-form discretisation, L = s(-lap+q)s
      TYPE(nldiel_params_type) :: dl
      INTEGER :: max_iter = 50
      INTEGER :: max_cg = 200
      REAL(KIND=dp) :: dx_max = 2.0_dp        ! Newton step limiter: max change of x = beta z phi where ions exist [kT]
      INTEGER :: output_unit = -1
      INTEGER :: io_unit = -1          ! default io unit: guard / line-search / breakdown lines print regardless of the print key
   END TYPE nlpb_params_type

   ! work fields of the nonlinear dielectric (module scope to keep the interfaces short)
   TYPE(pw_r3d_rs_type), SAVE :: chi_perp, chi_par, ukern
   TYPE(pw_r3d_rs_type), DIMENSION(3), SAVE :: ghat, pvec, work3b
   TYPE(pw_c1d_gs_type), SAVE :: ukern_g   ! FFT of ukern, fused smoothing
   TYPE(pw_r3d_rs_type), SAVE :: psi_s, psi_q   ! s = sqrt(eps), q = lap(s)/s (psi-form operator)
   LOGICAL, SAVE :: fields_ready = .FALSE.
   REAL(KIND=dp), SAVE :: rms_scf_min = -1.0_dp   ! running minimum of the SCF-induced entry residual
   REAL(KIND=dp), SAVE :: psi_q_min_last = 0.0_dp ! min of q = lap(s)/s of the last PSI_FORM solve (all ranks)
   ! convention of the direct (FFT-free) coarse-space transforms, determined at the first solve:
   ! 0 untested, 1..4 = (sign, origin) pair that reproduces the FFT coefficients, -1 = none did (FFT path)
   INTEGER, SAVE :: cc_dir_mode = 0

   ! history of the converged potentials of the last solves of this SCF run (index 1 = most recent) for the
   ! initial guess phi + sum_i c_i D_i over the differences D_i of consecutive entries (HISTORY_GUESS). A plain array (not
   ! a pw pool object): it survives a rebuilt pw environment and is dropped when the grid key changes.
   INTEGER, PARAMETER :: mh = 4
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :), SAVE :: hx
   INTEGER, SAVE :: n_hx = 0
   INTEGER, DIMENSION(10), SAVE :: h_key = 0
   REAL(KIND=dp), SAVE :: h_vol = 0.0_dp
   ! final potential of the last solve (previous SCF run), start of the first solve of the next run (COLD_START_FROM_LAST)
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), SAVE :: hy
   INTEGER, DIMENSION(10), SAVE :: hy_key = 0
   REAL(KIND=dp), SAVE :: hy_vol = 0.0_dp
   LOGICAL, SAVE :: hy_ok = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief forget the running minimum of the SCF-induced entry residual. Called by the SCF driver at the
!>        start of every SCF run and after every N_e update of the grand-canonical SCF; the next solve then sets
!>        the tolerance from its own entry residual instead of the smallest one seen in an earlier run.
! **************************************************************************************************
   SUBROUTINE nlpb_reset_scf_memory()
      rms_scf_min = -1.0_dp
   END SUBROUTINE nlpb_reset_scf_memory

! **************************************************************************************************
!> \brief a new SCF run (next MD / geometry step) starts; the potential history of the previous run is dropped.
!>        Called once per SCF run by the SCF driver.
! **************************************************************************************************
   SUBROUTINE nlpb_new_scf_run()
      n_hx = 0
   END SUBROUTINE nlpb_new_scf_run

! **************************************************************************************************
!> \brief diagnostics: global minimum of q = lap(s)/s of the last PSI_FORM solve (0 if there was none)
!> \return ...
! **************************************************************************************************
   FUNCTION nlpb_last_psi_q_min() RESULT(qmin)
      REAL(KIND=dp)                                      :: qmin
      qmin = psi_q_min_last
   END FUNCTION nlpb_last_psi_q_min

! **************************************************************************************************
!> \brief a_ion(phi) density, g_ion, and a_ion' = d(x g)/dx at one point (all in Hartree a.u.)
! **************************************************************************************************
   PURE SUBROUTINE ion_response(p, phi, gfac, afac, aion_dens)
      TYPE(nlpb_params_type), INTENT(IN)                 :: p
      REAL(KIND=dp), INTENT(IN)                          :: phi
      REAL(KIND=dp), INTENT(OUT)                         :: gfac, afac, aion_dens

      REAL(KIND=dp)                                      :: c2, ch, d, sh, th, x, xx

      IF (.NOT. p%ions) THEN
         gfac = 0.0_dp; afac = 0.0_dp; aion_dens = 0.0_dp
         RETURN
      END IF
      x = p%beta*p%zion*phi
      th = p%theta_b
      c2 = p%eps_b*p%kappa0sq/(fourpi*p%beta*p%zion**2)   ! 2 c0 [bohr^-3]
      IF (p%linear_ions) THEN
         ! linear PB: g = a' = 1 and a_ion = -(c2/beta) x^2/2 = -eps_b kappa0^2 phi^2/(8 pi)
         gfac = 1.0_dp; afac = 1.0_dp
         aion_dens = -0.5_dp*c2/p%beta*x*x
         RETURN
      END IF
      IF (ABS(x) > 300.0_dp) x = SIGN(300.0_dp, x)        ! guard overflow of cosh/sinh (exp(300) ~ 1e130)
      ch = COSH(x); sh = SINH(x)
      d = 1.0_dp - th + th*ch
      IF (ABS(x) < 1.0E-6_dp) THEN
         gfac = 1.0_dp/d
      ELSE
         gfac = sh/(x*d)
      END IF
      afac = (ch*(1.0_dp - th) + th)/(d*d)
      IF (th > 1.0E-12_dp) THEN
         xx = LOG(1.0_dp + th*(ch - 1.0_dp))/th
      ELSE
         xx = ch - 1.0_dp
      END IF
      aion_dens = -c2/p%beta*xx
   END SUBROUTINE ion_response

! **************************************************************************************************
!> \brief gradient of a real-space field via FFT: df(i) = d f/dx_i
! **************************************************************************************************
   SUBROUTINE gradient(f, df, wg1, wg2)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: f
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(INOUT)  :: df
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: wg1, wg2

      INTEGER                                            :: i
      INTEGER, DIMENSION(3)                              :: n

      CALL pw_transfer(f, wg1)
      DO i = 1, 3
         n(:) = 0; n(i) = 1
         CALL pw_copy(wg1, wg2)
         CALL pw_derive(wg2, n)
         CALL pw_transfer(wg2, df(i))
      END DO
   END SUBROUTINE gradient

! **************************************************************************************************
!> \brief divergence of a real-space vector field via FFT: div = sum_i d v(i)/dx_i
! **************************************************************************************************
   SUBROUTINE divergence(v, div, wg1, wr)
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(IN)     :: v
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: div, wr
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: wg1

      INTEGER                                            :: i
      INTEGER, DIMENSION(3)                              :: n

      CALL pw_zero(div)
      DO i = 1, 3
         n(:) = 0; n(i) = 1
         CALL pw_transfer(v(i), wg1)
         CALL pw_derive(wg1, n)
         CALL pw_transfer(wg1, wr)
         CALL pw_axpy(wr, div)
      END DO
   END SUBROUTINE divergence

! **************************************************************************************************
!> \brief In-place Gaussian smoothing v(l) <- w_b * v(l) of a vector field
! **************************************************************************************************
   SUBROUTINE smooth3(v, pool, wr)
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(INOUT)  :: v
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: wr

      INTEGER                                            :: l

      DO l = 1, 3
         CALL sol_convolution(pool=pool, func=v(l), convolution=wr, u=ukern)
         CALL pw_copy(wr, v(l))
      END DO
   END SUBROUTINE smooth3

! **************************************************************************************************
!> \brief df(l) = F^-1[ i G_l k(G) f(G) ] from the G-space field fg; k = optional smoothing
!>        kernel. 3 inverse FFTs; the separate convolution pass of smooth3 is folded into the multiplier.
! **************************************************************************************************
   SUBROUTINE grad_from_g(fg, df, wg2, kern)
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: fg
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(INOUT)  :: df
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: wg2
      TYPE(pw_c1d_gs_type), INTENT(IN), OPTIONAL         :: kern

      INTEGER                                            :: l
      INTEGER, DIMENSION(3)                              :: n

      DO l = 1, 3
         n(:) = 0; n(l) = 1
         CALL pw_copy(fg, wg2)
         IF (PRESENT(kern)) wg2%array(:) = wg2%array(:)*kern%array(:)
         CALL pw_derive(wg2, n)
         CALL pw_transfer(wg2, df(l))
      END DO
   END SUBROUTINE grad_from_g

! **************************************************************************************************
!> \brief acc(G) += sum_l i G_l k(G) F[v(l)]  (3 forward FFTs). The caller adds any other
!>        G-space term (e.g. -G^2 phi) and performs the single inverse FFT.
! **************************************************************************************************
   SUBROUTINE div_to_g(v, acc, wg2, kern)
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(IN)     :: v
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: acc, wg2
      TYPE(pw_c1d_gs_type), INTENT(IN), OPTIONAL         :: kern

      INTEGER                                            :: l
      INTEGER, DIMENSION(3)                              :: n

      DO l = 1, 3
         n(:) = 0; n(l) = 1
         CALL pw_transfer(v(l), wg2)
         IF (PRESENT(kern)) wg2%array(:) = wg2%array(:)*kern%array(:)
         CALL pw_derive(wg2, n)
         CALL pw_axpy(wg2, acc)
      END DO
   END SUBROUTINE div_to_g

! **************************************************************************************************
!> \brief Normalised Gaussian kernel exp(-r^2/2a^2)/((2pi)^(3/2) a^3), scaled by the cell volume for
!>        pw_func_u_convolution; minimum-image distance from the origin (any cell shape)
! **************************************************************************************************
   SUBROUTINE build_gauss_kernel(w, a)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: w
      REAL(KIND=dp), INTENT(IN)                          :: a

      INTEGER                                            :: i, j, k
      INTEGER, DIMENSION(3)                              :: lb, npts, ub
      INTEGER, DIMENSION(3)                              :: lbg
      INTEGER                                            :: s1, s2, s3
      REAL(KIND=dp)                                      :: f, r2, r2min, wint
      REAL(KIND=dp), DIMENSION(3)                        :: rv, x
      REAL(KIND=dp), DIMENSION(3, 3)                     :: dh

      lb(1:3) = w%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = w%pw_grid%bounds_local(2, 1:3)
      npts(1:3) = w%pw_grid%npts(1:3)
      lbg(1:3) = w%pw_grid%bounds(1, 1:3)
      dh(1:3, 1:3) = w%pw_grid%dh(1:3, 1:3)
      f = w%pw_grid%vol/((twopi)**1.5_dp*a**3)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,r2,x,rv,r2min,s1,s2,s3) SHARED(lbg,a,dh,f,lb,npts,ub,w)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               ! minimum-image distance from the cell origin = array index bounds(1,:) (index 0 is L/2)
               x(1) = REAL(MODULO(i - lbg(1), npts(1)), dp); IF (x(1) > 0.5_dp*npts(1)) x(1) = x(1) - npts(1)
               x(2) = REAL(MODULO(j - lbg(2), npts(2)), dp); IF (x(2) > 0.5_dp*npts(2)) x(2) = x(2) - npts(2)
               x(3) = REAL(MODULO(k - lbg(3), npts(3)), dp); IF (x(3) > 0.5_dp*npts(3)) x(3) = x(3) - npts(3)
               ! true distance for any cell shape: r = sum_i x_i dh(:,i); minimum over the 27 neighbouring images
               r2min = HUGE(1.0_dp)
               DO s3 = -1, 1
                  DO s2 = -1, 1
                     DO s1 = -1, 1
                        rv(1:3) = (x(1) + s1*npts(1))*dh(1:3, 1) + (x(2) + s2*npts(2))*dh(1:3, 2) + &
                                  (x(3) + s3*npts(3))*dh(1:3, 3)
                        r2 = rv(1)*rv(1) + rv(2)*rv(2) + rv(3)*rv(3)
                        IF (r2 < r2min) r2min = r2
                     END DO
                  END DO
               END DO
               w%array(i, j, k) = f*EXP(-0.5_dp*r2min/(a*a))
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
      ! normalisation check: analytic int w d^3r = 1 for the normalised Gaussian.
      wint = pw_integrate_function(w)/w%pw_grid%vol
      IF (ABS(wint - 1.0_dp) > 1.0E-2_dp) &
         CPABORT("SCCS smoothing kernel: int w /= 1 (a too small for the grid?)")
   END SUBROUTINE build_gauss_kernel

! **************************************************************************************************
!> \brief Apply the linearised operator J dphi = -(1/4pi)[nabla.(eps_r nabla dphi) - sa dphi], sa = eps_b kappa0^2 S_ion a_ion
! **************************************************************************************************
   SUBROUTINE apply_j(p, dphi, eps, sa, jdphi, work3, wr, wg1, wg2, wg3, pool, lb, ub)
      TYPE(nlpb_params_type), INTENT(IN)                 :: p
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: dphi, eps, sa
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: jdphi, wr
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(INOUT)  :: work3
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: wg1, wg2, wg3
      TYPE(pw_pool_type), POINTER                        :: pool
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'apply_j'

      INTEGER                                            :: handle, i, ig, j, k, l
      REAL(KIND=dp)                                      :: f, gu

      CALL timeset(routineN, handle)
      IF (p%psi_form) THEN
         ! -L p = s[lap(s p) - q s p]  (2 FFTs); the common tail below turns it into (1/4pi)[L p + sa p]
         wr%array(:, :, :) = psi_s%array(:, :, :)*dphi%array(:, :, :)
         CALL pw_transfer(wr, wg1)
         DO ig = 1, SIZE(wg1%array)
            wg1%array(ig) = -wg1%pw_grid%gsq(ig)*wg1%array(ig)
         END DO
         CALL pw_transfer(wg1, jdphi)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(jdphi,lb,psi_q,psi_s,ub,wr)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  jdphi%array(i, j, k) = psi_s%array(i, j, k)*(jdphi%array(i, j, k) - psi_q%array(i, j, k)*wr%array(i, j, k))
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
      ELSE IF (p%dev_l(5)) THEN
         ! unfused operators (DEV_L(5)): separate convolution passes
         CALL gradient(dphi, work3, wg1, wg2)
         IF (p%nldiel) THEN
            DO l = 1, 3
               CALL pw_copy(work3(l), work3b(l))
            END DO
            IF (p%a_smooth > 0.0_dp) CALL smooth3(work3, pool, wr)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,gu) SHARED(chi_par,chi_perp,ghat,lb,ub,work3)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     gu = ghat(1)%array(i, j, k)*work3(1)%array(i, j, k) + ghat(2)%array(i, j, k)*work3(2)%array(i, j, k) + &
                          ghat(3)%array(i, j, k)*work3(3)%array(i, j, k)
                     gu = (chi_par%array(i, j, k) - chi_perp%array(i, j, k))*gu
                     work3(1)%array(i, j, k) = chi_perp%array(i, j, k)*work3(1)%array(i, j, k) + gu*ghat(1)%array(i, j, k)
                     work3(2)%array(i, j, k) = chi_perp%array(i, j, k)*work3(2)%array(i, j, k) + gu*ghat(2)%array(i, j, k)
                     work3(3)%array(i, j, k) = chi_perp%array(i, j, k)*work3(3)%array(i, j, k) + gu*ghat(3)%array(i, j, k)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
            IF (p%a_smooth > 0.0_dp) CALL smooth3(work3, pool, wr)
            DO l = 1, 3
               CALL pw_axpy(work3b(l), work3(l))
            END DO
         ELSE
            DO l = 1, 3
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(eps,l,lb,ub,work3)
               DO k = lb(3), ub(3)
                  DO j = lb(2), ub(2)
                     DO i = lb(1), ub(1)
                        work3(l)%array(i, j, k) = work3(l)%array(i, j, k)*eps%array(i, j, k)
                     END DO
                  END DO
               END DO
!$OMP END PARALLEL DO
            END DO
         END IF
         CALL divergence(work3, jdphi, wg1, wr)
      ELSE
         ! fused FFT operators: 8 FFTs per matrix-vector product for every dielectric model
         ! (the unfused path needs 10 for scalar eps, 22 for nonlinear dielectric + smoothing)
         CALL pw_transfer(dphi, wg1)
         IF (p%nldiel) THEN
            IF (p%a_smooth > 0.0_dp) THEN
               CALL grad_from_g(wg1, work3, wg2, ukern_g)     ! w_b * grad dphi
            ELSE
               CALL grad_from_g(wg1, work3, wg2)
            END IF
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,gu) SHARED(chi_par,chi_perp,ghat,lb,ub,work3)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  gu = ghat(1)%array(i, j, k)*work3(1)%array(i, j, k) + ghat(2)%array(i, j, k)*work3(2)%array(i, j, k) + &
                       ghat(3)%array(i, j, k)*work3(3)%array(i, j, k)
                  gu = (chi_par%array(i, j, k) - chi_perp%array(i, j, k))*gu
                  work3(1)%array(i, j, k) = chi_perp%array(i, j, k)*work3(1)%array(i, j, k) + gu*ghat(1)%array(i, j, k)
                  work3(2)%array(i, j, k) = chi_perp%array(i, j, k)*work3(2)%array(i, j, k) + gu*ghat(2)%array(i, j, k)
                  work3(3)%array(i, j, k) = chi_perp%array(i, j, k)*work3(3)%array(i, j, k) + gu*ghat(3)%array(i, j, k)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
            ! J dphi (G) = -G^2 dphi(G)  [vacuum part]  +  sum_l i G_l w(G) F[chi . (w * grad dphi)]_l
            DO ig = 1, SIZE(wg3%array)
               wg3%array(ig) = -wg1%pw_grid%gsq(ig)*wg1%array(ig)
            END DO
            IF (p%a_smooth > 0.0_dp) THEN
               CALL div_to_g(work3, wg3, wg2, ukern_g)
            ELSE
               CALL div_to_g(work3, wg3, wg2)
            END IF
         ELSE
            CALL grad_from_g(wg1, work3, wg2)
         DO l = 1, 3
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(eps,l,lb,ub,work3)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     work3(l)%array(i, j, k) = work3(l)%array(i, j, k)*eps%array(i, j, k)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
         END DO
            CALL pw_zero(wg3)
            CALL div_to_g(work3, wg3, wg2)
         END IF
         CALL pw_transfer(wg3, jdphi)
      END IF
      f = -1.0_dp/fourpi
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(dphi,f,jdphi,lb,sa,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               jdphi%array(i, j, k) = f*(jdphi%array(i, j, k) - sa%array(i, j, k)*dphi%array(i, j, k))
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
      CALL timestop(handle)
   END SUBROUTINE apply_j

! **************************************************************************************************
!> \brief Residual R = rho_sol + rho_b + rho_ion + nabla^2 phi/4pi, free energy A[phi], and the field
!>        quantities needed later (rho_ion, aion = S_ion a_ion, adiel = n_mol S_diel lambda_diel, chi, ghat)
! **************************************************************************************************
   SUBROUTINE residual(p, phi, rho_sol, eps, s_ion, s_diel, resid, rho_ion, aion, adiel, a_tot, q_ion, &
                       e_ion, e_diel, work3, wr, wg1, wg2, wg3, pool, lb, ub, lap_out, lap_in)
      TYPE(nlpb_params_type), INTENT(IN)                 :: p
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: phi, rho_sol, eps, s_ion, s_diel
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: resid, rho_ion, aion, adiel, wr
      REAL(KIND=dp), INTENT(OUT)                         :: a_tot, q_ion, e_ion, e_diel
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(INOUT)  :: work3
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: wg1, wg2, wg3
      TYPE(pw_pool_type), POINTER                        :: pool
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub
      TYPE(pw_r3d_rs_type), INTENT(INOUT), OPTIONAL      :: lap_out   ! psi-form only, receives lap(s phi)
      TYPE(pw_r3d_rs_type), INTENT(IN), OPTIONAL         :: lap_in    ! psi-form only, lap(s phi) already known (no FFT)

      REAL(KIND=dp), PARAMETER                           :: gtol = 1.0E-14_dp

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'residual'

      INTEGER                                            :: handle, i, ig, j, k, l
      REAL(KIND=dp)                                      :: ad, afac, cpar, cperp, e_grad, ek, emag, &
                                                            gfac, gmag, gsq, lam, pmag, rmean

      CALL timeset(routineN, handle)
      IF (p%psi_form) THEN
         ! L phi = s A(s phi), A = -lap + q; e_grad = (1/8pi) int phi L phi; R = rho_sol + rho_ion - L phi/4pi
         ek = p%eps_b*p%kappa0sq/fourpi
         wr%array(:, :, :) = psi_s%array(:, :, :)*phi%array(:, :, :)
         IF (PRESENT(lap_in)) THEN
            resid%array(:, :, :) = lap_in%array(:, :, :)                   ! lap(s phi) carried over from the CG
         ELSE
            CALL pw_transfer(wr, wg1)
            DO ig = 1, SIZE(wg1%array)
               wg1%array(ig) = -wg1%pw_grid%gsq(ig)*wg1%array(ig)
            END DO
            CALL pw_transfer(wg1, resid)                                   ! lap(s phi)
         END IF
         IF (PRESENT(lap_out)) lap_out%array(:, :, :) = resid%array(:, :, :)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,gfac,afac,ad) SHARED(aion,ek,lb,p,phi,psi_q,psi_s,resid,rho_ion,s_ion,ub,work3,wr)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  resid%array(i, j, k) = psi_s%array(i, j, k)*(-resid%array(i, j, k) + psi_q%array(i, j, k)*wr%array(i, j, k))
                  work3(1)%array(i, j, k) = phi%array(i, j, k)*resid%array(i, j, k)
                  CALL ion_response(p, phi%array(i, j, k), gfac, afac, ad)
                  rho_ion%array(i, j, k) = -ek*s_ion%array(i, j, k)*gfac*phi%array(i, j, k)
                  aion%array(i, j, k) = s_ion%array(i, j, k)*ad
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         e_grad = 0.5_dp/fourpi*pw_integrate_function(work3(1))
         e_diel = 0.0_dp
         CALL pw_zero(adiel)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(lb,resid,rho_ion,rho_sol,ub)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  resid%array(i, j, k) = rho_sol%array(i, j, k) + rho_ion%array(i, j, k) - resid%array(i, j, k)/fourpi
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
      ELSE IF (p%dev_l(5)) THEN
         ! unfused operators (DEV_L(5)): separate convolution passes
         CALL gradient(phi, work3, wg1, wg2)
         ek = p%eps_b*p%kappa0sq/fourpi
         IF (p%nldiel) THEN
            ! vacuum self-energy from the raw gradient; keep the raw gradient in work3b for nabla^2 phi
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(lb,ub,wr,work3,work3b)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     wr%array(i, j, k) = work3(1)%array(i, j, k)**2 + work3(2)%array(i, j, k)**2 + work3(3)%array(i, j, k)**2
                     work3b(1)%array(i, j, k) = work3(1)%array(i, j, k)
                     work3b(2)%array(i, j, k) = work3(2)%array(i, j, k)
                     work3b(3)%array(i, j, k) = work3(3)%array(i, j, k)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
            e_grad = 0.5_dp/fourpi*pw_integrate_function(wr)
            IF (p%a_smooth > 0.0_dp) CALL smooth3(work3, pool, wr)    ! G = w_b * grad phi
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,gmag,emag,pmag,lam,cperp,cpar,gfac,afac,ad) &
!$OMP             SHARED(adiel,aion,chi_par,chi_perp,ek,ghat,lb,p,phi,pvec,rho_ion,s_diel,s_ion,ub,work3)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     gmag = SQRT(work3(1)%array(i, j, k)**2 + work3(2)%array(i, j, k)**2 + work3(3)%array(i, j, k)**2)
                     CALL nldiel_point(p%dl, gmag, s_diel%array(i, j, k), emag, pmag, lam, cperp, cpar)
                     IF (gmag > gtol) THEN
                        ghat(1)%array(i, j, k) = work3(1)%array(i, j, k)/gmag
                        ghat(2)%array(i, j, k) = work3(2)%array(i, j, k)/gmag
                        ghat(3)%array(i, j, k) = work3(3)%array(i, j, k)/gmag
                     ELSE
                        ghat(1)%array(i, j, k) = 0.0_dp; ghat(2)%array(i, j, k) = 0.0_dp; ghat(3)%array(i, j, k) = 0.0_dp
                     END IF
                     ! macroscopic polarisation P = n_mol S_diel p, antiparallel to G
                     pmag = -p%dl%n_mol*s_diel%array(i, j, k)*pmag
                     pvec(1)%array(i, j, k) = pmag*ghat(1)%array(i, j, k)
                     pvec(2)%array(i, j, k) = pmag*ghat(2)%array(i, j, k)
                     pvec(3)%array(i, j, k) = pmag*ghat(3)%array(i, j, k)
                     chi_perp%array(i, j, k) = cperp
                     chi_par%array(i, j, k) = cpar
                     adiel%array(i, j, k) = p%dl%n_mol*s_diel%array(i, j, k)*lam
                     CALL ion_response(p, phi%array(i, j, k), gfac, afac, ad)
                     rho_ion%array(i, j, k) = -ek*s_ion%array(i, j, k)*gfac*phi%array(i, j, k)
                     aion%array(i, j, k) = s_ion%array(i, j, k)*ad
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
            e_diel = pw_integrate_function(adiel)
            ! bound charge rho_b = -w_b * nabla.P  (stored in resid temporarily)
            CALL divergence(pvec, resid, wg1, wr)
            IF (p%a_smooth > 0.0_dp) THEN
               CALL sol_convolution(pool=pool, func=resid, convolution=wr, u=ukern)
               CALL pw_copy(wr, resid)
            END IF
            CALL pw_scale(resid, -1.0_dp)
            ! + nabla^2 phi/4pi from the raw gradient
            CALL divergence(work3b, wr, wg1, aion)   ! aion used as scratch, restored below
            CALL pw_axpy(wr, resid, 1.0_dp/fourpi)
            ! restore aion scratch (recompute pointwise ion term; cheap)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,gfac,afac,ad) SHARED(aion,lb,p,phi,resid,rho_ion,rho_sol,s_ion,ub)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     CALL ion_response(p, phi%array(i, j, k), gfac, afac, ad)
                     aion%array(i, j, k) = s_ion%array(i, j, k)*ad
                     resid%array(i, j, k) = resid%array(i, j, k) + rho_sol%array(i, j, k) + rho_ion%array(i, j, k)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
         ELSE
            ! linear local dielectric
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,gsq,gfac,afac,ad) &
!$OMP             SHARED(aion,ek,eps,lb,p,phi,rho_ion,s_ion,ub,wr,work3)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     gsq = work3(1)%array(i, j, k)**2 + work3(2)%array(i, j, k)**2 + work3(3)%array(i, j, k)**2
                     wr%array(i, j, k) = eps%array(i, j, k)*gsq
                     CALL ion_response(p, phi%array(i, j, k), gfac, afac, ad)
                     rho_ion%array(i, j, k) = -ek*s_ion%array(i, j, k)*gfac*phi%array(i, j, k)
                     aion%array(i, j, k) = s_ion%array(i, j, k)*ad
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
            e_grad = 0.5_dp/fourpi*pw_integrate_function(wr)
            e_diel = 0.0_dp
            CALL pw_zero(adiel)
            DO l = 1, 3
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(eps,l,lb,ub,work3)
               DO k = lb(3), ub(3)
                  DO j = lb(2), ub(2)
                     DO i = lb(1), ub(1)
                        work3(l)%array(i, j, k) = work3(l)%array(i, j, k)*eps%array(i, j, k)
                     END DO
                  END DO
               END DO
!$OMP END PARALLEL DO
            END DO
            CALL divergence(work3, resid, wg1, wr)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(lb,resid,rho_ion,rho_sol,ub)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     resid%array(i, j, k) = rho_sol%array(i, j, k) + rho_ion%array(i, j, k) + resid%array(i, j, k)/fourpi
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
         END IF
      ELSE IF (p%nldiel) THEN
         ! fused FFT operators (11 FFTs with smoothing; the unfused path needs 24)
         ek = p%eps_b*p%kappa0sq/fourpi   ! eps_b kappa0^2/(4 pi)
         CALL pw_transfer(phi, wg1)
         CALL grad_from_g(wg1, work3b, wg2)                  ! raw gradient -> vacuum self-energy
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(lb,ub,wr,work3b)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  wr%array(i, j, k) = work3b(1)%array(i, j, k)**2 + work3b(2)%array(i, j, k)**2 + work3b(3)%array(i, j, k)**2
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         e_grad = 0.5_dp/fourpi*pw_integrate_function(wr)
         IF (p%a_smooth > 0.0_dp) THEN
            CALL grad_from_g(wg1, work3, wg2, ukern_g)        ! G = w_b * grad phi
         ELSE
            DO l = 1, 3
               CALL pw_copy(work3b(l), work3(l))
            END DO
         END IF
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,gmag,emag,pmag,lam,cperp,cpar,gfac,afac,ad) &
!$OMP             SHARED(adiel,aion,chi_par,chi_perp,ek,ghat,lb,p,phi,pvec,rho_ion,s_diel,s_ion,ub,work3)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  gmag = SQRT(work3(1)%array(i, j, k)**2 + work3(2)%array(i, j, k)**2 + work3(3)%array(i, j, k)**2)
                  CALL nldiel_point(p%dl, gmag, s_diel%array(i, j, k), emag, pmag, lam, cperp, cpar)
                  IF (gmag > gtol) THEN
                     ghat(1)%array(i, j, k) = work3(1)%array(i, j, k)/gmag
                     ghat(2)%array(i, j, k) = work3(2)%array(i, j, k)/gmag
                     ghat(3)%array(i, j, k) = work3(3)%array(i, j, k)/gmag
                  ELSE
                     ghat(1)%array(i, j, k) = 0.0_dp; ghat(2)%array(i, j, k) = 0.0_dp; ghat(3)%array(i, j, k) = 0.0_dp
                  END IF
                  ! macroscopic polarisation P = n_mol S_diel p, antiparallel to G
                  pmag = -p%dl%n_mol*s_diel%array(i, j, k)*pmag
                  pvec(1)%array(i, j, k) = pmag*ghat(1)%array(i, j, k)
                  pvec(2)%array(i, j, k) = pmag*ghat(2)%array(i, j, k)
                  pvec(3)%array(i, j, k) = pmag*ghat(3)%array(i, j, k)
                  chi_perp%array(i, j, k) = cperp
                  chi_par%array(i, j, k) = cpar
                  adiel%array(i, j, k) = p%dl%n_mol*s_diel%array(i, j, k)*lam
                  CALL ion_response(p, phi%array(i, j, k), gfac, afac, ad)
                  rho_ion%array(i, j, k) = -ek*s_ion%array(i, j, k)*gfac*phi%array(i, j, k)
                  aion%array(i, j, k) = s_ion%array(i, j, k)*ad
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         e_diel = pw_integrate_function(adiel)
         ! resid(G) = -w(G) sum_l i G_l P_l(G)  [rho_b]  -  G^2 phi(G)/4pi  [nabla^2 phi/4pi]
         CALL pw_zero(wg3)
         IF (p%a_smooth > 0.0_dp) THEN
            CALL div_to_g(pvec, wg3, wg2, ukern_g)
         ELSE
            CALL div_to_g(pvec, wg3, wg2)
         END IF
         DO ig = 1, SIZE(wg3%array)
            wg3%array(ig) = -wg3%array(ig) - wg1%pw_grid%gsq(ig)*wg1%array(ig)/fourpi
         END DO
         CALL pw_transfer(wg3, resid)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(lb,resid,rho_ion,rho_sol,ub)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  resid%array(i, j, k) = resid%array(i, j, k) + rho_sol%array(i, j, k) + rho_ion%array(i, j, k)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
      ELSE
         ! fused FFT operators, scalar dielectric (8 FFTs; the unfused path needs 10)
         ek = p%eps_b*p%kappa0sq/fourpi
         CALL pw_transfer(phi, wg1)
         CALL grad_from_g(wg1, work3, wg2)
         ! linear local dielectric
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,gsq,gfac,afac,ad) &
!$OMP             SHARED(aion,ek,eps,lb,p,phi,rho_ion,s_ion,ub,wr,work3)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  gsq = work3(1)%array(i, j, k)**2 + work3(2)%array(i, j, k)**2 + work3(3)%array(i, j, k)**2
                  wr%array(i, j, k) = eps%array(i, j, k)*gsq
                  CALL ion_response(p, phi%array(i, j, k), gfac, afac, ad)
                  rho_ion%array(i, j, k) = -ek*s_ion%array(i, j, k)*gfac*phi%array(i, j, k)
                  aion%array(i, j, k) = s_ion%array(i, j, k)*ad
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         e_grad = 0.5_dp/fourpi*pw_integrate_function(wr)
         e_diel = 0.0_dp
         CALL pw_zero(adiel)
         DO l = 1, 3
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(eps,l,lb,ub,work3)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     work3(l)%array(i, j, k) = work3(l)%array(i, j, k)*eps%array(i, j, k)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
         END DO
         CALL pw_zero(wg3)
         CALL div_to_g(work3, wg3, wg2)
         CALL pw_transfer(wg3, resid)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(lb,resid,rho_ion,rho_sol,ub)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  resid%array(i, j, k) = rho_sol%array(i, j, k) + rho_ion%array(i, j, k) + resid%array(i, j, k)/fourpi
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
      END IF
      e_ion = pw_integrate_function(aion)
      q_ion = pw_integrate_function(rho_ion)
      a_tot = pw_integral_ab(rho_sol, phi) - e_grad + e_diel + e_ion
      IF (.NOT. p%ions) THEN
         ! no ions: the periodic Poisson operator has no G = 0 component, so the cell-average residual
         ! (int rho_sol + rho_b, ~1e-7 from the FFT divergence) is not solvable and would leave a residual
         ! floor (about 1e-7 RMS) at which Newton/CG stall. Project it out.
         rmean = pw_integrate_function(resid)/resid%pw_grid%vol
         resid%array(:, :, :) = resid%array(:, :, :) - rmean
      END IF
      CALL timestop(handle)
   END SUBROUTINE residual

! **************************************************************************************************
!> \brief Newton-CG solution of the NLPB equation.
!> \param phi        in: starting guess (restart), out: converged potential [Hartree]
!> \param rho_sol    solute charge density (electrons positive, CP2K convention)
!> \param eps        local dielectric function (linear dielectric path)
!> \param s_ion      ionic cavity; s_diel dielectric cavity (1 in bulk, 0 in the solute)
!> \param greenfn    screened Green function 4pi/(G^2+kappa0^2) (preconditioner)
!> \param a_tot      free energy A[phi] at convergence; q_ion, e_ion, e_ion_phi = int rho_ion phi, e_diel
!> \param aion_unit  a_ion(phi) per unit S_ion (for the KS potential); adiel_unit n_mol lambda_diel per unit S_diel
! **************************************************************************************************
!> \brief preconditioner z = M^-1 r: screened Green function, optionally symmetrically scaled by eps^-1/2
!>        (DEV_L(3)): M^-1 = eps^-1/2 G eps^-1/2 approximates the inverse of -(1/4pi) div eps grad
! **************************************************************************************************
   SUBROUTINE precond(r, z, seps, tmp, use_seps, pool, greenfn)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: r
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: z, tmp
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: seps
      LOGICAL, INTENT(IN)                                :: use_seps
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: greenfn

      IF (use_seps) THEN
         tmp%array(:, :, :) = r%array(:, :, :)*seps%array(:, :, :)
         CALL sol_poisson(pool=pool, density=tmp, vhartree=z, greenfn=greenfn)
         z%array(:, :, :) = z%array(:, :, :)*seps%array(:, :, :)
      ELSE
         CALL sol_poisson(pool=pool, density=r, vhartree=z, greenfn=greenfn)
      END IF
   END SUBROUTINE precond

! **************************************************************************************************
   SUBROUTINE nlpb_solve(p, phi, rho_sol, eps, s_ion, s_diel, greenfn, pool, para_env, &
                         a_tot, q_ion, e_ion, e_ion_phi, e_diel, aion_unit, adiel_unit, geff_unit, rms_out, niter, &
                         converged, ncg_tot)
      TYPE(nlpb_params_type), INTENT(IN)                 :: p
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: phi
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_sol, eps, s_ion, s_diel
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: greenfn
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), INTENT(OUT)                         :: a_tot, q_ion, e_ion, e_ion_phi, e_diel, rms_out
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: aion_unit, adiel_unit, geff_unit
      INTEGER, INTENT(OUT)                               :: niter, ncg_tot
      LOGICAL, INTENT(OUT)                               :: converged

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'nlpb_solve'

      INTEGER                                            :: handle, i, icg, ils, iter, j, k, l, ncg, nstall
      INTEGER, DIMENSION(3)                              :: lb, ub
      LOGICAL                                            :: cold_start, ls_ok, warm_cold, warm_path, cc_ions_here
      REAL(KIND=dp)                                      :: rms_prev
      REAL(KIND=dp)                                      :: a_new, a_old, ad, afac, alpha, beta_cg, &
                                                            cg_tol, cpar, cperp, ek, emag, gfac, &
                                                            gmag, lam, pap, pmag, rms, rz, rz_new, &
                                                            step, tol, vol, xmax
      REAL(KIND=dp), DIMENSION(4)                        :: chk, dmax
      REAL(KIND=dp)                                      :: qmax, qsol
      REAL(KIND=dp)                                      :: rdphi   ! R(phi).dphi before the line search
      REAL(KIND=dp)                                      :: vion, phi_need   ! counter-charge capacity guard
      LOGICAL                                            :: ls_retried
      REAL(KIND=dp), DIMENSION(2)                        :: fmin
      TYPE(pw_c1d_gs_type)                               :: wg1, wg2, wg3
      TYPE(pw_r3d_rs_type)                               :: adiel, aion, cgd, cgp, cgr, cgz, dphi, &
                                                            phi_try, resid, rho_ion, sa, seps, wr
      LOGICAL                                            :: linear_path, use_fis, use_seps
      LOGICAL                                            :: cg_estop, etol_active   ! energy-gain stop of the CG
      REAL(KIND=dp)                                      :: cg_gain_sum, rms_entry_cg, cg_rms_it
      REAL(KIND=dp), DIMENSION(20)                       :: cg_gain_ring
      REAL(KIND=dp), DIMENSION(6)                        :: cg_hist
      INTEGER                                            :: max_newton
      TYPE(pw_r3d_rs_type), DIMENSION(3)                 :: work3
      TYPE(pw_r3d_rs_type)                               :: qfis   ! q = lap(sqrt eps)/sqrt eps
      TYPE(pw_r3d_rs_type)                               :: lap_keep   ! psi-form: lap(s phi) of the last residual()
      LOGICAL                                            :: res_fresh   ! resid/rho_ion/aion/... are those of the current phi
      INTEGER                                            :: h_cg, h_fin
      ! initial guess from the potential history of this SCF run (HISTORY_GUESS), ion-free linear problem only
      LOGICAL                                            :: use_guess, ion_free_lin, a_analytic, rf_final
      REAL(KIND=dp)                                      :: rms_plain
      TYPE(pw_r3d_rs_type), DIMENSION(3)                 :: dv, jd
      ! coarse-space correction (low plane waves s*exp(iGr)) for the psi-form CG; the coarse operator is built
      ! from the Fourier coefficients of eps = s^2 (<e_i, J e_j> = V/(4pi) G_i.G_j eps^(G_i-G_j))
      LOGICAL                                            :: use_cc, cc_ok
      INTEGER                                            :: cc_mxy, cc_mz, cc_n
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: cc_m
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cc_map
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cc_cart
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: cc_l
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: cc_rhs
      ! FFT-free coarse transforms (separable phase sums over the local grid block) and the
      ! lap(s phi) carried from the CG (the closing residual() then needs no FFT)
      LOGICAL                                            :: cc_dir, lapu_ok, use_lapre, lap_stale, lap_new
      INTEGER                                            :: h_ccs, h_ccc
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: cc_ez
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: cc_wt, cc_tt, cc_ss
      TYPE(pw_r3d_rs_type)                               :: lapu, lapg, lapn   ! sum alpha F^-1[G^2 p^], per-iteration part, scratch

      CALL timeset(routineN, handle)

      lb(1:3) = phi%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = phi%pw_grid%bounds_local(2, 1:3)
      vol = phi%pw_grid%vol

      CALL pool%create_pw(wg1); CALL pool%create_pw(wg2); CALL pool%create_pw(wg3)
      CALL pool%create_pw(wr); CALL pool%create_pw(resid); CALL pool%create_pw(rho_ion)
      IF (p%psi_form) CALL pool%create_pw(lap_keep)
      use_lapre = p%psi_form .AND. p%lap_carry
      lapu_ok = .FALSE.
      lap_stale = .FALSE.   ! set when the history guess changes phi (lap_keep then belongs to the old phi)
      lap_new = .FALSE.     ! lapn holds lap(s phi) of the trial point of the last res_try
      cc_dir = .FALSE.
      IF (use_lapre) THEN
         CALL pool%create_pw(lapu); CALL pool%create_pw(lapg); CALL pool%create_pw(lapn)
      END IF
      res_fresh = .FALSE.
      ! ion-free linear dielectric: J is a constant operator. The Newton loop is then a sequence of CG runs on the same
      ! system, which is where the potential history and the CG target floor pay off; with Debye ions the single linear CG
      ! converges in few iterations and the guess is cost-neutral.
      ion_free_lin = (.NOT. p%nldiel) .AND. (.NOT. p%ions)
      use_guess = ion_free_lin .AND. (p%history_guess > 0) .AND. (.NOT. p%selfcheck)
      cc_ions_here = .FALSE.
      a_analytic = .FALSE.
      IF (use_guess) THEN
         CALL reuse_check_key(phi)
         DO l = 1, 3
            CALL pool%create_pw(dv(l)); CALL pool%create_pw(jd(l))
         END DO
      END IF
      CALL pool%create_pw(aion); CALL pool%create_pw(adiel); CALL pool%create_pw(sa); CALL pool%create_pw(dphi)
      CALL pool%create_pw(phi_try); CALL pool%create_pw(cgr); CALL pool%create_pw(cgz)
      CALL pool%create_pw(cgp); CALL pool%create_pw(cgd)
      ! symmetric eps^-1/2 scaling of the preconditioner (DEV_L(3)); Newton cap per call (DEV_I(3))
      ! default on for a scalar eps(r), off with the anisotropic nonlinear dielectric
      use_seps = (p%dev_l(3) .OR. (.NOT. p%nldiel)) .AND. (.NOT. p%dev_l(4))
      ! chi_perp-based symmetric preconditioner for the nonlinear dielectric (DEV_L(8))
      IF (p%nldiel .AND. p%dev_l(8)) use_seps = .TRUE.
      IF (p%psi_form) use_seps = .TRUE.   ! s = 1/seps is the operator itself
      CALL pool%create_pw(seps)
      IF (use_seps) THEN
         seps%array(:, :, :) = 1.0_dp/SQRT(MAX(eps%array(:, :, :), 1.0_dp))
      ELSE
         CALL pw_zero(seps)
      END IF
      ! transformed inner operator K = -(1/4pi)[lap - (q + sa/eps)] for scalar dielectrics
      use_fis = (use_seps .AND. (.NOT. p%nldiel) .AND. p%dev_l(7)) .OR. p%psi_form   ! DEV_L(7) or psi-form
      CALL pool%create_pw(qfis)
      IF (use_fis) THEN
         wr%array(:, :, :) = 1.0_dp/seps%array(:, :, :)             ! s = sqrt(eps)
         CALL pw_transfer(wr, wg1)
         DO i = 1, SIZE(wg1%array)
            wg1%array(i) = -wg1%pw_grid%gsq(i)*wg1%array(i)          ! lap(s) in G space
         END DO
         CALL pw_transfer(wg1, qfis)
         qfis%array(:, :, :) = qfis%array(:, :, :)*seps%array(:, :, :)   ! q = lap(s)/s
         IF (p%psi_form) THEN
            ! the transformed operator is the discretisation itself (residual / Jacobian / energy use s, q)
            CALL pool%create_pw(psi_s); CALL pool%create_pw(psi_q)
            psi_s%array(:, :, :) = wr%array(:, :, :)
            CALL pw_copy(qfis, psi_q)
            psi_q_min_last = MINVAL(psi_q%array)   ! diagnostics (collective, all ranks take this branch)
            CALL para_env%min(psi_q_min_last)
            IF (p%output_unit > 0) WRITE (UNIT=p%output_unit, FMT="(T3,A)") &
               "NLPB| electrostatics: psi-form discretisation, 2 FFTs per matvec"
         ELSE IF (p%output_unit > 0) THEN
            WRITE (UNIT=p%output_unit, FMT="(T3,A)") "NLPB| inner solver: transformed operator, 2 FFTs per matvec"
         END IF
      ELSE
         CALL pw_zero(qfis)
      END IF
      max_newton = p%max_iter
      IF (p%dev_i(3) > 0) max_newton = MIN(p%max_iter, p%dev_i(3))
      DO l = 1, 3
         CALL pool%create_pw(work3(l))
      END DO
      IF (p%nldiel) THEN
         CALL pool%create_pw(chi_perp); CALL pool%create_pw(chi_par)
         DO l = 1, 3
            CALL pool%create_pw(ghat(l)); CALL pool%create_pw(pvec(l)); CALL pool%create_pw(work3b(l))
         END DO
         IF (p%a_smooth > 0.0_dp) THEN
            CALL pool%create_pw(ukern)
            CALL build_gauss_kernel(ukern, p%a_smooth)
            CALL pool%create_pw(ukern_g); CALL pw_transfer(ukern, ukern_g)   ! fused smoothing
         END IF
         fields_ready = .TRUE.
      END IF

      tol = p%eps_resid
      IF (p%ions .AND. p%theta_b > 0.0_dp) THEN
         ! lattice gas: the electrolyte region cannot hold more counter charge than z * (sites) = z * int S_ion / R_ion^3
         ! (R_ion^3 = 3 beta z^2 theta_b/(eps_b kappa0^2))
         qmax = p%zion*pw_integrate_function(s_ion)*p%eps_b*p%kappa0sq/(3.0_dp*p%beta*p%zion**2*p%theta_b)
         qsol = pw_integrate_function(rho_sol)
         IF (ABS(qsol) > 0.99_dp*qmax) THEN
            CPABORT("SCCS NLPB: q_solute exceeds cell ion capacity")
         END IF
      ELSE IF (p%ions .AND. .NOT. p%linear_ions .AND. p%beta > 10.0_dp) THEN
         CPWARN("SCCS NLPB: point ions unbounded in the cavity tail")
      END IF
      ! counter-charge capacity of the electrolyte region. Linear PB holds q_ion = -(eps_b kappa0^2/4pi) int S_ion phi,
      ! so neutralising the solute needs phi ~ |q_sol| 4pi/(eps_b kappa0^2 V_ion). A cell whose implicit region is thinner
      ! than 2 (density edge + R_ion) has V_ion -> 0 and no finite potential compensates the charge: the SCF then runs
      ! away. Report V_ion every solve and stop with the reason instead of diverging.
      IF (p%ions) THEN
         vion = pw_integrate_function(s_ion)
         qsol = pw_integrate_function(rho_sol)
         phi_need = ABS(qsol)*fourpi/(p%eps_b*p%kappa0sq*MAX(vion, 1.0E-30_dp))
         IF (p%io_unit > 0) WRITE (UNIT=p%io_unit, FMT="(T3,A,ES12.4,A,ES12.4,A,ES12.4)") &
            "NLPB| V_ion [bohr^3] ", vion, "  q_solute [e] ", qsol, "  phi needed to compensate (LPB) [Ha] ", phi_need
         ! with V_ion ~ 0 even a neutral cell is ill-posed (the counter-charge potential is pinned by a
         ! vanishing volume): single points look normal but MD does not conserve energy.
         ! Healthy cells keep V_ion >= 10 bohr^3 at every solve; < 1 bohr^3 is empty or runaway.
         IF (vion < 1.0_dp) THEN
            CPABORT("SCCS NLPB: empty electrolyte (V_ion < 1 bohr^3), drop DEBYE_LENGTH")
         ELSE IF (phi_need > 1.0_dp) THEN
            CPABORT("SCCS NLPB: no ion volume for the counter charge, lengthen the cell")
         ELSE IF (phi_need > 0.2_dp) THEN
            CPWARN("SCCS NLPB: counter charge needs > 5 V, ion volume small")
         END IF
      END IF
      converged = .FALSE.
      niter = 0
      nstall = 0
      ls_retried = .FALSE.
      rms_prev = HUGE(1.0_dp)
      ek = p%eps_b*p%kappa0sq

      ! previous-run potential as start of a cold solve (COLD_START_FROM_LAST): the tolerance stays that of a cold start
      cold_start = (pw_integral_ab(phi, phi) <= 0.0_dp)
      warm_cold = .FALSE.
      ! COLD_START_FROM_LAST acts on the linearised-ion problem only (the others keep their own start)
      warm_path = p%cold_start_from_last .AND. (.NOT. p%nldiel) .AND. p%ions .AND. p%linear_ions
      IF (cold_start .AND. warm_path .AND. hy_ok .AND. ALLOCATED(hy)) THEN
         warm_cold = ALL(hy_key(1:3) == phi%pw_grid%bounds_local(1, 1:3)) .AND. &
                     ALL(hy_key(4:6) == phi%pw_grid%bounds_local(2, 1:3)) .AND. ALL(hy_key(7:9) == phi%pw_grid%npts(1:3)) .AND. &
                     (hy_key(10) == phi%pw_grid%id_nr) .AND. ABS(phi%pw_grid%vol - hy_vol) <= 1.0E-10_dp*hy_vol
         IF (warm_cold) phi%array(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)) = hy
      END IF
      CALL residual(p, phi, rho_sol, eps, s_ion, s_diel, resid, rho_ion, aion, adiel, a_old, q_ion, e_ion, e_diel, &
                    work3, wr, wg1, wg2, wg3, pool, lb, ub, lap_keep)
      res_fresh = .TRUE.
      rms = SQRT(pw_integral_ab(resid, resid)/vol)
      ! SCF coupling (tol = SCF_TOL_FACTOR x SCF residual): on a restart R0 = rho_sol - rho_sol(previous SCF step), so the
      ! entry residual is the SCF-induced perturbation; the first call (phi = 0) converges to eps_resid
      ! cold start = the restart field is exactly zero (first activated SCF step, or a rebuilt grid).
      ! pw_integral_ab is collective, so it must be evaluated on every rank unconditionally.
      ! Tolerance schedule (tol = SCF_TOL_FACTOR x SCF residual): the entry residual is the SCF-induced
      ! perturbation; its running minimum never lets the tolerance loosen again; no absolute floor except
      ! eps_floor, so the solver follows the SCF down to EPS_SCF. Cold start: SOLVER_TOLERANCE.
      IF (cold_start) THEN
         rms_scf_min = -1.0_dp
         tol = p%eps_resid
      ELSE IF (p%scf_factor > 0.0_dp) THEN
         IF (rms_scf_min < 0.0_dp) THEN
            rms_scf_min = rms
         ELSE IF (p%scf_tol_jump > 0.0_dp .AND. rms > p%scf_tol_jump*rms_scf_min) THEN
            ! SCF_TOL_JUMP: the entry residual jumped by more than this factor above the running minimum (an N_e update or
            ! another abrupt change of the solute density reached the solver after the memory was reset): follow it
            IF (p%io_unit > 0) WRITE (UNIT=p%io_unit, FMT="(T3,A,ES12.4,A,ES12.4)") &
               "NLPB| entry RMS jumped, tolerance memory restarted: ", rms_scf_min, " -> ", rms
            rms_scf_min = rms
         ELSE
            rms_scf_min = MIN(rms_scf_min, rms)
         END IF
         tol = MAX(p%eps_floor, p%scf_factor*rms_scf_min)
         IF (p%dev_r(2) > 0.0_dp) tol = MAX(p%dev_r(2), p%scf_factor*rms_scf_min)
      ELSE
         tol = p%eps_resid
      END IF
      IF (p%output_unit > 0) WRITE (UNIT=p%output_unit, FMT="(T3,A,ES12.4,A,ES12.4,A,ES12.4)") &
         "NLPB| initial RMS residual ", rms, "  tol ", tol, "  rms_scf_min ", rms_scf_min
      rms_plain = rms

      ! with a linear dielectric and linearised ions J is exact and constant, so one
      ! preconditioned CG straight to the coupled tolerance replaces the Newton loop (no restarts, no line search:
      ! A is a concave quadratic and CG increases it monotonically, no step limiter). DEV_L(6) selects the
      ! Newton loop; NLPB_SELFCHECK needs it for the finite-difference checks.
      ! Ions are required: the ion-free problem has a residual floor that only the Newton loop handles.
      linear_path = (.NOT. p%nldiel) .AND. (p%ions .AND. p%linear_ions) .AND. (.NOT. p%dev_l(6)) .AND. &
                    (.NOT. p%selfcheck) .AND. ((.NOT. use_fis) .OR. p%psi_form) .AND. (.NOT. p%coarse)
      ! the transformed operator and a coarser solver grid use the Newton loop; the psi-form operator is the
      ! discretisation itself, so it may take the linear path
      ncg_tot = 0
      cg_estop = .FALSE.; etol_active = .FALSE.; cg_hist = -1.0_dp; cg_gain_sum = 0.0_dp
      IF (p%skip) THEN
         ! SOLVER_STRIDE: intermittent relaxation - keep phi, report the entry residual
         iter = 0
         IF (p%output_unit > 0) WRITE (UNIT=p%output_unit, FMT="(T3,A,ES12.4)") &
            "NLPB| solve skipped (SOLVER_STRIDE), phi frozen, entry RMS ", rms
         converged = .TRUE.
      ELSE IF (linear_path) THEN
         iter = 1
         IF (rms > tol) THEN
            niter = 1
         ! sa = eps_b kappa0^2 S_ion a_ion(phi)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,gfac,afac,ad) SHARED(ek,lb,p,phi,s_ion,sa,ub)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  CALL ion_response(p, phi%array(i, j, k), gfac, afac, ad)
                  sa%array(i, j, k) = ek*s_ion%array(i, j, k)*afac
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
            ! COARSE_SPACE_IONS: two-level preconditioner for the linearised-ion problem (this path only)
            cc_ions_here = p%coarse_space_ions
            rms_entry_cg = rms
            etol_active = (p%e_tol > 0.0_dp)
            res_fresh = .FALSE.
            CALL timeset('nlpb_cg', h_cg)
            CALL run_cg(tol, 5*p%max_cg)
            CALL timestop(h_cg)
            etol_active = .FALSE.
            CALL pw_axpy(dphi, phi)
            CALL res_try(phi, 1.0_dp, a_new)
            IF (lap_new) lap_keep%array(:, :, :) = lapn%array(:, :, :)
            res_fresh = .TRUE.
            a_old = a_new
            rms = SQRT(pw_integral_ab(resid, resid)/vol)
            IF (p%output_unit > 0) WRITE (UNIT=p%output_unit, FMT="(T3,A,I5,A,ES12.4,A,F20.10)") &
               "NLPB| linear CG ", MIN(icg, 5*p%max_cg), "  RMS ", rms, "  A ", a_new
            ! the linear problem is solved when the CG energy gain has run dry, even if the RMS residual has
            ! stalled above the SCF-coupled tolerance (round-off / grid-scale modes at the cavity edge)
            IF (rms > tol .AND. cg_estop) converged = .TRUE.
            IF (p%cg_diag .AND. p%io_unit > 0) &
               WRITE (UNIT=p%io_unit, FMT="(T3,A,I5,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,A,A,6ES9.2)") &
               "NLPB| linear CG n ", MIN(icg, 5*p%max_cg), "  entry RMS ", rms_entry_cg, "  final RMS ", rms, "  tol ", tol, &
               "  last-20 energy gain [Ha] ", cg_gain_sum, "  stop ", &
               MERGE("energy", MERGE("rms   ", "cap   ", rms <= tol), cg_estop), &
               "  RMS@50,100,200,400,800,1000 ", cg_hist
            IF (rms > tol .AND. (.NOT. p%ions) .AND. (.NOT. p%dev_l(1)) .AND. &
                rms < MERGE(p%dev_r(1), 1.0E-6_dp, p%dev_r(1) > 0.0_dp)) THEN
               IF (p%output_unit > 0) WRITE (UNIT=p%output_unit, FMT="(T3,A,ES12.4)") &
                  "NLPB| residual floor reached, accepting RMS ", rms
               converged = .TRUE.
            END IF
         END IF
      ELSE
      IF (use_guess .AND. rms > tol) THEN
         ! initial guess from the potential history (the Jacobian of the linear problem is the same operator in every
         ! solve; sa is evaluated here because the guess applies J)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,gfac,afac,ad) SHARED(ek,lb,p,phi,s_ion,sa,ub)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  CALL ion_response(p, phi%array(i, j, k), gfac, afac, ad)
                  sa%array(i, j, k) = ek*s_ion%array(i, j, k)*afac
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         CALL reuse_guess()
      END IF
      newton: DO iter = 1, max_newton
         ! Make at least one Newton step per SCF call (the density changed); skip only if the entry
         ! residual is already far below the tolerance
         IF (rms <= tol .AND. (iter > 1 .OR. rms <= 0.01_dp*tol .OR. p%dev_l(2))) THEN
            converged = .TRUE.
            EXIT newton
         END IF
         niter = iter
         res_fresh = .FALSE.
         ! sa = eps_b kappa0^2 S_ion a_ion(phi)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,gfac,afac,ad) SHARED(ek,lb,p,phi,s_ion,sa,ub)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  CALL ion_response(p, phi%array(i, j, k), gfac, afac, ad)
                  sa%array(i, j, k) = ek*s_ion%array(i, j, k)*afac
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         ! ---- self-check of the Jacobian against the residual (first Newton step only): (R(phi+h p)-R(phi))/h = -J p
         ! NOTE: every operation here is collective; only the WRITEs are rank-0
         IF (iter == 1 .AND. p%selfcheck) THEN
            CALL pw_copy(resid, cgr)
            CALL sol_poisson(pool=pool, density=cgr, vhartree=cgz, greenfn=greenfn)   ! p = M^-1 R
            CALL apply_j(p, cgz, eps, sa, cgd, work3, wr, wg1, wg2, wg3, pool, lb, ub)               ! J p
            step = 1.0E-4_dp
            CALL pw_copy(phi, phi_try)
            CALL pw_axpy(cgz, phi_try, step)
            CALL residual(p, phi_try, rho_sol, eps, s_ion, s_diel, cgp, rho_ion, aion, adiel, a_new, q_ion, e_ion, &
                          e_diel, work3, wr, wg1, wg2, wg3, pool, lb, ub, lap_keep)
            CALL pw_axpy(resid, cgp, -1.0_dp)          ! cgp = R(phi+hp) - R(phi)
            CALL pw_scale(cgp, 1.0_dp/step)            ! finite-difference derivative
            CALL pw_axpy(cgd, cgp, 1.0_dp)             ! + J p  (should vanish)
            chk(1) = SQRT(pw_integral_ab(cgp, cgp)/vol); chk(2) = SQRT(pw_integral_ab(cgd, cgd)/vol)
            chk(3) = SQRT(pw_integral_ab(cgz, cgz)/vol); chk(4) = pw_integral_ab(resid, cgz)
            dmax(1) = MAXVAL(ABS(cgz%array)); dmax(2) = MAXVAL(ABS(cgd%array)); dmax(3) = MAXVAL(ABS(resid%array))
            dmax(4) = MAXVAL(ABS(cgp%array)); CALL para_env%max(dmax)
            IF (p%output_unit > 0) THEN
               WRITE (UNIT=p%output_unit, FMT="(T3,A,4ES14.5)") "NLPB| check max|p|, max|Jp|, max|R|, max|dR/dh+Jp|: ", dmax
               WRITE (UNIT=p%output_unit, FMT="(T3,A,3ES14.5)") "NLPB| check |dR/dh + Jp|, |Jp|, |p| RMS: ", chk(1:3)
               ! also: (dA/dh)/(R.p) should be 1  (R is the gradient of A)
               WRITE (UNIT=p%output_unit, FMT="(T3,A,2ES14.5)") "NLPB| check dA/dh vs R.p: ", (a_new - a_old)/step, chk(4)
            END IF
            ! residual() above left the module-level nonlinear-dielectric fields (chi_perp, chi_par, ghat,
            ! pvec, work3b) evaluated at phi + h*p. Re-evaluate at phi so that the production Jacobian below
            ! is the one at phi: a debug switch must not move the numbers.
            CALL residual(p, phi, rho_sol, eps, s_ion, s_diel, cgp, rho_ion, aion, adiel, a_new, q_ion, e_ion, &
                          e_diel, work3, wr, wg1, wg2, wg3, pool, lb, ub, lap_keep)
         END IF
         ! ---- preconditioned CG for J dphi = resid (inexact Newton: 0.1*rms, DEV_R(3)) ----
         cg_tol = MERGE(p%dev_r(3), 0.1_dp, p%dev_r(3) > 0.0_dp)*rms
         ! LINEAR_CG_FLOOR: a linear ion-free problem needs no tighter solve than the SCF-coupled tolerance itself
         IF (ion_free_lin .AND. p%linear_cg_floor > 0.0_dp) cg_tol = MAX(cg_tol, p%linear_cg_floor*tol)
         ! NEWTON_CG_FORCING: nonlinear path, the Newton loop checks the residual after the step, so a tighter CG
         ! target cannot change what it accepts; cg_tol = MAX(0.1 rms, f tol)
         IF (p%newton_cg_forcing > 0.0_dp .AND. (p%nldiel .OR. (p%ions .AND. (.NOT. p%linear_ions)))) &
            cg_tol = MAX(cg_tol, p%newton_cg_forcing*tol)
         ! a cold start solves its first (linear) step accurately before the nonlinear steps
         IF (cold_start .AND. iter == 1) cg_tol = MIN(cg_tol, 1.0E-3_dp*rms)
         ! J is exact and constant for a linear dielectric without ions, so one CG run goes straight to the final
         ! tolerance instead of 10x per Newton step (LINEAR_SINGLE_CG)
         IF (p%linear_single_cg .AND. p%psi_form .AND. ion_free_lin) cg_tol = MIN(cg_tol, tol)
         ! near the residual floor (RMS < 1e-5, no ions) the CG cannot converge anyway: cap it at 60 iterations
         ! so a floor exit costs 2 x 60 instead of 2 x 200 CG steps
         ncg = p%max_cg
         IF (rms < 1.0E-5_dp .AND. .NOT. p%ions) ncg = MIN(p%max_cg, MERGE(p%dev_i(1), 60, p%dev_i(1) > 0))
         ! chi_perp-based symmetric preconditioner scaling (DEV_L(8))
         IF (p%nldiel .AND. p%dev_l(8)) seps%array(:, :, :) = 1.0_dp/SQRT(1.0_dp + MAX(chi_perp%array(:, :, :), 0.0_dp))
         CALL timeset('nlpb_cg', h_cg)
         CALL run_cg(cg_tol, ncg)
         CALL timestop(h_cg)
         IF (p%cg_diag .AND. p%io_unit > 0) WRITE (UNIT=p%io_unit, FMT="(T3,A,I3,A,I5,A,2ES11.3)") &
            "NLPB| Newton CG diag iter ", iter, "  CG ", icg, "  entry RMS, cg_tol ", rms, cg_tol
         IF (p%selfcheck .AND. iter == 1) THEN
            CALL apply_j(p, dphi, eps, sa, cgd, work3, wr, wg1, wg2, wg3, pool, lb, ub)
            CALL pw_axpy(resid, cgd, -1.0_dp)
            chk(1) = SQRT(pw_integral_ab(cgd, cgd)/vol); chk(2) = SQRT(pw_integral_ab(dphi, dphi)/vol)
            dmax(1) = MAXVAL(ABS(dphi%array)); dmax(2) = MAXVAL(ABS(s_ion%array*dphi%array)); CALL para_env%max(dmax(1:2))
            IF (p%output_unit > 0) WRITE (UNIT=p%output_unit, FMT="(T3,A,4ES14.5)") &
               "NLPB| check |J dphi - R| RMS, |dphi| RMS, max|dphi|, max|S dphi|: ", chk(1:2), dmax(1:2)
         END IF
         ! ---- step limiter: the point-ion sinh is stiff, keep |beta z dphi| <= dx_max wherever S_ion is non-negligible
         xmax = 0.0_dp
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(dphi,lb,s_ion,ub) REDUCTION(MAX:xmax)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  IF (s_ion%array(i, j, k) > 1.0E-8_dp) xmax = MAX(xmax, ABS(dphi%array(i, j, k)))
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         CALL para_env%max(xmax)
         xmax = p%beta*p%zion*xmax
         ! ---- backtracking line search on A (concave: maximise) ----
         ! Cold start (phi == 0 on entry): the first Newton step is the bounded linear-PB solution, so the
         ! dx_max limiter would only slow it down. The Armijo backtracking below still protects.
         ! From iteration 2 the limiter applies.
         step = 1.0_dp
         ! limiter only for point ions; the lattice-gas response is bounded, Armijo suffices
         IF (p%ions .AND. .NOT. p%linear_ions .AND. xmax > p%dx_max .AND. &
             .NOT. (cold_start .AND. iter == 1)) step = p%dx_max/xmax
         rdphi = pw_integral_ab(resid, dphi)   ! > 0 for an ascent direction of the concave A
         ls_ok = .FALSE.
         DO ils = 1, 12
            CALL pw_copy(phi, phi_try)
            CALL pw_axpy(dphi, phi_try, step)
            CALL res_try(phi_try, step, a_new)
            IF (p%selfcheck .AND. iter == 1) THEN
               chk(1) = pw_integral_ab(rho_sol, phi_try); chk(2) = SQRT(pw_integral_ab(resid, resid)/vol)
               IF (p%output_unit > 0) WRITE (UNIT=p%output_unit, FMT="(T3,A,F8.4,5ES13.4)") &
                  "NLPB|   ls step, A, rho.phi, e_grad, e_ion, |R|: ", step, a_new, chk(1), chk(1) - a_new + e_diel + e_ion, &
                  e_ion, chk(2)
            END IF
            ! Armijo sufficient increase on the concave A (as in the nonlocal/nonlinear electrolyte model of
            ! https://doi.org/10.1063/5.0176308 : c = 0.1, step halved): A(phi + s dphi) >= A(phi) + c s R.dphi. When
            ! R.dphi <= 0 (dphi not an ascent direction) only the monotone test applies and the restart below follows.
            IF (a_new >= a_old + 0.1_dp*step*MAX(rdphi, 0.0_dp) - 1.0E-12_dp*ABS(a_old)) THEN
               ls_ok = .TRUE.
               EXIT
            END IF
            step = 0.5_dp*step
         END DO
         IF (.NOT. ls_ok .AND. .NOT. ls_retried) THEN
            ! one restart from phi = 0 (the next Newton step is then the bounded linear-PB solution of the current
            ! solute density, as on a cold start) before giving up; the second failure keeps the old phi below.
            ls_retried = .TRUE.
            CALL pw_zero(phi)
            CALL residual(p, phi, rho_sol, eps, s_ion, s_diel, resid, rho_ion, aion, adiel, a_old, q_ion, e_ion, e_diel, &
                          work3, wr, wg1, wg2, wg3, pool, lb, ub, lap_keep)
            lap_stale = .FALSE.
            res_fresh = .TRUE.
            rms = SQRT(pw_integral_ab(resid, resid)/vol)
            rms_prev = HUGE(1.0_dp); nstall = 0
            IF (p%io_unit > 0) WRITE (UNIT=p%io_unit, FMT="(T3,A,I4,A,ES12.4)") &
               "NLPB| line search failed at Newton ", iter, "  restarting from phi = 0, RMS ", rms
            CYCLE newton
         END IF
         IF (.NOT. ls_ok) THEN
            ! all 12 backtracks decreased the concave A: keep the OLD phi, do not accept a downhill step.
            ! rms is still the (> tol) value of the previous iteration, so converged stays .FALSE.
            ! diagnostics (collective min, every rank reaches this branch): the last step tried is 2*step
            fmin(1) = MINVAL(eps%array); fmin(2) = MINVAL(s_ion%array); CALL para_env%min(fmin)
            IF (p%io_unit > 0) WRITE (UNIT=p%io_unit, FMT="(T3,A,5ES12.4)") &
               "NLPB| diag R.dphi, A(last try)-A0, last step, min eps, min S_ion: ", rdphi, a_new - a_old, &
               2.0_dp*step, fmin(1), fmin(2)
            step = 0.0_dp
            IF (p%io_unit > 0) WRITE (UNIT=p%io_unit, FMT="(T3,A,I4,A,F6.3,A,ES12.4)") &
               "NLPB| line search failed at Newton ", iter, "  step ", step, "  RMS kept ", rms
            CPWARN("SCCS NLPB: line search failed, phi kept, Newton stopped")
            EXIT newton
         END IF
         CALL pw_copy(phi_try, phi)
         IF (lap_new) lap_keep%array(:, :, :) = lapn%array(:, :, :)
         IF (lap_new) lap_stale = .FALSE.
         res_fresh = .TRUE.
         a_old = a_new
         rms = SQRT(pw_integral_ab(resid, resid)/vol)
         IF (p%output_unit > 0) WRITE (UNIT=p%output_unit, FMT="(T3,A,I4,A,I5,A,F6.3,A,ES12.4,A,F20.10,A,ES9.2)") &
            "NLPB| Newton ", iter, "  CG ", icg, "  step ", step, "  RMS ", rms, "  A ", a_new, "  xmax[kT] ", xmax
         ! residual floor: without ions the FFT-derivative noise leaves RMS ~1e-7 that neither Newton nor the
         ! preconditioned CG can reduce, while A is already converged. Three Newton steps without halving the RMS
         ! at RMS < 1e-6 (DEV_R(1) overrides) -> accept.
         IF (rms > 0.5_dp*rms_prev) THEN
            nstall = nstall + 1
         ELSE
            nstall = 0
         END IF
         rms_prev = rms
         ! only the ion-free problem has a residual floor (with ions the solver reaches 1e-9); accepting RMS up to
         ! 1e-5 with electrolyte lets phi drift between SCF steps
         IF ((.NOT. p%ions) .AND. (.NOT. p%dev_l(1)) .AND. &
             ((nstall >= MERGE(p%dev_i(2), 3, p%dev_i(2) > 0)) .OR. (nstall >= 1 .AND. icg > ncg)) .AND. &
             rms < MERGE(p%dev_r(1), 1.0E-6_dp, p%dev_r(1) > 0.0_dp)) THEN
            IF (p%output_unit > 0) WRITE (UNIT=p%output_unit, FMT="(T3,A,ES12.4)") &
               "NLPB| residual floor reached, accepting RMS ", rms
            converged = .TRUE.
            EXIT newton
         END IF
      END DO newton
      END IF ! linear_path
      IF (rms <= tol) converged = .TRUE.
      a_tot = a_old
      rms_out = rms
      IF (warm_path) THEN
         IF (.NOT. ALLOCATED(hy)) THEN
            ALLOCATE (hy(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))
         ELSE IF (ANY(SHAPE(hy) /= (/ub(1) - lb(1) + 1, ub(2) - lb(2) + 1, ub(3) - lb(3) + 1/))) THEN
            DEALLOCATE (hy)
            ALLOCATE (hy(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))
         END IF
         hy(:, :, :) = phi%array(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3))
         hy_key(1:3) = phi%pw_grid%bounds_local(1, 1:3); hy_key(4:6) = phi%pw_grid%bounds_local(2, 1:3)
         hy_key(7:9) = phi%pw_grid%npts(1:3); hy_key(10) = phi%pw_grid%id_nr; hy_vol = phi%pw_grid%vol
         hy_ok = .NOT. p%skip
      END IF
      IF (use_guess) THEN
         CALL reuse_store(phi)
         DO l = 1, 3
            CALL pool%give_back_pw(jd(l)); CALL pool%give_back_pw(dv(l))
         END DO
      END IF
      ! final state quantities for the energy bookkeeping and the Kohn-Sham potential. When the last residual()
      ! call was made at the final phi (res_fresh: entry with no step, end of the linear CG, accepted line search) its
      ! outputs are exactly these quantities, so the evaluation is not repeated.
      CALL timeset('nlpb_final', h_fin)
      rf_final = res_fresh
      IF (.NOT. res_fresh) &
         CALL residual(p, phi, rho_sol, eps, s_ion, s_diel, resid, rho_ion, aion, adiel, a_new, q_ion, e_ion, e_diel, &
                       work3, wr, wg1, wg2, wg3, pool, lb, ub, lap_keep)
      IF (a_analytic .AND. .NOT. rf_final) a_tot = a_new   ! A was carried analytically by the history guess: take the evaluated one
      e_ion_phi = pw_integral_ab(rho_ion, phi)
      CALL pw_zero(adiel_unit)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,gfac,afac,ad,gmag,emag,pmag,lam,cperp,cpar) &
!$OMP             SHARED(adiel_unit,aion_unit,lb,p,phi,s_diel,ub,work3)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               CALL ion_response(p, phi%array(i, j, k), gfac, afac, ad)
               aion_unit%array(i, j, k) = ad
               IF (p%nldiel) THEN
                  ! work3 holds G = w_b*grad phi after residual(); n_mol lambda_diel per unit S_diel
                  gmag = SQRT(work3(1)%array(i, j, k)**2 + work3(2)%array(i, j, k)**2 + work3(3)%array(i, j, k)**2)
                  CALL nldiel_point(p%dl, gmag, s_diel%array(i, j, k), emag, pmag, lam, cperp, cpar)
                  adiel_unit%array(i, j, k) = p%dl%n_mol*lam
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      IF (p%nldiel) THEN
         IF (p%a_smooth > 0.0_dp) THEN
            CALL pool%give_back_pw(ukern_g); CALL pool%give_back_pw(ukern)
         END IF
         DO l = 1, 3
            CALL pool%give_back_pw(work3b(l)); CALL pool%give_back_pw(pvec(l)); CALL pool%give_back_pw(ghat(l))
         END DO
         CALL pool%give_back_pw(chi_par); CALL pool%give_back_pw(chi_perp)
         fields_ready = .FALSE.
      END IF
      DO l = 1, 3
         CALL pool%give_back_pw(work3(l))
      END DO
      IF (p%psi_form) THEN
         ! g_eff = [2 phi A(s phi) + lap(s phi^2) - q s phi^2]/(2 s) replaces |grad phi|^2 in the KS term
         wr%array(:, :, :) = psi_s%array(:, :, :)*phi%array(:, :, :)                 ! s phi
         cgd%array(:, :, :) = lap_keep%array(:, :, :)                               ! lap(s phi), saved by residual()
         phi_try%array(:, :, :) = wr%array(:, :, :)*phi%array(:, :, :)               ! s phi^2
         CALL pw_transfer(phi_try, wg1)
         DO i = 1, SIZE(wg1%array)
            wg1%array(i) = -wg1%pw_grid%gsq(i)*wg1%array(i)
         END DO
         CALL pw_transfer(wg1, cgp)                                                 ! lap(s phi^2)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(cgd,cgp,geff_unit,lb,phi,psi_q,seps,ub,wr)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  geff_unit%array(i, j, k) = seps%array(i, j, k)*(phi%array(i, j, k)*(-cgd%array(i, j, k) + &
                                             psi_q%array(i, j, k)*wr%array(i, j, k)) + 0.5_dp*cgp%array(i, j, k) - &
                                             0.5_dp*psi_q%array(i, j, k)*wr%array(i, j, k)*phi%array(i, j, k))
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         CALL pool%give_back_pw(psi_q); CALL pool%give_back_pw(psi_s)
      ELSE
         CALL pw_zero(geff_unit)
      END IF
      CALL timestop(h_fin)
      IF (p%psi_form) CALL pool%give_back_pw(lap_keep)
      IF (use_lapre) THEN
         CALL pool%give_back_pw(lapn); CALL pool%give_back_pw(lapg); CALL pool%give_back_pw(lapu)
      END IF
      CALL pool%give_back_pw(qfis)
      CALL pool%give_back_pw(seps)
      CALL pool%give_back_pw(cgd); CALL pool%give_back_pw(cgp); CALL pool%give_back_pw(cgz)
      CALL pool%give_back_pw(cgr); CALL pool%give_back_pw(phi_try); CALL pool%give_back_pw(dphi)
      CALL pool%give_back_pw(sa); CALL pool%give_back_pw(adiel); CALL pool%give_back_pw(aion)
      CALL pool%give_back_pw(rho_ion); CALL pool%give_back_pw(resid); CALL pool%give_back_pw(wr)
      CALL pool%give_back_pw(wg3); CALL pool%give_back_pw(wg2); CALL pool%give_back_pw(wg1)

      CALL timestop(handle)

   CONTAINS

! **************************************************************************************************
!> \brief residual() at phi_arg = phi + fac*dphi. After a psi-form CG the Laplacian of s*phi_arg is
!>        lap(s phi) - fac*lapu (lapu = sum alpha F^-1[G^2 p^] accumulated by the CG), so the closing residual needs
!>        no FFT; the result is left in lapn (the caller copies it to lap_keep on acceptance). When the history guess has
!>        changed phi, lap_keep is not lap(s phi): the Laplacian of the trial point is transformed (into lapn) and
!>        lap_stale stays set until a trial point is accepted.
! **************************************************************************************************
      SUBROUTINE res_try(phi_arg, fac, a_out)
         TYPE(pw_r3d_rs_type), INTENT(IN)                   :: phi_arg
         REAL(KIND=dp), INTENT(IN)                          :: fac
         REAL(KIND=dp), INTENT(OUT)                         :: a_out

         lap_new = .FALSE.
         IF (lapu_ok .AND. (.NOT. lap_stale)) THEN
            lapn%array(:, :, :) = lap_keep%array(:, :, :) - fac*lapu%array(:, :, :)
            CALL residual(p, phi_arg, rho_sol, eps, s_ion, s_diel, resid, rho_ion, aion, adiel, a_out, q_ion, e_ion, &
                          e_diel, work3, wr, wg1, wg2, wg3, pool, lb, ub, lap_in=lapn)
            lap_new = .TRUE.
         ELSE IF (lapu_ok) THEN
            CALL residual(p, phi_arg, rho_sol, eps, s_ion, s_diel, resid, rho_ion, aion, adiel, a_out, q_ion, e_ion, &
                          e_diel, work3, wr, wg1, wg2, wg3, pool, lb, ub, lap_out=lapn)
            lap_new = .TRUE.
         ELSE
            CALL residual(p, phi_arg, rho_sol, eps, s_ion, s_diel, resid, rho_ion, aion, adiel, a_out, q_ion, e_ion, &
                          e_diel, work3, wr, wg1, wg2, wg3, pool, lb, ub, lap_keep)
         END IF
      END SUBROUTINE res_try

! **************************************************************************************************
!> \brief drop the potential history when the grid / cell changed (host association: phi only used for its grid)
! **************************************************************************************************
      SUBROUTINE reuse_check_key(ph)
         TYPE(pw_r3d_rs_type), INTENT(IN)                   :: ph

         INTEGER, DIMENSION(10)                             :: key

         key(1:3) = ph%pw_grid%bounds_local(1, 1:3)
         key(4:6) = ph%pw_grid%bounds_local(2, 1:3)
         key(7:9) = ph%pw_grid%npts(1:3)
         key(10) = ph%pw_grid%id_nr
         IF (ALLOCATED(hx)) THEN
            IF (ALL(key == h_key) .AND. ABS(ph%pw_grid%vol - h_vol) <= 1.0E-10_dp*h_vol) RETURN
            DEALLOCATE (hx)
         END IF
         n_hx = 0
         h_key = key; h_vol = ph%pw_grid%vol
         ALLOCATE (hx(key(1):key(4), key(2):key(5), key(3):key(6), mh))
      END SUBROUTINE reuse_check_key

! **************************************************************************************************
!> \brief end of a solve, the potential joins the within-run history
! **************************************************************************************************
      SUBROUTINE reuse_store(ph)
         TYPE(pw_r3d_rs_type), INTENT(IN)                   :: ph

         IF (.NOT. ALLOCATED(hx)) RETURN
         IF (n_hx >= 1) hx(:, :, :, 2:mh) = hx(:, :, :, 1:mh - 1)
         hx(:, :, :, 1) = ph%array(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3))
         n_hx = MIN(n_hx + 1, mh)
      END SUBROUTINE reuse_store

! **************************************************************************************************
!> \brief HISTORY_GUESS: initial guess phi + sum_i c_i D_i for the linear ion-free problem, c = argmin
!>        |resid - sum_i c_i J D_i|_2 (exact for a constant J). D_i = differences of consecutive potentials of the
!>        history of this SCF run. The guess is taken when it lowers the entry RMS residual below 0.9 of the plain
!>        start; resid, phi and the reference value A are updated analytically (A is quadratic in phi).
! **************************************************************************************************
      SUBROUTINE reuse_guess()
         INTEGER                                            :: i, j, m
         REAL(KIND=dp)                                      :: r00, rms_pred
         REAL(KIND=dp), DIMENSION(3)                        :: bv, cv, rdv
         REAL(KIND=dp), DIMENSION(3, 3)                     :: gm, djm

         IF (.NOT. ALLOCATED(hx) .OR. cold_start .OR. p%skip) RETURN
         m = MIN(n_hx - 1, p%history_guess, 3)
         IF (m < 1) RETURN
         DO i = 1, m
            dv(i)%array(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)) = hx(:, :, :, i) - hx(:, :, :, i + 1)
            CALL apply_j(p, dv(i), eps, sa, jd(i), work3, wr, wg1, wg2, wg3, pool, lb, ub)
         END DO
         r00 = pw_integral_ab(resid, resid)
         DO i = 1, m
            bv(i) = pw_integral_ab(resid, jd(i))
            rdv(i) = pw_integral_ab(resid, dv(i))
            DO j = 1, m
               djm(i, j) = pw_integral_ab(dv(i), jd(j))
               gm(i, j) = pw_integral_ab(jd(i), jd(j))
            END DO
         END DO
         CALL small_ls(m, gm, bv, cv)
         ! |resid - sum_i c_i J D_i|^2
         rms_pred = SQRT(MAX(r00 - 2.0_dp*DOT_PRODUCT(cv(1:m), bv(1:m)) + &
                             DOT_PRODUCT(cv(1:m), MATMUL(gm(1:m, 1:m), cv(1:m))), 0.0_dp)/vol)
         IF (.NOT. (rms_pred < 0.9_dp*rms_plain)) RETURN
         DO i = 1, m
            CALL pw_axpy(dv(i), phi, cv(i))
            CALL pw_axpy(jd(i), resid, -cv(i))
         END DO
         lap_stale = .TRUE.   ! phi changed, lap_keep still belongs to the old phi
         rms = SQRT(pw_integral_ab(resid, resid)/vol)
         res_fresh = .FALSE.
         ! A(phi + d) = A + R.d - d.J d/2 with d = sum c_i D_i; the residual-derived fields are refreshed by the loop or the
         ! final residual call
         a_old = a_old + DOT_PRODUCT(cv(1:m), rdv(1:m)) - 0.5_dp*DOT_PRODUCT(cv(1:m), MATMUL(djm(1:m, 1:m), cv(1:m)))
         a_analytic = .TRUE.
      END SUBROUTINE reuse_guess

! **************************************************************************************************
!> \brief tiny SPD solve (Gaussian elimination, partial pivoting, relative ridge 1e-10) of the leading n x n block
! **************************************************************************************************
      SUBROUTINE small_ls(n, gm, bv, cv)
         INTEGER, INTENT(IN)                                :: n
         REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: gm
         REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: bv
         REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: cv

         INTEGER                                            :: i, j, k, piv
         REAL(KIND=dp)                                      :: f, tmp
         REAL(KIND=dp), DIMENSION(3, 4)                     :: aug

         cv(:) = 0.0_dp
         aug(:, :) = 0.0_dp
         DO i = 1, n
            DO j = 1, n
               aug(i, j) = gm(i, j)
            END DO
            aug(i, i) = aug(i, i)*(1.0_dp + 1.0E-10_dp) + TINY(1.0_dp)
            aug(i, 4) = bv(i)
         END DO
         DO k = 1, n
            piv = k
            DO i = k + 1, n
               IF (ABS(aug(i, k)) > ABS(aug(piv, k))) piv = i
            END DO
            IF (ABS(aug(piv, k)) <= 0.0_dp) RETURN
            IF (piv /= k) THEN
               DO j = 1, 4
                  tmp = aug(k, j); aug(k, j) = aug(piv, j); aug(piv, j) = tmp
               END DO
            END IF
            DO i = k + 1, n
               f = aug(i, k)/aug(k, k)
               DO j = k, 4
                  aug(i, j) = aug(i, j) - f*aug(k, j)
               END DO
            END DO
         END DO
         DO i = n, 1, -1
            tmp = aug(i, 4)
            DO j = i + 1, n
               tmp = tmp - aug(i, j)*cv(j)
            END DO
            cv(i) = tmp/aug(i, i)
         END DO
      END SUBROUTINE small_ls

! **************************************************************************************************
!> \brief preconditioned CG on J dphi = resid (host association: shared Newton / linear path)
! **************************************************************************************************
      SUBROUTINE run_cg(cg_tol_in, ncg_max)
         REAL(KIND=dp), INTENT(IN)                          :: cg_tol_in
         INTEGER, INTENT(IN)                                :: ncg_max

         lapu_ok = .FALSE.
         IF (use_fis) THEN
            ! two-level preconditioner for the ion-free psi-form problem (no Debye screening: the low-G modes of
            ! the layered dielectric converge slowly). DEV_I(5): 0 = on when ion-free, > 0 = on with M_z = DEV_I(5), < 0 = off.
            ! DEV_R(6) = bit mask: 1 FFT coarse transforms, 2 no lap(s phi) carry-over, 4 carry-over check,
            ! 8 CG only 10x per Newton step (ion-free linear dielectric)
            use_cc = p%psi_form .AND. ((p%dev_i(5) > 0) .OR. (p%dev_i(5) == 0 .AND. ((.NOT. p%ions) .OR. cc_ions_here)))
            IF (use_cc) CALL cc_setup()
            IF (use_cc .AND. cc_ok) THEN
               CALL run_cg_cc(cg_tol_in, ncg_max)
            ELSE
               CALL run_cg_fis(cg_tol_in, ncg_max)
            END IF
            RETURN
         END IF
         CALL pw_zero(dphi)
         CALL pw_copy(resid, cgr)
         CALL precond(cgr, cgz, seps, phi_try, use_seps, pool, greenfn)
         CALL pw_copy(cgz, cgp)
         rz = pw_integral_ab(cgr, cgz)
         icg = 0
         DO icg = 1, ncg_max
            CALL apply_j(p, cgp, eps, sa, cgd, work3, wr, wg1, wg2, wg3, pool, lb, ub)
            pap = pw_integral_ab(cgp, cgd)
            IF (p%selfcheck .AND. iter == 1) THEN
               chk(1) = SQRT(pw_integral_ab(cgr, cgr)/vol)
               IF (p%output_unit > 0) WRITE (UNIT=p%output_unit, FMT="(T3,A,I4,4ES13.4)") &
                  "NLPB|   cg ", icg, rz, pap, rz/MAX(pap, 1.0E-300_dp), chk(1)
            END IF
            IF (pap <= 0.0_dp) THEN
               ! J is SPD for eps >= 1, S_ion >= 0; a non-positive curvature means an indefinite operator
               IF (p%io_unit > 0) WRITE (UNIT=p%io_unit, FMT="(T3,A,I5,A,ES12.4)") &
                  "NLPB| CG breakdown at CG ", icg, "  p.Jp ", pap
               EXIT
            END IF
            alpha = rz/pap
            CALL pw_axpy(cgp, dphi, alpha)
            CALL pw_axpy(cgd, cgr, -alpha)
            cg_rms_it = SQRT(pw_integral_ab(cgr, cgr)/vol)
            CALL cg_track(icg, 0.5_dp*alpha*rz, cg_rms_it)
            IF (cg_rms_it <= cg_tol_in) EXIT
            IF (cg_estop) EXIT
            CALL precond(cgr, cgz, seps, phi_try, use_seps, pool, greenfn)
            rz_new = pw_integral_ab(cgr, cgz)
            beta_cg = rz_new/rz
            rz = rz_new
            CALL pw_scale(cgp, beta_cg)
            CALL pw_axpy(cgz, cgp)
         END DO
         ncg_tot = ncg_tot + MIN(icg, ncg_max)
      END SUBROUTINE run_cg

! **************************************************************************************************
!> \brief index of Miller triple (m1,m2,m3) in a box |m1|,|m2| <= a, |m3| <= b (1-based, 0 outside)
! **************************************************************************************************
      PURE FUNCTION box_index(m1, m2, m3, a, b) RESULT(ix)
         INTEGER, INTENT(IN)                                :: m1, m2, m3, a, b
         INTEGER                                            :: ix

         IF (ABS(m1) > a .OR. ABS(m2) > a .OR. ABS(m3) > b) THEN
            ix = 0
         ELSE
            ix = (m1 + a) + (2*a + 1)*((m2 + a) + (2*a + 1)*(m3 + b)) + 1
         END IF
      END FUNCTION box_index

! **************************************************************************************************
!> \brief coarse operator A0_ij = V/(4pi) (G_i.G_j) eps^(G_i - G_j) over the low plane waves
!>        0 < |m_xy| <= cc_mxy, |m_z| <= cc_mz (G = 0 left out: the constant potential is the null space), Cholesky factor.
! **************************************************************************************************
      SUBROUTINE cc_setup()
         INTEGER                                            :: i, ig, ix, j, k, m1, m2, m3, mode, n, nE, ntr
         LOGICAL                                            :: bad, ion
         REAL(KIND=dp)                                      :: dg, errmax, refmax, s
         COMPLEX(KIND=dp)                                   :: acc
         COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: e3, sa3, ref
         COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: a0

         CALL timeset('cc_setup', h_ccs)
         IF (ALLOCATED(cc_m)) DEALLOCATE (cc_m, cc_map, cc_cart, cc_l, cc_rhs)
         IF (ALLOCATED(cc_ez)) DEALLOCATE (cc_ez, cc_wt, cc_tt, cc_ss)
         cc_ok = .FALSE.
         cc_dir = .FALSE.
         ion = p%ions   ! with ions the constant mode is not a null vector of the operator, the coarse space keeps it
         cc_mxy = MERGE(0, MERGE(p%dev_i(4), 1, p%dev_i(4) > 0), p%dev_i(4) < 0)   ! DEV_I(4) < 0: 1-D (m_xy = 0) space only
         cc_mz = MERGE(p%dev_i(5), 10, p%dev_i(5) > 0)
         cc_n = (2*cc_mxy + 1)**2*(2*cc_mz + 1) - MERGE(0, 1, ion)
         ALLOCATE (cc_m(3, cc_n), cc_map((2*cc_mxy + 1)**2*(2*cc_mz + 1)), cc_cart(3, cc_n), cc_l(cc_n, cc_n), cc_rhs(cc_n))
         cc_map = 0
         n = 0
         DO m3 = -cc_mz, cc_mz
            DO m2 = -cc_mxy, cc_mxy
               DO m1 = -cc_mxy, cc_mxy
                  IF (m1 == 0 .AND. m2 == 0 .AND. m3 == 0 .AND. (.NOT. ion)) CYCLE
                  n = n + 1
                  cc_m(:, n) = (/m1, m2, m3/)
                  cc_map(box_index(m1, m2, m3, cc_mxy, cc_mz)) = n
                  DO k = 1, 3
                     cc_cart(k, n) = twopi*(REAL(m1, dp)/REAL(wg1%pw_grid%npts(1), dp)*wg1%pw_grid%dh_inv(1, k) + &
                                            REAL(m2, dp)/REAL(wg1%pw_grid%npts(2), dp)*wg1%pw_grid%dh_inv(2, k) + &
                                            REAL(m3, dp)/REAL(wg1%pw_grid%npts(3), dp)*wg1%pw_grid%dh_inv(3, k))
                  END DO
               END DO
            END DO
         END DO
         ! the reciprocal vectors must agree with the grid's own G list
         bad = .FALSE.
         DO ig = 1, wg1%pw_grid%ngpts_cut_local
            ix = box_index(wg1%pw_grid%g_hat(1, ig), wg1%pw_grid%g_hat(2, ig), wg1%pw_grid%g_hat(3, ig), cc_mxy, cc_mz)
            IF (ix > 0) THEN
               IF (cc_map(ix) > 0) THEN
                  dg = MAXVAL(ABS(cc_cart(:, cc_map(ix)) - wg1%pw_grid%g(:, ig)))
                  IF (dg > 1.0E-8_dp) THEN
                     bad = .TRUE.
                     IF (p%io_unit > 0) WRITE (UNIT=p%io_unit, FMT="(T3,A,3I4,A,3ES12.4,A,3ES12.4)") &
                        "NLPB| coarse G check m ", wg1%pw_grid%g_hat(:, ig), " grid g ", wg1%pw_grid%g(:, ig), &
                        " mine ", cc_cart(:, cc_map(ix))
                  END IF
               END IF
            END IF
         END DO
         s = MERGE(1.0_dp, 0.0_dp, bad)
         CALL para_env%max(s)
         IF (s > 0.5_dp) THEN
            IF (p%io_unit > 0) WRITE (UNIT=p%io_unit, FMT="(T3,A)") "NLPB| coarse space: G vectors disagree, switched off"
            CALL timestop(h_ccs)
            RETURN
         END IF
         ! eps = s^2 in G space, stored for the differences |m_i - m_j| <= 2 box
         wr%array(:, :, :) = 1.0_dp/(seps%array(:, :, :)*seps%array(:, :, :))
         CALL pw_transfer(wr, wg3)
         nE = (4*cc_mxy + 1)**2*(4*cc_mz + 1)
         ALLOCATE (e3(nE), sa3(nE))
         e3 = (0.0_dp, 0.0_dp)
         sa3 = (0.0_dp, 0.0_dp)
         DO ig = 1, wg1%pw_grid%ngpts_cut_local
            m1 = wg1%pw_grid%g_hat(1, ig); m2 = wg1%pw_grid%g_hat(2, ig); m3 = wg1%pw_grid%g_hat(3, ig)
            ix = box_index(m1, m2, m3, 2*cc_mxy, 2*cc_mz)
            IF (ix > 0) THEN
               e3(ix) = wg3%array(ig)
               IF (wg1%pw_grid%grid_span == HALFSPACE .AND. (m1 /= 0 .OR. m2 /= 0 .OR. m3 /= 0)) &
                  e3(box_index(-m1, -m2, -m3, 2*cc_mxy, 2*cc_mz)) = CONJG(wg3%array(ig))
            END IF
         END DO
         CALL para_env%sum(e3)
         IF (ion) THEN
            ! ion term of the Jacobian, <e_i, sa e_j> = V sa^(G_i - G_j) (K carries sa s^2/4pi on u = s dphi)
            wr%array(:, :, :) = sa%array(:, :, :)
            CALL pw_transfer(wr, wg3)
            DO ig = 1, wg1%pw_grid%ngpts_cut_local
               m1 = wg1%pw_grid%g_hat(1, ig); m2 = wg1%pw_grid%g_hat(2, ig); m3 = wg1%pw_grid%g_hat(3, ig)
               ix = box_index(m1, m2, m3, 2*cc_mxy, 2*cc_mz)
               IF (ix > 0) THEN
                  sa3(ix) = wg3%array(ig)
                  IF (wg1%pw_grid%grid_span == HALFSPACE .AND. (m1 /= 0 .OR. m2 /= 0 .OR. m3 /= 0)) &
                     sa3(box_index(-m1, -m2, -m3, 2*cc_mxy, 2*cc_mz)) = CONJG(wg3%array(ig))
               END IF
            END DO
            CALL para_env%sum(sa3)
         END IF
         ! FFT-free coarse transforms. The first solve finds the (sign, origin) convention under which the direct
         ! phase sums reproduce the FFT coefficients of eps (e3 above, no extra FFT); the result is kept for the run.
         IF (p%coarse_direct_sum .AND. cc_dir_mode > 0) THEN
            CALL cc_tables(cc_dir_mode)
            cc_dir = .TRUE.
         ELSE IF (p%coarse_direct_sum .AND. cc_dir_mode == 0) THEN
            ntr = 4
            ALLOCATE (ref(cc_n))
            DO n = 1, cc_n
               ref(n) = vol*e3(box_index(cc_m(1, n), cc_m(2, n), cc_m(3, n), 2*cc_mxy, 2*cc_mz))
            END DO
            refmax = MAXVAL(ABS(ref))
            DO k = 1, ntr
               mode = k
               CALL cc_tables(mode)
               wr%array(:, :, :) = 1.0_dp/(seps%array(:, :, :)*seps%array(:, :, :))
               CALL cc_proj(wr, cc_rhs)
               errmax = MAXVAL(ABS(cc_rhs - ref))
               IF (p%cg_diag .AND. p%io_unit > 0) WRITE (UNIT=p%io_unit, FMT="(T3,A,I2,A,2ES12.4)") &
                  "NLPB| coarse direct transform, convention ", mode, "  max err, max ref ", errmax, refmax
               IF (errmax <= 1.0E-9_dp*refmax) THEN
                  cc_dir_mode = mode
                  cc_dir = .TRUE.
                  EXIT
               END IF
               DEALLOCATE (cc_ez, cc_wt, cc_tt, cc_ss)
            END DO
            IF (.NOT. cc_dir) THEN
               cc_dir_mode = -1
               IF (p%io_unit > 0) WRITE (UNIT=p%io_unit, FMT="(T3,A)") &
                  "NLPB| coarse direct transform does not match the FFT, FFT path kept"
            END IF
            DEALLOCATE (ref)
         END IF
         ALLOCATE (a0(cc_n, cc_n))
         DO j = 1, cc_n
            DO i = 1, cc_n
               ix = box_index(cc_m(1, i) - cc_m(1, j), cc_m(2, i) - cc_m(2, j), cc_m(3, i) - cc_m(3, j), 2*cc_mxy, 2*cc_mz)
               a0(i, j) = vol/fourpi*((cc_cart(1, i)*cc_cart(1, j) + cc_cart(2, i)*cc_cart(2, j) + &
                                       cc_cart(3, i)*cc_cart(3, j))*e3(ix) + sa3(ix))
            END DO
         END DO
         ! Cholesky A0 = L L^H
         cc_l = (0.0_dp, 0.0_dp)
         DO j = 1, cc_n
            s = REAL(a0(j, j), dp)
            DO k = 1, j - 1
               s = s - REAL(cc_l(j, k)*CONJG(cc_l(j, k)), dp)
            END DO
            IF (s <= 1.0E-12_dp*MAX(REAL(a0(j, j), dp), 1.0E-300_dp)) THEN
               IF (p%io_unit > 0) WRITE (UNIT=p%io_unit, FMT="(T3,A,I6)") "NLPB| coarse operator not positive at ", j
               DEALLOCATE (e3, sa3, a0)
               cc_dir = .FALSE.
               CALL timestop(h_ccs)
               RETURN
            END IF
            cc_l(j, j) = CMPLX(SQRT(s), 0.0_dp, dp)
            DO i = j + 1, cc_n
               acc = a0(i, j)
               DO k = 1, j - 1
                  acc = acc - cc_l(i, k)*CONJG(cc_l(j, k))
               END DO
               cc_l(i, j) = acc/cc_l(j, j)
            END DO
         END DO
         DEALLOCATE (e3, sa3, a0)
         cc_ok = .TRUE.
         IF (p%cg_diag .AND. p%io_unit > 0) WRITE (UNIT=p%io_unit, FMT="(T3,A,I6,A,2I4,A,L2)") &
            "NLPB| coarse space dim ", cc_n, "  Mxy, Mz ", cc_mxy, cc_mz, "  direct ", cc_dir
         CALL timestop(h_ccs)
      END SUBROUTINE cc_setup

! **************************************************************************************************
!> \brief phase tables of the direct coarse transforms for the local grid block. mode 1..4 = (sign, origin):
!>        forward phase exp(-i sgn 2pi m n/N), n = global index measured from bounds(1, :) (origin 1) or as given (2).
!>        The local block (i, j) is treated as one contiguous column index ic = (i - lb1 + 1) + nxl (j - lb2); only
!>        m3 >= 0 is tabulated (a real field has T(-m3) = conj T(m3)).
! **************************************************************************************************
      SUBROUTINE cc_tables(mode)
         INTEGER, INTENT(IN)                                :: mode

         INTEGER                                            :: i, ic, j, k, m1, m2, m3, ixy, nxl, nyl, o1, o2, o3, sgn
         REAL(KIND=dp)                                      :: ang
         COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: ex, ey

         sgn = MERGE(1, -1, mode <= 2)
         IF (MOD(mode, 2) == 1) THEN
            o1 = wg1%pw_grid%bounds(1, 1); o2 = wg1%pw_grid%bounds(1, 2); o3 = wg1%pw_grid%bounds(1, 3)
         ELSE
            o1 = 0; o2 = 0; o3 = 0
         END IF
         nxl = ub(1) - lb(1) + 1; nyl = ub(2) - lb(2) + 1
         ALLOCATE (ex(-cc_mxy:cc_mxy, lb(1):ub(1)), ey(-cc_mxy:cc_mxy, lb(2):ub(2)), cc_ez(0:cc_mz, lb(3):ub(3)))
         ALLOCATE (cc_wt(nxl*nyl, (2*cc_mxy + 1)**2), cc_tt(nxl*nyl, 0:cc_mz), cc_ss(nxl*nyl, 0:cc_mz))
         DO i = lb(1), ub(1)
            DO m1 = -cc_mxy, cc_mxy
               ang = twopi*REAL(MODULO(m1*(i - o1), wg1%pw_grid%npts(1)), dp)/REAL(wg1%pw_grid%npts(1), dp)
               ex(m1, i) = CMPLX(COS(ang), -sgn*SIN(ang), dp)
            END DO
         END DO
         DO j = lb(2), ub(2)
            DO m2 = -cc_mxy, cc_mxy
               ang = twopi*REAL(MODULO(m2*(j - o2), wg1%pw_grid%npts(2)), dp)/REAL(wg1%pw_grid%npts(2), dp)
               ey(m2, j) = CMPLX(COS(ang), -sgn*SIN(ang), dp)
            END DO
         END DO
         DO k = lb(3), ub(3)
            DO m3 = 0, cc_mz
               ang = twopi*REAL(MODULO(m3*(k - o3), wg1%pw_grid%npts(3)), dp)/REAL(wg1%pw_grid%npts(3), dp)
               cc_ez(m3, k) = CMPLX(COS(ang), -sgn*SIN(ang), dp)
            END DO
         END DO
         DO m2 = -cc_mxy, cc_mxy
            DO m1 = -cc_mxy, cc_mxy
               ixy = (m1 + cc_mxy) + (2*cc_mxy + 1)*(m2 + cc_mxy) + 1
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     ic = (i - lb(1) + 1) + nxl*(j - lb(2))
                     cc_wt(ic, ixy) = ex(m1, i)*ey(m2, j)
                  END DO
               END DO
            END DO
         END DO
         DEALLOCATE (ex, ey)
      END SUBROUTINE cc_tables

! **************************************************************************************************
!> \brief out(n) = V * c(G_n) of the real grid function f, c = FFT coefficient (1/N) sum_r f e^(-iGr),
!>        by separable phase sums over the local block (11 + 9 phases per point instead of a distributed 3D FFT)
! **************************************************************************************************
      SUBROUTINE cc_proj(f, out)
         TYPE(pw_r3d_rs_type), INTENT(IN)                   :: f
         COMPLEX(KIND=dp), DIMENSION(:), INTENT(INOUT)      :: out

         INTEGER                                            :: h_pr, ic, ixy, kk, m3, n, ncol, nk
         COMPLEX(KIND=dp)                                   :: acc, c
         REAL(KIND=dp)                                      :: fac
         REAL(KIND=dp), DIMENSION(:, :), POINTER            :: f2

         CALL timeset('cc_proj', h_pr)
         ncol = (ub(1) - lb(1) + 1)*(ub(2) - lb(2) + 1); nk = ub(3) - lb(3) + 1
         CALL C_F_POINTER(C_LOC(f%array(lb(1), lb(2), lb(3))), f2, (/ncol, nk/))
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(m3,kk,ic,c) SHARED(cc_ez,cc_mz,cc_tt,f2,lb,ncol,nk)
         DO m3 = 0, cc_mz
            cc_tt(:, m3) = (0.0_dp, 0.0_dp)
            DO kk = 1, nk
               c = cc_ez(m3, lb(3) + kk - 1)
               DO ic = 1, ncol
                  cc_tt(ic, m3) = cc_tt(ic, m3) + c*f2(ic, kk)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         fac = vol/(REAL(wg1%pw_grid%npts(1), dp)*REAL(wg1%pw_grid%npts(2), dp)*REAL(wg1%pw_grid%npts(3), dp))
         DO n = 1, cc_n
            ixy = (cc_m(1, n) + cc_mxy) + (2*cc_mxy + 1)*(cc_m(2, n) + cc_mxy) + 1
            m3 = cc_m(3, n)
            acc = (0.0_dp, 0.0_dp)
            IF (m3 >= 0) THEN
               DO ic = 1, ncol
                  acc = acc + cc_wt(ic, ixy)*cc_tt(ic, m3)
               END DO
            ELSE
               DO ic = 1, ncol
                  acc = acc + cc_wt(ic, ixy)*CONJG(cc_tt(ic, -m3))
               END DO
            END IF
            out(n) = fac*acc
         END DO
         CALL para_env%sum(out)
         CALL timestop(h_pr)
      END SUBROUTINE cc_proj

! **************************************************************************************************
!> \brief f(r) = Re sum_n coef(n) e^(+iG_n r) on the local block (the phases of cc_proj, conjugated)
! **************************************************************************************************
      SUBROUTINE cc_expand(coef, f)
         COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: coef
         TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: f

         INTEGER                                            :: h_ex, ic, ixy, kk, m3, n, ncol, nk
         COMPLEX(KIND=dp)                                   :: e, w
         COMPLEX(KIND=dp), DIMENSION((2*cc_mxy + 1)**2, -cc_mz:cc_mz) :: cm
         REAL(KIND=dp), DIMENSION(:, :), POINTER            :: f2

         CALL timeset('cc_expand', h_ex)
         ncol = (ub(1) - lb(1) + 1)*(ub(2) - lb(2) + 1); nk = ub(3) - lb(3) + 1
         cm = (0.0_dp, 0.0_dp)
         DO n = 1, cc_n
            ixy = (cc_m(1, n) + cc_mxy) + (2*cc_mxy + 1)*(cc_m(2, n) + cc_mxy) + 1
            cm(ixy, cc_m(3, n)) = coef(n)
         END DO
         ! S'(m3) = S(m3) + conj S(-m3), S(m3) = sum_xy conj(w_xy) c(xy, m3), so that f = Re sum_{m3>=0} S'(m3) conj(e(m3))
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(m3,ic,ixy,w) SHARED(cc_mxy,cc_mz,cc_ss,cc_wt,cm,ncol)
         DO m3 = 0, cc_mz
            cc_ss(:, m3) = (0.0_dp, 0.0_dp)
            DO ixy = 1, (2*cc_mxy + 1)**2
               DO ic = 1, ncol
                  w = cc_wt(ic, ixy)
                  cc_ss(ic, m3) = cc_ss(ic, m3) + CONJG(w)*cm(ixy, m3)
                  IF (m3 > 0) cc_ss(ic, m3) = cc_ss(ic, m3) + w*CONJG(cm(ixy, -m3))
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         CALL C_F_POINTER(C_LOC(f%array(lb(1), lb(2), lb(3))), f2, (/ncol, nk/))
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(kk,m3,ic,e) SHARED(cc_ez,cc_mz,cc_ss,f2,lb,ncol,nk)
         DO kk = 1, nk
            f2(:, kk) = 0.0_dp
            DO m3 = 0, cc_mz
               e = cc_ez(m3, lb(3) + kk - 1)
               DO ic = 1, ncol
                  f2(ic, kk) = f2(ic, kk) + REAL(cc_ss(ic, m3), dp)*REAL(e, dp) + AIMAG(cc_ss(ic, m3))*AIMAG(e)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         CALL timestop(h_ex)
      END SUBROUTINE cc_expand

! **************************************************************************************************
!> \brief z = G r' + s * dphi_c, dphi_c = coarse solve of the low-G part of s r' (transformed system K u = r')
! **************************************************************************************************
      SUBROUTINE cc_apply(rr, zz)
         TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rr
         TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: zz

         CALL pw_transfer(rr, wg2)
         wg2%array(:) = wg2%array(:)*greenfn%array(:)
         CALL pw_transfer(wg2, zz)
         CALL cc_coarse(rr, zz)
      END SUBROUTINE cc_apply

! **************************************************************************************************
!> \brief zz = zz + s * dphi_c, dphi_c = coarse solve of the low-G part of s r' (K u = r' with u = s dphi)
! **************************************************************************************************
      SUBROUTINE cc_coarse(rr, zz)
         TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rr
         TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: zz

         INTEGER                                            :: i, ig, ix, k, m1, m2, m3
         COMPLEX(KIND=dp)                                   :: acc, v

         CALL timeset('cc_coarse', h_ccc)
         phi_try%array(:, :, :) = rr%array(:, :, :)/seps%array(:, :, :)
         IF (cc_dir) THEN
            CALL cc_proj(phi_try, cc_rhs)
         ELSE
         CALL pw_transfer(phi_try, wg3)
         cc_rhs = (0.0_dp, 0.0_dp)
         DO ig = 1, wg3%pw_grid%ngpts_cut_local
            m1 = wg3%pw_grid%g_hat(1, ig); m2 = wg3%pw_grid%g_hat(2, ig); m3 = wg3%pw_grid%g_hat(3, ig)
            ix = box_index(m1, m2, m3, cc_mxy, cc_mz)
            IF (ix > 0) THEN
               IF (cc_map(ix) > 0) THEN
                  v = vol*wg3%array(ig)
                  cc_rhs(cc_map(ix)) = v
                  IF (wg3%pw_grid%grid_span == HALFSPACE) THEN
                     IF (m1 /= 0 .OR. m2 /= 0 .OR. m3 /= 0) &
                        cc_rhs(cc_map(box_index(-m1, -m2, -m3, cc_mxy, cc_mz))) = CONJG(v)
                  END IF
               END IF
            END IF
         END DO
         CALL para_env%sum(cc_rhs)
         END IF
         DO i = 1, cc_n
            acc = cc_rhs(i)
            DO k = 1, i - 1
               acc = acc - cc_l(i, k)*cc_rhs(k)
            END DO
            cc_rhs(i) = acc/cc_l(i, i)
         END DO
         DO i = cc_n, 1, -1
            acc = cc_rhs(i)
            DO k = i + 1, cc_n
               acc = acc - CONJG(cc_l(k, i))*cc_rhs(k)
            END DO
            cc_rhs(i) = acc/CONJG(cc_l(i, i))
         END DO
         IF (cc_dir) THEN
            CALL cc_expand(cc_rhs, phi_try)
         ELSE
         CALL pw_zero(wg3)
         DO ig = 1, wg3%pw_grid%ngpts_cut_local
            ix = box_index(wg3%pw_grid%g_hat(1, ig), wg3%pw_grid%g_hat(2, ig), wg3%pw_grid%g_hat(3, ig), cc_mxy, cc_mz)
            IF (ix > 0) THEN
               IF (cc_map(ix) > 0) wg3%array(ig) = cc_rhs(cc_map(ix))
            END IF
         END DO
         CALL pw_transfer(wg3, phi_try)
         END IF
         zz%array(:, :, :) = zz%array(:, :, :) + phi_try%array(:, :, :)/seps%array(:, :, :)
         CALL timestop(h_ccc)
      END SUBROUTINE cc_coarse

! **************************************************************************************************
!> \brief cgd = K cgp, K = (G^2 + q + sa/eps)/4pi on the transformed unknown (2 FFTs)
! **************************************************************************************************
      SUBROUTINE apply_k_cc()
         INTEGER                                            :: ig

         CALL pw_transfer(cgp, wg1)
         DO ig = 1, SIZE(wg1%array)
            wg1%array(ig) = wg1%pw_grid%gsq(ig)*wg1%array(ig)
         END DO
         CALL pw_transfer(wg1, cgd)
         IF (use_lapre) lapg%array(:, :, :) = cgd%array(:, :, :)   ! G^2 part, for the carried lap(s phi)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(cgd,cgp,lb,qfis,sa,seps,ub)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  cgd%array(i, j, k) = (cgd%array(i, j, k) + (qfis%array(i, j, k) + &
                                        sa%array(i, j, k)*seps%array(i, j, k)**2)*cgp%array(i, j, k))/fourpi
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
      END SUBROUTINE apply_k_cc

! **************************************************************************************************
!> \brief run_cg_fis with the two-level preconditioner (Green function + coarse correction)
! **************************************************************************************************
      SUBROUTINE run_cg_cc(cg_tol_in, ncg_max)
         REAL(KIND=dp), INTENT(IN)                          :: cg_tol_in
         INTEGER, INTENT(IN)                                :: ncg_max

         LOGICAL                                            :: once
         REAL(KIND=dp)                                      :: rres

         ! DEV_I(6) = 1: coarse correction only once, as the initial guess; afterwards Green-function PCG
         once = (p%dev_i(6) == 1)
         CALL pw_zero(dphi)
         IF (use_lapre) CALL pw_zero(lapu)
         cgr%array(:, :, :) = resid%array(:, :, :)*seps%array(:, :, :)
         IF (once) THEN
            CALL pw_zero(cgz)
            CALL cc_coarse(cgr, cgz)
            CALL pw_copy(cgz, dphi)
            CALL pw_copy(cgz, cgp)
            CALL apply_k_cc()
            IF (use_lapre) CALL pw_axpy(lapg, lapu, 1.0_dp)
            CALL pw_axpy(cgd, cgr, -1.0_dp)
            CALL pw_transfer(cgr, wg2)
            wg2%array(:) = wg2%array(:)*greenfn%array(:)
            CALL pw_transfer(wg2, cgz)
         ELSE
            CALL cc_apply(cgr, cgz)
         END IF
         CALL pw_copy(cgz, cgp)
         rz = pw_integral_ab(cgr, cgz)
         icg = 0
         DO icg = 1, ncg_max
            CALL apply_k_cc()
            pap = pw_integral_ab(cgp, cgd)
            IF (pap <= 0.0_dp) THEN
               IF (p%io_unit > 0) WRITE (UNIT=p%io_unit, FMT="(T3,A,I5,A,ES12.4)") &
                  "NLPB| CG breakdown at CG ", icg, "  p.Jp ", pap
               EXIT
            END IF
            alpha = rz/pap
            CALL pw_axpy(cgp, dphi, alpha)
            IF (use_lapre) CALL pw_axpy(lapg, lapu, alpha)
            CALL pw_axpy(cgd, cgr, -alpha)
            wr%array(:, :, :) = cgr%array(:, :, :)/seps%array(:, :, :)
            rres = SQRT(pw_integral_ab(wr, wr)/vol)
            CALL cg_track(icg, 0.5_dp*alpha*rz, rres)
            IF (rres <= cg_tol_in) EXIT
            IF (cg_estop) EXIT
            IF (once) THEN
               CALL pw_transfer(cgr, wg2)
               wg2%array(:) = wg2%array(:)*greenfn%array(:)
               CALL pw_transfer(wg2, cgz)
            ELSE
               CALL cc_apply(cgr, cgz)
            END IF
            rz_new = pw_integral_ab(cgr, cgz)
            beta_cg = rz_new/rz
            rz = rz_new
            CALL pw_scale(cgp, beta_cg)
            CALL pw_axpy(cgz, cgp)
         END DO
         dphi%array(:, :, :) = dphi%array(:, :, :)*seps%array(:, :, :)
         lapu_ok = use_lapre
         ncg_tot = ncg_tot + MIN(icg, ncg_max)
         IF (ALLOCATED(cc_m)) DEALLOCATE (cc_m, cc_map, cc_cart, cc_l, cc_rhs)
         IF (cc_dir) DEALLOCATE (cc_ez, cc_wt, cc_tt, cc_ss)
      END SUBROUTINE run_cg_cc

! **************************************************************************************************
!> \brief preconditioned CG on the transformed system K u = R/s, then dphi = u/s (host association)
! **************************************************************************************************
      SUBROUTINE run_cg_fis(cg_tol_in, ncg_max)
         REAL(KIND=dp), INTENT(IN)                          :: cg_tol_in
         INTEGER, INTENT(IN)                                :: ncg_max

         INTEGER                                            :: ig
         REAL(KIND=dp)                                      :: rres
         TYPE(pw_c1d_gs_type)                               :: pg

         ! the search direction is kept in G space (pg = F[cgp]) next to its real-space copy. The Green function
         ! is applied in G space, so F[z] is available without a transform and z = F^-1[F[z]] needs one inverse FFT; the
         ! operator then costs 1 inverse FFT (G^2 part), the preconditioner 1 forward + 1 inverse: 3 FFTs per CG
         ! iteration. F[F^-1[g]] = g for g inside the G sphere, so this equals the real-space iteration up to round-off.
         CALL pool%create_pw(pg)
         CALL pw_zero(dphi)                                          ! dphi holds u during the iteration
         IF (use_lapre) CALL pw_zero(lapu)
         cgr%array(:, :, :) = resid%array(:, :, :)*seps%array(:, :, :)   ! r' = R/s
         CALL pw_transfer(cgr, wg2)
         wg2%array(:) = wg2%array(:)*greenfn%array(:)                ! F[z]
         CALL pw_copy(wg2, pg)
         CALL pw_transfer(wg2, cgz)
         CALL pw_copy(cgz, cgp)
         rz = pw_integral_ab(cgr, cgz)
         icg = 0
         DO icg = 1, ncg_max
            ! K p = (1/4pi)[F^-1(G^2 F[p]) + (q + sa/eps) p]  (1 FFT)
            DO ig = 1, SIZE(wg1%array)
               wg1%array(ig) = wg1%pw_grid%gsq(ig)*pg%array(ig)
            END DO
            CALL pw_transfer(wg1, cgd)
            IF (use_lapre) lapg%array(:, :, :) = cgd%array(:, :, :)   ! G^2 part, for the carried lap(s phi)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(cgd,cgp,lb,qfis,sa,seps,ub)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     cgd%array(i, j, k) = (cgd%array(i, j, k) + (qfis%array(i, j, k) + &
                                           sa%array(i, j, k)*seps%array(i, j, k)**2)*cgp%array(i, j, k))/fourpi
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
            pap = pw_integral_ab(cgp, cgd)
            IF (p%selfcheck .AND. iter == 1) THEN
               chk(1) = SQRT(pw_integral_ab(cgr, cgr)/vol)
               IF (p%output_unit > 0) WRITE (UNIT=p%output_unit, FMT="(T3,A,I4,4ES13.4)") &
                  "NLPB|   cg(fis) ", icg, rz, pap, rz/MAX(pap, 1.0E-300_dp), chk(1)
            END IF
            IF (pap <= 0.0_dp) THEN
               ! J is SPD for eps >= 1, S_ion >= 0; a non-positive curvature means an indefinite operator
               IF (p%io_unit > 0) WRITE (UNIT=p%io_unit, FMT="(T3,A,I5,A,ES12.4)") &
                  "NLPB| CG breakdown at CG ", icg, "  p.Jp ", pap
               EXIT
            END IF
            alpha = rz/pap
            CALL pw_axpy(cgp, dphi, alpha)
            IF (use_lapre) CALL pw_axpy(lapg, lapu, alpha)
            CALL pw_axpy(cgd, cgr, -alpha)
            ! convergence on the residual of the ORIGINAL system, s*r' = R - J dphi
            wr%array(:, :, :) = cgr%array(:, :, :)/seps%array(:, :, :)
            rres = SQRT(pw_integral_ab(wr, wr)/vol)
            CALL cg_track(icg, 0.5_dp*alpha*rz, rres)
            IF (rres <= cg_tol_in) EXIT
            IF (cg_estop) EXIT
            CALL pw_transfer(cgr, wg2)
            wg2%array(:) = wg2%array(:)*greenfn%array(:)
            CALL pw_transfer(wg2, cgz)
            rz_new = pw_integral_ab(cgr, cgz)
            beta_cg = rz_new/rz
            rz = rz_new
            CALL pw_scale(cgp, beta_cg)
            CALL pw_axpy(cgz, cgp)
            CALL pw_scale(pg, beta_cg)
            CALL pw_axpy(wg2, pg)
         END DO
         CALL pool%give_back_pw(pg)
         dphi%array(:, :, :) = dphi%array(:, :, :)*seps%array(:, :, :)   ! p = u/s
         lapu_ok = use_lapre
         ncg_tot = ncg_tot + MIN(icg, ncg_max)
      END SUBROUTINE run_cg_fis

! **************************************************************************************************
!> \brief bookkeeping of one CG iteration: ring buffer of the last 20 energy gains of the concave A
!>        (gain_j = alpha_j r_j.z_j / 2 is the exact increase of A in iteration j), RMS history samples, and the
!>        energy-gain stop. The sum of the last 20 gains bounds from below the error energy 20 iterations earlier and
!>        approximates it from above once the CG converges linearly; when it is below E_TOL the remaining error energy
!>        is far below E_TOL.
!> \param it ...
!> \param gain ...
!> \param rmsit ...
! **************************************************************************************************
      SUBROUTINE cg_track(it, gain, rmsit)
         INTEGER, INTENT(IN)                                :: it
         REAL(KIND=dp), INTENT(IN)                          :: gain, rmsit

         INTEGER                                            :: ih

         IF (it == 1) THEN
            cg_gain_ring = 0.0_dp
            cg_estop = .FALSE.
            cg_hist = -1.0_dp
         END IF
         cg_gain_ring(MOD(it - 1, 20) + 1) = gain
         cg_gain_sum = SUM(cg_gain_ring)
         ih = 0
         IF (it == 50) ih = 1
         IF (it == 100) ih = 2
         IF (it == 200) ih = 3
         IF (it == 400) ih = 4
         IF (it == 800) ih = 5
         IF (it == 1000) ih = 6
         IF (ih > 0) cg_hist(ih) = rmsit
         IF (etol_active .AND. it >= 20 .AND. cg_gain_sum < p%e_tol) cg_estop = .TRUE.
      END SUBROUTINE cg_track

   END SUBROUTINE nlpb_solve

! **************************************************************************************************
!> \brief ionic charge density rho_ion(r) = -(eps_b kappa0^2/4pi) S_ion g(x) phi for a given potential
!>        (same expression as in residual(); used for the IONIC_CHARGE_DENSITY cube print)
! **************************************************************************************************
   SUBROUTINE nlpb_ion_density(p, phi, s_ion, rho_ion)
      TYPE(nlpb_params_type), INTENT(IN)                 :: p
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: phi, s_ion
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho_ion

      INTEGER                                            :: i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: ad, afac, ek, gfac

      lb(1:3) = phi%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = phi%pw_grid%bounds_local(2, 1:3)
      ek = p%eps_b*p%kappa0sq/fourpi
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,gfac,afac,ad) SHARED(ek,lb,p,phi,rho_ion,s_ion,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               CALL ion_response(p, phi%array(i, j, k), gfac, afac, ad)
               rho_ion%array(i, j, k) = -ek*s_ion%array(i, j, k)*gfac*phi%array(i, j, k)
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
   END SUBROUTINE nlpb_ion_density

END MODULE qs_sccs_nlpb
