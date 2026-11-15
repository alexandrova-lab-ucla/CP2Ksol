!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Nonlocal solvent-excluded solute cavities for the SCCS implicit electrolyte.
!>          Implementation of Islam et al., J. Chem. Phys. 159, 234117 (2023),
!>          Eq. (1)-(13) (cavity chain) and Eq. (62)-(66) (chained functional derivatives).
!>
!>            x = ln(n/n_c),  S{n} = 1/2 erfc(x/(sigma sqrt 2)),  S'{n} = dS/dx = -exp(-x^2/2sigma^2)/(sqrt(2pi) sigma)
!>            S_vdW  = S{n_e}
!>            n_solv = n_c w_solv * (1 - S_vdW),   S_solv = S{n_solv}
!>            n_ion  = n_c w_ion  * (1 - S_vdW),   S_ion  = S{n_ion}
!>            n_diel = n_c w_diel * S_solv,        S_diel = 1 - S{n_diel}
!>            n_cav  = n_c w_cav  * S_solv,        S_cav  = 1 - S{n_cav}
!>            w_i(r) = [2b/(2b+R_i)] (1/4 pi b^3) exp(-(r - R_i)/b),  r = |r|,  b = a/sigma   (Eq. 7)
!>
!>          Kernel prefactor: Eq. (7) as printed in the paper is a monotonically decaying exponential of r (no
!>          absolute value): "This and other density kernels are taken to be exponential functions decaying over
!>          the length b = a/sigma", with the prefactor "chosen so that the effective electron density attains a
!>          value of n_c a distance R_i away from the S = 0.5 isosurface of a planar cavity function". That
!>          planar condition fixes the prefactor uniquely: for a planar step the half-space integral is
!>          H(z) = 2 pi C b^2 (z + 2b) exp((R_i - z)/b), and H(R_i) = 1 gives
!>          C = 1/(2 pi b^2 (2b + R_i)) = [2b/(2b + R_i)]/(4 pi b^3), i.e. exactly the printed prefactor.
!>          The kernel therefore saturates the interior (it does not vanish at r = 0). The growth toward
!>          r -> 0 is limited by a clamp on the exponent (xclamp), which never binds for R_i/b <= 200.
!>          Analytic normalisation of the model:  int w d^3r = 4 b exp(R_i/b)/(R_i + 2 b).
!>          Convolutions use pw_func_u_convolution (FFT); the kernel is pre-scaled by the cell volume so that
!>          the result is the true convolution (checked at runtime against the analytic int w d^3r, 2 %).
!>          Zero radius: the convolution is bypassed and the cavity equals its parent (S_solv <- S_vdW,
!>          S_ion <- S_vdW, S_diel <- S_solv, S_cav <- S_solv).
!>          Chained potentials: v_S{A,n}(r) = S'{n}(r) A(r)/n(r) (Eq. 62-63) applied through the chain with the
!>          (self-adjoint) convolution kernels (Eq. 64-66).
!> \author  ai-cp2k-lpb (2026-09-10)
! **************************************************************************************************
MODULE qs_sccs_cavity

   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi,&
                                              twopi
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integrate_function,&
                                              pw_transfer,&
                                              pw_zero
   USE qs_sccs_pw,                      ONLY: sol_convolution_g
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_sccs_cavity'

   INTEGER, PARAMETER, PUBLIC :: cav_solv = 1, cav_ion = 2, cav_diel = 3, cav_cav = 4

   PUBLIC :: cavity_params_type, cavity_state_type, cavity_build, cavity_release, cavity_potential_sum

   TYPE cavity_params_type
      REAL(KIND=dp) :: nc = 0.0_dp        ! cutoff density [bohr^-3]
      REAL(KIND=dp) :: sigma = 0.6_dp     ! width in ln(n)
      REAL(KIND=dp) :: a = 0.0_dp         ! kernel decay parameter a (b = a/sigma) [bohr]
      REAL(KIND=dp), DIMENSION(4) :: radius = 0.0_dp   ! R_solv, R_ion, R_diel, R_cav [bohr]
   END TYPE cavity_params_type

   TYPE cavity_state_type
      LOGICAL :: built = .FALSE.
      TYPE(pw_r3d_rs_type) :: n_e, n_solv, n_ion, n_diel, n_cav          ! effective densities of the chain
      TYPE(pw_r3d_rs_type) :: s_vdw, s_solv, s_ion, s_diel, s_cav        ! cavity functions
      LOGICAL, DIMENSION(4) :: conv = .FALSE.                            ! radius > 0: kernel built and convolved
      ! extra dE/dS_solv set by the caller (flexible region mask built from 1 - S_solv); owned by the caller
      LOGICAL :: has_a_solv = .FALSE.
      TYPE(pw_r3d_rs_type) :: a_solv
      TYPE(cavity_params_type) :: p
   END TYPE cavity_state_type

   ! The probe kernels depend only on the grid (cell, G-sphere, distribution), the kernel width b and the radii,
   ! never on the density, positions or SCF step. They are built once and kept in G space (volume-scaled, as
   ! sol_convolution would transform them); any change of the key below rebuilds them, so a cell change in
   ! GEO_OPT / CELL_OPT / NpT MD invalidates the cache exactly. Plain Fortran arrays (no pool object) so that a
   ! rebuilt pw environment cannot leave a dangling pool reference.
   TYPE kernel_cache_type
      LOGICAL :: valid = .FALSE.
      INTEGER :: id_nr = -1
      INTEGER, DIMENSION(3) :: npts = 0
      INTEGER, DIMENSION(2, 3) :: bounds = 0, bounds_local = 0
      REAL(KIND=dp), DIMENSION(3, 3) :: dh = 0.0_dp
      REAL(KIND=dp) :: b = 0.0_dp, vol = 0.0_dp
      REAL(KIND=dp), DIMENSION(4) :: radius = 0.0_dp
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: g
   END TYPE kernel_cache_type

   TYPE(kernel_cache_type), SAVE :: kc

CONTAINS

! **************************************************************************************************
!> \brief S{n} and S'{n} = dS/dx pointwise (n <= tiny -> S = 1, S' = 0)
! **************************************************************************************************
   ELEMENTAL SUBROUTINE shape(n, nc, sigma, s, sp)
      REAL(KIND=dp), INTENT(IN)                          :: n, nc, sigma
      REAL(KIND=dp), INTENT(OUT)                         :: s, sp

      REAL(KIND=dp), PARAMETER                           :: ntol = 1.0E-14_dp

      REAL(KIND=dp)                                      :: x

      IF (n < ntol) THEN
         s = 1.0_dp; sp = 0.0_dp
      ELSE
         x = LOG(n/nc)
         s = 0.5_dp*ERFC(x/(sigma*SQRT(2.0_dp)))
         sp = -EXP(-0.5_dp*(x/sigma)**2)/(SQRT(twopi)*sigma)
      END IF
   END SUBROUTINE shape

! **************************************************************************************************
!> \brief S{n} on a grid
! **************************************************************************************************
   SUBROUTINE shape_grid(n, nc, sigma, s, one_minus)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: n
      REAL(KIND=dp), INTENT(IN)                          :: nc, sigma
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: s
      LOGICAL, INTENT(IN)                                :: one_minus

      INTEGER                                            :: i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: sp, sv

      lb(1:3) = n%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = n%pw_grid%bounds_local(2, 1:3)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,sv,sp) SHARED(lb,n,nc,one_minus,s,sigma,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               CALL shape(n%array(i, j, k), nc, sigma, sv, sp)
               IF (one_minus) THEN
                  s%array(i, j, k) = 1.0_dp - sv
               ELSE
                  s%array(i, j, k) = sv
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
   END SUBROUTINE shape_grid

! **************************************************************************************************
!> \brief Kernel w_i on the grid (min-image distance from the cell origin), scaled by the cell volume so that
!>        pw_func_u_convolution(f, w) returns the true convolution f * w
! **************************************************************************************************
   SUBROUTINE build_kernel(w, b, radius)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: w
      REAL(KIND=dp), INTENT(IN)                          :: b, radius

      REAL(KIND=dp), PARAMETER                           :: xclamp = 2.0E2_dp

      INTEGER                                            :: i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      INTEGER, DIMENSION(3)                              :: npts
      INTEGER, DIMENSION(3)                              :: lbg
      INTEGER                                            :: s1, s2, s3
      REAL(KIND=dp)                                      :: f, r, r2, r2min, wint_ana, wint_grid
      REAL(KIND=dp), DIMENSION(3)                        :: rv, x
      REAL(KIND=dp), DIMENSION(3, 3)                     :: dh

      IF (radius > xclamp*b) CPABORT("SCCS cavity kernel: R_i/(a/sigma) too large")
      lb(1:3) = w%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = w%pw_grid%bounds_local(2, 1:3)
      npts(1:3) = w%pw_grid%npts(1:3)
      lbg(1:3) = w%pw_grid%bounds(1, 1:3)
      dh(1:3, 1:3) = w%pw_grid%dh(1:3, 1:3)
      f = 2.0_dp*b/(2.0_dp*b + radius)/(fourpi*b**3)*w%pw_grid%vol
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,r,x,rv,r2,r2min,s1,s2,s3) SHARED(lbg,b,dh,f,lb,npts,radius,ub,w)
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
               r = SQRT(r2min)
               ! Eq. (7): exp(-(r - R)/b), clamped at exp(+xclamp) so that the r -> 0 growth cannot overflow
               w%array(i, j, k) = f*EXP(-MAX((r - radius)/b, -xclamp))
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
      ! normalisation check: the grid sum must reproduce the analytic
      ! int w d^3r = 4 b exp(R/b)/(R + 2b) of Eq. (7). This gates the prefactor, the radial shape and the
      ! grid resolution at once (the periodic min-image tail beyond L/2 is < 1e-9 of the total for b << L).
      wint_grid = pw_integrate_function(w)/w%pw_grid%vol
      wint_ana = 4.0_dp*b*EXP(radius/b)/(radius + 2.0_dp*b)
      ! The analytic value is that of an ISOLATED kernel; on a periodic grid the images overlap once the kernel
      ! radius approaches half the cell (e.g. ION_RADIUS 4 A in an 11 A slab cell), which is physical for a slab.
      ! Warn, do not abort; a coarse grid shows up as a deviation at SMALL radii too.
      IF (ABS(wint_grid - wint_ana) > 2.0E-2_dp*wint_ana) &
         CPWARN("SCCS cavity kernel: periodic int w differs from the isolated analytic value by >2%")
   END SUBROUTINE build_kernel

! **************************************************************************************************
!> \brief Make kc hold the G-space kernels of this grid / b / radii (rebuild only when the key changed)
! **************************************************************************************************
   SUBROUTINE ensure_kernels(n_e, b, radius, pool)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: n_e
      REAL(KIND=dp), INTENT(IN)                          :: b
      REAL(KIND=dp), DIMENSION(4), INTENT(IN)            :: radius
      TYPE(pw_pool_type), POINTER                        :: pool

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cavity_kernels'

      INTEGER                                            :: handle, ic
      LOGICAL                                            :: key_ok
      TYPE(pw_c1d_gs_type)                               :: ug
      TYPE(pw_r3d_rs_type)                               :: w

      key_ok = kc%valid
      IF (key_ok) key_ok = (kc%id_nr == n_e%pw_grid%id_nr) .AND. ALL(kc%npts == n_e%pw_grid%npts) .AND. &
                           ALL(kc%bounds == n_e%pw_grid%bounds) .AND. ALL(kc%bounds_local == n_e%pw_grid%bounds_local) .AND. &
                           ALL(kc%dh == n_e%pw_grid%dh) .AND. (kc%vol == n_e%pw_grid%vol) .AND. (kc%b == b) .AND. &
                           ALL(kc%radius == radius)
      IF (key_ok) RETURN
      CALL timeset(routineN, handle)
      IF (ALLOCATED(kc%g)) DEALLOCATE (kc%g)
      kc%valid = .FALSE.
      DO ic = 1, 4
         IF (radius(ic) > 0.0_dp) THEN
            CALL pool%create_pw(w)
            CALL build_kernel(w, b, radius(ic))
            CALL pool%create_pw(ug)
            CALL pw_transfer(w, ug)
            IF (.NOT. ALLOCATED(kc%g)) ALLOCATE (kc%g(SIZE(ug%array), 4))
            kc%g(:, ic) = ug%array(:)
            CALL pool%give_back_pw(ug)
            CALL pool%give_back_pw(w)
         END IF
      END DO
      IF (.NOT. ALLOCATED(kc%g)) ALLOCATE (kc%g(1, 4))
      kc%id_nr = n_e%pw_grid%id_nr; kc%npts = n_e%pw_grid%npts; kc%bounds = n_e%pw_grid%bounds
      kc%bounds_local = n_e%pw_grid%bounds_local; kc%dh = n_e%pw_grid%dh; kc%vol = n_e%pw_grid%vol
      kc%b = b; kc%radius = radius
      kc%valid = .TRUE.
      CALL timestop(handle)
   END SUBROUTINE ensure_kernels

! **************************************************************************************************
!> \brief Build the five cavities from the solute electron density n_e
! **************************************************************************************************
   SUBROUTINE cavity_build(st, p, n_e, pool)
      TYPE(cavity_state_type), INTENT(INOUT)             :: st
      TYPE(cavity_params_type), INTENT(IN)               :: p
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: n_e
      TYPE(pw_pool_type), POINTER                        :: pool

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cavity_build'

      INTEGER                                            :: handle, ic
      REAL(KIND=dp)                                      :: b
      TYPE(pw_c1d_gs_type)                               :: fg, hg
      TYPE(pw_r3d_rs_type)                               :: tmp

      CALL timeset(routineN, handle)
      IF (st%built) CALL cavity_release(st, pool)
      st%p = p
      b = p%a/p%sigma
      IF (b <= 0.0_dp) CPABORT("SCCS cavity: kernel parameter a must be > 0")

      CALL pool%create_pw(st%n_e); CALL pool%create_pw(st%n_solv); CALL pool%create_pw(st%n_ion)
      CALL pool%create_pw(st%n_diel); CALL pool%create_pw(st%n_cav)
      CALL pool%create_pw(st%s_vdw); CALL pool%create_pw(st%s_solv); CALL pool%create_pw(st%s_ion)
      CALL pool%create_pw(st%s_diel); CALL pool%create_pw(st%s_cav)
      ! radius = 0 means "no displacement": the convolution is bypassed and the cavity equals its parent.
      ! Nothing else in the chain may use the kernel then.
      DO ic = 1, 4
         st%conv(ic) = (p%radius(ic) > 0.0_dp)
      END DO
      CALL ensure_kernels(n_e, b, p%radius, pool)
      CALL pool%create_pw(tmp)

      CALL pw_copy(n_e, st%n_e)
      CALL shape_grid(st%n_e, p%nc, p%sigma, st%s_vdw, .FALSE.)
      ! tmp = 1 - S_vdW
      tmp%array(:, :, :) = 1.0_dp - st%s_vdw%array(:, :, :)
      ! the solvent and ion kernels act on the same field: one forward transform
      IF (st%conv(cav_solv) .OR. st%conv(cav_ion)) THEN
         CALL pool%create_pw(fg)
         CALL pool%create_pw(hg)
         CALL pw_transfer(tmp, fg)
      END IF
      IF (st%conv(cav_solv)) THEN
         hg%array(:) = fg%array(:)*kc%g(:, cav_solv)
         CALL pw_transfer(hg, st%n_solv)
         st%n_solv%array(:, :, :) = p%nc*st%n_solv%array(:, :, :)
         CALL shape_grid(st%n_solv, p%nc, p%sigma, st%s_solv, .FALSE.)
      ELSE
         CALL pw_zero(st%n_solv); CALL pw_copy(st%s_vdw, st%s_solv)
      END IF
      IF (st%conv(cav_ion)) THEN
         hg%array(:) = fg%array(:)*kc%g(:, cav_ion)
         CALL pw_transfer(hg, st%n_ion)
         st%n_ion%array(:, :, :) = p%nc*st%n_ion%array(:, :, :)
         CALL shape_grid(st%n_ion, p%nc, p%sigma, st%s_ion, .FALSE.)
      ELSE
         CALL pw_zero(st%n_ion); CALL pw_copy(st%s_vdw, st%s_ion)
      END IF
      IF (st%conv(cav_solv) .OR. st%conv(cav_ion)) THEN
         CALL pool%give_back_pw(hg)
         CALL pool%give_back_pw(fg)
      END IF
      IF (st%conv(cav_diel)) THEN
         CALL sol_convolution_g(pool=pool, func=st%s_solv, convolution=st%n_diel, ug=kc%g(:, cav_diel))
         st%n_diel%array(:, :, :) = p%nc*st%n_diel%array(:, :, :)
         CALL shape_grid(st%n_diel, p%nc, p%sigma, st%s_diel, .TRUE.)
      ELSE
         CALL pw_zero(st%n_diel); CALL pw_copy(st%s_solv, st%s_diel)
      END IF
      IF (st%conv(cav_cav)) THEN
         CALL sol_convolution_g(pool=pool, func=st%s_solv, convolution=st%n_cav, ug=kc%g(:, cav_cav))
         st%n_cav%array(:, :, :) = p%nc*st%n_cav%array(:, :, :)
         CALL shape_grid(st%n_cav, p%nc, p%sigma, st%s_cav, .TRUE.)
      ELSE
         CALL pw_zero(st%n_cav); CALL pw_copy(st%s_solv, st%s_cav)
      END IF

      CALL pool%give_back_pw(tmp)
      st%built = .TRUE.
      CALL timestop(handle)
   END SUBROUTINE cavity_build

! **************************************************************************************************
!> \brief Release all fields
! **************************************************************************************************
   SUBROUTINE cavity_release(st, pool)
      TYPE(cavity_state_type), INTENT(INOUT)             :: st
      TYPE(pw_pool_type), POINTER                        :: pool

      IF (.NOT. st%built) RETURN
      CALL pool%give_back_pw(st%n_e); CALL pool%give_back_pw(st%n_solv); CALL pool%give_back_pw(st%n_ion)
      CALL pool%give_back_pw(st%n_diel); CALL pool%give_back_pw(st%n_cav)
      CALL pool%give_back_pw(st%s_vdw); CALL pool%give_back_pw(st%s_solv); CALL pool%give_back_pw(st%s_ion)
      CALL pool%give_back_pw(st%s_diel); CALL pool%give_back_pw(st%s_cav)
      st%conv(:) = .FALSE.
      st%built = .FALSE.
   END SUBROUTINE cavity_release

! **************************************************************************************************
!> \brief v(r) = S'{n}(r) A(r)/n(r)  (Eq. 62-63), in place on A
! **************************************************************************************************
   SUBROUTINE vs_apply(a, n, nc, sigma)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: a
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: n
      REAL(KIND=dp), INTENT(IN)                          :: nc, sigma

      REAL(KIND=dp), PARAMETER                           :: ntol = 1.0E-14_dp

      INTEGER                                            :: i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: sp, sv

      lb(1:3) = n%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = n%pw_grid%bounds_local(2, 1:3)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,sv,sp) SHARED(a,lb,n,nc,sigma,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               IF (n%array(i, j, k) < ntol) THEN
                  a%array(i, j, k) = 0.0_dp
               ELSE
                  CALL shape(n%array(i, j, k), nc, sigma, sv, sp)
                  a%array(i, j, k) = sp*a%array(i, j, k)/n%array(i, j, k)
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
   END SUBROUTINE vs_apply

! **************************************************************************************************
!> \brief Summed functional derivative  v(r) = delta/delta n_e(r) [ int S_diel A_d + S_ion A_i + S_cav A_c ]  through the
!>        chain (Eq. 64-66). The three chains are linear in their A fields and meet at shared links, so they are summed
!>        there: the dielectric and cavitation fields share the S_solv link; that link and the ion link end at the
!>        1 - S_vdW level, where they are added in G space and leave with one inverse transform.
!>        a_diel (dielectric), a_ion and a_cav (optional) are overwritten with work data; the sum is returned in v.
!>        S_ion  = S{n_c w_ion*(1-S_vdW)}
!>        S_diel = 1 - S{n_c w_diel * S{n_c w_solv*(1-S_vdW)}}   (S_cav analogous with w_cav)
! **************************************************************************************************
   SUBROUTINE cavity_potential_sum(st, a_diel, v, pool, a_ion, a_cav)
      TYPE(cavity_state_type), INTENT(IN)                :: st
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: a_diel, v
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_r3d_rs_type), INTENT(INOUT), OPTIONAL      :: a_ion, a_cav

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cavity_potential'

      INTEGER                                            :: handle
      LOGICAL                                            :: have_g
      REAL(KIND=dp)                                      :: nc, sigma
      TYPE(pw_c1d_gs_type)                               :: acc, gt
      TYPE(pw_r3d_rs_type)                               :: tmp

      CALL timeset(routineN, handle)
      nc = st%p%nc; sigma = st%p%sigma
      CALL pool%create_pw(tmp)
      CALL pw_zero(v)
      ! dielectric: v'' = vS{-A, n_diel}; v' = vS{n_c (w_diel * v''), n_solv}
      IF (st%conv(cav_diel)) THEN
         a_diel%array(:, :, :) = -a_diel%array(:, :, :)
         CALL vs_apply(a_diel, st%n_diel, nc, sigma)
         CALL sol_convolution_g(pool=pool, func=a_diel, convolution=tmp, ug=kc%g(:, cav_diel))
         a_diel%array(:, :, :) = nc*tmp%array(:, :, :)
      END IF
      ! cavitation joins at the S_solv level (a bypassed kernel makes S_cav its parent)
      IF (PRESENT(a_cav)) THEN
         IF (st%conv(cav_cav)) THEN
            a_cav%array(:, :, :) = -a_cav%array(:, :, :)
            CALL vs_apply(a_cav, st%n_cav, nc, sigma)
            CALL sol_convolution_g(pool=pool, func=a_cav, convolution=tmp, ug=kc%g(:, cav_cav))
            a_cav%array(:, :, :) = nc*tmp%array(:, :, :)
         END IF
         a_diel%array(:, :, :) = a_diel%array(:, :, :) + a_cav%array(:, :, :)
      END IF
      ! a_diel is dE/dS_solv here: add what depends on S_solv directly
      IF (st%has_a_solv) a_diel%array(:, :, :) = a_diel%array(:, :, :) + st%a_solv%array(:, :, :)
      have_g = .FALSE.
      ! S_solv link: v = vS{-n_c (w_solv * vS{a, n_solv}), n_e}; a bypassed link leaves a at the S_vdW level
      IF (st%conv(cav_solv)) THEN
         CALL vs_apply(a_diel, st%n_solv, nc, sigma)
         CALL pool%create_pw(acc); CALL pool%create_pw(gt)
         CALL pw_transfer(a_diel, gt)
         acc%array(:) = (-nc)*(gt%array(:)*kc%g(:, cav_solv))
         have_g = .TRUE.
      ELSE
         CALL pw_axpy(a_diel, v)
      END IF
      ! ion: v = vS{-n_c (w_ion * vS{A, n_ion}), n_e}; bypassed kernel: v = vS{A, n_e}
      IF (PRESENT(a_ion)) THEN
         IF (st%conv(cav_ion)) THEN
            CALL vs_apply(a_ion, st%n_ion, nc, sigma)
            IF (.NOT. have_g) THEN
               CALL pool%create_pw(acc); CALL pool%create_pw(gt)
               acc%array(:) = (0.0_dp, 0.0_dp)
               have_g = .TRUE.
            END IF
            CALL pw_transfer(a_ion, gt)
            acc%array(:) = acc%array(:) + (-nc)*(gt%array(:)*kc%g(:, cav_ion))
         ELSE
            CALL pw_axpy(a_ion, v)
         END IF
      END IF
      IF (have_g) THEN
         CALL pw_transfer(acc, tmp)
         CALL pw_axpy(tmp, v)
         CALL pool%give_back_pw(gt); CALL pool%give_back_pw(acc)
      END IF
      CALL vs_apply(v, st%n_e, nc, sigma)
      CALL pool%give_back_pw(tmp)
      CALL timestop(handle)
   END SUBROUTINE cavity_potential_sum

END MODULE qs_sccs_cavity
