!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Nonlinear (rotationally saturating) dielectric response of the implicit solvent, pointwise part
!>          (after Islam et al., J. Chem. Phys. 159, 234117 (2023),
!>          Eq. (27)-(38), (44)-(45), (74)-(78); model of Gunceler et al., MSMSE 21, 074005 (2013)).
!>          Hartree atomic units, eps0 -> 1:
!>            alpha0_rot = (4 pi/3) beta p_mol^2,  alpha_pol = alpha0_rot (eps_inf-1)/(eps_b-eps_inf)
!>            1/alpha_sic = [(eps_b-eps_inf)/alpha0_rot - n_mol]/(eps_b-1)
!>            local field  E = -G + (4 pi/alpha_sic) p(E),  G = w_b * grad(phi)
!>            p_rot = (alpha0_rot/4 pi) g_rot(beta p_mol |E|) E,  p_pol = (alpha_pol/4 pi) E
!>            g_rot(x) = 3 (x coth x - 1)/x^2
!>            lambda_diel = lambda_rot + lambda_pol + lambda_sic,
!>              lambda_rot = -(1/beta) ln(sinh(x)/x),  lambda_pol = -(1/2)(4 pi/alpha_pol) p_pol^2,
!>              lambda_sic = +(1/2)(4 pi/alpha_sic) p^2
!>            chi_perp = n_mol S (1/alpha_perp - 1/alpha_sic)^-1,  alpha_perp = alpha_pol + alpha0_rot g_rot
!>            chi_par  = n_mol S (1/alpha_par  - 1/alpha_sic)^-1,  alpha_par  = alpha_perp + alpha0_rot g_rot'(x) x
!>          The collinear local-field equation |E| (1 - (4 pi/alpha_sic) chi_loc(|E|)) = |G| is solved per grid point
!>          by safeguarded Newton (chi_loc = (alpha_pol + alpha0_rot g_rot)/4 pi), E antiparallel to G.
!> \author  ai-cp2k-lpb (2026-09-10)
! **************************************************************************************************
MODULE qs_sccs_nldiel

   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_sccs_nldiel'

   PUBLIC :: nldiel_params_type, nldiel_setup, nldiel_point, g_rot

   TYPE nldiel_params_type
      REAL(KIND=dp) :: eps_b = 78.4_dp, eps_inf = 1.78_dp
      REAL(KIND=dp) :: n_mol = 0.0_dp      ! bulk solvent molecular density [bohr^-3]
      REAL(KIND=dp) :: p_mol = 0.0_dp      ! molecular dipole [e bohr]
      REAL(KIND=dp) :: beta = 0.0_dp       ! 1/(k_B T) [1/Hartree]
      REAL(KIND=dp) :: alpha0_rot = 0.0_dp, alpha_pol = 0.0_dp, inv_alpha_sic = 0.0_dp
      LOGICAL :: linear = .FALSE.          ! g_rot = 1 (linear response)
   END TYPE nldiel_params_type

CONTAINS

! **************************************************************************************************
!> \brief g_rot(x) = 3 (x coth x - 1)/x^2 and its derivative (series for small x)
! **************************************************************************************************
   ELEMENTAL SUBROUTINE g_rot(x, g, dg)
      REAL(KIND=dp), INTENT(IN)                          :: x
      REAL(KIND=dp), INTENT(OUT)                         :: g, dg

      REAL(KIND=dp)                                      :: c, s, xx

      xx = ABS(x)
      IF (xx < 1.0E-2_dp) THEN
         ! 3(x coth x - 1)/x^2 = 1 - x^2/15 + 2 x^4/315 - ...
         g = 1.0_dp - xx*xx/15.0_dp + 2.0_dp*xx**4/315.0_dp
         dg = -2.0_dp*xx/15.0_dp + 8.0_dp*xx**3/315.0_dp
      ELSE
         IF (xx > 30.0_dp) THEN
            c = 1.0_dp   ! coth x -> 1
            s = 0.0_dp   ! 1/sinh^2 -> 0
         ELSE
            c = COSH(xx)/SINH(xx)
            s = 1.0_dp/SINH(xx)**2
         END IF
         g = 3.0_dp*(xx*c - 1.0_dp)/(xx*xx)
         ! d/dx [3(x coth x - 1)/x^2] = 3[(coth x - x/sinh^2 x) x^2 - 2 x (x coth x - 1)]/x^4
         dg = 3.0_dp*((c - xx*s)*xx*xx - 2.0_dp*xx*(xx*c - 1.0_dp))/xx**4
      END IF
      IF (x < 0.0_dp) dg = -dg
   END SUBROUTINE g_rot

! **************************************************************************************************
!> \brief Derive alpha0_rot, alpha_pol, alpha_sic from eps_b, eps_inf, n_mol, p_mol, beta (Eq. 33, 44, 45)
! **************************************************************************************************
   SUBROUTINE nldiel_setup(p)
      TYPE(nldiel_params_type), INTENT(INOUT)            :: p

      IF (p%beta <= 0.0_dp .OR. p%n_mol <= 0.0_dp .OR. p%p_mol <= 0.0_dp) &
         CPABORT("SCCS nonlinear dielectric: beta, n_mol and p_mol must be > 0")
      IF (p%eps_b <= p%eps_inf .OR. p%eps_inf <= 1.0_dp) &
         CPABORT("SCCS nonlinear dielectric: need eps_b > eps_inf > 1")
      p%alpha0_rot = fourpi/3.0_dp*p%beta*p%p_mol**2
      p%alpha_pol = p%alpha0_rot*(p%eps_inf - 1.0_dp)/(p%eps_b - p%eps_inf)
      p%inv_alpha_sic = ((p%eps_b - p%eps_inf)/p%alpha0_rot - p%n_mol)/(p%eps_b - 1.0_dp)
   END SUBROUTINE nldiel_setup

! **************************************************************************************************
!> \brief Pointwise nonlinear dielectric response.
!> \param gmag   |G| = |w_b * grad phi| at the point
!> \param s      S_diel at the point
!> \param emag   out: |E| local field (E antiparallel to G)
!> \param pmag   out: |p| total molecular polarisation (p parallel to E)
!> \param lam    out: lambda_diel (ground-state molecular free energy density per molecule)
!> \param chi_perp, chi_par  out: local susceptibilities (Eq. 75-76), already including n_mol S
! **************************************************************************************************
   ELEMENTAL SUBROUTINE nldiel_point(p, gmag, s, emag, pmag, lam, chi_perp, chi_par)
      TYPE(nldiel_params_type), INTENT(IN)               :: p
      REAL(KIND=dp), INTENT(IN)                          :: gmag, s
      REAL(KIND=dp), INTENT(OUT)                         :: emag, pmag, lam, chi_perp, chi_par

      INTEGER                                            :: it
      REAL(KIND=dp)                                      :: a_par, a_perp, chi_loc, csic, dg, e, f, fp, g, &
                                                            lam_pol, lam_rot, lam_sic, ppol, prot, x

      csic = fourpi*p%inv_alpha_sic
      ! solve f(e) = e (1 - csic chi_loc(e)) - gmag = 0, chi_loc = (alpha_pol + alpha0_rot g_rot(beta p_mol e))/4pi
      ! low-field solution as start: e = gmag/(1 - csic (alpha_pol + alpha0_rot)/4pi)
      e = gmag/(1.0_dp - csic*(p%alpha_pol + p%alpha0_rot)/fourpi)
      IF (e < 0.0_dp) e = gmag
      DO it = 1, 50
         x = p%beta*p%p_mol*e
         IF (p%linear) THEN
            g = 1.0_dp; dg = 0.0_dp
         ELSE
            CALL g_rot(x, g, dg)
         END IF
         chi_loc = (p%alpha_pol + p%alpha0_rot*g)/fourpi
         f = e*(1.0_dp - csic*chi_loc) - gmag
         fp = 1.0_dp - csic*chi_loc - e*csic*p%alpha0_rot*dg*p%beta*p%p_mol/fourpi
         IF (ABS(fp) < 1.0E-30_dp) EXIT
         e = e - f/fp
         IF (e < 0.0_dp) e = 0.5_dp*ABS(e)
         IF (ABS(f) <= 1.0E-12_dp*MAX(1.0_dp, gmag)) EXIT
      END DO
      emag = e
      x = p%beta*p%p_mol*e
      IF (p%linear) THEN
         g = 1.0_dp; dg = 0.0_dp
      ELSE
         CALL g_rot(x, g, dg)
      END IF
      prot = p%alpha0_rot/fourpi*g*e
      ppol = p%alpha_pol/fourpi*e
      pmag = prot + ppol
      ! free energies (Eq. 36-38)
      IF (p%linear .OR. ABS(x) < 1.0E-3_dp) THEN
         lam_rot = -x*x/(6.0_dp*p%beta)                 ! -(1/beta) ln(sinh x/x) ~ -x^2/6 beta
      ELSE
         IF (x > 30.0_dp) THEN
            lam_rot = -(x - LOG(2.0_dp*x))/p%beta
         ELSE
            lam_rot = -LOG(SINH(x)/x)/p%beta
         END IF
      END IF
      lam_pol = -0.5_dp*fourpi/p%alpha_pol*ppol*ppol
      lam_sic = 0.5_dp*fourpi*p%inv_alpha_sic*pmag*pmag
      lam = lam_rot + lam_pol + lam_sic
      ! susceptibilities (Eq. 75-78), multiplied by n_mol S
      a_perp = p%alpha_pol + p%alpha0_rot*g
      a_par = a_perp + p%alpha0_rot*dg*x
      chi_perp = p%n_mol*s/(1.0_dp/a_perp - p%inv_alpha_sic)
      chi_par = p%n_mol*s/(1.0_dp/a_par - p%inv_alpha_sic)
   END SUBROUTINE nldiel_point

END MODULE qs_sccs_nldiel
