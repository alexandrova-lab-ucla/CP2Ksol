!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief State and integrator of the electron-number dynamics of the grand-canonical electrode
!>        (&SCF &GCE UPDATE_SCHEME DYNAMICS).
!>
!>        The electron number n is a dynamical variable with the fictitious mass M_n and the momentum P_n,
!>        driven by the difference between the electrode potential U and its target U_T and, optionally,
!>        coupled to a Nose-Hoover chain of its own:
!>
!>           dn/dt    = P_n / M_n
!>           dP_n/dt  = (U - U_T) - P_n P_xi1 / Q_1
!>           dxi_1/dt = P_xi1 / Q_1
!>           dP_xi1/dt = P_n^2 / M_n - kT - P_xi1 P_xi2 / Q_2
!>           dxi_j/dt = P_xij / Q_j,   dP_xij/dt = P_xi(j-1)^2 / Q_(j-1) - kT - P_xij P_xi(j+1) / Q_(j+1)
!>
!>        (last link without the drag term). U - U_T = WF - WF_T = -dOmega/dn with the grand potential
!>        Omega = E + WF_T (n - n_input), so the conserved quantity of the chain-length-0 flow is
!>        K_n + Omega with K_n = P_n^2 / (2 M_n); with a chain it is K_n + Omega + E_th,
!>        E_th = sum_j P_xij^2 / (2 Q_j) + kT xi_1 + kT sum_(j>1) xi_j.
!>
!>        The chain is integrated with the symmetric Suzuki-Yoshida factorisation of Martyna, Tuckerman,
!>        Tobias and Klein (Mol. Phys. 87, 1117 (1996), https://doi.org/10.1080/00268979600100761; chain
!>        itself: Martyna, Klein, Tuckerman, J. Chem. Phys. 97, 2635 (1992), https://doi.org/10.1063/1.463940).
!>        Units are Hartree atomic units (energy Ha, charge e, time a.u.): U in Ha per e, n in e,
!>        P_n in Ha a.u./e, M_n in Ha a.u.^2/e^2.
!> \author  ai-cp2k-lpb (2026-11-08)
! **************************************************************************************************
MODULE gce_dyn_types

   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gce_dyn_types'

   INTEGER, PARAMETER, PUBLIC           :: gce_dyn_max_chain = 16
   ! Suzuki-Yoshida order (number of weights) and number of sub-steps of the chain propagator
   INTEGER, PARAMETER, PRIVATE          :: n_ys = 5, n_c = 1

   TYPE, PUBLIC :: gce_dyn_type
      LOGICAL                                :: active = .FALSE.
      ! ---- input (Hartree atomic units) ----
      REAL(KIND=dp)                          :: period = 0.0_dp, &              ! POTENTIOSTAT_PERIOD
                                                thermostat_period = 0.0_dp, &   ! POTENTIOSTAT_THERMOSTAT_PERIOD (0 = period)
                                                kt_input = -1.0_dp, &           ! POTENTIOSTAT_TEMPERATURE as kT (<0 = unset)
                                                capacitance_input = 0.0_dp, &   ! POTENTIOSTAT_CAPACITANCE [e/Ha] (0 = measured)
                                                initial_rate = 0.0_dp           ! POTENTIOSTAT_INITIAL_RATE [e/a.u.]
      INTEGER                                :: chain_length = 0              ! POTENTIOSTAT_CHAIN_LENGTH (0 = constant energy)
      ! ---- read from MOTION&MD by the hook in qs_init ----
      LOGICAL                                :: have_md = .FALSE., is_md = .FALSE.
      REAL(KIND=dp)                          :: dt = 0.0_dp, md_kt = -1.0_dp
      ! ---- derived at the end of the first SCF (or restored) ----
      LOGICAL                                :: masses_ready = .FALSE.
      REAL(KIND=dp)                          :: mass_n = 0.0_dp, kt_xi = 0.0_dp
      REAL(KIND=dp), DIMENSION(gce_dyn_max_chain) :: mass_xi = 0.0_dp
      ! ---- dynamical state: the named fields that go through the restart ----
      LOGICAL                                :: open_pending = .FALSE.        ! the opening half-step of the next MD step is due
      INTEGER                                :: next_len = 1, &               ! MD steps until the next potential measurement
                                                ion_step_saved = 0            ! gce%ion_step after the last completed step
      REAL(KIND=dp)                          :: n_e = 0.0_dp, &               ! electron number n [e]
                                                p_n = 0.0_dp, &               ! P_n
                                                n_input = 0.0_dp, &           ! electron number of the input charge state
                                                du_last = 0.0_dp, &           ! (U - U_T) of the last measurement [Ha/e]
                                                capacitance = 0.0_dp          ! C used for M_n [e/Ha]
      REAL(KIND=dp), DIMENSION(gce_dyn_max_chain) :: xi = 0.0_dp, p_xi = 0.0_dp
      ! ---- transient (not persisted) ----
      LOGICAL                                :: initialised = .FALSE., restored = .FALSE., &
                                                restart_valid = .FALSE., is_initial_call = .FALSE., &
                                                is_repeat = .FALSE., measure_now = .TRUE., have_wf_on = .FALSE., &
                                                dt_baseline = .FALSE.
      INTEGER                                :: last_sim_step = 0, n_steps = 0, max_outer_cap = 0
      REAL(KIND=dp)                          :: wf_on = 0.0_dp, e_on = 0.0_dp, dt_last = 0.0_dp, sim_time_last = 0.0_dp
   END TYPE gce_dyn_type

   PUBLIC :: gce_dyn_set_masses, gce_dyn_kick, gce_dyn_drift, gce_dyn_nhc_propagate, &
             gce_dyn_kinetic_n, gce_dyn_thermostat_energy, gce_dyn_next_len

CONTAINS

! **************************************************************************************************
!> \brief M_n = 1/(C omega^2) with omega = 2 pi / period; Q_j = kT / omega_th^2 for every link.
!> \param dyn ...
!> \param capacitance C [e/Ha], measured dn/dU (positive)
!> \param kt kT of the potentiostat thermostat [Ha]
! **************************************************************************************************
   SUBROUTINE gce_dyn_set_masses(dyn, capacitance, kt)
      TYPE(gce_dyn_type), INTENT(INOUT)                  :: dyn
      REAL(KIND=dp), INTENT(IN)                          :: capacitance, kt

      REAL(KIND=dp)                                      :: omega, omega_th, twopi_local

      twopi_local = 8.0_dp*ATAN(1.0_dp)
      CPASSERT(capacitance > 0.0_dp)
      CPASSERT(dyn%period > 0.0_dp)
      omega = twopi_local/dyn%period
      dyn%capacitance = capacitance
      dyn%mass_n = 1.0_dp/(capacitance*omega*omega)
      dyn%kt_xi = kt
      IF (dyn%thermostat_period > 0.0_dp) THEN
         omega_th = twopi_local/dyn%thermostat_period
      ELSE
         omega_th = omega
      END IF
      dyn%mass_xi = 0.0_dp
      IF (dyn%chain_length > 0) THEN
         CPASSERT(kt > 0.0_dp)
         dyn%mass_xi(1:dyn%chain_length) = kt/(omega_th*omega_th)
      END IF
      dyn%masses_ready = .TRUE.
   END SUBROUTINE gce_dyn_set_masses

! **************************************************************************************************
!> \brief P_n <- P_n + h (U - U_T); h is the duration of the kick
! **************************************************************************************************
   SUBROUTINE gce_dyn_kick(dyn, h, du)
      TYPE(gce_dyn_type), INTENT(INOUT)                  :: dyn
      REAL(KIND=dp), INTENT(IN)                          :: h, du

      dyn%p_n = dyn%p_n + h*du
   END SUBROUTINE gce_dyn_kick

! **************************************************************************************************
!> \brief n <- n + dt P_n / M_n
! **************************************************************************************************
   SUBROUTINE gce_dyn_drift(dyn, dt)
      TYPE(gce_dyn_type), INTENT(INOUT)                  :: dyn
      REAL(KIND=dp), INTENT(IN)                          :: dt

      dyn%n_e = dyn%n_e + dt*dyn%p_n/dyn%mass_n
   END SUBROUTINE gce_dyn_drift

! **************************************************************************************************
!> \brief Nose-Hoover chain acting on P_n for the duration h (a half step of the MD step when h = dt/2):
!>        symmetric Strang factorisation (force and drag sweeps down and up, particle scaling and chain
!>        positions in between), repeated with the Suzuki-Yoshida weights.
!> \param dyn ...
!> \param h duration [a.u.]
! **************************************************************************************************
   SUBROUTINE gce_dyn_nhc_propagate(dyn, h)
      TYPE(gce_dyn_type), INTENT(INOUT)                  :: dyn
      REAL(KIND=dp), INTENT(IN)                          :: h

      INTEGER                                            :: ic, iy, j, m
      REAL(KIND=dp)                                      :: aa, hs, w1
      REAL(KIND=dp), DIMENSION(n_ys)                     :: w

      m = dyn%chain_length
      IF (m < 1) RETURN
      ! Suzuki-Yoshida weights of order 4 (five factors)
      w1 = 1.0_dp/(4.0_dp - 4.0_dp**(1.0_dp/3.0_dp))
      w = [w1, w1, 1.0_dp - 4.0_dp*w1, w1, w1]

      DO iy = 1, n_ys
         hs = w(iy)*h/REAL(n_c, KIND=dp)
         DO ic = 1, n_c
            ! ---- down sweep: forces and drags, last link first ----
            CALL chain_force(dyn, m, 0.5_dp*hs)
            DO j = m - 1, 1, -1
               aa = EXP(-0.25_dp*hs*dyn%p_xi(j + 1)/dyn%mass_xi(j + 1))
               dyn%p_xi(j) = dyn%p_xi(j)*aa
               CALL chain_force(dyn, j, 0.5_dp*hs)
               dyn%p_xi(j) = dyn%p_xi(j)*aa
            END DO
            ! ---- scale P_n by the first link, advance the chain positions ----
            aa = EXP(-hs*dyn%p_xi(1)/dyn%mass_xi(1))
            dyn%p_n = dyn%p_n*aa
            DO j = 1, m
               dyn%xi(j) = dyn%xi(j) + hs*dyn%p_xi(j)/dyn%mass_xi(j)
            END DO
            ! ---- up sweep ----
            DO j = 1, m - 1
               aa = EXP(-0.25_dp*hs*dyn%p_xi(j + 1)/dyn%mass_xi(j + 1))
               dyn%p_xi(j) = dyn%p_xi(j)*aa
               CALL chain_force(dyn, j, 0.5_dp*hs)
               dyn%p_xi(j) = dyn%p_xi(j)*aa
            END DO
            CALL chain_force(dyn, m, 0.5_dp*hs)
         END DO
      END DO
   END SUBROUTINE gce_dyn_nhc_propagate

! **************************************************************************************************
!> \brief P_xij <- P_xij + h G_j with G_1 = P_n^2/M_n - kT and G_j = P_xi(j-1)^2/Q_(j-1) - kT
! **************************************************************************************************
   SUBROUTINE chain_force(dyn, j, h)
      TYPE(gce_dyn_type), INTENT(INOUT)                  :: dyn
      INTEGER, INTENT(IN)                                :: j
      REAL(KIND=dp), INTENT(IN)                          :: h

      REAL(KIND=dp)                                      :: g

      IF (j == 1) THEN
         g = dyn%p_n*dyn%p_n/dyn%mass_n - dyn%kt_xi
      ELSE
         g = dyn%p_xi(j - 1)*dyn%p_xi(j - 1)/dyn%mass_xi(j - 1) - dyn%kt_xi
      END IF
      dyn%p_xi(j) = dyn%p_xi(j) + h*g
   END SUBROUTINE chain_force

! **************************************************************************************************
!> \brief K_n = P_n^2 / (2 M_n)
! **************************************************************************************************
   FUNCTION gce_dyn_kinetic_n(dyn) RESULT(ek)
      TYPE(gce_dyn_type), INTENT(IN)                     :: dyn
      REAL(KIND=dp)                                      :: ek

      ek = 0.0_dp
      IF (dyn%mass_n > 0.0_dp) ek = 0.5_dp*dyn%p_n*dyn%p_n/dyn%mass_n
   END FUNCTION gce_dyn_kinetic_n

! **************************************************************************************************
!> \brief E_th = sum_j P_xij^2/(2 Q_j) + kT xi_1 + kT sum_(j>1) xi_j
! **************************************************************************************************
   FUNCTION gce_dyn_thermostat_energy(dyn) RESULT(eth)
      TYPE(gce_dyn_type), INTENT(IN)                     :: dyn
      REAL(KIND=dp)                                      :: eth

      INTEGER                                            :: j

      eth = 0.0_dp
      DO j = 1, dyn%chain_length
         eth = eth + 0.5_dp*dyn%p_xi(j)*dyn%p_xi(j)/dyn%mass_xi(j) + dyn%kt_xi*dyn%xi(j)
      END DO
   END FUNCTION gce_dyn_thermostat_energy

! **************************************************************************************************
!> \brief number of MD steps from force evaluation index ion_step to the next potential measurement.
!>        Without SOLVENT_ONLY_FOR_CHARGE the potential is measured every step (1). With it the solvent is
!>        on at ion_step = 0 and whenever MOD(ion_step + 1, interval) == 0.
! **************************************************************************************************
   FUNCTION gce_dyn_next_len(ion_step, interval, flash) RESULT(len_steps)
      INTEGER, INTENT(IN)                                :: ion_step, interval
      LOGICAL, INTENT(IN)                                :: flash
      INTEGER                                            :: len_steps

      INTEGER                                            :: m_next

      IF (.NOT. flash .OR. interval <= 1) THEN
         len_steps = 1
      ELSE
         m_next = ((ion_step + 1)/interval + 1)*interval - 1
         len_steps = m_next - ion_step
      END IF
   END FUNCTION gce_dyn_next_len

END MODULE gce_dyn_types
