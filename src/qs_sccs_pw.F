!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-3.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Pool-based plane-wave helpers for the SCCS solver: periodic Poisson solve with an
!>          explicit Green function, FFT convolution, and grid restriction / prolongation between MGRID levels.
!>          The solver uses them in place of pw_poisson_solve / pw_func_u_convolution so that every field may live
!>          on a coarser level than the Poisson environment's grid.
!> \author  ai-cp2k-lpb (2026-09-19)
! **************************************************************************************************
MODULE qs_sccs_pw

   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE pw_methods,                      ONLY: pw_copy,&
                                              pw_derive,&
                                              pw_integral_ab,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_sccs_pw'

   PUBLIC :: sol_poisson, sol_convolution, sol_convolution_g, sol_restrict, sol_prolong

CONTAINS

! **************************************************************************************************
!> \brief vhartree = F^-1[ g(G) F[density] ], g = greenfn or the periodic 4 pi/G^2 (G = 0 -> 0);
!>        optional Hartree energy 0.5 int rho v and gradient of the potential
! **************************************************************************************************
   SUBROUTINE sol_poisson(pool, density, vhartree, greenfn, ehartree, dvhartree)
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: density
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: vhartree
      TYPE(pw_c1d_gs_type), INTENT(IN), OPTIONAL         :: greenfn
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: ehartree
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(INOUT), &
         OPTIONAL                                        :: dvhartree

      INTEGER                                            :: i, ig
      INTEGER, DIMENSION(3)                              :: n
      TYPE(pw_c1d_gs_type)                               :: wg1, wg2

      CALL pool%create_pw(wg1)
      CALL pw_transfer(density, wg1)
      IF (PRESENT(greenfn)) THEN
         wg1%array(:) = wg1%array(:)*greenfn%array(:)
      ELSE
         DO ig = 1, SIZE(wg1%array)
            IF (wg1%pw_grid%gsq(ig) > 1.0E-12_dp) THEN
               wg1%array(ig) = fourpi*wg1%array(ig)/wg1%pw_grid%gsq(ig)
            ELSE
               wg1%array(ig) = 0.0_dp
            END IF
         END DO
      END IF
      CALL pw_transfer(wg1, vhartree)
      IF (PRESENT(ehartree)) ehartree = 0.5_dp*pw_integral_ab(density, vhartree)
      IF (PRESENT(dvhartree)) THEN
         CALL pool%create_pw(wg2)
         DO i = 1, 3
            n(:) = 0; n(i) = 1
            CALL pw_copy(wg1, wg2)
            CALL pw_derive(wg2, n)
            CALL pw_transfer(wg2, dvhartree(i))
         END DO
         CALL pool%give_back_pw(wg2)
      END IF
      CALL pool%give_back_pw(wg1)
   END SUBROUTINE sol_poisson

! **************************************************************************************************
!> \brief convolution = func * u through one G-space product (same semantics as pw_func_u_convolution)
! **************************************************************************************************
   SUBROUTINE sol_convolution(pool, func, convolution, u)
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: func
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: convolution
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: u

      TYPE(pw_c1d_gs_type)                               :: funcg, ug

      CALL pool%create_pw(funcg)
      CALL pool%create_pw(ug)
      CALL pw_transfer(func, funcg)
      CALL pw_transfer(u, ug)
      funcg%array(:) = funcg%array(:)*ug%array(:)
      CALL pw_transfer(funcg, convolution)
      CALL pool%give_back_pw(ug)
      CALL pool%give_back_pw(funcg)
   END SUBROUTINE sol_convolution

! **************************************************************************************************
!> \brief convolution = func * u with the kernel already in G space (ug = F[u] as sol_convolution would build it);
!>        same operations as sol_convolution minus the kernel transform
! **************************************************************************************************
   SUBROUTINE sol_convolution_g(pool, func, convolution, ug)
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: func
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: convolution
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: ug

      TYPE(pw_c1d_gs_type)                               :: funcg

      CALL pool%create_pw(funcg)
      CALL pw_transfer(func, funcg)
      funcg%array(:) = funcg%array(:)*ug(:)
      CALL pw_transfer(funcg, convolution)
      CALL pool%give_back_pw(funcg)
   END SUBROUTINE sol_convolution_g

! **************************************************************************************************
!> \brief restriction fine -> coarse level: keep the G vectors of the coarse sphere (orthogonal projector)
! **************************************************************************************************
   SUBROUTINE sol_restrict(fine_pool, coarse_pool, f_fine, f_coarse)
      TYPE(pw_pool_type), POINTER                        :: fine_pool, coarse_pool
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: f_fine
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: f_coarse

      TYPE(pw_c1d_gs_type)                               :: g_c, g_f

      IF (ASSOCIATED(fine_pool, coarse_pool)) THEN
         CALL pw_copy(f_fine, f_coarse)
         RETURN
      END IF
      CALL fine_pool%create_pw(g_f)
      CALL coarse_pool%create_pw(g_c)
      CALL pw_transfer(f_fine, g_f)
      CALL pw_zero(g_c)
      CALL pw_copy(g_f, g_c)
      CALL pw_transfer(g_c, f_coarse)
      CALL coarse_pool%give_back_pw(g_c)
      CALL fine_pool%give_back_pw(g_f)
   END SUBROUTINE sol_restrict

! **************************************************************************************************
!> \brief prolongation coarse -> fine level: zero padding in G space (transpose of the restriction)
! **************************************************************************************************
   SUBROUTINE sol_prolong(coarse_pool, fine_pool, f_coarse, f_fine)
      TYPE(pw_pool_type), POINTER                        :: coarse_pool, fine_pool
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: f_coarse
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: f_fine

      TYPE(pw_c1d_gs_type)                               :: g_c, g_f

      IF (ASSOCIATED(fine_pool, coarse_pool)) THEN
         CALL pw_copy(f_coarse, f_fine)
         RETURN
      END IF
      CALL coarse_pool%create_pw(g_c)
      CALL fine_pool%create_pw(g_f)
      CALL pw_transfer(f_coarse, g_c)
      CALL pw_zero(g_f)
      CALL pw_copy(g_c, g_f)
      CALL pw_transfer(g_f, f_fine)
      CALL fine_pool%give_back_pw(g_f)
      CALL coarse_pool%give_back_pw(g_c)
   END SUBROUTINE sol_prolong

END MODULE qs_sccs_pw
