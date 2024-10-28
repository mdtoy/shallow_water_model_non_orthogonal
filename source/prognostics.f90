module prognostics

!-----------------------------------------------------------------
!
!   Contains arrays related to the prognostic variables, and
!       those diagnostic variables derived from the prognostic
!       variables.  Variables are initialized to zero here.
!       Surface geopotential also declared here.
!
!-----------------------------------------------------------------


use kinds
use model_parameters
use physical_parameters


implicit none
save


!
! Declare prognostic variables
!
real (kind = dbl_kind), dimension(im,jm,ntprog) ::          &
                     u1_cov, u2_cov,   &  ! covariant velocity components (m/s)
                     sqrt_G_h_star        ! sqrt_G*thickness of fluid (m) sqrt_G*(h-hs)

!
! Declare tendencies of prognostic variables
!
real (kind = dbl_kind), dimension(im,jm,nttend) ::          &
                     u1_cov_f,        &   ! d/dt (u1_cov)   (m/s^2)
                     u2_cov_f,        &   ! d/dt (u2_cov)   (m/s^2)
                     sqrt_G_h_star_f      ! d/dt (sqrt_G*h_star)   (m/s)

!
! Declare diagnostic variables
!
real (kind = dbl_kind), dimension(im,jm) ::                 &
                     ke_horiz,  &     ! contribution to kinetic energy
                                      ! from horizontal velocity (J/kg)
                     zeta,      &     ! relative vorticity (s^-1)
                     pv,        &     ! vertical component of potential 
                                      ! vorticity (m^2 eta/kg/s)
                     h,         &     ! height of free surface (m)
                     h_star,    &     ! thickness of fluid (m)
                     u_init_u1, &     ! initial Cartesian u and v (m/s)
                     u_init_u2, &     ! at u1 and u2 points
                     v_init_u1, &
                     v_init_u2, &
                     u, v,      &     ! Cartesian-based velocities (m/s)
                     u1_cont,   &     ! contravariant velocity
                     u2_cont          ! components (m/s)



!
! Declare topographic height
!
real (kind = dbl_kind), dimension (im,jm) ::                &
                     hs             ! topographic height (m)






contains


!======================================================================
! BEGINNING OF INIT_PROGNOSTICS
!======================================================================

subroutine init_prognostics

implicit none

! initialize prognostic arrays
u1_cov = c0
u2_cov = c0
sqrt_G_h_star = c0

u1_cov_f = c0
u2_cov_f = c0
sqrt_G_h_star_f = c0


end subroutine init_prognostics

!======================================================================
! END OF INIT_PROGNOSTICS
!======================================================================


end module prognostics
