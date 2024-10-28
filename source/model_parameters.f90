module model_parameters

!-----------------------------------------------------------------
!    This module specifies the grid size, time step length, and
!    run duration.
!-----------------------------------------------------------------


use kinds
use physical_parameters


implicit none
save

!
! Set grid dimensions
!
integer, parameter :: im = 100,    &    ! number of x1-direction grid points
                      jm = 100          ! number of x2-direction grid points

real (kind=dbl_kind), parameter :: &
                      dx1 = c2*pi*a_rad/im,     &   ! x1-direction grid spacing (m)
                      dx2 = c2*pi*a_rad/jm          ! x2-direction grid spacing (m)


!
! Set time step length
!
real (kind=dbl_kind), parameter :: dt = 30._dbl_kind      ! time step length (s)


!
! Set termination time
!
real (kind=dbl_kind), parameter ::  &
                       tau_duration = 5*dt/3600   !  168._dbl_kind        ! run duration (hours)


!
! Set frequency of output
! (i.e. number of timesteps per output calls)
!
integer, parameter :: out_freq = 1   !  360


real (kind=dbl_kind), parameter :: &
                      invdx1  = 1.0_dbl_kind/dx1,     &   ! inverse dx1
                      invdx2  = 1.0_dbl_kind/dx2          ! inverse dx2


integer, parameter ::   &
           ntprog = 2,  &   ! no. of time slots needed for prognostic variables
           nttend = 3       ! no. of time slots needed for tendency variables



! Declare horizontal grid point "+n/-n" indices:
!  ip1 represents "i + 1", im1 represents "i - 1"
!  ip2 represents "i + 2", im2 represents "i - 2"
!  Periodic boundary conditions are accounted for in that grid point
!  index zero has the value im, and index im+1 gets the value 1, etc.
!  Ditto the above for "j".
!  Indices are set in subroutine initialize_model.
integer, dimension (1:im) :: ip1, im1, ip2, im2
integer, dimension (1:jm) :: jp1, jm1, jp2, jm2



!
! Specify weighting factors for time-stepping
!

! Euler forward factors
real (kind=dbl_kind), parameter ::    &
          w3_ef = 1.0_dbl_kind,       &
          w2_ef = 0.0_dbl_kind,       &
          w1_ef = 0.0_dbl_kind

! Adams-Bashworth 2nd. order factors
real (kind=dbl_kind), parameter ::    &
          w3_ab2 =  1.5_dbl_kind,     &
          w2_ab2 = -0.5_dbl_kind,     &
          w1_ab2 =  0.0_dbl_kind

! Adams-Bashworth 3rd. order factors
real (kind=dbl_kind), parameter ::    &
          w3_ab3 =  1.91666666666666666666666666666666_dbl_kind,     &
          w2_ab3 = -1.33333333333333333333333333333333_dbl_kind,     &
          w1_ab3 =  0.41666666666666666666666666666666_dbl_kind


real (kind = dbl_kind), dimension(im,jm) ::                       &
          f_cor       ! Coriolis parameter (s^-1)


! Domain coordinate (x1-x2) values
real (kind=dbl_kind), dimension (1:im) :: x1_h, x1_u1
real (kind=dbl_kind), dimension (1:jm) :: x2_h, x2_u2


! Domain Cartesian (x-y) values
real (kind=dbl_kind), dimension (1:im,1:jm) :: x_h, x_u, x_v, x_q
real (kind=dbl_kind), dimension (1:im,1:jm) :: y_h, y_u, y_v, y_q

! Non-orthogonl coordinate parameters
real (kind=dbl_kind), parameter ::    &
          C_fac = c0,  &  ! 0.5_dbl_kind,   &
          D_fac = c0,  &  ! 1.0_dbl_kind,   &
          alf   = 0.0_dbl_kind,   &  ! latitude offset of coordinate
          bet   = 0.0_dbl_kind       ! longitude offset of coordinate


! Coordinate tranformation metric terms
real (kind=dbl_kind), dimension (1:im,1:jm) ::     &
          sqrt_G_h, inv_sqrt_G_h, inv_sqrt_G_q,    &
          G_cont_11_u1, G_cont_22_u2, G_cont_12_u1, G_cont_12_u2

! Components of the transformation Jacobian
real (kind=dbl_kind), dimension (1:im,1:jm) ::     &
          dx1_dx_u1, dx1_dy_u2, dx2_dx_u1, dx2_dy_u2,  &
          dx_dx1_u1, dx_dx2_u2, dy_dx1_u1, dy_dx2_u2

! G_cont_12 weightings for averaging contravariant/covariant velocities
real (kind=dbl_kind), dimension (1:im,1:jm) ::     &
          mean_G_cont_12_u1_nw, mean_G_cont_12_u1_ne,   &
          mean_G_cont_12_u1_sw, mean_G_cont_12_u1_se,   &
          mean_G_cont_12_u2_nw, mean_G_cont_12_u2_ne,   &
          mean_G_cont_12_u2_sw, mean_G_cont_12_u2_se




end module model_parameters
