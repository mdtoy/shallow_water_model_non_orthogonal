module physical_parameters

!-----------------------------------------------------------------
!	this modules specifices physical parameters and the units
!	  of those parameters (MKS is standard)
!-----------------------------------------------------------------

use kinds

implicit none
save

real (kind=dbl_kind), parameter ::                                &
       pi           = 3.14159265358979323846_dbl_kind,   &
       grav         = 9.8100_dbl_kind,                   &
       a_rad        = 6371220._dbl_kind, &   ! Earth radius
       Omega        = 7.292E-05_dbl_kind  ! Earth rotation rate

real (kind=dbl_kind), parameter ::                                &
                c0        = 0.00000_dbl_kind,                     &
                c1        = 1.00000_dbl_kind,                     &
                c2        = 2.00000_dbl_kind,                     &
                c3        = 3.00000_dbl_kind,                     &
                c4        = 4.00000_dbl_kind,                     &
                p5        = 0.50000_dbl_kind,                     &
                p25       = 0.25000_dbl_kind,                     &
                inv6      = c1/6.0_dbl_kind,                      &
                inv8      = c1/8.0_dbl_kind,                      &
                inv12     = c1/12.0_dbl_kind,                     &
                invgrav   = c1/grav,                              &
                inv_a_rad = c1/a_rad


!-------------------------------------------------------------------
!   VARIABLE DEFINITION
!-------------------------------------------------------------------
! pi = pi
! grav = gravitational acceleration (m/s^2)
! f_cor = Coriolis parameter
!-------------------------------------------------------------------

end module physical_parameters

