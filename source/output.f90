module output

!-----------------------------------------------------------------
!   This module outputs the model data in ascii and/or binary 
!   form.
!-----------------------------------------------------------------

use kinds
use model_parameters
use physical_parameters
use prognostics
use step


implicit none


contains



subroutine initial_output (tau)

implicit none

real (kind = dbl_kind), intent(in) :: tau    ! Time in hours

integer :: nt,i,j
integer :: ntm

ntm = int(tau_duration*3600._dbl_kind/(dt*out_freq))+1

write (31) ntm, im, jm
write (31) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (31) (x1_h(i),                  i = 1,im)
write (31) (x1_u1(i),                 i = 1,im)
write (31) (x2_h(j),                  j = 1,jm)
write (31) (x2_u2(j),                 j = 1,jm)


do j = 1,jm
   write (31) (hs(i,j),        i = 1,im)
end do


! Write out Cartesian-coordinate data

do j = 1,jm
   write (31) (x_h(i,j),        i = 1,im)
end do

do j = 1,jm
   write (31) (y_h(i,j),        i = 1,im)
end do

do j = 1,jm
   write (31) (x_u(i,j),        i = 1,im)
end do

do j = 1,jm
   write (31) (y_u(i,j),        i = 1,im)
end do

do j = 1,jm
   write (31) (x_v(i,j),        i = 1,im)
end do

do j = 1,jm
   write (31) (y_v(i,j),        i = 1,im)
end do

do j = 1,jm
   write (31) (x_q(i,j),        i = 1,im)
end do

do j = 1,jm
   write (31) (y_q(i,j),        i = 1,im)
end do

! Output non-orthogonal coordinate field as function of x and y (Cartesian coords)
! NOTE:  X1 and X2 coords stand in for x and y
! X1-coord
do j = 1,jm
   write (31) (x1_h(i)+a_rad*C_fac*Sin(x2_h(j)*inv_a_rad-alf),   i = 1,im)
end do
! X2-coord
do j = 1,jm
   write (31) (x2_h(j)+a_rad*D_fac*Sin(x1_h(i)*inv_a_rad-bet),   i = 1,im)
end do


! Output metric terms
do j = 1,jm
   write (31) (sqrt_G_h(i,j),                                i = 1,im)
end do
do j = 1,jm
   write (31) (G_cont_11_u1(i,j),                            i = 1,im)
end do
do j = 1,jm
   write (31) (G_cont_22_u2(i,j),                            i = 1,im)
end do
do j = 1,jm
   write (31) (G_cont_12_u1(i,j),                            i = 1,im)
end do
do j = 1,jm
   write (31) (G_cont_12_u2(i,j),                            i = 1,im)
end do

! Output theoretical initial values
do j = 1,jm
   write (31) (u_init_u1(i,j),                               i = 1,im)
end do
do j = 1,jm
   write (31) (u_init_u2(i,j),                               i = 1,im)
end do
do j = 1,jm
   write (31) (v_init_u1(i,j),                               i = 1,im)
end do
do j = 1,jm
   write (31) (v_init_u2(i,j),                               i = 1,im)
end do


end subroutine initial_output



subroutine output_data(stp_cnt,tau)

implicit none


!----------------------------------------------------------------------------
! INTENT IN 
!----------------------------------------------------------------------------
integer (kind = int_kind) :: stp_cnt         ! Number of time step
real (kind = dbl_kind), intent(in) :: tau    ! Time in hours

!----------------------------------------------------------------------------
! LOCAL
!----------------------------------------------------------------------------
integer :: i,j
real (kind = dbl_kind), dimension(im,jm) ::            &
     u, v      ! Diagnosed Cartesian-coordinate velocities



! Diagnose Cartesian-coordinate velocities for output
u(:,:) = dx1_dx_u1(:,:)*u1_cov(:,:,n4) + dx2_dx_u1(:,:)*u2_cov_u1(:,:)
v(:,:) = dx1_dy_u2(:,:)*u1_cov_u2(:,:) + dx2_dy_u2(:,:)*u2_cov(:,:,n4)


do j = 1, jm
   write (31) (u1_cov(i,j,n4),      i = 1,im)
   write (31) (u2_cov(i,j,n4),      i = 1,im)
   write (31) (u1_cont(i,j),        i = 1,im)
   write (31) (u2_cont(i,j),        i = 1,im)
   write (31) (u(i,j),              i = 1,im)
   write (31) (v(i,j),              i = 1,im)
   write (31) (zeta(i,j),           i = 1,im)
   write (31) (pv(i,j),             i = 1,im)
   write (31) (h_star(i,j),         i = 1,im)
   write (31) (h(i,j),              i = 1,im)
end do





call calc_gmeans(tau)

print *, "Output data has been written."




end subroutine output_data





subroutine calc_gmeans (tau)
! Calculates and outputs global mass-weighted means of various variables

implicit none

!----------------------------------------------------------------------------
! INTENT IN 
!----------------------------------------------------------------------------
real (kind = dbl_kind), intent(in) :: tau    ! Time in hours

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: i,j,k

real (kind = dbl_kind) ::                                            &
      sqrt_G_h_star_bar, & ! area-weighted mean mass (m)
      s_G_h_star_temp,   &
      ke_bar,            & ! mass-weighted mean kinetic energy (J/kg)
      geop_bar,          & ! mass-weighted mean geopotential energy (J/kg)
      zeta_bar,          & ! mass-weighted mean relative vorticity (s^-1)
      pv_bar,            & ! mass-weighted mean potential vort. (Pa^-1 s^-1)
      pot_enstr_bar,     & ! mass-weighted mean potential enstrophy (Pa^-2 s^-2)
      total_energy_bar     ! mass-weighted mean total energy (J/kg)


! Initialize
sqrt_G_h_star_bar = c0
ke_bar = c0
geop_bar = c0
zeta_bar = c0
pv_bar = c0
pot_enstr_bar = c0
total_energy_bar = c0


! Calculate area-weighted mass
do j = 1,jm
   do i = 1,im
      sqrt_G_h_star_bar = sqrt_G_h_star_bar + sqrt_G_h_star(i,j,n4)
   end do
end do
sqrt_G_h_star_bar = sqrt_G_h_star_bar / (im*jm)



! Calculate global mass-weighted kinetic energy, geopotential energy,
! potential vorticity and enstrophy
do j = 1,jm
   do i = 1,im
      ke_bar = ke_bar + sqrt_G_h_star(i,j,n4)*ke_horiz(i,j)
      geop_bar = geop_bar + sqrt_G_h_star(i,j,n4)*grav*  &
         ( p5*h_star(i,j)+hs(i,j) )
      s_G_h_star_temp = p25*( sqrt_G_h_star(i,j,n4) +                   &
         sqrt_G_h_star(im1(i),j,n4) + sqrt_G_h_star(im1(i),jm1(j),n4) + &
            sqrt_G_h_star(i,jm1(j),n4) )
      zeta_bar = zeta_bar + s_G_h_star_temp*zeta(i,j)
      pv_bar = pv_bar + s_G_h_star_temp*pv(i,j)
      pot_enstr_bar = pot_enstr_bar + s_G_h_star_temp*p5*pv(i,j)**2
   end do
end do
ke_bar = ke_bar / ( sqrt_G_h_star_bar*im*jm )
geop_bar = geop_bar / ( sqrt_G_h_star_bar*im*jm )
zeta_bar = zeta_bar / ( sqrt_G_h_star_bar*im*jm )
pv_bar = pv_bar / ( sqrt_G_h_star_bar*im*jm )
pot_enstr_bar = pot_enstr_bar / ( sqrt_G_h_star_bar*im*jm )


! Calculate global mass-weighted mean total energy (i.e. geopotential
! and kinetic energy)
total_energy_bar = geop_bar + ke_bar




! Write global means to file
write (45, "(F24.10,F24.6,7(ES24.14E3))" )                                 &
             tau, tau*3600._dbl_kind, sqrt_G_h_star_bar, geop_bar, ke_bar, &
             total_energy_bar, zeta_bar, pv_bar, pot_enstr_bar



print "(A9,F24.13,A26,F24.13,A19,F24.13)", "ke_bar =", ke_bar,         &
      "total_energy_bar =", total_energy_bar, "h_star_bar =", sqrt_G_h_star_bar

end subroutine calc_gmeans





end module output
