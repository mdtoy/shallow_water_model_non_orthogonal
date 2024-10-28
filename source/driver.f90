!=====================================================================
! PROGRAM SHALLOW_WATER_1_LAYER
!=====================================================================

program SHALLOW_WATER_1_LAYER

!---------------------------------------------------------------------
! PURPOSE:  Driver for a typical one-layer shallow water model.  
!           Adams-Bashforth 3rd order time stepping scheme is used.
!
! AUTHOR:   Michael Toy,  January 2012
!---------------------------------------------------------------------


use kinds
use model_parameters
use physical_parameters
use initialize
use step
use output


implicit none


integer (kind=int_kind) :: step_count

real (kind=dbl_kind) ::   &
                tau,      &           ! time in hours
                tau_end               ! run termination time (hours)
                      

!
! Initialize step_count
!
step_count = 0


!
! Initialize the model (i.e. prognostic variables, initial and
!                     boundary conditions, eta_coordinate and time)
!
call initialize_model(step_count,tau)


print "(I10,A10,F18.10,A6)", step_count, "Time =", tau, "hours"


!
! Calculate run termination time
!
tau_end = tau + tau_duration


!
! Open output files to be written in subroutine output_data.
!
open (unit = 31, file = "./output/sh_water_model.1-layer.out",          &
                                             form = "unformatted")
open (unit = 45, file = "./output/gmeans.out")


call initial_output (tau)

write (45, "(A16,8(A24))" ) "Tau(hr)", "Tau(sec)", "h_star",            &
                     "geop", "ke", "total_energy", "zeta", "pv", "pot_enstrophy"



!
!  Take first timestep (Euler forward)
!
call update_diagnostics(tau,w1_ef,w2_ef,w3_ef)
call output_data(step_count,tau)        ! Output initial conditions
step_count = step_count + 1
tau = tau + dt/3600._dbl_kind
call step_dynamics(step_count,w1_ef,w2_ef,w3_ef)
print "(I10,A10,F18.10,A6)", step_count, "Time =", tau, "hours"
call update_diagnostics(tau,w1_ab2,w2_ab2,w3_ab2)
if (mod(step_count,out_freq)==0) then
   call output_data(step_count,tau)
end if

!
!  Take second (AB2 forward) timestep
!
step_count = step_count + 1
tau = tau + dt/3600._dbl_kind
call step_dynamics(step_count,w1_ab2,w2_ab2,w3_ab2)
print "(I10,A10,F18.10,A6)", step_count, "Time =", tau, "hours"
call update_diagnostics(tau,w1_ab3,w2_ab3,w3_ab3)
if (mod(step_count,out_freq)==0) then
   call output_data(step_count,tau)
end if


do while (tau .lt. tau_end)
   step_count = step_count + 1
   tau = tau + dt/3600._dbl_kind
   call step_dynamics(step_count,w1_ab3,w2_ab3,w3_ab3)
   print "(I10,A10,F18.10,A6)", step_count, "Time =", tau, "hours"
   call update_diagnostics(tau,w1_ab3,w2_ab3,w3_ab3)
   if (mod(step_count,out_freq)==0) then
      call output_data(step_count,tau)
   end if
end do


print *
print *, "*******  END PROGRAM  ******* "
print *


close (31)
close (45)


end program SHALLOW_WATER_1_LAYER

!=====================================================================
! END OF SHALLOW_WATER_1_LAYER
!=====================================================================
