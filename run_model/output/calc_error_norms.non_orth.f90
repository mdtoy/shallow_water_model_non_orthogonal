program calc_error_norms

!--------------------------------------------------------------------
! This program calculates L1, L2 and Linf norms of input variable.
! Data is read from netCDF files.
!--------------------------------------------------------------------



implicit none

include 'netcdf.inc'




integer, parameter :: dbl_kind = selected_real_kind(13)



integer :: i
integer :: im,jm

integer :: ncid_in,err
integer :: varid
integer :: dimid


real (kind = dbl_kind), allocatable ::                               &
           zeta(:), zeta_conv(:)

real (kind = dbl_kind) ::                                            &
           L1, L2, Linf, error, L1_denom, L2_denom, Linf_denom



! read in array limits
err = NF_OPEN("zetas_for_error_norm_calc.non_orth.dx_200km.nc",NF_NOWRITE,ncid_in)
err = NF_INQ_DIMID(ncid_in,'x1_u1',dimid)
err = NF_INQ_DIMLEN(ncid_in,dimid,im)
err = NF_INQ_DIMID(ncid_in,'x2_u2',dimid)
err = NF_INQ_DIMLEN(ncid_in,dimid,jm)




print *
print *, "   im =", im, "   jm =", jm
print *



! read in variables
allocate (zeta(im*jm))
allocate (zeta_conv(im*jm))
err = NF_INQ_VARID(ncid_in,'zeta_200km',varid)
err = NF_GET_VAR_DOUBLE(ncid_in,varid,zeta)
err = NF_INQ_VARID(ncid_in,'zeta_conv_interp',varid)
err = NF_GET_VAR_DOUBLE(ncid_in,varid,zeta_conv)





! Calculate error norms (Williamson et al. (1992))
L1   = 0._dbl_kind
L2   = 0._dbl_kind
Linf = 0._dbl_kind
L1_denom  = 0._dbl_kind
L2_denom  = 0._dbl_kind
Linf_denom = 0._dbl_kind
do i =  1, im*jm
   error = zeta(i) - zeta_conv(i)
   L1 = L1 + abs(error)
   L2 = L2 + error**2
   Linf = max(Linf,abs(error))
   L1_denom = L1_denom + abs(zeta_conv(i))
   L2_denom = L2_denom + zeta_conv(i)**2
   Linf_denom = max(Linf_denom,abs(zeta_conv(i)))
end do
L1 = L1 / L1_denom
L2 = L2**0.5 / L2_denom**0.5
Linf = Linf / Linf_denom

print *
print *, "L1 =", L1
print *, "L2 =", L2
print *, "Linf =", Linf
print *


err = NF_CLOSE(ncid_in)


end program calc_error_norms
