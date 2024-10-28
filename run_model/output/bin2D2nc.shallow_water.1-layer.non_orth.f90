program bin2D2nc_sh_water_1_layer

!--------------------------------------------------------------------
! This program reads the surface pressure field from the binary
! output file surface_pressure.out and writes it to a NetCDF file.
!--------------------------------------------------------------------



implicit none

include 'netcdf.inc'




integer, parameter :: dbl_kind = selected_real_kind(13)

integer :: nt,i,j
integer :: ntm,im,jm

integer :: ncid,err,x1id,x2id,y1id,y2id,timeid,varid
integer, dimension (3) :: dimids_h,dimids_u1,dimids_u2,dimids_q
integer, dimension (3) :: start_vector, count_vector
integer :: count

real (kind = dbl_kind), allocatable :: time_vector(:)
real (kind = dbl_kind), allocatable :: x_vector_1(:)
real (kind = dbl_kind), allocatable :: x_vector_2(:)
real (kind = dbl_kind), allocatable :: y_vector_1(:)
real (kind = dbl_kind), allocatable :: y_vector_2(:)

real (kind = dbl_kind), allocatable ::                               &
           fld1_2d(:), fld2_2d(:), fld3_2d(:), fld4_2d(:),           &
           fld5_2d(:), fld6_2d(:), fld7_2d(:), fld8_2d(:),           &
           fld9_2d(:), fld10_2d(:)



open (unit = 20, file = "./sh_water_model.1-layer.out", action = "read", &
         form = "unformatted")


read (20) ntm, im, jm


print *
print *, "ntm = ", ntm, "   im = ", im, "   jm = ", jm
print *

allocate (time_vector(ntm))
allocate (x_vector_1(im))
allocate (x_vector_2(im))
allocate (y_vector_1(jm))
allocate (y_vector_2(jm))


read (20) (time_vector(nt),     nt = 1,ntm)
read (20) (x_vector_1(i),       i = 1,im  )
read (20) (x_vector_2(i),       i = 1,im  )
read (20) (y_vector_1(j),       j = 1,jm  )
read (20) (y_vector_2(j),       j = 1,jm  )



! Create and open the netCDF file
err = NF_CREATE("sh_water_model.1-layer.non_orth.nc", NF_CLOBBER, ncid)

err = NF_REDEF(ncid)

! Define the three dimensions (time and 2D space)
err = NF_DEF_DIM(ncid, 'time', ntm, timeid)
err = NF_DEF_DIM(ncid, 'x1_h', im, x1id)
err = NF_DEF_DIM(ncid, 'x1_u1', im, x2id)
err = NF_DEF_DIM(ncid, 'x2_h', jm, y1id)
err = NF_DEF_DIM(ncid, 'x2_u2', jm, y2id)

! Define the "measurement" variables of the three dimensions and
! attach attributes to them.
err = NF_DEF_VAR(ncid, 'time', NF_DOUBLE, 1, timeid, varid)
err = NF_PUT_ATT_TEXT(ncid,varid,'units',5,'hours')
err = NF_PUT_ATT_TEXT(ncid,varid,'long_name',4,'time')
err = NF_PUT_ATT_TEXT(ncid,varid,'axis',1,'T')

err = NF_DEF_VAR(ncid, 'x1_h', NF_DOUBLE, 1, x1id, varid)
err = NF_PUT_ATT_TEXT(ncid,varid,'units',1,'m')
err = NF_PUT_ATT_TEXT(ncid,varid,'axis',1,'X')

err = NF_DEF_VAR(ncid, 'x1_u1', NF_DOUBLE, 1, x2id, varid)
err = NF_PUT_ATT_TEXT(ncid,varid,'units',1,'m')
err = NF_PUT_ATT_TEXT(ncid,varid,'axis',1,'X')

err = NF_DEF_VAR(ncid, 'x2_h', NF_DOUBLE, 1, y1id, varid)
err = NF_PUT_ATT_TEXT(ncid,varid,'units',1,'m')
err = NF_PUT_ATT_TEXT(ncid,varid,'axis',1,'Y')

err = NF_DEF_VAR(ncid, 'x2_u2', NF_DOUBLE, 1, y2id, varid)
err = NF_PUT_ATT_TEXT(ncid,varid,'units',1,'m')
err = NF_PUT_ATT_TEXT(ncid,varid,'axis',1,'Y')

err = NF_ENDDEF(ncid)


! Write the "measurement" dimensional variables to the netCDF file
err = NF_INQ_VARID(ncid, 'time', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, time_vector)

err = NF_INQ_VARID(ncid, 'x1_h', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, x_vector_1)

err = NF_INQ_VARID(ncid, 'x1_u1', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, x_vector_2)

err = NF_INQ_VARID(ncid, 'x2_h', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, y_vector_1)

err = NF_INQ_VARID(ncid, 'x2_u2', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, y_vector_2)


! Read in 'z_surf' as the vector 'fld1_2d'
allocate (fld1_2d(im*jm))
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do


! Define the 'dimensions vector' dimids to be used for writing
! the 2- and 3-dimensional variables to the netCDF file
dimids_h(1) = x1id
dimids_h(2) = y1id
dimids_h(3) = timeid

dimids_u1(1) = x2id
dimids_u1(2) = y1id
dimids_u1(3) = timeid

dimids_u2(1) = x1id
dimids_u2(2) = y2id
dimids_u2(3) = timeid

dimids_q(1) = x2id
dimids_q(2) = y2id
dimids_q(3) = timeid


! Define 'z_surf' in the netCDF file
err = NF_REDEF(ncid)

err = NF_DEF_VAR(ncid, 'z_surf', NF_DOUBLE, 2, dimids_h, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, 'm')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 14, 'surface height')

err = NF_ENDDEF(ncid)

! Write values of 'z_surf' to netCDF file
err = NF_INQ_VARID(ncid, 'z_surf', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)



! Write coordinate data

! x_h as function of X1,X2
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'x_h', NF_DOUBLE, 2, dimids_h, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, 'm')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 27, 'Cartesian x values at h-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'x_h', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! y_h as function of X1,X2
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'y_h', NF_DOUBLE, 2, dimids_h, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, 'm')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 27, 'Cartesian y values at h-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'y_h', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! x_u as function of X1,X2
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'x_u', NF_DOUBLE, 2, dimids_u1, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, 'm')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 28, 'Cartesian x values at u1-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'x_u', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! y_u as function of X1,X2
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'y_u', NF_DOUBLE, 2, dimids_u1, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, 'm')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 28, 'Cartesian y values at u1-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'y_u', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! x_v as function of X1,X2
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'x_v', NF_DOUBLE, 2, dimids_u2, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, 'm')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 28, 'Cartesian x values at u2-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'x_v', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! y_v as function of X1,X2
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'y_v', NF_DOUBLE, 2, dimids_u2, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, 'm')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 28, 'Cartesian y values at u2-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'y_v', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! x_q as function of X1,X2
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'x_q', NF_DOUBLE, 2, dimids_q, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, 'm')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 27, 'Cartesian x values at q-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'x_q', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! y_q as function of X1,X2
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'y_q', NF_DOUBLE, 2, dimids_q, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, 'm')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 27, 'Cartesian y values at q-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'y_q', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! x1_h as function of X,Y
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'x1_coord_h', NF_DOUBLE, 2, dimids_h, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, 'm')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 29, 'Non-orthog x1 values at h-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'x1_coord_h', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! x2_h as function of X,Y
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'x2_coord_h', NF_DOUBLE, 2, dimids_h, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, 'm')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 29, 'Non-orthog x2 values at h-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'x2_coord_h', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)


! Metric terms

! sqrt_G
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'sqrt_G_h', NF_DOUBLE, 2, dimids_h, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, '-')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 20, 'Metric term at h-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'sqrt_G_h', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! G_cont_11_u1
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'G_cont_11_u1', NF_DOUBLE, 2, dimids_u1, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, '-')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 35, 'Contrav. metric tensor 11 at u1-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'G_cont_11_u1', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! G_cont_22_u2
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'G_cont_22_u2', NF_DOUBLE, 2, dimids_u2, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, '-')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 35, 'Contrav. metric tensor 12 at u2-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'G_cont_22_u2', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! G_cont_12_u1
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'G_cont_12_u1', NF_DOUBLE, 2, dimids_u1, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, '-')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 35, 'Contrav. metric tensor 12 at u1-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'G_cont_12_u1', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! G_cont_12_u2
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'G_cont_12_u2', NF_DOUBLE, 2, dimids_u2, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, '-')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 35, 'Contrav. metric tensor 12 at u2-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'G_cont_12_u2', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! u_init_u1
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'u_init_u1', NF_DOUBLE, 2, dimids_u1, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 5, 'm s-1')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 35, 'Ideal initial Cartesian u at u1-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'u_init_u1', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! u_init_u2
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'u_init_u2', NF_DOUBLE, 2, dimids_u2, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 5, 'm s-1')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 35, 'Ideal initial Cartesian u at u2-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'u_init_u2', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! v_init_u1
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'v_init_u1', NF_DOUBLE, 2, dimids_u1, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 5, 'm s-1')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 35, 'Ideal initial Cartesian v at u1-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'v_init_u1', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)

! v_init_u2
count = 0
do j = 1,jm
   read (20) ( fld1_2d(count+i),        i = 1,im )
   count = count + im
end do
err = NF_REDEF(ncid)
err = NF_DEF_VAR(ncid, 'v_init_u2', NF_DOUBLE, 2, dimids_u2, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 5, 'm s-1')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 35, 'Ideal initial Cartesian v at u2-pts')
err = NF_ENDDEF(ncid)
err = NF_INQ_VARID(ncid, 'v_init_u2', varid)
err = NF_PUT_VAR_DOUBLE(ncid, varid, fld1_2d)


deallocate (fld1_2d)



! Define variables in the netCDF file
err = NF_REDEF(ncid)

err = NF_DEF_VAR(ncid, 'u1_cov', NF_DOUBLE, 3, dimids_u1, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 5, 'm s-1')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 34, 'x1-component of covariant velocity')

err = NF_DEF_VAR(ncid, 'u2_cov', NF_DOUBLE, 3, dimids_u2, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 5, 'm s-1')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 34, 'x2-component of covariant velocity')

err = NF_DEF_VAR(ncid, 'u1_cont', NF_DOUBLE, 3, dimids_u1, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 5, 'm s-1')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 38, 'x1-component of contravariant velocity')

err = NF_DEF_VAR(ncid, 'u2_cont', NF_DOUBLE, 3, dimids_u2, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 5, 'm s-1')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 38, 'x2-component of contravariant velocity')

err = NF_DEF_VAR(ncid, 'u', NF_DOUBLE, 3, dimids_u1, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 5, 'm s-1')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 23, 'x-component of velocity')

err = NF_DEF_VAR(ncid, 'v', NF_DOUBLE, 3, dimids_u2, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 5, 'm s-1')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 23, 'y-component of velocity')

err = NF_DEF_VAR(ncid, 'zeta', NF_DOUBLE, 3, dimids_q, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 3, 's-1')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 18, 'relative vorticity')

err = NF_DEF_VAR(ncid, 'pv', NF_DOUBLE, 3, dimids_q, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 7, 'm-1 s-1')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 19, 'potential vorticity')

err = NF_DEF_VAR(ncid, 'h_star', NF_DOUBLE, 3, dimids_h, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, 'm')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 22, 'column fluid thickness')

err = NF_DEF_VAR(ncid, 'h', NF_DOUBLE, 3, dimids_h, varid)
err = NF_PUT_ATT_TEXT(ncid, varid, 'units', 1, 'm')
err = NF_PUT_ATT_TEXT(ncid, varid, 'long_name', 22, 'height of free surface')

err = NF_ENDDEF(ncid)


! Allocate space for the variable fields
allocate (fld1_2d(im*jm))
allocate (fld2_2d(im*jm))
allocate (fld3_2d(im*jm))
allocate (fld4_2d(im*jm))
allocate (fld5_2d(im*jm))
allocate (fld6_2d(im*jm))
allocate (fld7_2d(im*jm))
allocate (fld8_2d(im*jm))
allocate (fld9_2d(im*jm))
allocate (fld10_2d(im*jm))


! Define count_vector to be used in the NF_PUT_VARA_DOUBLE function
count_vector(1) = im
count_vector(2) = jm
count_vector(3) = 1

! Loop through time and read in fields as the vectors 'fldx_2d'
! Write these fields to the netCDF file
do nt = 1, ntm
   count = 0
   do j = 1,jm
      read (20) ( fld1_2d(count+i),   i = 1,im )
      read (20) ( fld2_2d(count+i),   i = 1,im )
      read (20) ( fld3_2d(count+i),   i = 1,im )
      read (20) ( fld4_2d(count+i),   i = 1,im )
      read (20) ( fld5_2d(count+i),   i = 1,im )
      read (20) ( fld6_2d(count+i),   i = 1,im )
      read (20) ( fld7_2d(count+i),   i = 1,im )
      read (20) ( fld8_2d(count+i),   i = 1,im )
      read (20) ( fld9_2d(count+i),   i = 1,im )
      read (20) ( fld10_2d(count+i),  i = 1,im )
      count = count + im
   end do

   ! Define start_vector to be used in the  NF_PUT_VARA_DOUBLE function
   start_vector(1) = 1
   start_vector(2) = 1
   start_vector(3) = nt

   err = NF_INQ_VARID(ncid, 'u1_cov', varid)
   err = NF_PUT_VARA_DOUBLE(ncid, varid, start_vector, count_vector, &
                            fld1_2d)

   err = NF_INQ_VARID(ncid, 'u2_cov', varid)
   err = NF_PUT_VARA_DOUBLE(ncid, varid, start_vector, count_vector, &
                            fld2_2d)

   err = NF_INQ_VARID(ncid, 'u1_cont', varid)
   err = NF_PUT_VARA_DOUBLE(ncid, varid, start_vector, count_vector, &
                            fld3_2d)

   err = NF_INQ_VARID(ncid, 'u2_cont', varid)
   err = NF_PUT_VARA_DOUBLE(ncid, varid, start_vector, count_vector, &
                            fld4_2d)


   err = NF_INQ_VARID(ncid, 'u', varid)
   err = NF_PUT_VARA_DOUBLE(ncid, varid, start_vector, count_vector, &
                            fld5_2d)

   err = NF_INQ_VARID(ncid, 'v', varid)
   err = NF_PUT_VARA_DOUBLE(ncid, varid, start_vector, count_vector, &
                            fld6_2d)

   err = NF_INQ_VARID(ncid, 'zeta', varid)
   err = NF_PUT_VARA_DOUBLE(ncid, varid, start_vector, count_vector, &
                            fld7_2d)

   err = NF_INQ_VARID(ncid, 'pv', varid)
   err = NF_PUT_VARA_DOUBLE(ncid, varid, start_vector, count_vector, &
                            fld8_2d)

   err = NF_INQ_VARID(ncid, 'h_star', varid)
   err = NF_PUT_VARA_DOUBLE(ncid, varid, start_vector, count_vector, &
                            fld9_2d)

   err = NF_INQ_VARID(ncid, 'h', varid)
   err = NF_PUT_VARA_DOUBLE(ncid, varid, start_vector, count_vector, &
                            fld10_2d)

end do


close (20)


err = NF_CLOSE(ncid)




end program bin2D2nc_sh_water_1_layer