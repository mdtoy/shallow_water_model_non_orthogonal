load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"


; Parameters
a_rad = 6371220.
pi = 3.14159265358979323846
C_fac = 0.5
D_fac = 1.0
alf   = 0.
bet   = 0.


a=addfile("converged_file.nc","r")
b=addfile("sh_water_model.1-layer.non_orth.Williamson_TC_5.dt_1s.dx_25km.new.nc","r")


; read in true cartesian coordinate values
xcart = a->x_u
ycart = a->y_v


dimsa = getfilevardimsizes(a,"zeta")
ntma = dimsa(0)
im = dimsa(1)
jm = dimsa(2)

dimsb = getfilevardimsizes(b,"zeta")
ntmb = dimsb(0)

zeta_conv = a->zeta(ntma-1,:,:)


; Create 2x2 matrix of non-orthogonal coordinate values on 
; converged Cartesian grid
x1_conv = new( (/jm,im/), double )
x2_conv = new( (/jm,im/), double )
do i = 0,im-1
   do j = 0,jm-1
      x1_conv(j,i) = xcart(i) + a_rad*C_fac*sin(ycart(j)/a_rad-alf)
      x2_conv(j,i) = ycart(j) + a_rad*D_fac*sin(xcart(i)/a_rad-bet)
   end do
end do


; Rearrange x1 values to constrain values within regular domain,
; i.e., within -pi and pi times a_rad
do i = 0,im-1
   do j = 0,jm-1
      if ( abs(x1_conv(j,i))/a_rad .gt. pi ) then
         x1_conv(j,i) = x1_conv(j,i) - (abs(x1_conv(j,i))/x1_conv(j,i))*2*pi*a_rad
      end if
   end do
end do

; Sort resulting x1 dimension
ii_rgt = dim_pqsort_n(x1_conv,2,1)    ; permutation vector
; Sort zeta_conv and x2_conv accordingly
do j = 0,jm-1
   zeta_conv(j,:) = zeta_conv(j,ii_rgt(j,:))
   x2_conv(j,:) = x2_conv(j,ii_rgt(j,:))
end do
delete(ii_rgt)


; Rearrange x2 values to constrain values within regular domain,
; i.e., within -pi and pi times a_rad
do i = 0,im-1
   do j = 0,jm-1
      if ( abs(x2_conv(j,i))/a_rad .gt. pi ) then
         x2_conv(j,i) = x2_conv(j,i) - (abs(x2_conv(j,i))/x2_conv(j,i))*2*pi*a_rad
      end if
   end do
end do

; Sort resulting x2 dimension
ii_lft = dim_pqsort_n(x2_conv,2,0)    ; permutation vector
; Sort zeta_conv and x1_conv accordingly
do i = 0,im-1
   zeta_conv(:,i) = zeta_conv(ii_lft(:,i),i)
   x1_conv(:,i) = x1_conv(ii_lft(:,i),i)
end do
delete(ii_lft)


; Read in zeta from non-orthog run
zeta = b->zeta(ntmb-1,:,:)

; Form x1-x2 1D vectors for interpolation
x1_vec = zeta&x1_u1
x2_vec = zeta&x2_u2


; Interpolate Cartesian (converged) field onto non-orthogonal grid
zeta_conv_interp = rcm2rgrid_Wrap(x2_conv,x1_conv,zeta_conv,x2_vec,x1_vec,1)



outfile = addfile("zetas_for_error_norm_calc.non_orth.cartesian_conv.dx_25km.new.nc","c")
outfile ->zeta_conv_interp=zeta_conv_interp
outfile ->zeta_25km=zeta

