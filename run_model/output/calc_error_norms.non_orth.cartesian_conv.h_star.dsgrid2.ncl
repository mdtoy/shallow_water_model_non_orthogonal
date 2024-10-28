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
b=addfile("sh_water_model.1-layer.non_orth.TC_5_Gauss_mtn.dt_1s.dx_100km.nc","r")


; read in cartesian coordinate values
xcart = a->x_h
ycart = a->y_h


dimsa = getfilevardimsizes(a,"h_star")
ntma = dimsa(0)
jma = dimsa(1)
ima = dimsa(2)

dimsb = getfilevardimsizes(b,"h_star")
ntmb = dimsb(0)
jmb = dimsb(1)
imb = dimsb(2)

h_star_conv = a->h_star(ntma-1,:,:)


; Create 2x2 matrix of non-orthogonal coordinate values on 
; converged Cartesian grid
x1_conv = new( (/jma,ima/), double )
x2_conv = new( (/jma,ima/), double )
do i = 0,ima-1
   do j = 0,jma-1
      x1_conv(j,i) = xcart(i) + a_rad*C_fac*sin(ycart(j)/a_rad-alf)
      x2_conv(j,i) = ycart(j) + a_rad*D_fac*sin(xcart(i)/a_rad-bet)
   end do
end do


; Rearrange x1 values to constrain values within regular domain,
; i.e., within -pi and pi times a_rad
do i = 0,ima-1
   do j = 0,jma-1
      if ( abs(x1_conv(j,i))/a_rad .gt. pi ) then
         x1_conv(j,i) = x1_conv(j,i) - (abs(x1_conv(j,i))/x1_conv(j,i))*2*pi*a_rad
      end if
   end do
end do

; Sort resulting x1 dimension
ii_rgt = dim_pqsort_n(x1_conv,2,1)    ; permutation vector
; Sort h_star_conv and x2_conv accordingly
do j = 0,jma-1
   h_star_conv(j,:) = h_star_conv(j,ii_rgt(j,:))
   x2_conv(j,:) = x2_conv(j,ii_rgt(j,:))
end do
delete(ii_rgt)


; Rearrange x2 values to constrain values within regular domain,
; i.e., within -pi and pi times a_rad
do i = 0,ima-1
   do j = 0,jma-1
      if ( abs(x2_conv(j,i))/a_rad .gt. pi ) then
         x2_conv(j,i) = x2_conv(j,i) - (abs(x2_conv(j,i))/x2_conv(j,i))*2*pi*a_rad
      end if
   end do
end do

; Sort resulting x2 dimension
ii_lft = dim_pqsort_n(x2_conv,2,0)    ; permutation vector
; Sort h_star_conv and x1_conv accordingly
do i = 0,ima-1
   h_star_conv(:,i) = h_star_conv(ii_lft(:,i),i)
   x1_conv(:,i) = x1_conv(ii_lft(:,i),i)
end do
delete(ii_lft)


; Read in h_star from non-orthog run
h_star = b->h_star(ntmb-1,:,:)

; Form x1-x2 1D vectors for interpolation
x1_vec = h_star&x1_h
x2_vec = h_star&x2_h


; Convert to 1D for dsgrid2
x2_conv_1d = ndtooned(x2_conv)
x1_conv_1d = ndtooned(x1_conv)
h_star_conv_1d = ndtooned(h_star_conv)



; Interpolate Cartesian (converged) field onto non-orthogonal grid
print("Interpolating onto non-orthogonal grid....")
h_star_conv_interp_1d = dsgrid2(x2_conv_1d,x1_conv_1d,h_star_conv_1d,x2_vec,x1_vec)
print("....finished interpolating onto non-orthogonal grid....")


; Convert h_star_conv_interp back into 2D array
h_star_conv_interp = onedtond(h_star_conv_interp_1d,(/jmb,imb/))


; Add some dimensions
h_star_conv_interp!0 = "x2_h"
h_star_conv_interp!1 = "x1_h"
h_star_conv_interp&x2_h = x2_vec
h_star_conv_interp&x1_h = x1_vec


outfile = addfile("h_stars_for_error_norm_calc.non_orth.cartesian_conv.TC_5_Gauss_mtn.dx_100km.dsgrid2.nc","c")
outfile ->h_star_conv_interp=h_star_conv_interp
outfile ->h_star_100km=h_star

