load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"


; Parameters
a_rad = 6371220.
pi = 3.14159265358979323846


a=addfile("sh_water_model.1-layer.non_orth.Williamson_TC_5.dt_30s.new.nc","r")
dimsa = getfilevardimsizes(a,"zeta")
ntm = dimsa(0)
im = dimsa(1)
jm = dimsa(2)

; read in cartesian coordinate values
x = a->x_q
y = a->y_q


zeta = a->zeta(ntm-1,:,:)     ; last timestep only

; Form x-y 1D vectors for interpolation
xvec = zeta&x1_u1
yvec = zeta&x2_u2


; Rearrange x values to constrain values within Cartesian domain,
;  i.e., within -pi and pi times a_rad
do i = 0,im-1
   do j = 0,jm-1
      if ( abs(x(j,i))/a_rad .gt. pi ) then
         x(j,i) = x(j,i) - (abs(x(j,i))/x(j,i))*2*pi*a_rad
      end if
   end do
end do

; Sort resulting x dimension
ii_rgt = dim_pqsort_n(x,2,1)    ; permutation vector
; Sort zeta and y accordingly
do j = 0,jm-1
   zeta(j,:) = zeta(j,ii_rgt(j,:))
   y(j,:) = y(j,ii_rgt(j,:))
end do
delete(ii_rgt)

; Rearrange y values to constrain values within Cartesian domain,
;  i.e., within -pi and pi times a_rad
do i = 0,im-1
   do j = 0,jm-1
      if ( abs(y(j,i))/a_rad .gt. pi ) then
         y(j,i) = y(j,i) - (abs(y(j,i))/y(j,i))*2*pi*a_rad
      end if
   end do
end do

; Sort resulting y dimension
ii_lft = dim_pqsort_n(y,2,0)    ; permutation vector
; Sort zeta and x accordingly
do i = 0,im-1
   zeta(:,i) = zeta(ii_lft(:,i),i)
   x(:,i) = x(ii_lft(:,i),i)
end do
delete(ii_lft)


zeta_xy = rcm2rgrid_Wrap(y,x,zeta,yvec,xvec,1)

outfile = addfile("zeta_xy.new.nc","c")
outfile ->zeta_xy=zeta_xy

