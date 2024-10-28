load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"



begin

; Parameters
a_rad = 6371220.    ; Earth radius
pi = 3.14159265358979323846


a=addfile("zeta.non_orth.extracted.nc","r")
dimsa = getfilevardimsizes(a,"zeta")
ntm = dimsa(0)
jm = dimsa(1)
im = dimsa(2)

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


print("Sorting in x-dimension....")
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

print("Sorting in y-dimension....")
; Sort resulting y dimension
ii_lft = dim_pqsort_n(y,2,0)    ; permutation vector
; Sort zeta and x accordingly
do i = 0,im-1
   zeta(:,i) = zeta(ii_lft(:,i),i)
   x(:,i) = x(ii_lft(:,i),i)
end do
delete(ii_lft)


; Convert to 1D for natgrid
y1d = ndtooned(y)
x1d = ndtooned(x)
zeta1d = ndtooned(zeta)


print("Interpolating onto Cartesian grid....")
zeta_xy1d = natgrid(y1d,x1d,zeta1d,yvec,xvec)

; Convert zeta_xy back into 2D array
zeta_xy = onedtond(zeta_xy1d,(/jm,im/))

zeta_xy!0 = "y"
zeta_xy!1 = "x"
zeta_xy&y = yvec
zeta_xy&x = xvec


outfile = addfile("zeta_xy.troubleshoot.nc","c")
outfile ->zeta_xy=zeta_xy



;;;;;
; Plot results (zeta_xy)
 output_format = "X11"
 
 wks = gsn_open_wks(output_format,"zeta_plot_x_y_interp.troubleshoot")
 plot = new(2,graphic)
 
 gsn_define_colormap(wks,"WhViBlGrYeOrReWh")


;;;;;
;  Plot-settings for zeta

    res                     = True
    res@gsnSpreadColors     = True         ; Use full colormap, but start
    res@gsnSpreadColorStart =  10           ; at color index 24.
    res@gsnSpreadColorEnd   = 100


    res@cnFillOn       = True     ; Turn on contour fill
    res@cnLinesOn      = True    ; Turn off contour lines
    res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                ; use "CellFill" and "RasterFill"
    res@cnLineLabelsOn      = False        ; Turn off contour line labels
    res@cnLevelSelectionMode = "ManualLevels"
    res@lbLabelAutoStride   = True 

    res@tmYLLabelFont = 21    ;  21 = helvetica
    res@tmXBLabelFont = 21    ;  22 = helvetica-bold
    res@lbLabelStride       = 1       ;
    res@lbTitleFontThicknessF   = .5
    res@lbLabelFontHeightF    = .016
    res@tiYAxisString = "y (m)"
    res@tiXAxisString = "x (m)"
    res@tiXAxisOffsetYF = 0.19

;;;;;
; Create plot

    minv = -3.0E-05
    maxv =  3.0E-05
    res@cnMinLevelValF = minv           
    res@cnMaxLevelValF = maxv          
    res@cnLevelSpacingF = 0.5E-05
     

    print("Plotting results....")
    plot(0) = gsn_csm_contour(wks,zeta(:,:),res)
    plot(1) = gsn_csm_contour(wks,zeta_xy(:,:),res)

    gsn_panel(wks,plot,(/2,1/),False)

end