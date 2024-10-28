;************************************************************************
; Plot of coordinate systems
;************************************************************************
;	Include NCL scripts libraries				 
;************************************************************************
;************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
; 	Main Program					 
;************************************************************************
;************************************************************************
begin

  fname =  "sh_water_model.1-layer.nc"
  ff      = addfile(fname,"r")       ; Open netCDF files.	
  
 ;output_format = "X11"
 output_format = "pdf"


  wks = gsn_open_wks(output_format,"zeta_plot_x1_x2.non_orth.offset_coord")      ; open workstation
  print ("Writing file Ps"+"."+output_format)


;;;;;;;;
; Load variables
  zeta = ff->zeta
  
  

; gsn_define_colormap(wks,"Rainbow")
; gsn_define_colormap(wks,"gui_default")
gsn_define_colormap(wks,"WhViBlGrYeOrReWh")


;;;;;
;  Plot settings for zeta

    res                     = True
    res@gsnDraw             = False        ; don't draw yet
    res@gsnFrame            = False        ; don't advance frame yet
;   res@gsnMaximize         = True         ; Maximize size of plot in frame
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

;   res@trYLog = True
;   res@trYReverse               = True
    res@tmYLLabelFont = 21    ;  21 = helvetica
    res@tmXBLabelFont = 21    ;  22 = helvetica-bold
    ; res@tmYBMode = "Explicit"
;   res@lbLabelBarOn = True
;   res@lbLabelFont = 21
    res@lbLabelStride       = 1       ;
    res@lbTitleFontThicknessF   = .5
    res@lbLabelFontHeightF    = .016
    res@tiYAxisString = "x2 (m)"
    res@tiXAxisString = "x1 (m)"
    res@tiXAxisOffsetYF = 0.19

    ;res@vpXF      = 0.10   ; Change the size and location of the
    ;res@vpYF      = 0.85   ; plot on the viewport.
    ;res@vpWidthF  = .80          ; choose width of plot
    ;res@vpHeightF = .50          ; choose height of plot





;;;;;
; Create plot


    ;arad = 6371220.
    ;pi = 3.14159265358979323846
    ;res@trXMinF =  45. *(pi/180.)*arad
    ;res@trXMaxF = 180. *(pi/180.)*arad
    ;res@trYMinF = -30. *(pi/180.)*arad
    ;res@trYMaxF = 150.*(pi/180.)*arad


    minv = -3.0E-05
    maxv =  3.0E-05
    res@cnMinLevelValF = minv           
    res@cnMaxLevelValF = maxv          
    res@cnLevelSpacingF = 0.5E-05
     


    contour = gsn_csm_contour(wks,zeta(56,:,:),res)

    draw(contour)
    frame(wks)


end

