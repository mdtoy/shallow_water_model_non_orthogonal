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

  fname =  "sh_water_model.1-layer.non_orth.nc"
  ff      = addfile(fname,"r")       ; Open netCDF files.	
  
 ;output_format = "X11"
 output_format = "pdf"


  wks = gsn_open_wks(output_format,"Coord_plot_x_y.non-orth.offset_coord")        ; open workstation
  print ("Writing file Ps"+"."+output_format)


;;;;;;;;
; Load variables
  x_h = ff->x_h
  y_h = ff->y_h
  
  


;   gsn_define_colormap(wks,"gui_default")
; gsn_define_colormap(wks,"WhViBlGrYeOrRe")


;;;;;
;  Plot settings for x1

    res                     = True
    res@gsnDraw             = False        ; don't draw yet
    res@gsnFrame            = False        ; don't advance frame yet
;   res@gsnMaximize         = True         ; Maximize size of plot in frame
    res@gsnSpreadColors     = True         ; Use full colormap, but start


    res@cnFillOn       = False     ; Turn on contour fill
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

    ;res@vpXF      = 0.10   ; Change the size and location of the
    ;res@vpYF      = 0.85   ; plot on the viewport.
    ;res@vpWidthF  = .80          ; choose width of plot
    ;res@vpHeightF = .50          ; choose height of plot


;;;;;
;  Plot settings for x2

    resx2                  = True
    resx2@gsnDraw          = False        ; don't draw yet
    resx2@gsnFrame         = False        ; don't advance frame yet

    resx2@cnLineLabelsOn      = False        ; Turn off contour line labels
    resx2@cnLevelSelectionMode = "ManualLevels"
    resx2@lbLabelAutoStride   = True 

    resx2@cnLineColor      = "Black"
    resx2@cnLineThicknessF = 1.
    resx2@cnFillOn            = False
    ;minv = 0.  
    ;maxv = 12000.
    ;resx2@cnMinLevelValF = minv           
    ;resx2@cnMaxLevelValF = maxv          
    resx2@cnLinesOn       = True
    resx2@tiMainOn  = False      ; turn main title off
    resx2@tiYAxisOn = False      ; turn y-axis title off
    ; Turn borders off
    resx2@tmXBBorderOn = False
    resx2@tmXTBorderOn = False
    resx2@tmYLBorderOn = False
    resx2@tmYRBorderOn = False
    resx2@tmXBLabelsOn = False    ; turn x labels off
    resx2@tmYLLabelsOn = False    ; turn y labels off
    ; turn tickmarks off
    resx2@tmYROn = False                                        
    resx2@tmYLOn = False 
    resx2@tmXBOn = False 
    resx2@tmXTOn = False 
    resx2@cnInfoLabelOn = False
    ; resx2@gsnContourZeroLineThicknessF = 2

    ;resx2@vpXF      = 0.10   ; Change the size and location of the
    ;resx2@vpYF      = 0.85   ; plot on the viewport.
    ;resx2@vpWidthF  = .80          ; choose width of plot
    ;resx2@vpHeightF = .50          ; choose height of plot





;;;;;
; Create plot


    minv = -100000000 
    maxv =  100000000
    res@cnMinLevelValF = minv           
    res@cnMaxLevelValF = maxv          
    res@cnLevelSpacingF = 2000000
    
    resx2@cnMinLevelValF = minv           
    resx2@cnMaxLevelValF = maxv   
    resx2@cnLevelSpacingF = 2000000
     


    contourx1 = gsn_csm_contour(wks,x_h,res)
    contourx2 = gsn_csm_contour(wks,y_h,resx2)

    overlay(contourx1,contourx2)

    draw(contourx1)
    frame(wks)


end

