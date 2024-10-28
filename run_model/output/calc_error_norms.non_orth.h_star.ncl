load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"


a=addfile("converged_file.nc","r")
b=addfile("sh_water_model.1-layer.non_orth.Williamson_TC_5.dt_1s.dx_400km.new.nc","r")


xa = a->x1_h
ya = a->x2_h

xb = b->x1_h
yb = b->x2_h

dimsa = getfilevardimsizes(a,"h_star")
ntma = dimsa(0)

h_star_conv = a->h_star(ntma-1,:,:)

h_star_conv_interp = linint2_Wrap(xa,ya,h_star_conv,True,xb,yb,0)


dimsb = getfilevardimsizes(b,"h_star")
ntmb = dimsb(0)


h_star = b->h_star(ntmb-1,:,:)



outfile = addfile("h_stars_for_error_norm_calc.non_orth.dx_400km.new.nc","c")
outfile ->h_star_conv_interp=h_star_conv_interp
outfile ->h_star_400km=h_star

