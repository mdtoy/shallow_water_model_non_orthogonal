load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"


a=addfile("converged_file.nc","r")
b=addfile("sh_water_model.1-layer.non_orth.Williamson_TC_5.dt_1s.dx_200km.nc","r")


xa = a->x1_u1
ya = a->x2_u2

xb = b->x1_u1
yb = b->x2_u2

dimsa = getfilevardimsizes(a,"zeta")
ntma = dimsa(0)

zeta_conv = a->zeta(ntma-1,:,:)

zeta_conv_interp = linint2_Wrap(xa,ya,zeta_conv,True,xb,yb,0)


dimsb = getfilevardimsizes(b,"zeta")
ntmb = dimsb(0)


zeta = b->zeta(ntmb-1,:,:)



outfile = addfile("zetas_for_error_norm_calc.non_orth.dx_200km.nc","c")
outfile ->zeta_conv_interp=zeta_conv_interp
outfile ->zeta_200km=zeta

