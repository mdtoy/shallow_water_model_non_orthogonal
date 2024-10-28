load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"

a=addfile("sh_water_model.1-layer.non_orth.Williamson_TC_5.dt_30s.nc","r")

x = a->x_q
y = a->y_q

; xmin = min(x)
; xmax = max(x)

; ymin = min(y)
; ymax = max(y)

zeta = a->zeta

time = a->time

; xvec = a->x1_u1
; yvec = a->x2_u2


;xvec =  fspan(-24000000,24000000,961)
;yvec = fspan(-27000000,27000000,1081)

xvec =  fspan(todouble(-24000000),todouble(24000000),241)
yvec = fspan(todouble(-27000000),todouble(27000000),271)

zeta_xy = rcm2rgrid_Wrap(y,x,zeta(56,:,:),yvec,xvec,1)

outfile = addfile("zeta_xy.new.nc","c")
outfile ->zeta_xy=zeta_xy

