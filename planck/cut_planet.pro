result = command_line_args(count=count)
lat = float(result[0])
lon = float(result[1])
radius = float(result[2])
nside = float(result[3])

;glon_glat  = [lon 33.75, lat -40.33]   ; Jupiter
;radius     = 1.5d              ; deg
;nside      = 512

ang2vec, lat , lon ,vector,/astro

query_disc,nside,vector,radius,listpix,/nested,/deg,/inclusive

print, listpix

exit, status = 0
