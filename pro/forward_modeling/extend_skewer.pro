FUNCTION extend_skewer,f,vgrid,v_data,newvgrid=newvgrid
;extend periodically simulated sightlines until they reach
;the length of 'v_data'

;'f' is the Npix x n_sightlines matrix, vgrid is the initial 
;velocity space to be extended. 

;the final velocity space can be retrieved through the 
;'newvgrid' keyword. 
;
;The output is a new flux matrix with the extended sightlines

npix=n_elements(f[*,0])
nsk=n_elements(f[0,*])
data_vsize=max(v_data)-min(v_data)
skewer_vsize=max(vgrid)-min(vgrid)
prop=data_vsize/skewer_vsize
n=fix(npix*prop)+1
;calculate the new v space
newvgrid=min(vgrid)+findgen(n)/(npix-1)*skewer_vsize
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;go on
newsk=make_array(n,nsk)
nint=long(prop) 
nmod= n mod npix
for i=0,nint-1 do begin
    newsk[i*npix:(i+1)*npix-1,*]=f
endfor
newsk[nint*npix:nint*npix+nmod-1,*]=f[0:nmod-1,*]
return,newsk
END
