FUNCTION noise_skewers,f,sigma_v,Nrealization
;add noise to skewers according to data
;skewers must be already rebinned

;f is the npix x n_sightlines matrix of the flux, sigma is
;a vector with npix elements where you specify the amplitude 
;of the noise as a function of the pixel. 

;At the moment, the noise is flux-independent, this may 
;change in the future. 

;if you are dealing with small numbers, you can draw
;multiple realization of the noise by setting Nrealization
;to whatever you like. 

if not(keyword_Set(Nrealization)) then Nrealization=1

npix=n_elements(f[*,0])
nsk=n_elements(f[0,*])
fout=make_array(npix,nsk*Nrealization)
if n_elements(sigma_v) ne npix then begin
    message,"Skewers and noise don't have the same length, something is wrong!"
endif
for i=0,Nrealization*nsk-1 do begin
    fout[*,i]=f[*,i mod nsk]+randomn(undef,npix)*sigma_v
endfor
return,fout
END
