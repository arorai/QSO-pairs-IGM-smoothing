FUNCTION add_resolution_effect,flux,vgrid,FWHM,RSGC=RSGC
;convolve the flux field 'flux', on a REGULAR velocity space 'vgrid',
;with a gaussian kernel with specified FWHM. By default, this is done
;in Fourier space, but you can do it in real space by setting the 
;keyword RSGC (mainly for testing of for numerical reasons)


nsk=n_elements(flux[0,*])
npix=n_elements(flux[*,0])
sigma=FWHM/2.335
output=make_array(npix,nsk,/double)
ft=make_array(npix,nsk,/complex)

;convolution kernel
xvec=vgrid-vgrid(npix/2)
dv=vgrid[1]-vgrid[0]

IF KEYWORD_SET(RSGC) then begin
    for i=0,nsk-1 do begin
        output[*,i]=gaussian_convolve(vgrid,flux[*,i],sigma,nsigma=3.0)
    endfor
endif else begin

    kernel=shift(dv/sqrt(2*!pi*sigma^2)*exp(-xvec^2/(2d*sigma^2)),-npix/2)
    fker=fft(kernel,1)
    fker=real_part(fker)
    for i=0,nsk-1 do begin
               
        temp=fft(flux[*,i],1)*fker
        output[*,i]=shift(real_part(fft(temp,-1)),-npix/2)
  endfor
endelse  

return,output
END
