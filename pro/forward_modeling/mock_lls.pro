FUNCTION mock_lls,f,vgrid,zmin,scale_factor=scale_factor,lognhi_range=lognhi_range
lya=1215.67d
c=3.0e5
wave=lya*(1+zmin)*exp(vgrid/c)
nsk=n_elements(f[0,*])
flls=f

IF ~KEYWORD_SET(scale_factor) THEN scale_factor=0.85
IF ~KEYWORD_SET(lognhi_range) THEN lognhi_range=[17.2,19.0]

for i=0,nsk-1 do begin
    flls[*,i]=add_lls(wave,f[*,i],/bottomheavy,plls_frac=1.8,$
      scale_factor=scale_factor,lognhi_range=lognhi_range)    
endfor
return, flls
END
