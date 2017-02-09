FUNCTION mock_metals,f,vgrid,zmin,extra_metals=extra_metals
;add metal lines to a flux f using KG's method.
;if you want more metals, set the keyword 'extra_metals'
;to the fraction you want metal abundance to be increased by.
;example: extra_metals=0.3 adds 30% more metals
;extra_metals must be positive. Subtracting metals is not 
;implemented yet

;this is potentially problematic: you need a metal line list in order to
;do this.

fmet=f
;prepare list of spectra with metals
lya=1216.
c=3.0e5

;modify redshift below z=2.3, as there is no wavelength coverage
;for metals below that
zlow=max([zmin,2.3])

wave=lya*(1+zlow)*exp(vgrid/c)

catalog_file='/home/lee/lya/pdf/lls_metals/dr7_z.gt.1.7_snr.list'
readcol,catalog_file,plate,fiber,MJD,zq,SN,skipline=1,$
  format=['I,I,LL,F,F']
cond1=SN gt 15.0
cond2=zq lt min(wave)/1260-1 
cond3=zq gt max(wave)/1390-1 
qso_good=where(cond1 and cond2 and cond3)
nqso=n_elements(qso_good)
print,nqso


;tau_0=5.
;dwdv=lya*(1+zlow)/c ; approximate linear conversion between wavelength and 

fmetals=make_array(n_elements(vgrid),nqso,/float)
for i=0,nqso-1 do begin  
    ix=qso_good(i)
    temp=add_metals_segment(wave,plate[ix],fiber[ix],mjd=mjd[ix],waveline=ww,ewline=ee,bline=bb)    
    fmetals[*,i]=temp
endfor

nsk=n_elements(f[0,*])

;pick a random skewer from metal contaminated ones
;and substitute to a metal free one. Do this 
;extra_metals X number of metal contaminated skewers.
IF KEYWORD_SET(extra_metals) THEN BEGIN
    mf=mean(fmetals,dimension=1)
    i_met=where(mf lt 0.999)
    i_no_met=where(mf ge 0.999)
    n_met=n_elements(i_met)
    for j=0,n_met*extra_metals DO BEGIN
        fmetals[*,i_no_met[fix(randomu(seed)*(nqso-n_met))]]= $
          fmetals[*,i_met[ fix(randomu(seed)*n_met) ] ]
    ENDFOR
ENDIF

for i=0,nsk-1 do begin
    ;chose one of the skewers of fmetals
    randompick=long(randomu(sdf)*nqso)
    fmet[*,i]=f[*,i]*fmetals[*,randompick]
endfor


return, fmet
END
