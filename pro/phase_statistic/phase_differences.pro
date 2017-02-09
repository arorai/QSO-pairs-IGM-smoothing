FUNCTION sign,x
s=make_array(size(x,/dimension),/int)
S(where(x ge 0))=1
S(where(x lt 0))=-1
return, S
END


FUNCTION common_frequencies,v1,v2
v_side=min([max(v1)-min(v1),max(v2)-min(v2)]);-max([min(v1),min(v2)])
nk=min([n_elements(v1),n_elements(v2)])/2
min_k=2*!pi/v_side
return, min_k*(dindgen(nk)+1)
END


FUNCTION phase_differences,y1,y2,x1,x2,kk=kk,not_recursive=not_recursive,revert=revert
;function that computes phase differences between two functions (y1 & y2)
;sample at the vector x1 and x2. If the x vectors are regular and coincide, it applies
;a standard FFT, otherwise it employs a Least-square approximate method

;if kk is specified, it sets the wavenumbers at which the phases are calculated, 
;otherwise it finds them and they can be retrieved via the same keyword.

npt1=n_elements(x1)
npt2=n_elements(x2)
coinc=mean(abs([x1,npt1]-[x2,npt2]))
;is the space regular?
check=variance(x1[1:*]-x1[0:-1])
IF check+coinc lt 1e-4 and (not keyword_set(kk) )  THEN BEGIN
    kk=common_frequencies(x1,x2)
;regular space, fft can be used
    ft1=FFT(y1,-1,/double,dimension=1)
    ft2=FFT(y2,-1,/double,dimension=1)
    temp=(ft1[1:*,*]*reverse(ft2[1:*,*]))[0:npt1/2-1,*]
	
    ;where the radius is very small there could bad discretization effects
    ;phases are randomized there to avoid artifacts
    bad=where(sqrt(real_part(temp)^2+imaginary(temp)^2) lt 1e-7)
    ph=acos(real_part(temp)/sqrt(real_part(temp)^2+imaginary(temp)^2))*sign(imaginary(temp))
    ph[bad]=randomu(seed,n_elements(bad))*2*!pi-!pi

ENDIF ELSE BEGIN
;space is not regular, Lomb PD must be used
    if not keyword_set(kk) then kk=common_frequencies(x1,x2)
    ft1=recursive_Lomb_FT(y1,x1,kk,not_recursive=not_recursive)
    ft2=recursive_Lomb_FT(y2,x2,kk,not_recursive=not_recursive)
    cross_ft=ft1*conj(ft2)
    bad=where(sqrt(real_part(cross_ft*conj(cross_ft))) lt 1e-7)
    ph= acos(real_part(cross_ft)/sqrt(real_part(cross_ft*conj(cross_ft))))*sign(imaginary(cross_ft))
    ph[bad]=randomu(seed,n_elements(bad))*2*!pi-!pi
ENDELSE
return,ph
END 
