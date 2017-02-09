FUNCTION phase_list,y1,y2,x1,x2,ks=ks,k_in=k_in,r_in=r_in,r_out=r_out,$
                    not_recursive=not_recursive,revert=revert
;produce a list of phases given a list of function pairs.

;y1 and y2 are matrices, with each row being the companion
;function of the correspondent row in the other matrix.
;x1 and x2 are the domains for y1 and y2, which can 
;be regular or not. In the latter
;cases, the fourier transformed is calculated with a Lomb
;algorithm (see 'phase_differences.pro' for more info)

;K_IN can be used to specify the desired wavenumbers, otherwise it
;will calculate them based on x1 and x2. 
;
;R_IN can be used to specify the transverse separations between the elements
;of y1 and y2. It must have as many elements as the number of columns of y1
;(see 'npairs' below) 


;NOT_RECURSIVE and REVERT can be used to slightly modify 
;the phase calculation if LSSA are used  (see recursive_lomb_ft
;for more information)



;OUTPUT
;the function return a list of phases. the keywords 'ks' and 'r_out'
;can be used to retrieve the correspondent list of wavenumbers 
;and transerse separations

npairs=n_elements(y1[0,*])
ph_now=phase_differences(y1,y2,x1,x2,kk=k_in,not_recursive=not_recursive,revert=revert)
nk=n_elements(k_in)
nptot=long(npairs*nk)

ph=reform(ph_now,nptot)
ks=reform(rebin(k_in,nk,npairs),nptot)
if n_elements(r_in) gt 0 then begin
    r_out=reform(transpose(rebin([r_in],npairs,nk)),nptot)
endif

return,ph
END
