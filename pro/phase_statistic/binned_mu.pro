FUNCTION binned_mu,ph,k,kbins
;this function takes a list of phases
;with relative wavenumbers and returns
;an array of wrapped-Cauchy concentration parameters 'mu'.
;These correspond to the best wrapped
;Cauchy function that reproduces the
;phase distribution in each k-bin. 

;the vector kbins specifies the limits
;of each k-interval

nb=n_elements(kbins)-1
mu=fltarr(nb)
FOR i=0,nb-1 DO BEGIN
    good_i=where(k lt kbins[i+1] and k gt kbins[i])
    if good_i[0] ne -1 then begin
    mu[i]=optimize_mu(ph[good_i],0.5)
    endif
ENDFOR
return,mu
END
