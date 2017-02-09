FUNCTION phase_pdf, phases,nbins,x_vec=x_vec, err=err
;return the probability distribution function of the
;absolute values of a set of phases, in the interval
;[0,pi]. 
;
;takes an ensemble of phases as input, and the number of bins  
;if x_vec is specified, the number of bins is ignored.

;X_VEC: centers of the the PDF binds. If not specified,  
;it can be returned as an output

;ERR: Poisson errorbars on each bins. will be NaN if the PDF
; is zero (improve this?)




if keyword_set(x_vec) then begin
	nbins=n_elements(x_vec)
endif else begin
	x_vec=(findgen(nbins)+0.5)/(nbins)*!pi
endelse

freq=histogram(abs(phases),min=0,max=!pi*(1-1./(nbins-1) ),nbins=nbins )

renormalization=!pi/nbins*n_elements(phases) ;renormalize by N*dtheta
pdf=freq/renormalization
err=1/sqrt(freq)*pdf

return,pdf

END
