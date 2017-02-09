;calculate the Fourier coefficients for the 
;wavenumbers 'kk' of the function 'y' sampled 
;in a vector 'x', which could be irregularly 
;spaced in the domain. Uses the least-square
;approximation and the Lomb phase method.

;NOT_RECURSIVE: by default, the algorithm starts
;from the lowest wavenumber k, and after the 
;calculation of each component it subtracts it 
;from the function (to have orthogonal components). 
;Set this keyword to skip this subtraction.

;REVERT: set this keyword if you want to start the
;calculation from the highest k instead of the lowest
;

FUNCTION Lomb_shift,k,x
;this can be done once for all skewers
tot=total(exp(dcomplex(0,2*k*x) ) )
return, 1/(2)*atan(real_part(tot)/imaginary(tot));you may needto putk back here!
END

FUNCTION Lomb_FT,y,x,k
nf=n_elements(y[0,*])
theta=Lomb_shift(k,x)
c_f=1d#cos(k*x-theta)/total(cos(k*x-theta)^2)
s_f=1d#sin(k*x-theta)/total(sin(k*x-theta)^2)
c=y##c_f
s=y##s_f
return, dcomplex(c,s)*exp(dcomplex(0,-theta))
END

FUNCTION recursive_Lomb_FT,y,x,kk,not_recursive=not_recursive,revert=revert
nk=n_elements(kk)
nf=n_elements(y[0,*])
ft=make_array(nk,nf,/dcomplex)
func=y
for l=0,nf-1 do begin
;subtract the mean before calculating the 
;fourier components
func[*,l]=y[*,l]-mean(y[*,l])
endfor
if keyword_set(revert) then kk=reverse(kk)
for i=0,nk-1 do begin
    ft[i,*]=Lomb_FT(double(func),double(x),double(kk[i]))
    if not keyword_set(not_recursive) then begin
        func-=real_part(ft[i,*])##cos(x*kk[i])+imaginary(ft[i,*])##sin(x*kk[i])
    endif
endfor
return, ft
END

