function draw_llscolumn, seed, prob=prob, bottomheavy=bottomheavy
;+
; This draws LLS columns densities in the range 17.2 < log(N_HI) <
; 20.3. This distribution is relative, i.e. it is normalized within
; the column density range, and the path length density l(x) and
; redshift evolution has to be specified externally.
;
; Returns the neutral hydrogen column density in cm^-2
;
; INPUTS:
;   seed    - Random number seed
;   prob    - Probability drawn from RANDOMU
; 
; EXAMPLE:
;   nhi = draw_llscolumn(seed, prob=prob)
;
; Modification History:
;   KG Lee 24/01/2014 - Changed to lower-limit of range from 17.2 to
;                       17.5. Added option for bottom-heavy CDDF
;-

prob = randomu(seed)

if keyword_set(bottomheavy) then begin

   logk1 = 2.819
   logk2 = 7.039
   
   b1 = -1.2
   b2 = -1.4
   
   if prob LE 0.52 then begin
      term1 = prob * (1.+b1)/10.^logk1
      term2 = 10.^(17.5*(1.+b1))
      nhi = (term1 + term2)^(1./(1.+b1))
   endif else begin
      term1 = (prob-0.52) * (1.+b2)/10.^logk2
      term2 = 10.^(19.*(1.+b2))
      nhi = (term1 + term2)^(1./(1.+b2))
   endelse
   
endif else begin
   logk1 = -4.477
   logk2 = 3.123
   
   b1 = -0.8
   b2 = -1.2
   
   if prob LE 0.525 then begin
      term1 = prob * (1.+b1)/10.^logk1
      term2 = 10.^(17.5*(1.+b1))
      nhi = (term1 + term2)^(1./(1.+b1))
   endif else begin
      term1 = (prob-0.525) * (1.+b2)/10.^logk2
      term2 = 10.^(19.*(1.+b2))
      nhi = (term1 + term2)^(1./(1.+b2))
   endelse
endelse

return, nhi
end

