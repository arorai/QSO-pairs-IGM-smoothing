;+
; NAME:
;   ADD_LLS
;
; PURPOSE:
;   For an input continuum, add optically thick Lyman-alpha absorbers
;   corresponding to a path length incidence of Ribaudo 2011, and a broken power-law
;   column density distribution derived by Prochaska et al 2010. We
;   assume that the column density distribution is constant at all
;   redshifts. The code is hardwired to draw only log(N_HI) =
;   [17.5,20.3] systems, but pLLS with log(N_HI)=[16.5,17.5]
;
;   We assume that the input segment is short enough that just the
;   central redshift value of the segment is used to sample the column
;   density distribution. Also assume that the input segment is
;   continuous w.r.t. wavelength. 
;
;
; CATEGORY:
;   Function
;
; CALLING SEQUENCE:
;   output_spectrum = add_lls(wave_in, cont_in [, seed=seed,
;   z_eff=z_eff, linestr=linestr, scale_factor=scale_factor,
;   plls_frac=plls_frac, bottomheavy=bottomheavy]) 
;
; INPUTS:
;   wave_in       - Input observed wavelength vector 
;   cont_in       - Input continuum or flux vector to which the
;                   optically thick systems are to be added
;
; OPTIONAL INPUTS:
;   z_eff         - Redshift used to sample the path length density
;                   distribution. If not set, the redshift
;                   corresponding to the median wavelength is used. 
;   seed          - Random number seed 
;   scale_factor  - Factor by which to scale the number of absorbers
;                   drawn per path length.
;   plls_frac     - If set, will also add pLLS corresponding to some
;                   factor of the LLS incidence
;
; KEYWORD OPTIONS:
;   bottomheavy   - Flag to be passed to DRAW_LLSCOLUMN and
;                   DRAW_PLLSCOLUMN to get steeper column density
;                   distribution 
;
; OUTPUTS:
;   output_spec   - CONT_IN with optically thick Lya systems added
;
; OPTIONAL OUTPUTS:
;   linestr       - Structure of LLS within this spectrum, that is fed
;                   to x_allvoigt
;
; PROCEDURES USED:
;   draw_llscolumn(), x_allvoigt()
;
; MODIFICATION HISTORY:
;  Khee-Gan Lee (MPIA) 27-01-2013
;        Original code
;  Khee-Gan Lee (MPIA) 04-02-2014
;        Changed LLS incidence to Ribaudo et al 2011 values
;
;-
function add_lls, wave_in, cont_in, z_eff=z_eff,seed=seed, linestr=linestr, $
                  scale_factor=scale_factor, lognhi_range=lognhi_range, $
                  plls_frac=plls_frac, bottomheavy=bottomheavy
  
  if not keyword_set(z_eff) then z_eff = (median(wave_in)/1215.7)-1.

  if not keyword_set(scale_factor) then scale_factor = 1.

  delta_z = (max(wave_in) - min(wave_in))/1215.7

  ; Number of LLS expected per unit redshift, at this redshift.
  ;; Following power-law of Songaila and Cowie 2010, but normalization
  ;; is using Prochaska et al 2010's estimates at z=3.7. We assume this
  ;; spans the column density range 17.2 < log N_HI < 20.3 
  ;lz = 1.7 * ((1.+z_eff)/4.7)^2 * scale_factor
  lz = 0.1157 * (1.+z_eff)^1.83 * scale_factor

  ;; This is the expectation value for the number of LLS
  nlls_exp = delta_z * lz
  ;; and this is the number of LLS in this particular segment
  nlls = randomn(seed, poisson=nlls_exp)
  
  tau = replicate(0., n_elements(wave_in))
  if nlls NE 0 then begin
     nhi_vec = fltarr(nlls)
     zlls_vec = fltarr(nlls)
     ;; Now draw from the column density distribution, also random
     ;; redshifts within the input wavelenght range
     for il=0, nlls-1 do begin
        nhi_vec[il] = draw_llscolumn(seed, bottomheavy=bottomheavy)
        zlls_vec[il] = randomu(seed)*delta_z + min(wave_in)/1215.67-1.
     endfor
     
     linestr = {ion:'HI', wrest:1215.6701, f:0.41460, gamma:6.265e8, $
                n:0., b:45., zabs:0.}
     
     linestr = replicate(linestr, nlls)
     
     linestr.n = alog10(nhi_vec)
     linestr.zabs = zlls_vec
     
     tau = x_allvoigt(wave_in, linestr)
  endif

  if keyword_set(plls_frac) then begin
     ;; This is the expectation value for the number of pLLS
     nplls_exp = plls_frac * float(nlls_exp)
     ;; and this is the number of LLS in this particular segment
     nplls = randomn(seed, poisson=nplls_exp)

     if nplls NE 0 then begin
     nhi_vec = fltarr(nplls)
     zlls_vec = fltarr(nplls)
     ;; Now draw from the column density distribution, also random
     ;; redshifts within the input wavelenght range
     for il=0, nplls-1 do begin
        nhi_vec[il] = draw_pllscolumn(seed, bottomheavy=bottomheavy)
        zlls_vec[il] = randomu(seed)*delta_z + min(wave_in)/1215.67-1.
     endfor
     
     linestr = {ion:'HI', wrest:1215.6701, f:0.41460, gamma:6.265e8, $
                n:0., b:45., zabs:0.}
     
     linestr = replicate(linestr, nplls)
     
     linestr.n = alog10(nhi_vec)
     linestr.zabs = zlls_vec
     
     tau += x_allvoigt(wave_in, linestr)
  endif

  endif     

return, exp(-tau) * cont_in
end

