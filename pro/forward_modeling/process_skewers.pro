FUNCTION process_skewers,f,vo,vn,res,sig,metals=metals,LyLS=LyLS,zmin=zmin, zmax=zmax $
                         ,LLsk=LLsk,METsk=METsk,rebin_matrix=rebin_matrix,extra_metals=extra_metals,$
			scale_factor=scale_factor



;function responsible for the forward-modeling of the simulation,putting
;together the different steps. 

;INPUT: npix X n_sightlines flux matrix 'f' 

;	original velocity grid of the simulation 'vo'

;	velocity  grid of the data to which the simulation
;	must be adapted 'vn'
	
;	resolution FWHM 'res'

;	noise array 'sig'. It must have the same lenght as 'vn', for consistency

;KEYWORDS: 

;	METALS: if set, it adds metals to the sightlines. CAREFUL! this requires 
;	a metal catalog.

;	LyLS: if set, add LLSs to the sightlines. 

;	ZMIN and ZMAX: delimit the redshift range for the metals and lls distributions

;	EXTRA_METALS: scale up metal contamination (see 'mock_metals.pro' for more info)

;	SCALE_FACTOR: artificially increase the number of LLS

;	LLsk,METsk can be used to obtain the added contaminants (for inspection)
;	REBIN_MATRIX to retrieve the matrix used for rebinning (also for inspection)


;OUTPUT
;	fully post-processe sightlines.
;	note that the right velocity space for the output is 'vn'


time=systime(/seconds)
datapix=min(vn[1:*]-vn[0:-1])
datapix_over_simpix=2^fix(alog(datapix/mean(vo[1:*]-vo[0:-1]))/alog(2)) ;biggest power of 2
rebin_factor=max([1,datapix_over_simpix/4]) ;require at least 4 simulated pixels for
;each data pixel,as starting point
f=rebin(f,n_elements(vo)/rebin_factor,n_elements(f[0,*]))
vsm=rebin(vo,n_elements(vo)/rebin_factor) ;makes things faster
print,'done first rebinning,time',systime(/second)-time

time=systime(/seconds)
;periodic extension
chunksize=max(vo)
f=extend_skewer(temporary(f),vsm,vn,newvgrid=vext)
print,'done extension,time',systime(/second)-time
time=systime(/seconds)


;if required add lls  (zmin must be declared)
if keyword_set(LyLS) then begin
    if ~ (keyword_set(LLsk)) then begin
        continuum_matrix=rebin([1d],n_elements(f[*,0]),n_elements(f[0,*]))
        LLsk=mock_lls(temporary(continuum_matrix),vext,zmin,scale_factor=scale_factor)
    endif
    f*=LLsk
ENDIF
print,'done lls,time',systime(/second)-time

time=systime(/seconds)
;if required add metals (already smoothed by resolution)
if keyword_set(metals) then begin
    IF ~ (keyword_set(METsk)) then begin
        continuum_matrix=rebin([1d],n_elements(f[*,0]),n_elements(f[0,*]))
        METsk=mock_metals(temporary(continuum_matrix),vext,zmin,extra_metals=extra_metals)
    ENDIF
    f*=METsk
ENDIF


;COMMENTED THE FOLLOWING, BUT COULD BE USED
;now renormalize to the mean flux of this redshift
;if keyword_set(zmin) and keyword_set(zmax) then begin
;	meanflux=lya_fbar_faucher((zmin+zmax)/2)
;	f=makeflux(-alog(f),meanflux)
;endif

print,'done metals,time',systime(/second)-time
time=systime(/seconds)


;add resolution smoothing
nchunks=fix(max(vext)/chunksize)

for i=0,nchunks do begin
    indsnow=where(vext gt i*chunksize and vext lt (i+1)*chunksize)
    if indsnow[0] ne -1 then begin
    f[indsnow,*]=add_resolution_effect(temporary(f[indsnow,*]),vext[indsnow],res)
    endif
endfor
print,'done resolution,time',systime(/second)-time
time=systime(/seconds)

;print,'now rebin pairs as in data'
f=rebin_function(temporary(f),vext,vn,rebin_matrix=rebin_matrix)
print,'done second rebinning,time',systime(/second)-time
time=systime(/seconds)


;print,'add noise to skewers as in data'
Nrealization=5
ffinal=noise_skewers(temporary(f),sig,Nrealization)
print,'done noise,time',systime(/second)-time
time=systime(/seconds)
return,ffinal
END

