FUNCTION bg_noise,Id,zmin,zmax, continuum_norm=cont
common qso_gen,qso_bg,qso_fg,qsodir
IF (z_min(Id)-zmin) gt 1e-4  or z_max(Id)-zmax lt -1e-4 THEN BEGIN
    message,'the redshift range you asked for is not fully covered'
ENDIF
IF strcompress(qso_bg[Id].LYA_FILE, /rem) NE '' THEN BEGIN
    file_bg = qsodir + qso_bg[Id].lya_filepath + qso_bg[Id].LYA_FILE
    inflg_bg = lya_choose_flag(qso_bg[Id])
    IF inflg_bg EQ 2 OR inflg_bg EQ 11 THEN BEGIN
        flux_bg = x_readspec(file_bg, inflg = inflg_bg, sig = sig_bg $
                             , wav = wave_bg)
    ENDIF ELSE flux_bg = x_readspec(file_bg, /auto, sig = sig_bg, wav = wave_bg)
    zspace=wl_to_z(wave_bg)
    i_good=where(zspace lt zmax and zspace gt zmin)
    cfile_bg =repstr(repstr(file_bg, '.gz', ''),'.fits','_c.fits')
    IF file_test(cfile_bg) GT 0 THEN cont_bg  = x_readspec(cfile_bg) $
    ELSE BEGIN ;; Autofit
        message,'ops,there is no continuum here!!'
    ENDELSE

    IF djs_median(flux_bg) LE 1d-10 THEN BEGIN
        flux_bg = 1d17 * flux_bg
        sig_bg = 1d17 * sig_bg
        cont_bg = 1d17 * cont_bg         
    ENDIF
    
    IF KEYWORD_SET(CONT) THEN BEGIN         
        signorm_bg = sig_bg/(cont_Bg + (cont_bg EQ 0.0))
    ENDIF ELSE BEGIN
        signorm_bg = sig_bg
    ENDELSE
    return, signorm_bg[i_good]
ENDIF
END

