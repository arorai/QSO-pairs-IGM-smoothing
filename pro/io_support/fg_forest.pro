FUNCTION fg_forest,Id,zmin,zmax, continuum_norm=cont
common qso_gen,qso_bg,qso_fg,qsodir
IF (z_min(Id)-zmin) gt 1e-4  or z_max(Id)-zmax lt -1e-4 THEN BEGIN
    message,'the redshift range you asked for is not fully covered'
ENDIF
IF strcompress(qso_fg[Id].LYA_FILE, /rem) NE '' THEN BEGIN
    file_fg = qsodir + qso_fg[Id].lya_filepath + qso_fg[Id].LYA_FILE
    inflg_fg = lya_choose_flag(qso_fg[Id])
    IF inflg_fg EQ 2 OR inflg_fg EQ 11 THEN BEGIN
        flux_fg = x_readspec(file_fg, inflg = inflg_fg, sig = sig_fg $
                             , wav = wave_fg)
    ENDIF ELSE flux_fg = x_readspec(file_fg, /auto, sig = sig_fg, wav = wave_fg)
    zspace=wl_to_z(wave_fg)
    i_good=where(zspace lt zmax and zspace gt zmin)
    cfile_fg =repstr(repstr(file_fg, '.gz', ''),'.fits','_c.fits')
    IF file_test(cfile_fg) GT 0 THEN cont_fg  = x_readspec(cfile_fg) $
    ELSE BEGIN ;; Autofit
        message,'ops,there is no continuum here!!'
    ENDELSE

    IF djs_median(flux_fg) LE 1d-10 THEN BEGIN
        flux_fg = 1d17 * flux_fg
        cont_fg = 1d17 * cont_fg         
    ENDIF

    
    IF KEYWORD_SET(CONT) THEN BEGIN         
        fnorm_fg = flux_fg/(cont_fg + (cont_fg EQ 0.0))
    ENDIF ELSE BEGIN
        fnorm_fg = flux_fg
    ENDELSE
    return, fnorm_fg[i_good]
ENDIF
END

