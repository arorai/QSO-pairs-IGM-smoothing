FUNCTION vdata_fg,Id,zmin,zmax
common qso_gen,qso_bg,qso_fg,qsodir
IF (z_min(Id)-zmin) gt 1e-4  or z_max(Id)-zmax lt -1e-4 THEN BEGIN
    message,'the redshift range you asked for is not fully covered'
ENDIF
IF strcompress(qso_bg[Id].LYA_FILE, /rem) NE '' THEN BEGIN
    file_fg = qsodir + qso_fg[Id].lya_filepath + qso_fg[Id].LYA_FILE
    inflg_fg = lya_choose_flag(qso_fg[Id])
    IF inflg_fg EQ 2 OR inflg_fg EQ 11 THEN BEGIN
        flux_fg = x_readspec(file_fg, inflg = inflg_fg, sig = sig_fg $
                             , wav = wave_fg)
    ENDIF ELSE flux_fg = x_readspec(file_fg, /auto, sig = sig_fg, wav = wave_fg)
    zspace_fg=wl_to_z(wave_fg)
    i_good=where(zspace_fg lt zmax and zspace_fg gt zmin)
    good_z=zspace_fg(i_good)
    vspace_fg=Dz_to_Dv(zmin,good_z)
    return, vspace_fg
ENDIF

END

