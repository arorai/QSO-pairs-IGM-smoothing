FUNCTION vdata_bg,Id,zmin,zmax
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
    zspace_bg=wl_to_z(wave_bg)
    i_good=where(zspace_bg lt zmax and zspace_bg gt zmin)
    good_z=zspace_bg(i_good)
    vspace_bg=Dz_to_Dv(zmin,good_z)
    return,vspace_bg
ENDIF
END
