FUNCTION phase_wn_bins
;define kbins used in the analysis. 
;log spaced. 
;Redefine it if you like


mink=0.005
maxk=0.1
nbins=12
;spaced in log
bins=alog10(mink)+findgen(nbins+1)/nbins*(alog10(maxk/mink))
return,10^bins
END

