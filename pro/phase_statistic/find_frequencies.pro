FUNCTION find_frequencies,v
;choose k-array based on the pixels distribution
;of a velocity array v.

v_side=max(v)-min(v)
min_k=2*!pi/v_side
nk=n_elements(v)/2
return, min_k*(findgen(nk)+1)
END

