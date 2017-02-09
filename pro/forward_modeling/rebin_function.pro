FUNCTION trans_matrix,x_old,x_new
;transfer matrix between the  old and new pixel grid. 
;I tried to make the calculation O(n), but it is still 
;a bit slow. There are likely more efficient implementations around.

no=n_elements(x_old)
nn=n_elements(x_new)
A=make_array(nn,no)

seps_old=(x_old[1:*]+x_old[0:-1])*0.5
seps_old=[2*x_old[0]-seps_old[0],seps_old,2*x_old[no-1]-seps_old[no-2]]
seps_new=(x_new[1:*]+x_new[0:-1])*0.5
seps_new=[2*x_new[0]-seps_new[0],seps_new,2*x_new[nn-1]-seps_new[nn-2]]

scan=max([seps_new[0],seps_old[0]])

i=0
j=0
totrow=0

ipre=where(seps_new le seps_old[1])
ipost=where(seps_new ge seps_old[no-1])
if ipre[0] ne -1 then A[ipre,0]=1.
if ipost[0] ne -1 then A[ipost,no-1]=1.
totmat=nn-total(A)
i=max(ipre)+1
WHILE abs(totmat) gt 1e-4 do begin
    DX=min([seps_new[i+1],seps_old[no]])-max([seps_new[i],seps_old[0]])    
    IF seps_old[j+1] ge seps_new[i] and seps_new[i+1] ge seps_old[j] THEN BEGIN
        sup=min([seps_new[i+1],seps_old[j+1]])
        inf=max([seps_old[j],seps_new[i]])
        A[i,j]=(sup-inf )/DX  
        totrow+=A[i,j]
    ENDIF
    j+=1
    IF abs(totrow-1) le 1e-4  THEN BEGIN
        A[i,*]*=1.0D/totrow
        totmat-=1
        totrow=0
        i+=1
        j-=1
    ENDIF 
ENDWHILE
return,A

END

FUNCTION rebin_function,y,x_old,x_new,rebin_matrix=rebin_matrix
;function that takes a function y (could be multi valued,i.e. a matrix)
;defined on the x_old space, and rebins it on the x_new space.
;
;rebinning is done weighting bins accordingly to the overlapping fraction
;between the old ones and the new ones (e.g. if a new pixel has 60% overlap
; with an old one, the value of y in that old pixel will weight 60%).
;This should preserve the number of photons.

y_new=dblarr(n_elements(x_new))

If ~ (keyword_set(rebin_matrix)) then rebin_matrix=trans_matrix(x_old,x_new)
y_new=rebin_matrix#Y
return,y_new
END

