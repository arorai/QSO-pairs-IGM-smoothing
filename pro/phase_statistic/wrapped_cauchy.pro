FUNCTION wrapped_cauchy, x, mu_1 
;computes the probability density at x assuming
;a concentration parameter mu_1

;only for symmetric Cauchy
c=1/sqrt(1-mu_1^2)
return, 1/(!pi*c*(1-mu_1*cos(x)))
END
