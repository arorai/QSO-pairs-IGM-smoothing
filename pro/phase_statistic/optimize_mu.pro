;given a sample of phases {x} and an initial 
;guess 'mu_in', it calculates the wrapped-Cauchy 
;concentration parameter that otpimizes the likelihood
;of the phase distribution.
FUNCTION optimize_mu, x, mu_in
mu=double(mu_in)
ctrl=1d

new_x=abs(x)
while ctrl gt 1e-5 DO begin
mu_new=total(cos(new_x)*wrapped_cauchy(new_x,mu))/total(wrapped_cauchy(new_x,mu))
ctrl=abs(mu-mu_new)/mu
mu=mu_new
endwhile
return,mu
END
