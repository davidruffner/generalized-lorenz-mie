;quick script to test angular functions tau pi

n = 2
m = 1
xmax = 2*!pi
xmin = 0
npts = 1000

x = (xmax-xmin)*findgen(npts)/(npts-1)+xmin

costheta = cos(x)

norm = sqrt((2*n+1)*factorial(n-m)/(4*!pi*factorial(n+m)))

pi_mn = -dbr_pi_mn(costheta,n,m)/norm
tau_mn = -dbr_tau_mn(costheta,n,m)/norm

pitau_n = dgg_pi_n(costheta,n)
pi_n = pitau_n[*,0];*norm
tau_n = pitau_n[*,1];*norm

ntau_mn = Deriv(x,sin(x)*dbr_pi_mn(costheta,n,m)/m)


p3 = polarplot(pi_mn,x,color="green",symbol="X",aspect_ratio=1)
p2 = polarplot(tau_mn,x,color="red",symbol="X",/overplot)
p1 = polarplot(tau_n,x,color="blue",/overplot)
p0 = polarplot(pi_n,x,color="black",/overplot)

;; p3.xtitle = "$\theta$"
;; p3.ytitle = "f(cos($\theta$))"
p3.title = "Angular functions vs. $\theta$"

t = text(.25,.7,"$n=$"+strtrim(n,2)+" $m=$"+strtrim(m,2))
p3.save,"testangularfunctions4.png"
