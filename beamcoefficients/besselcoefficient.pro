;+
;NAME:
;    BESSELCOEFFICIENT
;
; PURPOSE:
;    Calculate the coefficient for the vector spherical harmonics of
;    a Bessel beam as a function of position in the field and m and n
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    [p_mn,q_mn] = besselcoefficient(pos,theta0,m,n)
;
;INPUTS:
;    pos:    [3,N] array of positions where you want the coefficients
;
;    theta0: Convergence angle of Bessel beam
;
;    m,n:    indicies of the coefficient
;
;OUTPUTS:
;    [p_mn,q_mn]: Complex coefficient for the vector spherical
;    harmonics
;
;REFERENCE:
;    Taylor, J. M. and G. D. Love. "Multipole expansion of
;    bessel and gaussian beams for mie scattering calculations"
;    J. Opt. Soc. Am. A. 26, 278(2009)
;
;MODIFICATION HISTORY:
; 03/06/2014 Written by David B. Ruffner, New York University
; 03/21/2014 Updated pi and tau functions
; 03/25/2014 Added a factor to match up with Jackson

function besselcoefficient,n,m,pos,theta0,k

;arranging imput
sz = size(pos)
if sz[0] eq 1 then begin
   x = pos[0]
   y = pos[1]
   z = pos[2]
endif else begin
   x = pos[0,*]
   y = pos[1,*]
   z = pos[2,*]
endelse
rho = sqrt(x^2 + y^2)
phi = atan(-y,x)-!pi/2

;Calculating factors
ci = complex(0,1)

jacksonfactor = sqrt(n*(double(n) + 1.d))

phasefactor_e = -complex(0,1.d)/(!pi*2.d);Determined by matching coefficients to the plane 
                           ;wave case
phasefactor_m = complex(1.d,0)/(!pi*2.d)

un = (4.d*!pi*ci^double(n))/(n*(double(n) + 1.d))

expftr = exp(ci*z*k*cos(theta0))

nmfactor = sqrt((2*n+1)*factorial(n-m)/(4*!pi*factorial(n+m)))

pi_mn = dbr_pi_mn(cos(theta0),n,m)

tau_mn = dbr_tau_mn(cos(theta0),n,m)

x = k*sin(theta0)*rho
i1 = !pi*exp(ci*(m-1)*phi)*beselj(x,1-m) 
i2 = !pi*exp(ci*(m+1)*phi)*beselj(x,-1-m)
i_plus = i1+i2
i_minus = i1-i2

ae_mn = jacksonfactor*un*expftr*(tau_mn*i_plus+pi_mn*i_minus)*phasefactor_e

am_mn = jacksonfactor*un*expftr*(pi_mn*i_plus+tau_mn*i_minus)*phasefactor_m

return,[am_mn,ae_mn]
end
















;Author: David Ruffner
