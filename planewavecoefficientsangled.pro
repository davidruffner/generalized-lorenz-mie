;+
;NAME:
;    planewaveecoefficientsangled
;
; PURPOSE:
;    Calculate the beam shape coefficients of a plane wave beam at an angle
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    E = planewaveecoefficientsangled(nmax,alpha,beta,gamma)
;
;INPUTS:
;
;    nmax:  max value of n in the vector spherical harmonic series
;expansion
;
;    alpha: azimuthal angle
;
;    beta: polar angle
;
;    gamma: rotation of polarization from x axis
;
;OUTPUTS:
;    AB:    A [nmax,2*nmax+1,2] matrix of the beam shape coefficients
;
;REFERENCE:
;
;MODIFICATION HISTORY:
; 03/21/2014 Written by David B. Ruffner, New York University

function planewavecoefficientsangled,nmax,alpha,beta,gamma

am_mn = fltarr(nmax + 1.d,2.d*nmax + 1.d)*complex(1.d,0)
ae_mn = fltarr(nmax + 1.d,2.d*nmax + 1.d)*complex(1.d,0)
ci = complex(0,1.d)

cosgamma = cos(double(gamma))
singamma = sin(double(gamma))
cosbeta = cos(double(beta))

for n=1,nmax do begin
   for m = -n,n do begin
      tau_mn = dbr_tau_mn(cosbeta,double(n),double(m))
      pi_mn = dbr_pi_mn(cosbeta,double(n),double(m))
      ;played around with factors to get it to match
      am_mn[n,m] = - ci^(double(n))*(tau_mn*singamma+ci*pi_mn*cosgamma)* $
                     exp(-ci*m*alpha);*sqrt(n*(double(n)+1))/(norm^2.
      ae_mn[n,m] = - ci^(n + 1.d)*(tau_mn*cosgamma-ci*pi_mn*singamma)* $
                     exp(-ci*m*alpha);*sqrt(n*(double(n)+1))/(norm^2.
      am_mn[n,m] *= ci*4.d*!pi/sqrt(n*(double(n) + 1.d));1/(4*!pi*sqrt(n*(double(n)+1)))
                         ;sqrt(n*(double(n)+1))/(e_mn*norm^2.)
      ae_mn[n,m] *= 4.d*!pi/sqrt(n*(double(n) + 1.d));1/(4*!pi*sqrt(n*(double(n)+1)))
                         ;sqrt(n*(double(n)+1))/(e_mn*norm^2.)
      ;tried this combination 3_24_14
      ;; e_mn = 4*!pi*(2*double(n) + 1.d)/(n*(double(n)+1))
      ;; ae_mn[n,m] = - e_mn*ci^(double(n) + 1.d) * $
      ;;              (tau_mn*cosgamma - ci*pi_mn*singamma)*exp(-ci*m*alpha)
      ;; am_mn[n,m] = - ci*e_mn*ci^(double(n)) * $
      ;;              (tau_mn*singamma + ci*pi_mn*cosgamma)*exp(-ci*m*alpha)
   endfor
endfor

b1 = where(finite([[[am_mn]],[[ae_mn]]],/nan),/null)
if b1 ne !NULL then begin $
      print, "We have infinities"
      stop
endif & $

return,[[[am_mn]],[[ae_mn]]]
end
