;+
;NAME:
;    dbr_riccatibesselprime
;
; PURPOSE:
;    Calculate the derivative of the riccati-bessel function
;    associated with the spherical bessel, jn(x)
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    out = dbr_riccatibesselprime(x,n)
;
;INPUTS:
;    n:    index
;
;    x:    input
;
;OUTPUTS:
;    out: value of the derivative of the riccati-bessel function of
;    order n at x
;
;DEPENDENCY:
;    Calls idl's BESELJ function.
;
;MODIFICATION HISTORY:
; 03/14/2010 Written by David B. Ruffner, New York University

function dbr_riccatibesselprime,x,n

psin = sqrt(!pi*x/2.)*beselj(x,n+1/2.)
b1 = where(finite(psin,/nan),/null)
if b1 ne !NULL then begin $
      psin(b1) = 0 & $
endif & $

psiprime = psin*(-n/x+beselj(x,n-1/2.
return,psiprime
end
