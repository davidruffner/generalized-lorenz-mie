;+
;NAME:
;    dbr_xi_n
;
; PURPOSE:
;    Calculate the riccati-bessel function associated with the
;    spherical bessel, jn(x)
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    out = dbr_xi_n(x,n)
;
;INPUTS:
;    n:    index
;
;    x:    input
;
;OUTPUTS:
;    out: value of the riccati-bessel function of order n at x
;
;DEPENDENCY:
;    Calls idl's BESELJ and BESELY functions.
;
;MODIFICATION HISTORY:
; 03/19/2010 Written by David B. Ruffner, New York University

function dbr_xi_n,x,n

xi_n = sqrt(!pi*x/2.d)*(beselj(x,n+1/2.d)+complex(0,1)*besely(x,n+1/2.d))

b1 = where(finite(real_part(xi_n),/nan),/null)
if b1 ne !NULL then begin $
      xi_n(b1) = complex(0,0) & $
endif & $

return,xi_n
end
