;+
;NAME:
;    dbr_j_n
;
; PURPOSE:
;    Calculate the first spherical hankel function j_n(x)
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    out = dbr_j_n(x,n)
;
;INPUTS:
;    n:    index
;
;    x:    input
;
;OUTPUTS:
;    out: value of the spherical bessel function of first kind of order n at x
;
;DEPENDENCY:
;    Calls idl's BESELJ  functions.
;
;MODIFICATION HISTORY:
; 10/01/2014 Written by David B. Ruffner, New York University

function dbr_j_n,x,n

j_n = sqrt(!pi/(2.d*x))*(beselj(x,n+1/2.d))

b1 = where(finite(real_part(j_n),/nan),/null)
if b1 ne !NULL then begin $
      j_n(b1) = 0. & $
endif & $

return,j_n
end
