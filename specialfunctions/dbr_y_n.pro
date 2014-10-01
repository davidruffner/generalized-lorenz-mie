;+
;NAME:
;    dbr_y_n
;
; PURPOSE:
;    Calculate the spherical bessel function of second kind y_n(x)
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    out = dbr_y_n(x,n)
;
;INPUTS:
;    n:    index
;
;    x:    input
;
;OUTPUTS:
;    out: value of the spherical bessel function of second kind of order n at x
;
;DEPENDENCY:
;    Calls idl's BESELY  functions.
;
;MODIFICATION HISTORY:
; 10/01/2014 Written by David B. Ruffner, New York University

function dbr_y_n,x,n

y_n = sqrt(!pi/(2.d*x))*(besely(x,n+1/2.d))

b1 = where(finite(real_part(y_n),/nan),/null)
if b1 ne !NULL then begin $
      y_n(b1) = 0. & $
endif & $

return,y_n
end
