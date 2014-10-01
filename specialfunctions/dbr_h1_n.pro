;+
;NAME:
;    dbr_h1_n
;
; PURPOSE:
;    Calculate the first spherical hankel function h1_n(x)
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    out = dbr_h1_n(x,n)
;
;INPUTS:
;    n:    index
;
;    x:    input
;
;OUTPUTS:
;    out: value of the first spherical hankel function of order n at x
;
;DEPENDENCY:
;    Calls idl's BESELJ and BESELY functions.
;
;MODIFICATION HISTORY:
; 10/01/2014 Written by David B. Ruffner, New York University

function dbr_h1_n,x,n

h1_n = sqrt(!pi/(2.d*x))*(beselj(x,n+1/2.d)+complex(0,1)*besely(x,n+1/2.d))

b1 = where(finite(real_part(h1_n),/nan),/null)
if b1 ne !NULL then begin $
      h1_n(b1) = complex(0,0) & $
endif & $

return,h1_n
end
