;+
;NAME:
;    dbr_xi_nprime
;
; PURPOSE:
;    Calculate the derivative of the riccati-bessel function
;    associated with the spherical bessel, jn(x)
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    out = dbr_xi_nprime(x,n)
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
; 03/19/2010 Written by David B. Ruffner, New York University

function dbr_xi_nprime,x,n

b = where(x eq 0,complement=bc,/null)
h1_n = beselj(x,n +  1/2.d)+complex(0,1)*besely(x,n +  1/2.d);temporary
if b ne !NULL then begin
   h1_n(b) = complex(0,0)
   h1_n(bc) = sqrt(!pi/(2.d*x(bc)))*h1_n
endif else h1_n = sqrt(!pi/(2.d*x))*h1_n


xi_np1 = sqrt(!pi*x/2.)*(beselj(x,n + 1.d + 1/2.d)+$
                             complex(0,1)*besely(x,n + 1.d + 1/2.d))
xi_nm1 = sqrt(!pi*x/2.)*(beselj(x,n - 1.d + 1/2.d)+$
                             complex(0,1)*besely(x,n - 1.d + 1/2.d))

;Fix bug which causes beselj to give NANs
b1 = where(finite(h1_n,/nan),/null)
b2 = where(finite(xi_np1,/nan),/null)
if b1 ne !NULL then begin 
   h1_n(b1) = complex(0,0) 
endif 
if b2 ne !NULL then begin 
   xi_np1(b2) = complex(0,0) 
endif 
b3 = where(finite(xi_nm1,/nan),/null)
if b3 ne !NULL then begin 
   xi_nm1(b3) = complex(0,0) 
endif 
;from dlmf.nist.gov/10.51
;j_nprime(x) = -j_n+1(x)+(n/x)j_n(x)
;Checked this recurrence relation with Numerical recipes
xi_nprime = -xi_np1+(1.d + n)*h1_n

xi_nprime2 = xi_nm1-(n)*h1_n

xi_nprimeavg = (xi_nprime + xi_nprime2)/2.d

return,xi_nprimeavg
end
