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

b = where(x eq 0,complement=bc,/null)
jn = beselj(x,n +  1/2.d);temporary
if b ne !NULL then begin
   jn(b) = 0
   jn(bc) = sqrt(!pi/(2.d*x(bc)))*jn
endif else jn = sqrt(!pi/(2.d*x))*jn


psinp1 = sqrt(!pi*x/2.)*beselj(x,n + 1.d + 1/2.d)
psinm1 = sqrt(!pi*x/2.)*beselj(x,n - 1.d + 1/2.d)

;Fix bug which causes beselj to give NANs
b1 = where(finite(jn,/nan),/null)
b2 = where(finite(psinp1,/nan),/null)
if b1 ne !NULL then begin 
   jn(b1) = 0 
endif 
if b2 ne !NULL then begin 
   psinp1(b2) = 0
endif 
b3 = where(finite(psinm1,/nan),/null)
if b3 ne !NULL then begin 
   psinm1(b3) = 0 
endif 
;from dlmf.nist.gov/10.51
;j_nprime(x) = -j_n+1(x)+(n/x)j_n(x)
;Checked this recurrence relation with Numerical recipes
psiprime = -psinp1+(1.d + n)*jn

psiprime2 = psinm1-(n)*jn

psiprimeavg = (psiprime + psiprime2)/2.d

return,psiprimeavg
end
