;+
;NAME:
;    dgg_pi_n
;
; PURPOSE:
;    Calculate the angular function pi_n from Lorentz Mie theory
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    out = dgg_pi_n(costheta,l,m)
;
;INPUTS:
;    costheta:  cosine of polar angle  
;
;    l:      total angular momentum number
;
;    m:      azimuthal angular momentum number
;
;OUTPUTS:
;    out: value of the angular function
;
;DEPENDENCY:
;
;REFERENCE:
;   1. W. J. Wiscombe,
;      Improved Mie scattering algorithms,
;      Applied Optics 19, 1505-1509 (1980).
;
;    Press W. H. et al. "Numerical Recipes 3rd ed."
;    Cambridge University Press (2007)
;
;MODIFICATION HISTORY:
; 03/19/2014 Adapted from sphericalfield.pro by David B. Ruffner,
;            New York University


function dgg_pi_n,costheta, nc

b = where(abs(costheta) gt 1.0,/null)
if b ne !NULL then begin
   message,"Bad arguments in routine dgg_pi_n (2)"
    return, -1
 endif

ci = dcomplex(0,1)

; ... angular functions (4.47), page 95
pi_nm1 = 0.d                    ; \pi_0(\cos\theta)
pi_n   = 1.d                    ; \pi_1(\cos\theta)

; Compute functions by summing multipole contributions
for n = 1.d, nc do begin

 ; upward recurrences ...
; ... Legendre factor (4.47)
; Method described by Wiscombe (1980)
    swisc = pi_n * costheta 
    twisc = swisc - pi_nm1
    tau_n = n * twisc - pi_nm1  ; \tau_n(\cos\theta)  

; upward recurrences ...
; ... angular functions (4.47)
; Method described by Wiscombe (1980)
    pi_nm1 = pi_n
    pi_n = swisc + (n + 1.d) * twisc / n
endfor
return,[[pi_nm1],[tau_n]]
end
  
