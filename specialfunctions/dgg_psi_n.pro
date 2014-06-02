;+
;NAME:
;    dgg_psi_n
;
; PURPOSE:
;    Calculate the riccati-bessel function associated with the
;    spherical bessel, j_n(x)
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    out = dgg_psi_n(x,n)
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
;    Calls idl's BESELJ function.
;
;MODIFICATION HISTORY:
; 03/19/2014  Adapted from sphericalfield by David B. Ruffner,
;             New York University

function dgg_psi_n,x,nc

; starting points for recursive function evaluation ...
; ... Riccati-Bessel radial functions, page 478
sinx = sin(x)
cosx = cos(x)
psi_nm2 = cosx ; \psi_{-1}(x)
psi_nm1 = sinx ; \psi_0(x)

; Compute field by summing multipole contributions
for n = 1.d, nc do begin
; upward recurrences ...

; ... Riccati-Bessel function, page 478
    psi_n   = (2.d*n - 1.d) * psi_nm1 / x - psi_nm2    ; \psi_n(x)
; ...and it's derivative
    dn = (n * psi_n)/x - psi_nm1

; ... Riccati-Bessel function
    psi_nm2 = psi_nm1
    psi_nm1 = psi_n
 endfor

return,[[psi_nm1],[dn]]
end
