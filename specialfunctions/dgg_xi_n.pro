;+
;NAME:
;    dgg_xi_n
;
; PURPOSE:
;    Calculate the riccati-bessel function associated with the
;    spherical bessel, h^(1)_n(x)
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    out = dgg_xi_n(x,n)
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

function dgg_xi_n,x,nc

; starting points for recursive function evaluation ...
; ... Riccati-Bessel radial functions, page 478
sinx = sin(x)
cosx = cos(x)
xi_nm2 = dcomplex(cosx, sinx) ; \xi_{-1}(x)
xi_nm1 = dcomplex(sinx,-cosx) ; \xi_0(x)

; Compute field by summing multipole contributions
for n = 1.d, nc do begin
; upward recurrences ...

; ... Riccati-Bessel function, page 478
    xi_n   = (2.d*n - 1.d) * xi_nm1 / x - xi_nm2    ; \xi_n(x)
; ...and it's derivative
    dn = (n * xi_n)/x - xi_nm1

; ... Riccati-Bessel function
    xi_nm2 = xi_nm1
    xi_nm1 = xi_n
 endfor

return,[[xi_nm1],[dn]]
end
