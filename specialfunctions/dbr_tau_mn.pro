;+
;NAME:
;    dbr_tau_mn
;
; PURPOSE:
;    Calculate the angular function tau_mn from Lorentz Mie theory
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    out = dbr_tau_mn(costheta,l,m)
;
;INPUTS:
;    costheta:  cosine of polar angle  
;
;    phi:    azimuthal angle
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
;
;REFERENCE:
;    Press W. H. et al. "Numerical Recipes 3rd ed."
;    Cambridge University Press (2007)
;
;    Wikipedia
;
;    Taylor and Love. "Multipole Expansion of Bessel and
;     Gaussian beams for Mie scattering
;     calculation". J. Opt. Soc. Am. A. 26, 278 (2009)
;
;MODIFICATION HISTORY:
; 03/18/2014 Written by David B. Ruffner, New York University
; 03/25/2014 Fixed a bug when m=0. Now it doesn't blow up
; 06/04/2014 DBR: commented out the low precision at poles message


function dbr_tau_mn,x, n, minput

tol = 1.e-12

parityfactor = 1.
if minput lt 0 then begin
   parityfactor = (-1)^(-minput)
   m = -minput
endif else m = minput

if m gt n then begin
    message,"Bad arguments in routine dbr_tau_mn"
    return, -1
 endif

b = where(abs(x) gt 1.0)
if b ne -1 then begin
   message,"Bad arguments in routine dbr_tau_mn"
    return, -1
endif

norm = sqrt((2*n + 1)*(n - m)/((2*n - 1)*(n + m)))

;; if m eq 0 then begin
;;    tau_mn = (n + m)*norm*dbr_plegendre(x,n-1,m)/sqrt(1 - x^2)-$
;;                  n*x*dbr_plegendre(x,n,m)/sqrt(1 - x^2)
;; endif else begin
;;    tau_mn = (n + m)*norm*dbr_pi_mn(x,n-1,m)/m - $
;;                  n*x*dbr_pi_mn(x,n,m)/m
;; endelse
normfactor = Sqrt((2.d*double(n)+ 3.d)*(n+ 1.d -m)/$
                          ((2.d*double(n) + 1.d)*(n+ 1.d +m)))
if m eq 0 then begin
   sintheta = sqrt(1-x^2)
   tau_mn = parityfactor*((n + 1.d - m)*dbr_plegendre(x,n+1,m)/normfactor - $
               (n + 1.d)*x*dbr_plegendre(x,n,m))/sintheta
   b = where(abs(sintheta) lt tol,/null)
   if ~(b eq !NULL) then begin
      ;print,"Low precision near pole, setting tau_mn to zero",n,m,sintheta
      tau_mn(b) = 0
   endif
endif else begin
   tau_mn = parityfactor*((n + 1.d - m)*dbr_pi_mn(x,n+1,m)/normfactor - $
               (n + 1.d)*x*dbr_pi_mn(x,n,m))/m
endelse

return,tau_mn

end
  
