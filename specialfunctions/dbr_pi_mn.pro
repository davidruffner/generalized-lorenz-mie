;+
;NAME:
;    dbr_pi_mn
;
; PURPOSE:
;    Calculate the angular function pi_mn from Lorentz Mie theory
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    out = dbr_pi_mn(costheta,l,m)
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
;    Calls dbr_plegendre.pro
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
; 03/24/2014 if m=0 then return zeros


function dbr_pi_mn,x, n, minput

parityfactor = 1.
if minput lt 0 then begin
   parityfactor = -(-1)^(-minput)
   m = -minput
endif else m = minput

if m gt n then begin
    ;message,"Bad arguments in routine dbr_pi_mn (1)"
    return, 0.d
 endif

b = where(abs(x) gt 1.0,/null)
if b ne !NULL then begin
   message,"Bad arguments in routine dbr_pi_mn (2)"
    return, -1
 endif

if m eq 0 then begin
   return,double(fltarr(n_elements(x)))
endif

;Get starting values for recursion
pi_mm=1.0d
if m gt 0 then begin
   omx2 = (1.0d - x)*(1.0d + x)
   fact = 1.0d + double( fltarr(n_elements(x)) )
   pi_mm *= fact/(fact+1.0d) ;This accounts for dividing by sin theta
   fact += 2.0d
   for i=1,m-1 do begin
      pi_mm *= omx2*fact/(fact+1.0d)
      fact += 2.0d
   endfor
endif else pi_mm=1/((1.0d - x)*(1.0d + x))

    
pi_mm = sqrt((2.d*m+1.d)*pi_mm/(4.0d*!pi))
pi_mm = pi_mm*m
  
if m mod 2 eq 1 then begin
    pi_mm=-pi_mm
endif

if n eq m then begin
    ;print,"Done return directly"
    return, pi_mm*parityfactor
 endif else begin
    pi_mmp1 = x*sqrt(2.0d*m+3.0d)*pi_mm
    if n eq (m+1) then begin
      return, pi_mmp1*parityfactor
    endif else begin
      oldfact = sqrt(2.0d*m+3.0d)
      for nn = m+2,n do begin
        fact = sqrt((4.0d*nn*double(nn)-1.0d)/(nn*double(nn)-m*double(m)))
	pi_nn = (x*pi_mmp1-pi_mm/oldfact)*fact
        ;print,"pnn",pnn," oldfact",oldfact," fact",fact," nn",nn," m",m
        oldfact = fact
	pi_mm = pi_mmp1
	pi_mmp1 = pi_nn
      endfor
      return, pi_nn*parityfactor
    endelse
	
 endelse
end
  
