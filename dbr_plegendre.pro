;+
;NAME:
;    dbr_plegendre
;
; PURPOSE:
;    Calculate the normalized associated legendre polynomial
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    pnm = dbr_plegendre(n,m,x)
;
;INPUTS:
;    n:    index
;
;    m:    index
;
;    x:    input
;
;OUTPUTS:
;    pnm: value of the nomalized associated legendre polynomial
;
;REFERENCE:
;    Press W. H. et al. "Numerical Recipes 3rd ed."
;    Cambridge University Press (2007)
;
;MODIFICATION HISTORY:
; 03/07/2010 Written by David B. Ruffner, New York University


function dbr_plegendre,x, n, minput
;computes associated legendre polynomial Plm(x) 
parityfactor = 1.
if minput lt 0 then begin
   parityfactor = (-1)^(-minput)
   m = -minput
endif else m = minput

if m gt n then begin
    ;message,"Bad arguments in routine dbr_plegendre"
    return, 0.d
 endif

b = where(abs(x) gt 1.0)
if b ne -1 then begin
   message,"Bad arguments in routine dbr_plegendre"
    return, -1
endif


pmm=1.0d
if m gt 0 then begin
   omx2 = (1.0d - x)*(1.0d + x)
   fact = 1.0d
   for i=0,m-1 do begin
      pmm *= omx2*fact/(fact+1.0d)
      fact += 2.0d
   endfor
endif
    
pmm = sqrt((2.d*m+1.d)*pmm/(4.0d*!pi))
  
if m mod 2 eq 1 then begin
    pmm=-pmm
endif
if n eq m then begin
    ;print,"Done return directly"
    return, pmm*parityfactor
 endif else begin
    pmmp1 = x*sqrt(2.0d*m+3.0d)*pmm
    if n eq (m+1) then begin
      return, pmmp1*parityfactor
    endif else begin
      oldfact = sqrt(2.0d*m+3.0d)
      for nn = m+2,n do begin
        fact = sqrt((4.0d*nn*double(nn)-1.0d)/(nn*double(nn)-m*double(m)))
	pnn = (x*pmmp1-pmm/oldfact)*fact
        ;print,"pnn",pnn," oldfact",oldfact," fact",fact," nn",nn," m",m
        oldfact = fact
	pmm = pmmp1
	pmmp1 = pnn
      endfor
      return, pnn*parityfactor
    endelse
	
 endelse
end
  
