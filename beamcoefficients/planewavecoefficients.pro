;+
;NAME:
;    planewaveecoefficients
;
; PURPOSE:
;    Calculate the beam shape coefficients of a plane wave beam
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    E = planewaveecoefficients(rvec,nmax)
;
;INPUTS:
;
;    nmax:  max value of n in the vector spherical harmonic series expansion
;
;OUTPUTS:
;    AB:    A [nmax,2*nmax+1,2] matrix of the beam shape coefficients
;
;REFERENCE:
;
;MODIFICATION HISTORY:
; 03/06/2010 Written by David B. Ruffner, New York University

function planewavecoefficients,nmax

am_mn = fltarr(nmax+1,2*nmax+1)*complex(1,0)
ae_mn = fltarr(nmax+1,2*nmax+1)*complex(1,0)
ci = complex(0,1)
for n=1,nmax do begin
   for m = -n,n do begin
      ;print,n,m
      if m eq 1 then begin
            am_mn[n,m] = -(ci^n)*sqrt(4.d*!pi*(2.d*n + 1.d))/2.d
            ae_mn[n,m] = -am_mn[n,m]*ci 
      endif
      if m eq -1 then begin 
            am_mn[n,m] = -(ci^n)*sqrt(4.d*!pi*(2.d*n + 1.d))/2.d
            ae_mn[n,m] = am_mn[n,m]*ci 
      endif     
   endfor
endfor

b1 = where(finite([[[am_mn]],[[ae_mn]]],/nan),/null)
if b1 ne !NULL then begin $
      print, "We have infinities"
      stop
endif & $

return,[[[am_mn]],[[ae_mn]]]
end
