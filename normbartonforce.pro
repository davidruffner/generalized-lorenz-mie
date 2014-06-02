;+
; NAME:
;       bartonforce
;
; PURPOSE:
;       Calculates the force on a sphere from a general beam
;
; CATEGORY:
;       Lorentz Mie theory
;
; CALLING SEQUENCE:
;       force = normbartonforce(bscs,a,np,nm, lambda)
;
; INPUTS:
;
;       bscs: beam shape coefficients (jackson conventions)
;
;       a: radius of sphere [micrometers]
;
;       np: (complex) refractive index of sphere
;
;       nm: (complex) refractive index of medium
;
;       lambda: wavelength of light [micrometers]
;
; KEYWORD FLAGS:
;
; OUTPUTS:
;       force: force on sphere
;
; EXAMPLE:
;
;
;
; REFERENCE:
;   1. Adapted from Barton 1989
;   2. C. F. Bohren and D. R. Huffman, 
;      Absorption and Scattering of Light by Small Particles,
;      (New York, Wiley, 1983).
;   3. Jackson
;
;  Mathematica Notebooks:
;        BartonForceSI_Jackson.nb
;        ConvertingBartontoJackson.nb
;
; MODIFICATION HISTORY:
; Written by David B. Ruffner, New York University, 2014/04/04
; 2014/04/16 Fixed a bug with nmax when only dipole term. There was a
; problem with nmax-1<nc-1 instead of (nmax-1)<(nc-1), David
; B. Ruffner
  
;
; Copyright (c) 2014-2019 David B. Ruffner and David G. Grier
;-

function normbartonforce,bscs,a,np,nm, lambda,ab=ab

;particle coeffiecients
if n_elements(ab) eq 0 then begin
   ab = sphere_coefficients(a,np,nm,lambda)
endif
an = ab[0,*]
bn = ab[1,*]
nc = n_elements(an)-1

;Beam shape coefficients
szb = size(bscs)
nmax = szb[1]-1
;print,"nmax0",nmax
;Pad the am_mn and ae_mn matrices with zeros to fix bug
if nmax eq 1 then begin
   ;print,"here0"
   bscs2 = complex(fltarr(szb[1]+1,szb[2],szb[3]))
   bscs2[0:szb[1]-1,*,*] = bscs 
   bscs = bscs2
endif

am_mn = bscs[*,*,0]
ae_mn = bscs[*,*,1]

am_mnsca = am_mn

;speed of light
c =  299792458.d;m/s
;wavevector
k = 2*!pi*nm/lambda

szp = size(am_mn)
nmax = szp[1]-1


   

;XY components of the force. Barton Eq. 5. Converted to notation of 
;Jackson and SI units using BartonFxySI_Jackson.nb
Fxytotal = complex(0.d,0.d)
;print,"Testing",nc,nmax,(nmax-1) < (nc-1)
;print, "limit n",(nmax-1) < (nc-1)
for n = 1.d, (nmax-1) < (nc-1) do begin
   for m = -n,n do begin
      ;Get factor 
      lmfactor1pm = sqrt(n*(double(n) + 2.d))*sqrt((double(n) + m + 1.d)* $
                      (double(n) + m + 2.d)/( (2.d*double(n) + 3.d)* $
                      (2.d*double(n) + 1.d)))/(double(n) + 1.d)
      lmfactor1mm = sqrt(n*(double(n)+2.))*sqrt((double(n) - m + 1.d)* $
                      (double(n) - m + 2.d)/( (2.d*double(n) + 3.d)* $
                      (2.d*double(n) + 1.d)))/(double(n) + 1.d)
      lmfactor2 = -sqrt((double(n)+m+1)*(double(n)-m))/(n*(double(n)+1))

      ;Extinction part of force
      Fxyinner1eA = lmfactor1pm*( $
                       (-an[n]*ae_mn[n,m])*conj(ae_mn[n + 1, m + 1]) + $
                       (-bn[n]*am_mn[n,m])*conj(am_mn[n + 1, m + 1]))+ $
                    lmfactor1mm*( $
                       conj(-an[n]*ae_mn[n,m])*ae_mn[n + 1, m - 1] + $
                       conj(-bn[n]*am_mn[n,m])*am_mn[n + 1, m - 1])
      Fxyinner1eB = lmfactor1pm*( $
                       (ae_mn[n,m])*conj(-an[n+1]*ae_mn[n + 1, m + 1]) + $
                       (am_mn[n,m])*conj(-bn[n+1]*am_mn[n + 1, m + 1]))+ $
                    lmfactor1mm*( $
                       conj(ae_mn[n,m])*(-an[n+1]*ae_mn[n + 1, m - 1]) + $
                       conj(am_mn[n,m])*(-bn[n+1]*am_mn[n + 1, m - 1]))
      Fxyinner2e = -lmfactor2*( $
                       -(-an[n]*ae_mn[n,m])*conj(am_mn[n,m+1])+ $
                       conj(-an[n]*ae_mn[n,m+1])*am_mn[n,m]+ $
                       (-bn[n]*am_mn[n,m])*conj(ae_mn[n,m+1])- $
                       conj(-bn[n]*am_mn[n,m+1])*ae_mn[n,m])

      ;Scattering part of the force
      Fxyinner1s = lmfactor1pm*( $
                    2*(-an[n]*ae_mn[n,m])*conj(-an[n+1]*ae_mn[n+1,m+1]) + $
                    2*(-bn[n]*am_mn[n,m])*conj(-bn[n+1]*am_mn[n+1,m+1])) + $
                   lmfactor1mm*( $
                    2*conj((-an[n]*ae_mn[n,m]))*(-an[n+1]*ae_mn[n+1,m-1]) + $
                    2*conj(-bn[n]*am_mn[n,m])*(-bn[n+1]*am_mn[n+1,m-1]))

      Fxyinner2s = -lmfactor2*( $
                       -2*(-an[n]*ae_mn[n,m])*conj(-bn[n]*am_mn[n,m+1]) $
                       +2*(-bn[n]*am_mn[n,m])*conj(-an[n]*ae_mn[n,m+1]))

      Fxyext = (Fxyinner1eA+Fxyinner1eB+Fxyinner2e)/(2*!pi*(k*a)^2) 
      Fxysca = (Fxyinner1s+Fxyinner2s)/(2*!pi*(k*a)^2) 
      ;; print,"Fxyext",Fxyext,"lmfactor1",lmfactor1pm
      ;; print,"an ",an[n]," bn ",bn[n]
      ;; print,"an+1 ",an[n+1]," bn+1 ",bn[n+1]
      ;; print,"n,m ",string(n,format='(I01)'),", "+string(m,format='(I01)'),$
      ;;    " fxpart ",real_part(Fxyext+Fxysca),$
      ;;    " fypart ",imaginary(Fxyext+Fxysca)
      
      Fxytotal = Fxytotal + complex(0.,-1.d)*(Fxyext+Fxysca);/(4*!pi*1.051)

   endfor

endfor
;print,"here"
Fxtotal = real_part(Fxytotal)
Fytotal = imaginary(Fxytotal)


;Z component of the force. Barton Eq. 6. Converted to notation of 
;Jackson and SI units using BartonForceSI_Jackson.nb
Fztotal = 0.d
for n = 1.d, (nmax-1) < (nc-1) do begin
   for m = -n,n do begin
      ;Get factor 
      lmfactor1 = sqrt(n*(double(n)+2.))*sqrt((double(n) - m + 1.d)* $
                      (double(n) + m + 1.d)/( (2.d*double(n) + 3.d)* $
                      (2.d*double(n) + 1.d)))/(double(n) + 1.d)
      lmfactor2 = -m/(n*(double(n)+1))

      ;Extinction part of force
      Fzinner1ab = lmfactor1*(ae_mn[n+1,m]*conj(-an[n]*ae_mn[n,m]) + $
                             am_mn[n+1,m]*conj(-bn[n]*am_mn[n,m]) + $
                             (-an[n+1]*ae_mn[n+1,m])*conj(ae_mn[n,m]) + $
                             (-bn[n+1]*am_mn[n+1,m])*conj(am_mn[n,m]))
      Fzinner2a = lmfactor2*((-an[n]*ae_mn[n,m])*conj(am_mn[n,m])+ $
                             ae_mn[n,m]*conj(-bn[n]*am_mn[n,m]))
      Fzext = imaginary(Fzinner1ab+Fzinner2a)/(!pi*(k*a)^2) 

      ;Scattering part of the force
      Fzinner1c = lmfactor1*( $
                     2*(-an[n+1]*ae_mn[n+1,m])*conj(-an[n]*ae_mn[n,m])+ $
                     2*(-bn[n+1]*am_mn[n+1,m])*conj(-bn[n]*am_mn[n,m]))

      Fzinner2b = lmfactor2*(2*(-an[n]*ae_mn[n,m])*conj(-bn[n]*am_mn[n,m]) )
      
      Fzsca = imaginary(Fzinner1c + Fzinner2b)/(!pi*(k*a)^2) 

      Fztotal = Fztotal-(Fzsca+Fzext)

   endfor

endfor



return,[Fxtotal,Fytotal,Fztotal]
end
