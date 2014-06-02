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
;       force = bartonforce(bscs,a,np,nm, lambda)
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
;
; Copyright (c) 2014-2019 David B. Ruffner and David G. Grier
;-

function normbartonforce,bscs,a,np,nm, lambda, 

;particle coeffiecients
ab = sphere_coefficients(a,np,nm,lambda)
an = ab[0,*]
bn = ab[1,*]
nc = n_elements(an)-1

;Beam shape coefficients
am_mn = bscs[*,*,0]
ae_mn = bscs[*,*,1]

am_mnsca = am_mn

;speed of light
c =  299792458.d;m/s
;wavevector
k = 2*!pi*nm/lambda

szp = size(am_mn)
nmax = szp[1]-1

;Z component of the force. Barton Eq. 6. Converted to notation of 
;Jackson and SI units using BartonForceSI_Jackson.nb
Fztotal = 0.d
for n = 1.d, nmax-1 < nc-1 do begin
   for m = -n,n do begin
      ;Get factor 
      lmfactor1 = sqrt(n*(double(n)+1))*sqrt((double(n) - m + 1.d)* $
                      (double(n) + m + 1.d)/( (2.d*double(n) + 3.d)* $
                      (2.d*double(n) + 1.d)))/(double(n) + 1.d)
      lmfactor2 = -m/(n*(double(n)+1))

      ;Extinction part of force
      Fzinner1ab = lmfactor1*(ae_mn[n+1,m]*conj(-an[n]*ae_mn[n,m]) + $
                             am_mn[n+1,m]*conj(-bn[n]*am_mn[n,m]) + $
                             (-an[n+1]*ae_mn[n+1,m])*conj(ae_mn[n,m]) + $
                             (-bn[n+1]*am_mn[n+1,m])*conj(am_mn[n,m])
      Fzinner2a = lmfactor2*((-an[n]*ae_mn[n,m])*conj(am_mn[n,m])+ $
                             ae_mn[n,m]*conj(-bn[n]*am_mn[n,m]))
      Fzext = imaginary(Fzinner1ab+Fzinner2a)/(!pi*(k*a)^2) 

      ;Scattering part of the force
      Fzinner1c = lmfactor1*(
                     2*(-an[n+1]*ae_mn[n,m])*conj(-an[n]*ae_mn[n,m])+ $
                     2*(-bn[n+1]*am_mn[n,m])*conj(-bn[n]*am_mn[n,m]))

      Fzinner2b = lmfactor2*(2*(-an[n]*ae_mn[n,m])*conj(-bn[n]*am_mn[n,m]) )
      
      Fzsca = imaginary(Fzinner1c + Fzinner2b)/(!pi*(k*a)^2) 

      Fztotal = Fztotal+Fzsca+Fzext

   endfor

endfor



return,Fztotal
end
