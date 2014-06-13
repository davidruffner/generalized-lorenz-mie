;+
;NAME:
;    GAUSSIANTRAPCOEFFICIENTS
;
; PURPOSE:
;    Calculate the coefficients for the vector spherical harmonics of
;    a gaussian beam as a function of position in the field
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    {p_mn,q_mn} = gaussiantrapcoefficients(pos,nmax,k,gamma,thetaG,NT)
;
;INPUTS:
;    pos:    array of positions where you want the coefficients
;
;    nmax:   maximum value of n in the sum
;
;    k  : Wave number of light
;
;    gamma : Ratio between focal lenght of objective and width of
;            gaussian beam in the objective aperture.
;
;    thetaG : maximum convergence angle of objective
;
;    NT : Number of terms to use in integration (100 is in general a
;         good number)
;
;OUTPUTS:
;    {p_mn,q_mn}: Complex coefficients for the vector spherical
;    harmonics
;
;REFERENCE:
;    Taylor, J. M. and G. D. Love. "Multipole expansion of
;    bessel and gaussian beams for mie scattering calculations"
;    J. Opt. Soc. Am. A. 26, 278(2009)
;
;MODIFICATION HISTORY:
; 06/11/2014 Written by Henrique W. Moyses and David B. Ruffner, New York University

function gaussiantrapcoefficients,pos,nmax,k,gamma,thetaG,NT,verbose=verbose
tol = 0.0001d
if n_elements(verbose) eq 0 then verbose = 0

am_mn = dcomplexarr(nmax+1,2*nmax+1)
ae_mn = dcomplexarr(nmax+1,2*nmax+1)

for n=1,nmax do begin
   for m = -n,n do begin
      ;first find the gaussian coefficient at required points
      am_ae = gaussiantrapcoefficient(n,m,pos,k,gamma,thetaG,NT)
      am_mn[n,m] = am_ae[0]
      ae_mn[n,m] = am_ae[1]
      if abs(am_ae[0]) gt tol and verbose then begin
         print,"m,n",m,n
         print,am_ae
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

