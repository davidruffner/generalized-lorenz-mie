;+
;NAME:
;    efieldmn
;
; PURPOSE:
;    Calculate the field with vector spherical harmonic coefficient m and n
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    Enm = efieldmn(r,pmn,qmn,m,n)
;
;INPUTS:
;    r:    [3,N] array of positions where you want the field
;
;    pmn:    coefficient of vector spherical harmonic
;
;    qmn:    other coefficient of vector spherical harmonic
;
;    m,n:    indicies of vector spherical harmonic
;
;OUTPUTS:
;    E:   [3,N] Electric field vector at each point.
;
;REFERENCE:
;    Barton, J. P. et al. "Internal and near-surface electromagnetic 
;     fields for a spherical particle irradiated by a focused
;     laser beam". J. Appl. Phys. 64, 1632 (1988)
;
;    Jackson

;
;MODIFICATION HISTORY:
; 03/06/2014 Written by David B. Ruffner, New York University
; 03/14/2014 Calculate psi psi' Ylm and Ylm' with custum functions
; 03/18/2014 Now use angular functions 

function efieldmn2, rvec, pmn,qmn,m,n,k,a,nm

sz = size(rvec)
if sz[0] eq 1 then begin
   x = rvec[0]
   y = rvec[1]
   z = rvec[2]
endif else begin
   x = rvec[0,*]
   y = rvec[1,*]
   z = rvec[2,*]
endelse
rho   = sqrt(x^2 + y^2)
r = sqrt(x^2+y^2+z^2)
theta = atan(rho, z)
phi   = atan(y, x)
costheta = cos(theta)
sintheta = sin(theta)
cosphi = cos(phi)
sinphi = sin(phi)


im = complex(0,1)
expimphi = exp(im*m*phi)

kr = k*r

;Calculate special functions

;Ricatti bessel functions
psin = dbr_riccatibessel(kr,n)
psiprime = dbr_riccatibesselprime(kr,n)
bm2 = where(finite(psin,/nan),/null)
if bm2 ne !NULL then begin 
   print,"we have nansss! psin in efieldmn2"
   stop
endif

bm1 = where(finite(psiprime,/nan),/null)
if bm1 ne !NULL then begin 
   print,"we have nansss! psiprime in efieldmn2"
   stop
endif

;spherical harmonic functions
ylm = dbr_sphericalharmonic(theta,phi,n,m)
b0 = where(finite(ylm,/nan),/null)
if b0 ne !NULL then begin 
   print,"we have nansss! ylm in efieldmn2"
   stop
endif

;angular functions

pi_mn = dbr_pi_mn(costheta,n,m)
tau_mn = dbr_tau_mn(costheta,n,m)


;Convert to Barton notation for coefficients

;; Amn = im*pmn/(2*!pi*k*a)

;; Bmn = nm*im*qmn/(2*!pi*k*a)
Amn = pmn
Bmn = qmn

;Calculate the field

;using field from Jackson

Ermn = -Amn*sqrt(n*(double(n)+1))*psin*ylm/(kr^2)

Ethmn = -Amn*psiprime*tau_mn*expimphi/(kr*sqrt(n*(double(n)+1))) $
        -Bmn*psin*pi_mn*expimphi/(kr*sqrt(n*(double(n)+1)))

Ephimn = -Amn*im*psiprime*pi_mn*expimphi/(kr*sqrt(n*(double(n)+1))) $
         -Bmn*im*psin*tau_mn*expimphi/(kr*sqrt(n*(double(n)+1)))
;; Ermn = (1/r^2)*n*(n+1)*Amn*psin*ylm
;; ;Ermn = (1/r^2)*Amn*psin*ylm

;; Ethmn = (1/r)*(k*Amn*psiprime*dylmdtheta+$
;;                 im*m*k*Bmn*psin*ylm/sintheta)
;; Ephimn = (1/r)*(im*m*k*Amn*psiprime*ylm/sintheta-$
;;                 k*Bmn*psin*dylmdtheta)
b1 = where(finite(Ermn,/nan),/null)
if b1 ne !NULL then begin 
         Ermn(b1) = 0
endif
b2 = where(finite(Ethmn,/nan),/null)
if b2 ne !NULL then begin 
         Ethmn(b2) = 0
endif
b3 = where(finite(Ephimn,/nan),/null)
if b3 ne !NULL then begin 
         Ephimn(b3) = 0
endif

E = [Ermn,Ethmn,Ephimn]

return,E
end
