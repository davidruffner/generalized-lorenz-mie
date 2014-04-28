;+
; NAME:
;       planewaveforce
;
; PURPOSE:
;       Calculates the normalized force on a sphere from a plane wave
;       in the direction of the wavevector
;
; CATEGORY:
;       Lorentz Mie theory
;
; CALLING SEQUENCE:
;       force = planewaveforce(a,np,nm, lambda)
;
; INPUTS:
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
;       force: normalized force on sphere from plane wave
;
; EXAMPLE:
;
; I: Intensity of incident wave (I=Power/area)
; a: radius of sphere
; c: speed of light
; 
; Force = planewaveforce(a,np,nm,lambda)*!pi*a^2*Power/(c*area)
;
; Say P = 20mW = 2*10^(10) pN*m/s, and area = (.2mm)^2 = 4*10^4.um^2,
; a=1.5um, and of course c = 3.*10^(8)m/s. Assume the particle
; refractive index is np=1.45, and the medium's refractive index
; nm = 1.33, and the light wavelength is lambda = .447um
;
;IDL> power = 2.*10.^(10.) & area = 40000. & a = 1.5 & c = 3*10.^8.
;IDL> np = 1.45 & nm = 1.33 & lambda = .447
;IDL> force = planewaveforce(a,np,nm,lambda)*!pi*a^2*Power/(c*area);
;                    ; in pN 
;
;
; REFERENCE:
;   1. Adapted from Chapter 4 in
;      C. F. Bohren and D. R. Huffman, 
;      Absorption and Scattering of Light by Small Particles,
;      (New York, Wiley, 1983).
;
; MODIFICATION HISTORY:
; Written by David B. Ruffner, New York University, 2014/03/31
;
; Copyright (c) 2014-2019 David B. Ruffner and David G. Grier
;-

function planewaveforce,a,np,nm, lambda,ab=ab
;a=8.9 & np =1.45 & nm=1.33 & lambda=.447 

if n_elements(ab) eq 0 then begin
   ab = sphere_coefficients(a,np,nm,lambda)
endif

nc = n_elements(ab[0,*])

k = 2.d*!pi*nm/lambda
x = k*a

n = dindgen(nc-2) + 1.d ;Don't use zero to avoid nan
;print,abs(ab)
;print,[transpose(n),abs(ab[0,1:nc-2]),abs(ab[0,2:nc-1])]

;Calculate the extinction efficiency
Qext_n1 = (2*n + 1.d)*real_part(ab[0,1:nc-2]+ab[1,1:nc-2]) ;Bohren Eq.~(4.62)
Qext = (2.d/x^2)*total(Qext_n1)

;Calculate the scattering efficiency
Qsca_n1 = (2*n+1)*(abs(ab[0,1:nc-2])^2.+abs(ab[1,1:nc-2])^2.) ;Bohren Eq.~(4.61)
Qsca = (2.d/x^2)*total(Qsca_n1)
print,"Qext",Qext
print,"Qsca",Qsca
;Print,"They are the same so I guess there is no absorption"

;Calculate the asymmetry parameter p. 120 Bohren
QscaCosTh_n1 = n*(n + 2.d)*real_part(ab[0,1:nc-2]*conj(ab[0,2:nc-1]) $
                                   + ab[1,1:nc-2]*conj(ab[1,2:nc-1]))/(n + 1.d)

QscaCosTh_n2 = (2.d*n + 1.d)*real_part(ab[0,1:nc-2]*conj(ab[1,1:nc-2]))/$
                                                          (n*(n+1))
;help,QscaCosTh_n2
;help,QscaCosTh_n1

QscaCosTh = (total(QscaCosTh_n1+QscaCosTh_n2))*(4.d/x^2)

;Calculate the normalized force
;print,"Qext",Qext
;print,"QscaCosTh",QscaCosTh
force =  Qext - QscaCosTh
;print,"asymmetry parameter",QscaCosTh/Qsca



return,force
end
