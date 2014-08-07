;+
;NAME:
;    GAUSSIANTRAPCOEFFICIENT
;
; PURPOSE:
;    Calculate the coefficient for the vector spherical harmonics of
;    a gaussian beam as a function of position in the field and m and
;    n. Second part of the code. This part performs the integration
;    over the objective  aperture angles using a tabulated method. 

;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    [p_mn,q_mn] = gaussiantrapcoefficient(n,m,pos,k,gamma,thetaG,NT)
;
;INPUTS:
;    pos:    [3,N] array of positions where you want the coefficients
;
;    m,n:    indicies of the coefficient
;
;    k  : Wave number of light
;
;    gamma : Ratio between focal length of objective and width of
;            gaussian beam in the objective aperture.
; 
;    thetaG : maximum convergence angle of objective
;
;    NT:      number of terms to use in integration (100 is usually a
;             good number)
;
;OUTPUTS:
;    [p_mn,q_mn]: Complex coefficient for the vector spherical
;    harmonics
;
;REFERENCE:
;    Taylor, J. M. and G. D. Love. "Multipole expansion of
;    bessel and gaussian beams for mie scattering calculations"
;    J. Opt. Soc. Am. A. 26, 278(2009)
;
;MODIFICATION HISTORY:
; 06/11/2014 Written by Henrique W. Moyses, New York University
; 06/13/2014 Edited by David B. Ruffner, New York University: combined
;             gaussiantrapcoefficient1.pro and
;gaussiantrapcoefficient2.pro 
; 06/23/2014 DBR: beam waist is also limited by gamma

function gaussiantrapcoefficientint,n,m,pos,k,gamma,thetaG,NT

;define integration step angles
delta = dindgen(NT) * thetaG / (NT-1)

;define arrays for storage of step coefficients
am_mn_temp = dcomplexarr(NT) 
ae_mn_temp = dcomplexarr(NT)

;calculate w0
w0 = SQRT(8.d * !Pi) / (k * sin(thetaG))
;Also gamma can limit the beam waist
if gamma ne 0 then begin
   ;sinthetamax = 1/(sqrt(2)*gamma) 
   sinthetamax = 1/(gamma) 
   w0b = SQRT(8.d * !Pi) / (k * sinthetamax);Not sure about this factor
   w0 = w0 > w0b
endif

;fill up arrays of integrand coefficients
for i=0L, NT-1L do begin
 theta = delta[i]
 ;calculate multiplicative factor for gaussian beam
 fac  =  k^2.d * w0^2.d / (4.d * !dpi)
 fac *= SQRT(COS(theta))
 fac *= EXP(-(gamma * SIN(theta))^2.d)
 fac *= SIN(theta)
 ;Cal besselcoefficient
 am_ae = fac*besselcoefficient(n,m,pos,theta,k)
 am_mn_temp[i] = am_ae[0]
 ae_mn_temp[i] = am_ae[1]
endfor 

;perform integration
am_mn_real = INT_TABULATED(delta,REAL_PART(am_mn_temp))
am_mn_imaginary = INT_TABULATED(delta,IMAGINARY(am_mn_temp))
ae_mn_real = INT_TABULATED(delta,REAL_PART(ae_mn_temp))
ae_mn_imaginary = INT_TABULATED(delta,IMAGINARY(ae_mn_temp))

ci = dcomplex(0.d,1.d)
am_mn = am_mn_real + ci * am_mn_imaginary
ae_mn = ae_mn_real + ci * ae_mn_imaginary

return,[am_mn,ae_mn]
end



