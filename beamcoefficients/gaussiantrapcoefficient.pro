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
;             gaussiantrapcoefficient1.pro and gaussiantrapcoefficient2.pro 


function am_integrandreal,theta
common share1,params
k = params.k
w0 = params.w0
n = params.n
m = params.m
pos = params.pos
gamma = params.gamma

;calculate multiplicative factor for gaussian beam
fac  =  k^2.d * w0^2.d / (4.d * !dpi)
fac *= SQRT(COS(theta))
fac *= EXP(-(gamma * SIN(theta))^2.d)
fac *= SIN(theta)
;Call besselcoefficient
am_aepart = fac*besselcoefficient(n,m,pos,theta,k)

return,real_part(am_aepart[0])
end

function ae_integrandreal,theta
common share1,params
k = params.k
w0 = params.w0
n = params.n
m = params.m
pos = params.pos
gamma = params.gamma

;calculate multiplicative factor for gaussian beam
fac  =  k^2.d * w0^2.d / (4.d * !dpi)
fac *= SQRT(COS(theta))
fac *= EXP(-(gamma * SIN(theta))^2.d)
fac *= SIN(theta)
;Call besselcoefficient
am_aepart = fac*besselcoefficient(n,m,pos,theta,k)


return,real_part(am_aepart[1])
end

function am_integrandimaginary,theta
common share1,params
k = params.k
w0 = params.w0
n = params.n
m = params.m
pos = params.pos
gamma = params.gamma

;calculate multiplicative factor for gaussian beam
fac  =  k^2.d * w0^2.d / (4.d * !dpi)
fac *= SQRT(COS(theta))
fac *= EXP(-(gamma * SIN(theta))^2.d)
fac *= SIN(theta)
;Call besselcoefficient
am_aepart = fac*besselcoefficient(n,m,pos,theta,k)

return,imaginary(am_aepart[0])
end

function ae_integrandimaginary,theta
common share1,params
k = params.k
w0 = params.w0
n = params.n
m = params.m
pos = params.pos
gamma = params.gamma

;calculate multiplicative factor for gaussian beam
fac  =  k^2.d * w0^2.d / (4.d * !dpi)
fac *= SQRT(COS(theta))
fac *= EXP(-(gamma * SIN(theta))^2.d)
fac *= SIN(theta)
;Call besselcoefficient
am_aepart = fac*besselcoefficient(n,m,pos,theta,k)

return,imaginary(am_aepart[1])
end



function gaussiantrapcoefficient,n,m,pos,k,gamma,thetaG,NT

common share1, params

;calculate w0
w0 = SQRT(8.d * !Pi) / (k * sin(thetaG))

params = {w0:w0,$
          n:n,$
          m:m,$
          pos:pos,$
          gamma:gamma,$
          k:k}

eps = 10^(-3.)
;print,"getting trap coefficient",m,n
am_mn_real = qsimp('am_integrandreal',0,thetaG,eps=eps)
ae_mn_real = qsimp('ae_integrandreal',0,thetaG,eps=eps)
am_mn_imaginary = qsimp('am_integrandimaginary',0,thetaG,eps=eps)
ae_mn_imaginary = qsimp('ae_integrandimaginary',0,thetaG,eps=eps)
;print,"done with integrals"
ci = dcomplex(0.d,1.d)
am_mn = am_mn_real + ci * am_mn_imaginary
ae_mn = ae_mn_real + ci * ae_mn_imaginary

return,[am_mn,ae_mn]
end



