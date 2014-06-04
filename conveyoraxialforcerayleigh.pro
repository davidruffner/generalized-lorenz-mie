;+
;NAME:
;    conveyoraxialforcerayleigh.pro
;
; PURPOSE:
;    Calculate the axial force along a conveyor optical axis using
;    rayleigh approximation
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    force = conveyoraxialforcerayleigh(ap,np,nm,eta1,eta2,intensity,npts)
;
;INPUTS:
;    ap:     radius of particle  
;
;    np:     index of refraction of partice
;
;    nm:     index of refraction of medium
;
;    lambda: vacuum wavelength of trapping light
;
;    eta1:   axial wavevector of first bessel beam component
;
;    eta2:   axial wavevector of second bessel beam component
;
;    int: in watts/um^2
;
;    npts:   number of points to calculate force along
;
;OUTPUTS:
;    force: [3,npts] force vector at each point along period of conveyor 
;
;DEPENDENCY:
;
;MODIFICATION HISTORY:
; 2014/04/14 Written by David B. Ruffner, New York University
; 2014/06/14 DBR: Note, I need to normalize the conveyor but not sure how

function conveyoraxialforcerayleigh, ap,np,nm,lambda,eta1,eta2,$
                                     int=int,npts=npts

if n_elements(npts) eq 0 then npts = 100
if n_elements(int) eq 0 then int = 1.

ci = complex(0,1.d)

;speed of light
c = 299792458.d; m/s

deta = eta1-eta2
meta = (eta1+eta2)/2

theta1 = acos(eta1)
theta2 = acos(eta2)

;print,"theta1",theta1
;print,"theta2",theta2


k = 2*!pi*nm/lambda

;Calculate the sphere coefficients
ab = sphere_coefficients(ap,np,nm,lambda)
a1 = ab[0,1]
b1 = ab[1,1]

;Calculate the force constant
;polint = -12*!pi*a1*int/(ci*c*nm*k^2)
polint = 6*!pi*ci*a1*int/(c*k^2)
fintmN = real_part(polint) & $;mN
fint = fintmN*10.^6. & $;pN

fphsmN = imaginary(polint) & $;mN
fphs = fphsmN*10.^6. & $;pN


;zvalues to calculate force over
zt = lambda/(nm*abs(eta1-eta2))
zmax = zt
zmin = 0
zs = (zmax-zmin)*findgen(npts)/(npts-1.d)+zmin

print,"fint",fint," fphs",fphs

fz = -fint*deta*sin(deta*k*zs)+4*fphs*meta*(cos(deta*k*zs/2))^2

forces = fltarr(3,npts)
forces[2,*] = fz

return,[transpose(zs),forces]

end



