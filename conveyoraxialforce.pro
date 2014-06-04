;+
;NAME:
;    conveyoraxialforce.pro
;
; PURPOSE:
;    Calculate the axial force along a conveyor optical axis
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    force = conveyoraxialforce(ap,np,nm,eta1,eta2,intensity,npts)
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
;KEYWORDS:
;
;    int: in mW/um^2
;
;    npts:   number of points to calculate force along
;
;    norm: When set the program outputs additional info
;
;OUTPUTS:
;    force: [3,npts] force vector at each point along period of conveyor 
;
;DEPENDENCY:
;
;MODIFICATION HISTORY:
; 2014/04/14 Written by David B. Ruffner, New York University
; 2014/06/02 DBR:Fixed bug with f0, updated documentation 
; 2014/06/04 DBR:Made the conveyor beam normalized at intensity maxima
;                by dividing the bessel beam coefficients by two.

function conveyoraxialforce, ap,np,nm,lambda,eta1,eta2,int=int,npts=npts,$
                             norm=norm

if n_elements(npts) eq 0 then npts = 100
if n_elements(int) eq 0 then int = 1.

;speed of light
c = 299792458.d; m/s

theta1 = acos(eta1)
theta2 = acos(eta2)

;print,"theta1",theta1
;print,"theta2",theta2


k = 2*!pi*nm/lambda

;Calculate the force constant
f0mN = !pi*(ap^2.)*int/c & $;mN
f0 = f0mN*10.^9. & $;pN
;print,f0

;Calculate the sphere coefficients
ab = sphere_coefficients(ap,np,nm,lambda)
an = ab[0,*]
bn = ab[1,*]
nc = n_elements(an)-1
;print,nc

;zvalues to calculate force over
zt = lambda/(nm*abs(eta1-eta2))
zmax = zt
zmin = 0
zs = (zmax-zmin)*findgen(npts)/(npts-1.d)+zmin


forces = fltarr(3,npts)
for i=0,npts-1 do begin $
   print,string(13b),"getting beam coefficients",i,format='(A,A,I,$)' & $
   pos=[0,0,zs[i]] & $
   bscs1 = besselcoefficients(pos,theta1,nc,k) & $
   bscs2 = besselcoefficients(pos,theta2,nc,k) & $
   bscs = (bscs1+bscs2)/2. & $
   ;Calculate the force
   if norm eq 1 then begin $
       forces[*,i] = normbartonforce(bscs,ap,np,nm,lambda) & $
   endif else forces[*,i] = f0*normbartonforce(bscs,ap,np,nm,lambda) & $
endfor

return,[transpose(zs),forces]

end



