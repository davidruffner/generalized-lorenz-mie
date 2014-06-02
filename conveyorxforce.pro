;+
;NAME:
;    conveyorxforce.pro
;
; PURPOSE:
;    Calculate the radial force at the axial fixed pt
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    force = conveyorxforce(z_fixed,ap,np,nm,eta1,eta2,intensity,npts)
;
;INPUTS:
;    z_fixed:the axial stable fixed point in um
;
;    ap:     radius of particle in um
;
;    np:     index of refraction of partice
;
;    nm:     index of refraction of medium
;
;    lambda: vacuum wavelength of trapping light in um
;
;    eta1:   norm axial wavevector of first bessel beam component
;
;    eta2:   norm axial wavevector of second bessel beam component
;
;KEYWORDS:
;
;    int: in mW/um^2
;
;    npts:   number of points to calculate force along
;
;    norm:   If set gives the force efficiency, Q = F/(!pi*ap^2*int/c)
;                where c is the speed of light
;
;    verbose: when set the code prints out additional info
;
;OUTPUTS:
;    force: [3,npts] force vector at each point along period of conveyor 
;
;DEPENDENCY:
;
;MODIFICATION HISTORY:
; 2014/04/14 Written by David B. Ruffner, New York University
; 2014/05/14 Added norm,verbose keywords DBR
; 2014/06/02 DBR:Fixed bug with f0, updated documentation 


function conveyorxforce, z_fixed,ap,np,nm,lambda,eta1,eta2,$
                              int=int,npts=npts,norm=norm,verbose=verbose

if n_elements(npts) eq 0 then npts = 100
if n_elements(int) eq 0 then int = 1.
if n_elements(norm) eq 0 then norm = 0
if n_elements(verbose) eq 0 then verbose = 0

;speed of light
c = 299792458.d; m/s

theta1 = acos(eta1)
theta2 = acos(eta2)

k = 2*!pi*nm/lambda

;Calculate the force constant
f0mN = !pi*(ap^2)*int/c & $;mN
f0 = f0mN*10.^9. & $;pN
if norm then f0 = 1


;Calculate the sphere coefficients
ab = sphere_coefficients(ap,np,nm,lambda)
an = ab[0,*]
bn = ab[1,*]
nc = n_elements(an)-1

;; ;xvalues to calculate force over
xmax = ap
xmin = -ap
xs = (xmax-xmin)*findgen(npts)/(npts-1)+xmin
;; ;yvalues to calculate force over
;; ymax = 2*ap
;; ymin = -2*ap
;; ys = (ymax-ymin)*findgen(npts)/(npts-1)+ymin


forces = fltarr(3,npts)
for i=0,npts-1 do begin $
   if verbose then $
      print,string(13b),"getting beam coefficients",i,format='(A,A,I,$)' & $
   pos=[xs[i],0,z_fixed] & $
   bscs1 = besselcoefficients(pos,theta1,nc,k) & $
   bscs2 = besselcoefficients(pos,theta2,nc,k) & $
   bscs = bscs1+bscs2 & $
   ;Calculate the force
   forces[*,i] = f0*normbartonforce(bscs,ap,np,nm,lambda) & $
endfor
if verbose then plot,xs,forces[0,*]
return,[transpose(xs),forces]

end



