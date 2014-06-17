;+
;NAME:
;    gaussianaxialforce.pro
;
; PURPOSE:
;    Calculate the axial force along a conveyor optical axis
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    force = gaussianaxialforce(ap,np,nm,lambda,gamma,thetaG,$
;                               int=int,NT=NT,npts=npts)
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
;    gamma : Ratio between focal length of objective and width of
;            gaussian beam in the objective aperture.
;
;    thetaG : maximum convergence angle of objective
;
;
;KEYWORDS:
;
;    int: in mW/um^2
;
;    NT : Number of terms to use in integration (100 is in general a
;         good number)
;
;    npts:   number of points to calculate force along
;
;    norm: When set the force is normalized
;
;OUTPUTS:
;    force: [3,npts] force vector at each point along period of conveyor 
;
;DEPENDENCY:
;
;MODIFICATION HISTORY:
; 2014/06/13 Written by David B. Ruffner, New York University


function gaussianaxialforce, ap,np,nm,lambda,gamma,thetaG,$
                             int=int,NT=NT,npts=npts,norm=norm,zrange=zrange

if n_elements(npts) eq 0 then npts = 100
if n_elements(int) eq 0 then int = 1.
if n_elements(NT) eq 0 then NT = 10.
if n_elements(norm) eq 0 then norm = 0

if n_elements(zrange) eq 0 then begin
   w0 = sqrt(8*!pi)*lambda/(2*!pi*sin(thetaG))
   zr = !pi*w0^2/lambda
   zrange = [-zr,zr]
endif

;speed of light
c = 299792458.d; m/s

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
zmax = zrange[1]
zmin = zrange[0]
zs = (zmax-zmin)*findgen(npts)/(npts-1.d)+zmin


forces = fltarr(3,npts)
for i=0,npts-1 do begin $
                                ;print,string(13b),"getting
                                ;beam
                                ;coefficients",i,format='(A,A,I,$)' & $' 
   print,"getting beam coefficients",i & $
   pos=[0,0,zs[i]] & $
   bscs = gaussiantrapcoefficients(pos,nc,k,gamma,thetaG,NT) & $
   ;Calculate the force
   if norm eq 1 then begin $
       forces[*,i] = normbartonforce(bscs,ap,np,nm,lambda) & $
   endif else forces[*,i] = f0*normbartonforce(bscs,ap,np,nm,lambda) & $
endfor

return,[transpose(zs),forces]

end



