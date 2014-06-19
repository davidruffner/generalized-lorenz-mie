;+
;NAME:
;    gaussianyforce.pro
;
; PURPOSE:
;    Calculate the radial force at the axial fixed pt
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    force = gaussianyforce(z_fixed,ap,np,nm,eta1,eta2,intensity,npts)
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
;    thetaG: Convergence angle of the Gaussian beam
;
;    gamma: Ratio between focal length of objective and width of
;            gaussian beam in the objective aperture.
;
;KEYWORDS:
;
;    nt:   number of points to integrate over angles
;
;    norm:   when set outputs the force efficiency or normalized force
;
;    int: in mW/um^2
;    verbose: when set the code prints out additional info
;
;OUTPUTS:
;    force: [3,npts] force vector at each point along period of conveyor 
;
;DEPENDENCY:
;
;MODIFICATION HISTORY:
; 2014/06/18 Written by David B. Ruffner, New York University


function gaussianyforce, z_fixed,ap,np,nm,lambda,thetaG,gamma,$
                         nt=nt,int=int,npts=npts,norm=norm,verbose=verbose

if n_elements(int) eq 0 then int = 1.
if n_elements(NT) eq 0 then NT = 10.
if n_elements(norm) eq 0 then norm = 0
if n_elements(verbose) eq 0 then verbose = 0

;speed of light
c = 299792458.d; m/s

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

;; ;; ;xvalues to calculate force over
;; xmax = ap
;; xmin = -ap
;; xs = (xmax-xmin)*findgen(npts)/(npts-1)+xmin
;yvalues to calculate force over
ymax = ap
ymin = -ap
ys = (ymax-ymin)*findgen(npts)/(npts-1)+ymin

forces = fltarr(3,npts)
for i=0,npts-1 do begin $
   if verbose then $
      print,string(13b),"getting beam coefficients",i,format='(A,A,I,$)' & $
   pos=[0,ys[i],z_fixed] & $
   bscs = gaussiantrapcoefficientsint(pos,nc,k,gamma,thetaG,nt) & $
   ;Calculate the force
   forces[*,i] = f0*normbartonforce(bscs,ap,np,nm,lambda) & $
endfor
if verbose then plot,ys,forces[1,*]
return,[transpose(ys),forces]

end



