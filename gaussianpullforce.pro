;+
;NAME:
;    gaussianpullforce.pro
;
; PURPOSE:
;    Calculate the max pulling force along a gaussian beam 
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    force = gaussianpullforce(ap,np,nm,thetaG,gamma,int=int,nt=nt,norm=norm)
;
;INPUTS:
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
;
;OUTPUTS:
;    pullforce: strongest pulling force of gaussian beam 
;
;DEPENDENCY:
;
;MODIFICATION HISTORY:
; 2014/06/17 Written by David B. Ruffner, New York University

function gaussianpullforce, ap,np,nm,lambda,thetaG,gamma,int=int,norm=norm,nt=nt

if n_elements(int) eq 0 then int = 1.
if n_elements(NT) eq 0 then NT = 10.
if n_elements(norm) eq 0 then norm = 0

;speed of light
c = 299792458.d; m/s

k = 2*!pi*nm/lambda

;Calculate the force constant
f0mN = !pi*(ap^2.)*int/c & $;mN
f0 = f0mN*10.^9. & $;pN
if norm eq 1 then f0 = 1
;print,f0

;Calculate the sphere coefficients
ab = sphere_coefficients(ap,np,nm,lambda)
an = ab[0,*]
bn = ab[1,*]
nc = n_elements(an)-1
;print,nc

;Calculate a few points over the range to find where to start looking
;for the minima. This assumes the axial force profile is smooth

;zvalues to calculate force over
w0 = sqrt(8*!pi)*lambda/(2*!pi*sin(thetaG))
zr = !pi*w0^2/lambda
zmax = 1.5*zr
zmin = 0
npts = 10
zs = (zmax-zmin)*findgen(npts)/(npts-1.d)+zmin


forces = fltarr(3,npts)
for i=0,npts-1 do begin $
   print,string(13b),"getting beam coefficients",i,format='(A,A,I,$)' & $
   pos=[0,0,zs[i]] & $
   bscs = gaussiantrapcoefficientsint(pos,nc,k,gamma,thetaG,nt) & $
   ;Calculate the force
   forces[*,i] = f0*normbartonforce(bscs,ap,np,nm,lambda) & $
endfor

order = sort(forces[2,*])
minforce = min(forces[2,*],minind)
if minind eq 0 or minind eq n_elements(forces[2,*])-1 then begin $
   print,"No bracket!!" & $
   plot,zs,forces[2,*] & $
   return,1 & $
   endif
abc = zs[minind-1:minind+1]
fabc = forces[2,minind-1:minind+1]

;Now use Golden Section Search method to find minima- 10.2 Numerical
;                                                     Recipes
tol = 10.^(-6.)
;; abc = zs[order[0:2]]
;; fabc = forces[2,order[0:2]]
gr = .38197
s = sort(abc)
abc = abc(s) 
fabc = fabc(s) 

count = 0
;print,"searching for pull force"
while abc[1]-abc[0] gt tol do begin $
   ;sort by position
   s = sort(abc) & $
   abc = abc(s) & $
   fabc = fabc(s) & $
   ;pick new point in bigger half
   if abc[1] -abc[0] gt abc[2]-abc[1] then begin $
      x = abc[1] - (abc[1] -abc[0])*gr & $
   endif else begin $
      x = abc[1] + (abc[2] -abc[1])*gr & $
   endelse & $
   ;evaluate the force at this new point
   pos = [0,0,x] & $
   bscs = gaussiantrapcoefficientsint(pos,nc,k,gamma,thetaG,nt) & $
   ;Calculate the force
   newf = f0*normbartonforce(bscs,ap,np,nm,lambda) & $
   fx = newf[2] & $
   ;form new bracket
   fs = [fabc,fx] & $
   abcx = [abc,x] & $
   s1 = sort(abcx) & $
   abcx = abcx(s1) & $
   fs = fs(s1) & $
   minf = min(fs,ind) & $
   abc = abcx[ind-1:ind+1] & $
   fabc = fs[ind-1:ind+1] & $

   
   count+=1 & $
endwhile
;print,"new min ",count,abc[1],fabc[1] 

return,fabc[1]

end



