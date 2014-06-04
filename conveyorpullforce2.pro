;+
;NAME:
;    conveyorpullforce.pro
;
; PURPOSE:
;    Calculate the axial force along a conveyor optical axis
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    force = conveyorpullforce2(ap,np,nm,eta1,eta2,intensity,npts)
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
;    eta1:   norm axial wavevector of first bessel beam component
;
;    eta2:   norm axial wavevector of second bessel beam component
;
;    int: in mW/um^2
;
;KEYWORDS:
;
;    npts:   number of points to calculate force along
;
;    norm:   when set outputs the force efficiency or normalized force
;
;OUTPUTS:
;    force: [3,npts] force vector at each point along period of conveyor 
;
;DEPENDENCY:
;
;MODIFICATION HISTORY:
; 2014/05/13 Written by David B. Ruffner, New York University
; 2014/06/04 DBR:Made the conveyor beam normalized at intensity maxima
;                by dividing the bessel beam coefficients by two.

function conveyorpullforce2, ap,np,nm,lambda,eta1,eta2,int=int,norm=norm

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
zt = lambda/(nm*abs(eta1-eta2))
zmax = zt
zmin = 0
npts = 10
zs = (zmax-zmin)*findgen(npts)/(npts-1.d)+zmin


forces = fltarr(3,npts)
for i=0,npts-1 do begin $
   print,string(13b),"getting beam coefficients",i,format='(A,A,I,$)' & $
   pos=[0,0,zs[i]] & $
   bscs1 = besselcoefficients(pos,theta1,nc,k) & $
   bscs2 = besselcoefficients(pos,theta2,nc,k) & $
   bscs = (bscs1+bscs2)/2. & $
   ;Calculate the force
   forces[*,i] = f0*normbartonforce(bscs,ap,np,nm,lambda) & $
   stable = teststableforce(zs[i],ap,np,nm,lambda,eta1,eta2) & $
   if not stable then forces[*,i] = [0.,0.,1.+i*.1] & $
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
;First check to see if it is bracketing a minima
;; if not (( (zs[order[0]] lt zs[order[1]]) and $
;;    (zs[order[0]] gt zs[order[2]])) or $
;;    ((zs[order[0]] gt zs[order[1]]) and $
;;    (zs[order[0]] lt zs[order[2]]))) then begin $
;;    ;print,"" & $   
;;    ;print,"we're bracketed" & $
;;    print,"No bracket!!" & $
;;    plot,zs,forces[2,*] & $
;;    stop & $
;;    return,1 & $

;; endif

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
   bscs1 = besselcoefficients(pos,theta1,nc,k) & $
   bscs2 = besselcoefficients(pos,theta2,nc,k) & $
   bscs = (bscs1+bscs2)/2. & $
   newf= f0*normbartonforce(bscs,ap,np,nm,lambda) & $
   stable = teststableforce(x,ap,np,nm,lambda,eta1,eta2) & $
   if not stable then newf = [0.,0.,1.+.001*count] & $
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

stableroot = abc[1]
;Check if it's transversly stable
nptsx = 10
forcesx = conveyorxforce(stableroot,ap,np,nm,lambda,eta1,eta2,$
                                               norm=1,int=1,npts=nptsx)
forcesy = conveyoryforce(stableroot,ap,np,nm,lambda,eta1,eta2,$
                                               norm=1,int=1,npts=nptsx)
forcesx[1,*] = -forcesx[1,*];For some reason transverse force is opposite what
                            ;it should be. FIX ME
forcesy[2,*] = -forcesy[2,*];For some reason transverse force is opposite what
                            ;it should be. FIX ME

if n_elements(fr_filename) ne 0 then write_gdf,forcesx,fr_filename

dforcesxdx = deriv(forcesx[0,*],forcesx[1,*])
dforcesydy = deriv(forcesy[0,*],forcesy[2,*])

xstiffness = -dforcesxdx[nptsx/2.+1]
ystiffness = -dforcesydy[nptsx/2.+1]



if xstiffness lt 0 or ystiffness lt 0 then begin 
  print,"unstable!" 
  print,"stable root", stableroot
  print,"x stiffness  ",xstiffness
  print,"y stiffness  ",ystiffness
  return,1.
endif

return,fabc[1]

end



