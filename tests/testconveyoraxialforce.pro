;Script for testing the code for calculating the force on a sphere in
;lorenz mie theory for a general beam. 
;Parameters for beam and particle 
;based on paper by Cizmar et. al. (Cizmar2006.pdf)

;Makes the graph in the top left hand corner of Figure 1 in Cizmar2006.pdf

;Author David Ruffner
;2014_06_02, Last edit 

;Input parameters
int = 10.^(3.); mW/um^2
c = 299792458.d; m/s

theta1 = 0.357924
theta2 = 0.132825

nm = 1.33
lambda = 0.532/nm
k = 2*!pi/lambda
a = .005;um
np = 1.59

;Calculate the force constant
f0mN = !pi*(a^2)*int/c & $;mN
f0 = f0mN*10.^9. & $;pN

;Calculate the sphere coefficients
ab = sphere_coefficients(a,np,nm,lambda)
an = ab[0,*]
bn = ab[1,*]
nc = n_elements(an)-1

;zvalues to calculate force over
zmax = 5.3
zmin = -2.
npts = 200
zs = (zmax-zmin)*findgen(npts)/(npts-1)+zmin


forces = fltarr(3,npts)
for i=0,npts-1 do begin $
   print,"getting coefficients" & $
   pos=[0,0,zs[i]] & $
   bscs1 = besselcoefficients(pos,theta1,nc,k) & $
   bscs2 = besselcoefficients(pos,theta2,nc,k) & $
   bscs = bscs1+bscs2 & $
   print,i & $
   ;Calculate the force
   forces[*,i] = f0*normbartonforce(bscs,a,np,nm,lambda) & $
endfor

p=plot(zs,forces[2,*],symbol="square",color="red")
p=plot(zs,fltarr(npts),linestyle=2,/overplot)
p.xtitle = "z ($\mu$m)"
p.ytitle = "F$_z$ (pN)"

p.title ="Axial force on Mie sphere from conveyor"

p.save,"conveyorFzvsZ.png"



