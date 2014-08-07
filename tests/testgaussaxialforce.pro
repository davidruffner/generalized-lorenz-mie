;Script for testing the code for calculating the force on a sphere in
;lorenz mie theory for a gaussian beam.
;Author Henrique Moyses

;Modification history
; 2014_08_07: David Ruffner edited 
 
;Parameters for beam and particle 
;Input parameters
int = 100; mW/um^2 (10)
c = 299792458.d; m/s

nm = 1.33
lambda = 0.532
k = 2*!pi*nm/lambda
a = 0.5;um
np = 1.45
thetaG= !dpi / 3.
;Calculate the force constant
f0mN = nm*!pi*(a^2)*int/c & $;mN
f0 = f0mN*10.^9. & $;pN
;; f0=1.
;Calculate the sphere coefficients
ab = sphere_coefficients(a,np,nm,lambda)
an = ab[0,*]
bn = ab[1,*]
nc = n_elements(an)-1

;zvalues to calculate force over
zmax = a * 2.
zmin = -a * 2.
npts = 10
zs = (zmax-zmin)*findgen(npts)/(npts-1)+zmin


forces = fltarr(3,npts)
for i=0,npts-1 do begin $
   print,"getting coefficients" & $
   pos=[0,0,zs[i]] & $
   bscs = gaussiantrapcoefficients(pos,nc,k,0.d,thetaG,100.) & $
   print,i & $
   ;Calculate the force
   forces[*,i] = f0*normbartonforce(bscs,a,np,nm,lambda) & $
endfor
;write_gdf,forces,'axialforce.gdf'
p=plot(zs,forces[2,*],symbol="square",color="red")
p=plot(zs,fltarr(npts),linestyle=2,/overplot)
p.xtitle = "z ($\mu$m)"
p.ytitle = "F$_z$ (pN)"

p.title ="Axial force on Mie sphere from gaussian trap 1.5 um Poly"

p.save,"GaussFzvsZ.png"



