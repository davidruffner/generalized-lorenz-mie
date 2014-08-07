;Script for testing the code for calculating the force on a sphere in
;lorenz mie theory for a general beam. 
;Parameters for beam and particle 
;based on paper by Cizmar et. al. (Cizmar2006.pdf)

;Makes the graph in the top left hand corner of Figure 1 in Cizmar2006.pdf

;Author David Ruffner
;2014_06_02, Last edit 
;Edited by Henrique Moyses

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

;Calculate the sphere coefficients
ab = sphere_coefficients(a,np,nm,lambda)
an = ab[0,*]
bn = ab[1,*]
nc = n_elements(an)-1

;xvalues to calculate force over
xmax = a * 2.
xmin = -a * 2.
npts = 10
xs = (xmax-xmin)*findgen(npts)/(npts-1)+xmin

forces = fltarr(3,npts)

for i=0,npts-1 do begin $
   print,"getting coefficients" & $
   print,xs[i] & $
   pos=[xs[i],0,0] & $
   bscs = gaussiantrapcoefficients(pos,nc,k,0.d,thetaG,100.) & $
   print,i & $
   ;Calculate the force
   forces[*,i] = f0*normbartonforce(bscs,a,np,nm,lambda) & $
   endfor

   p=plot(xs,-forces[0,*],symbol="square",color="red") ;There is a bug and
                                ;The force in xy is opposite; FIXME
p=plot(xs,fltarr(npts),linestyle=2,/overplot)
p.xtitle = "x ($\mu$m)"
p.ytitle = "F$_x$ (pN)"

p.title ="Radial force on Mie sphere from gaussian trap 1.5 um Poly"
p.save,"GaussFXvsX.png"



