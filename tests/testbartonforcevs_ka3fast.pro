;Script for testing the code for calculating the force on a sphere in
;lorenz mie theory for a general beam.

; For this test we'll use
;a plane wave and compare the force we get using our code
;to the results in Bohren and Huffman.

;The parameters used are in Irvine1964.pdf in this directory and we 
;see good agreement with the bottom curve in Figure 3.

;Additionally it is a good test of the special functions in the
;code since it uses plane waves at an angle

;Author David Ruffner
;Last edited 2014_06_02


np = 1.5
nm = 1.
lambda = .532
k = 2*!pi*nm/lambda

amin = .5/k
amax = 30/k
npts = 50
as = (amax-amin)*findgen(npts)/(npts-1)+amin

ab = sphere_coefficients(amax,np,nm,lambda)
an = ab[0,*]
nc = n_elements(an)

;Calculate the beam shape coefficients for a plane wave
forces = fltarr(npts)
forcesbh = fltarr(npts)

;angled plane wave
alpha = !pi/4.
beta = !pi/2
gamma = !pi/4.
for i=0,npts-1 do begin $
   bscs = planewavecoefficientsangled(nc,alpha,beta,gamma) & $

   ;Calculate the norm force
   force = normbartonforce(bscs,as[i],np,nm,lambda) & $
   forces[i] = sqrt(total(force^2)) & $
   forcesbh[i] = planewaveforce(as[i],np,nm,lambda) & $
   print,i & $
endfor 

p = plot(k*as,forces,symbol=".",color="red",$
                              name="glmt")
p = plot(k*as,forcesbh,linestyle=2,color="blue",$
                              name="bh",/overplot)
p.xtitle = "x"
p.ytitle = "abs(norm force)"
;; p.title = "Difference between Barton and Bohren force vs plane wave tilt angle"

p.save,"planewaveforcevs_ka3fast.png"
