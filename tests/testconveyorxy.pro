;;Script for testing the code for calculating the field of the optical
;;conveyor
;Parameters for beam based on paper by Cizmar et. al. (Cizmar2006.pdf)

;Author David Ruffner
;2014_06_02, Last edit 

;Set up parameters of the beam
theta1 = 0.357924;Convergence angle of two Bessel beam components
theta2 = 0.132825

axis1 = 'x'
axis2 = 'y'
plane  = axis1+axis2
dv0 = [1,0,0]
dv1 = [0,1,0]
nmax = 100
n0 = 30.
n1 = 30.


nm = 1.33
lambda = 0.532
k = 2*!pi*nm/lambda



;Setting up points to calculate field
r0max = 5.2
r0min = -5.2
r1max = 5.2
r1min = -5.2
dv0 = dv0*(r0max-r0min)/n0
dv1 = dv1*(r1max-r1min)/n1
r0 = (r0max-r0min)*findgen(n0)/(n0-1)+r0min
r1 = (r1max-r1min)*findgen(n1)/(n1-1)+r1min
rcorner = -dv0*n0/2-dv1*n1/2

pts = fieldcrossectionpoints(rcorner,dv0,dv1,n0,n1)
ptslist = reform(pts,3,n0*n1)
x = ptslist[0,*]
y = ptslist[1,*]
z = ptslist[2,*]


print,"getting coefficients"
pos=[0,0,0]
bscs1 = besselcoefficients(pos,theta1,nmax,k)
bscs2 = besselcoefficients(pos,theta2,nmax,k)
bscs = bscs1+bscs2

print,"calculating the field"
e = efield_vsh_sum(x,y,z,bscs,lambda,nm,nmax,/cartesian)

;write_gdf,e,"conveyorefield"+plane+".gdf"
help,e

intx = abs(e[0,*])^2
inty = abs(e[1,*])^2
intz = abs(e[2,*])^2

intx2d = reform(intx,n0,n1)
inty2d = reform(inty,n0,n1)
intz2d = reform(intz,n0,n1)
;plotimage,bytscl(intx2d+inty2d+intz2d),/iso
inttot = intx2d+inty2d+intz2d

index = bytscl(findgen(256))   
img = contour(inttot,r0,r1,rgb_table=34,/fill,$
                n_levels=256,c_color=index,aspect_ratio=1,$
                xrange=[r0min,r0max],yrange=[r1min,r1max])
img.xtitle = ""+axis1+" ($\mu$m)"
img.ytitle = ""+axis2+" ($\mu$m)"
img.title = "Conveyor intensity crossection "+plane
img.xtickinterval = floor(.3*r0max)
img.ytickinterval = floor(.3*r1max)
img.save,"conveyorint"+plane+".png",border=5


