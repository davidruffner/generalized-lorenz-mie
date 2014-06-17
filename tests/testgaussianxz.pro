;;Script for testing the code for calculating the field of the
;;gaussian beam
;;
;

;Author Henrique Moyses
;2014_06_03, Last edit 

;Set up parameters of the beam

axis1 = 'x'
axis2 = 'z'
plane  = axis1+axis2
nmax = 40
n0 = 40L
n1 = 40L

nm = 1.33
lambda = 0.532
k = 2*!pi*nm/lambda
thetaG= !dpi/3.
w0 = 2. / (k * sin(thetaG))
NT = 10.
gamma = 0.


;Setting up points to calculate field
r0max = 2.
r0min = -2.
r1max = 2.
r1min = -2.

r0 = findgen(n0)/(n0-1.)  * (r0max - r0min) + r0min
r1 = findgen(n1)/(n1-1.)  * (r1max - r1min) + r1min

y0 = 0.
npts = n0 * n1
x = ((findgen(npts) mod n0) / (n0-1.)) * (r0max - r0min) + r0min
z = ((float(floor(findgen(npts) / n0)) ) / (n1-1.)) * (r1max - r1min) + r1min
y = fltarr(npts) + y0

print,min(x)
print,min(z)

print,"getting coefficients"
pos=[0.,0.,0.d]
bscs = gaussiantrapcoefficientsint(pos,nmax,k,0.d,thetaG,nt)


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
help,inttot
print,max(inttot)
s=surface(inttot)
index = bytscl(findgen(256))+70   
img = contour(alog(inttot),r0,r1,rgb_table=34,/fill,$
                n_levels=256,c_color=index,aspect_ratio=1,$
                xrange=[r0min,r0max],yrange=[r1min,r1max])
img.xtitle = ""+axis1+" ($\mu$m)"
img.ytitle = ""+axis2+" ($\mu$m)"
img.title = "Gaussian intensity crossection "+plane
img.xtickinterval = floor(1.3*r0max)
img.ytickinterval = floor(.3*r1max)
img.save,"gaussianint"+plane+".png",border=5


