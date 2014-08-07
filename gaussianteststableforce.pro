;program to test if the gaussian beam is transversely stable at a given
;axial position.

function gaussianteststableforce,stableroot,ap,np,nm,lambda,thetaG,gamma,$
                                     nt=nt,kx=kx,ky=ky
if n_elements(nt) eq 0 then nt=30

;Check if it's transversly stable
nptsx = 5 
forcesx = gaussianxforce(stableroot,ap,np,nm,lambda,thetaG,gamma,$
                                               nt=nt,norm=1,npts=nptsx)
forcesy = gaussianyforce(stableroot,ap,np,nm,lambda,thetaG,gamma,$
                                               nt=nt,norm=1,npts=nptsx)
forcesx[1,*] = -forcesx[1,*];For some reason transverse force is opposite what
                            ;it should be. FIX ME
forcesy[2,*] = -forcesy[2,*];For some reason transverse force is opposite what
                            ;it should be. FIX ME

;if n_elements(fr_filename) ne 0 then write_gdf,forcesx,fr_filename

dforcesxdx = deriv(forcesx[0,*],forcesx[1,*])
dforcesydy = deriv(forcesy[0,*],forcesy[2,*])

xstiffness = -dforcesxdx[nptsx/2.+1]
ystiffness = -dforcesydy[nptsx/2.+1]

;print,"stable root", stableroot
;print,"x stiffness  ",xstiffness
;print,"y stiffness  ",ystiffness
ky = ystiffness
kx = xstiffness

if xstiffness lt 0 or ystiffness lt 0 then begin
   return, 0
endif

return,1
end
