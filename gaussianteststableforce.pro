;program to test if the gaussian beam is transversely stable at a given
;axial position.

function gaussianteststableforce,stableroot,ap,np,nm,lambda,eta1,eta2

;Check if it's transversly stable
nptsx = 5
forcesx = gaussianxforce(stableroot,ap,np,nm,lambda,eta1,eta2,$
                                               norm=1,int=1,npts=nptsx)
forcesy = gaussianyforce(stableroot,ap,np,nm,lambda,eta1,eta2,$
                                               norm=1,int=1,npts=nptsx)
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

if xstiffness lt 0 or ystiffness lt 0 then begin
   return, 0
endif

return,1
end
