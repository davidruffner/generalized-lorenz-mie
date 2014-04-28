;+
;NAME:
;    conveyorfixedpoint.pro
;
; PURPOSE:
;    Calculate the axial fixed point and the axial and transvers
;stiffness at that point
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    fixedpt = conveyorfixedpoint(ap,np,nm,eta1,eta2,intensity,npts,$
;                                  fr_filename=fr_filename,$
;                                  fz_filename=fz_filename)
;
;INPUTS:
;    ap:     radius of particle  
;
;    np:     index of refraction of partice
;
;    nm:     index of refraction of medium
;
;    lambda: vacuum wavelength of trapping light
;
;    eta1:   axial wavevector of first bessel beam component
;
;    eta2:   axial wavevector of second bessel beam component
;
;    int: in watts/um^2
;
;    npts:   number of points to calculate force along
;
;KEYWORDS:
;   fr_filename:  name of the file name to write the radial force
;profile
;
;   fz_filename: name of the filename to write the axial force profile.
;
;OUTPUTS:
;    fixedpt: [zroot,zstiff,rstiff] defining parameters of fixed point 
;
;DEPENDENCY:
;
;MODIFICATION HISTORY:
; 2014/04/14 Written by David B. Ruffner, New York University

function conveyorfixedpoint, ap,np,nm,lambda,eta1,eta2,int=int,npts=npts,$
                             fr_filename=fr_filename,$
                             fz_filename=fz_filename

if n_elements(npts) eq 0 then npts = 100
if n_elements(int) eq 0 then int = 1.
print,eta1,eta2
;Calculate the axial forces on axis
forcesz = conveyoraxialforce(ap,np,nm,lambda,eta1,eta2,int=int,npts=npts)

if n_elements(fz_filename) ne 0 then write_gdf,forcesz,fz_filename

dforceszdz = deriv(forcesz[0,*],forcesz[3,*])

;Find the axial fixed point
sign1 = forcesz[3,0]
sign2 = forcesz[3,1]
print,sign1,sign2

count = 0
found = 1
while ~(sign1 gt 0 and sign2 le 0) and count lt npts-2 do begin $
   print,string(13b),"searching for root...",count,format='(A,A,I,$)' & $
   count+=1 & $
   sign1 = forcesz[3,count] & $
   sign2 = forcesz[3,count+1] & $
   found = 0 & $
   if sign1 gt 0 and sign2 le 0 then found=1
endwhile

;If you don't find any axial stable point then exit
if ~found then begin
   return, [-1,0,0]
endif

stableroot2 = forcesz[0,count]-forcesz[3,count]/dforceszdz[count]
stableroot3 = forcesz[0,count+1]-forcesz[3,count+1]/dforceszdz[count+1]
stableroot = mean([stableroot2,stableroot2])


zstiffness = -mean([dforceszdz[count],dforceszdz[count+1]])

;Now Calculate the radial stiffness
;Let's go the size of a particle radius on either side

forcesx = conveyorxforce(stableroot,ap,np,nm,lambda,eta1,eta2,$
                                                   int=int,npts=npts)
forcesx[1,*] = -forcesx[1,*];For some reason transverse force is opposite what
                            ;it should be. FIX ME
print,"printing key parameters"
print,stableroot,ap,np,nm,lambda,eta1,eta2,npts

if n_elements(fr_filename) ne 0 then write_gdf,forcesx,fr_filename

dforcesxdx = deriv(forcesx[0,*],forcesx[1,*])

rstiffness = -dforcesxdx[npts/2.+1]

print,"stable root", stableroot
print,"z stiffness  ",zstiffness
print,"r stiffness  ",rstiffness

return,[stableroot,zstiffness,rstiffness]

end


