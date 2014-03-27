;+
;NAME:
;    efield_vsh_sum
;
; PURPOSE:
;    Calculate the field from its vector spherical
;    harmonic coefficients pmn and qmn
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    E = efield_vsh_sum(rvec,nmax)
;
;INPUTS:
;    rvec:    [3,N] array of positions where you want the field
;
;    pnm:     Beam shape coefficients 1
;
;    qnm:     Beam shape coefficients 2
;
;OUTPUTS:
;    E:   [3,N] Electric field vector at each point.
;
;REFERENCE:
;
;MODIFICATION HISTORY:
; 03/14/2014 Written by David B. Ruffner, New York University

function efield_vsh_sum,rvec,k,nm,a,abmn,cartesian=cartesian

sz = size(rvec)
E = fltarr(sz[1],sz[2])
im = complex(0,1)

Amn = abmn[*,*,0]
Bmn = abmn[*,*,1]
szp = size(Amn)
nmax = szp[1]-1

for n=1,nmax do begin
   for m = -n,n do begin
      if abs(Amn[n,m]) eq 0 and abs(Bmn[n,m]) eq 0 then begin
         continue
      endif
      print,"n,m ",n,m
      emn = efieldmn2(rvec,Amn[n,m],Bmn[n,m],m,n,k,a,nm)
      b1 = where(finite(emn,/nan),/null)
      if b1 ne !NULL then begin 
         print,"we have nansss!"
         stop
      endif
      E = E+emn
   endfor
endfor

x = rvec[0,*]
y = rvec[1,*]
z = rvec[2,*]
rho   = sqrt(x^2 + y^2)
r     = sqrt(rho^2 + z^2)
theta = atan(rho, z)
phi   = atan(y, x)
costheta = cos(theta)
sintheta = sin(theta)
cosphi = cos(phi)
sinphi = sin(phi)

if keyword_set(cartesian) then begin
    Ec = E
    Ec[0,*] =  E[0,*] * sintheta * cosphi
    Ec[0,*] += E[1,*] * costheta * cosphi
    Ec[0,*] -= E[2,*] * sinphi

    Ec[1,*] =  E[0,*] * sintheta * sinphi
    Ec[1,*] += E[1,*] * costheta * sinphi
    Ec[1,*] += E[2,*] * cosphi

    Ec[2,*] =  E[0,*] * costheta - E[1,*] * sintheta

    return, Ec
endif


return,E
end
