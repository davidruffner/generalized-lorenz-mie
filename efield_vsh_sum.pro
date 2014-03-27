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
; 03/21/2014 Overhaul based on dbr_sphericalfield

function efield_vsh_sum,x_,y_,z_,bscs,lambda,nc,$
                        cartesian=cartesian ; project to cartesian coordinates
tol = 1.e-12
npts = n_elements(x_)

k = 2.d * !dpi / lambda         ; wavenumber in medium [pixel^-1]

ci = dcomplex(0,1)

; convert to spherical coordinates centered on the sphere.
; (r, theta, phi) is the spherical coordinate of the pixel
; at (x,y) in the imaging plane at distance z from the
; center of the sphere.
rho   = sqrt(x_^2 + y_^2)
r     = sqrt(rho^2 + z_^2)
theta = atan(rho, z_)
phi   = atan(y_, x_)
costheta = cos(theta)
sintheta = sin(theta)
cosphi = cos(phi)
sinphi = sin(phi)

kr = k*r                        ; reduced radial coordinate

; storage for vector spherical harmonics: [r,theta,phi]
Mmn = dcomplexarr(3,npts)
Nmn = dcomplexarr(3,npts)

; storage for field
E = dcomplexarr(3,npts)

;Beam shape coefficients
am_mn = bscs[*,*,0]
ae_mn = bscs[*,*,1]

szp = size(am_mn)
nmax = szp[1]-1

; Compute field by summing multipole contributions
for n = 1.d, nc do begin
   for m = -n,n do begin
      if abs(am_mn[n,m]) lt tol and abs(ae_mn[n,m]) lt tol then continue
      ;if m eq 0 then continue
      print,"m,n",m,n,abs(am_mn[n,m]),abs(ae_mn[n,m])
; Calculate the special functions
      pi_mn = -dbr_pi_mn(costheta,n,m);normalized
      tau_mn = -dbr_tau_mn(costheta,n,m);normalized
      psi_n = dbr_psi_n(kr,n)
      dn = -dbr_psi_nprime(kr,n) ;derivative of psi_n
      p_mn = dbr_plegendre(costheta,n,m);normalized
    
;Vector spherical harmonics Jackson (10.55)
;    M1n[0,*] = 0.d
      Mmn[1,*] = -psi_n*pi_mn*exp(ci*m*phi)/sqrt(n*(double(n)+1))
; ... divided by 1/kr
      Mmn[2,*] = -ci*psi_n*tau_mn*exp(ci*m*phi)/sqrt(n*(double(n)+1)) 
; ... divided by 1/kr

      Nmn[0,*] = -psi_n*p_mn*exp(ci*m*phi)*sqrt(n*(double(n)+1))
; ... divided by sintheta/(kr)^2
      Nmn[1,*] = -dn*tau_mn*exp(ci*m*phi)/sqrt(n*(double(n)+1))
; ... divided by 1/kr
      Nmn[2,*] = -ci*dn*pi_mn*exp(ci*m*phi)/sqrt(n*(double(n)+1))
; ... divided by 1/kr
      
;Debugging
      b1 = where(finite(tau_mn,/nan),/null)
      if b1 ne !NULL then begin 
         print, "We have infinities in taumn"
      endif 
      b1 = where(finite(pi_mn,/nan),/null)
      if b1 ne !NULL then begin 
         print, "We have infinities in pimn"
      endif 
      ;; if m eq 0 then begin
      ;;    print,"Mmn",max(abs(Mmn),indM,dimension=2)
      ;;    print,"indM",indM
      ;;    print,"r,theta,phi",[transpose(r[indM]),transpose(theta[indM]),$
      ;;                   transpose(phi[indM])]
      ;;    print,"Nmn",max(abs(Nmn),indN,dimension=2)
      ;;    help,r
      ;; endif

; the scattered field in spherical coordinates (4.45)
      E += (ae_mn[n,m]* Nmn - am_mn[n,m]* Mmn)
   endfor

endfor

b1 = where(finite(E[0,*],/nan),/null)
if b1 ne !NULL then begin $
      print, "We have infinities in Er"
      stop
endif & $
b1 = where(finite(E[1,*],/nan),/null)
if b1 ne !NULL then begin $
      print, "We have infinities in Eth"
      stop
endif & $
b1 = where(finite(E[2,*],/nan),/null)
if b1 ne !NULL then begin $
      print, "We have infinities in Ephi"
      stop
endif & $

; geometric factors were divided out of the vector
; spherical harmonics for accuracy and efficiency ...
; ... put them back at the end.
E[0,*] *= 1 / kr^2
E[1,*] *= 1/ kr
E[2,*] *= 1 / kr


; By default, the scattered wave is returned in spherical
; coordinates.  Project components onto Cartesian coordinates.
; Assumes that the incident wave propagates along z and 
; is linearly polarized along x
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

return, E
end
;; for n=1,nmax do begin
;;    for m = -n,n do begin
;;       if abs(Amn[n,m]) eq 0 and abs(Bmn[n,m]) eq 0 then begin
;;          continue
;;       endif
;;       print,"n,m ",n,m
;;       emn = efieldmn2(rvec,Amn[n,m],Bmn[n,m],m,n,k,a,nm)
;;       b1 = where(finite(emn,/nan),/null)
;;       if b1 ne !NULL then begin 
;;          print,"we have nansss!"
;;          stop
;;       endif
;;       E = E+emn
;;    endfor
;; endfor

;; x = rvec[0,*]
;; y = rvec[1,*]
;; z = rvec[2,*]
;; rho   = sqrt(x^2 + y^2)
;; r     = sqrt(rho^2 + z^2)
;; theta = atan(rho, z)
;; phi   = atan(y, x)
;; costheta = cos(theta)
;; sintheta = sin(theta)
;; cosphi = cos(phi)
;; sinphi = sin(phi)

;; if keyword_set(cartesian) then begin
;;     Ec = E
;;     Ec[0,*] =  E[0,*] * sintheta * cosphi
;;     Ec[0,*] += E[1,*] * costheta * cosphi
;;     Ec[0,*] -= E[2,*] * sinphi

;;     Ec[1,*] =  E[0,*] * sintheta * sinphi
;;     Ec[1,*] += E[1,*] * costheta * sinphi
;;     Ec[1,*] += E[2,*] * cosphi

;;     Ec[2,*] =  E[0,*] * costheta - E[1,*] * sintheta

;;     return, Ec
;; endif


