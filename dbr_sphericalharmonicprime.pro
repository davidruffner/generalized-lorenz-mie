;+
;NAME:
;    dbr_sphericalharmonicprime
;
; PURPOSE:
;    Calculate the derivative of spherical harmonic Y_l^m with respect
;    to theta
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    out = dbr_sphericalharmonicprime(theta,phi,l,m)
;
;INPUTS:
;    theta:  polar angle  
;
;    phi:    azimuthal angle
;
;    l:      total angular momentum number
;
;    m:      azimuthal angular momentum number
;
;OUTPUTS:
;    out: value of the derivative of the spherical harmonic
;
;DEPENDENCY:
;    Calls dbr_plegendre.pro
;
;REFERENCE:
;    Press W. H. et al. "Numerical Recipes 3rd ed."
;    Cambridge University Press (2007)
;
;MODIFICATION HISTORY:
; 03/14/2010 Written by David B. Ruffner, New York University


function dbr_sphericalharmonicprime,theta,phi, l, m

theta = double(theta)
phi = double(phi)
costheta = cos(theta)
sintheta = sin(theta)
plm = dbr_plegendre(cos(theta),l,m)

ylm = plm*exp(complex(0,1)*m*phi)

im = complex(0,1.d)

if m le 0 then begin
   mp = -m
   parityfactor = (-1.d)^mp
   ;print,"parity factor",parityfactor
   normfactor = Sqrt((2.d*double(l)+ 3.d)*(l+ 1.d -mp)/$
                          ((2.d*double(l) + 1.d)*(l+ 1.d +mp))) 
   dylmdtheta = parityfactor* $
                exp(im*m*phi)*((l-mp+ 1.d)*$
                dbr_plegendre(costheta,l+ 1.d,mp)/normfactor -$
                           (l+ 1.d)*costheta*dbr_plegendre(costheta,l,mp))/ $
                           sintheta
endif else begin
   normfactor = Sqrt((2.d*double(l)+ 3.d)*(l+ 1.d -m)/$
                             ((2.d*double(l)+ 1.d)*(l+ 1.d +m)))
   dylmdtheta = exp(im*m*phi)*((l-m+ 1.d)*$
                dbr_plegendre(costheta,l+ 1.d,m)/normfactor -$
                           (l+ 1.d)*costheta*dbr_plegendre(costheta,l,m))/ $
                           sintheta
endelse
b1 = where(finite(dylmdtheta,/nan),/null)
if b1 ne !NULL then begin 
   dylmdtheta(b1) = 0
endif
return,dylmdtheta
end
  
