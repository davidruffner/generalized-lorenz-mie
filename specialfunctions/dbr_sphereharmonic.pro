;+
;NAME:
;    dbr_sphericalharmonic
;
; PURPOSE:
;    Calculate the spherical harmonic Y_l^m
;
;CATEGORY:
;    Mathematics
;
;CALLING SEQUENCE:
;    out = dbr_sphericalharmonic(theta,phi,l,m)
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
;    out: value of the spherical harmonic
;
;DEPENDENCY:
;    Calls dbr_plegendre.pro
;
;REFERENCE:
;    Press W. H. et al. "Numerical Recipes 3rd ed."
;    Cambridge University Press (2007)
;
;MODIFICATION HISTORY:
; 03/07/2010 Written by David B. Ruffner, New York University


function dbr_sphericalharmonic,theta,phi, l, m

plm = dbr_plegendre(l,m,cos(double(theta)))

ylm = plm*exp(complex(0.d,1.d)*m*double(phi))

return,ylm
end
  
