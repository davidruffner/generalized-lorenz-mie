;+
; NAME:
;       dbr_sphericalfield
;
; PURPOSE:
;       Calculates the complex electric field defined by an array of scattering
;       coefficients.
;
; CATEGORY:
;       Holography, light scattering, microscopy
;
; CALLING SEQUENCE:
;       field = dbr_sphericalfield(x, y, z, a, lambda, $
;                           mpp = mpp)
;
; INPUTS:
;       x: [npts] array of pixel coordinates [pixels]
;       y: [npts] array of pixel coordinates [pixels]
;       z: If field is required in a single plane, then
;          z is the plane's distance from the sphere's center
;          [pixels].
;          Otherwise, z is an [npts] array of coordinates.
;
;       NOTE: Ideally, x, y and z should be double precision.
;             This is left to the calling program for efficiency.
;
;       a: [2,nc] array of a and b scattering coefficients, where
;          nc is the number of terms required for convergence.
;
;       lambda: wavelength of light in medium [pixels]
;
; KEYWORD FLAGS:
;       cartesian: If set, return field components in Cartesian
;           coordinates.  Default: Spherical polar coordinates
;
; OUTPUTS:
;       field: [3,npts] complex values of field at the positions r.
;              [0,*]: r component
;              [1,*]: theta component
;              [2,*]: phi component
;
;              If CARTESIAN is set:
;              [0,*]: x component (incident polarization)
;              [1,*]: y component (transverse component)
;              [2,*]: z component (axial component, relative to
;              incident beam).
;
; REFERENCE:
;   1. Adapted from Chapter 4 in
;      C. F. Bohren and D. R. Huffman, 
;      Absorption and Scattering of Light by Small Particles,
;      (New York, Wiley, 1983).
;   2. W. J. Wiscombe,
;      Improved Mie scattering algorithms,
;      Applied Optics 19, 1505-1509 (1980).
;
; MODIFICATION HISTORY:
; Written by David G. Grier, New York University, 5/2007
; 6/9/2007: DGG finally read Section 4.8 in Bohren and Huffman about
;    numerical stability of the recursions used to compute the scattering
;    coefficients.  Feh.  Result is a total rewrite.
; 6/20/2007: DGG Calculate \tau_n(\cos\theta) and \pi_n(\cos\theta)
;    according to recurrence relations in 
;    W. J. Wiscombe, Appl. Opt. 19, 1505-1509 (1980).
;    This is supposed to improve numerical accuracy.
; 2/8/2008: DGG. Replaced single [3,npts] array of input coordinates
;    with two [npts] arrays for x and y, and a separate input for z.
;    Eliminated double() call for coordinates.  Z may have 1 element or
;    npts elements. Small documentation fixes.
; 4/3/2008: Bo Sun (Sephiroth), NYU: Calculate Lorenz-Mie a and b
;    coefficients using continued fractions rather than recursion.
;    Osman Akcakir from Arryx pointed out that the results are
;    more accurate in extreme cases.  Method described in
;    William J. Lentz, "Generating Bessel functions in Mie scattering
;    calculations using continued fractions," Appl. Opt. 15, 668-671
;    (1976).
; 4/4/2008: DGG small code clean-ups and documentation.  Added
;    RECURSIVE keyword for backward compatibility in computing a and b
;    coefficients.
; 4/11/2008: Sephiroth: Corrected small error in jump code for
;    repeated fractions in Mie coefficients.
; 6/25/2008: DGG Don't clobber x coordinate input values.
; 10/9/2008: DGG adapted from SPHEREFIELD by separating out
;    calculation of scattering coefficients, a_n and b_n.  This
;    is therefore more general, and can be replaced more
;    readily with a GPU-accelerated version.
; 10/13/2008: DGG eliminated RECURSIVE keyword.
; 3/19/2014: Adapted from sphericalfield.pro by David B. Ruffner,
;            New York University
;
; Copyright (c) 2007-2010 Bo Sun and David G. Grier
;-

function dbr_sphericalfield, x_, y_, z_, ab, lambda, $
                         cartesian = cartesian ; project to cartesian coordinates

npts = n_elements(x_)
nc = n_elements(ab[0,*])-1      ; number of terms required for convergence

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

; storage for scattered field
Es = dcomplexarr(3,npts)

; Compute field by summing multipole contributions
for n = 1.d, nc do begin
   for m = -n,n do begin
      if m ne 1 and m ne -1 then continue
;Calculate the coefficients 
      if m eq 1 then begin
         am_mn = - ci^n*sqrt(4.d*!pi*(2.d*n + 1.d ))
         ae_mn = -ci*am_mn
      endif else begin
         am_mn = -ci^n*sqrt(4.d*!pi*(2.d*n + 1.d ))
         ae_mn = ci*am_mn
      endelse
      print,"m,n",m,n
; Calculate the special functions
      norm = sqrt((2*n+1)*factorial(n-m)/(4*!pi*factorial(n+m)))
      pi_mn = -dbr_pi_mn(costheta,n,m)
      tau_mn = -dbr_tau_mn(costheta,n,m)
      xi_n = dbr_xi_n(kr,n)
      dn = -dbr_xi_nprime(kr,n) ;derivative of xi_n
      p_mn = dbr_plegendre(costheta,n,m)
    
;Vector spherical harmonics Jackson (10.55)
;    M1n[0,*] = 0.d
      Mmn[1,*] = -xi_n*pi_mn/sqrt(n*(double(n)+1))
; ... divided by exp(ci*m*phi)/kr
      Mmn[2,*] = -ci*xi_n*tau_mn/sqrt(n*(double(n)+1)) 
; ... divided by exp(ci*m*phi)/kr

      Nmn[0,*] = -xi_n*pi_mn*sqrt(n*(double(n)+1))/m
; ... divided by exp(ci*m*phi)/(kr)^2
      Nmn[1,*] = -dn*tau_mn/sqrt(n*(double(n)+1))
; ... divided by exp(ci*m*phi)/kr
      Nmn[2,*] = -ci*dn*pi_mn/sqrt(n*(double(n)+1))
; ... divided by exp(ci*m*phi)/kr

; the scattered field in spherical coordinates (4.45)
      Es += (ae_mn*ab[0,n] * Nmn - am_mn*ab[1,n] * Mmn)
   endfor

endfor

; geometric factors were divided out of the vector
; spherical harmonics for accuracy and efficiency ...
; ... put them back at the end.
Es[0,*] *= cosphi*sintheta / kr^2
Es[1,*] *= cosphi/ kr
Es[2,*] *= sinphi / kr

; By default, the scattered wave is returned in spherical
; coordinates.  Project components onto Cartesian coordinates.
; Assumes that the incident wave propagates along z and 
; is linearly polarized along x
if keyword_set(cartesian) then begin
    Ec = Es
    Ec[0,*] =  Es[0,*] * sintheta * cosphi
    Ec[0,*] += Es[1,*] * costheta * cosphi
    Ec[0,*] -= Es[2,*] * sinphi

    Ec[1,*] =  Es[0,*] * sintheta * sinphi
    Ec[1,*] += Es[1,*] * costheta * sinphi
    Ec[1,*] += Es[2,*] * cosphi

    Ec[2,*] =  Es[0,*] * costheta - Es[1,*] * sintheta

    return, Ec
endif

return, Es
end
