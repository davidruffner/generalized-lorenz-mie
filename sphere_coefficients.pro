;+
; NAME:
;       sphere_coefficients
;
; PURPOSE:
;       Calculates the Mie scattering coefficients for a homogeneous
;       isotropic sphere illuminated
;       by a coherent plane wave linearly polarized in the x direction.
;
; CATEGORY:
;       Holography, light scattering, microscopy
;
; CALLING SEQUENCE:
;       ab = sphere_coefficients(a, np, nm, lambda)
;
; INPUTS:
;       a: radius of sphere [micrometers]
;
;       np: (complex) refractive index of sphere
;
;       nm: (complex) refractive index of medium
;
;       lambda: wavelength of light [micrometers]
;
; OUTPUTS:
;       ab: [2,nc] complex a and b scattering coefficients.
;
; REFERENCES:
;   1. Adapted from Chapter 4 in
;      C. F. Bohren and D. R. Huffman, 
;      Absorption and Scattering of Light by Small Particles,
;      (New York, Wiley, 1983).
;   2. W. J. Lentz, 
;      Generating Bessel functions in Mie scattering
;      calculations using continued fractions,
;      Applied Optics 15, 668-671 (1976).
;   3. M. I. Mischenko, L. D. Travis and A. A. Lacis,
;      Scattering, Absorption and Emission of Light by
;      Small Particles,
;      (Cambridge, Cambridge University Press, 2001).
;
; MODIFICATION HISTORY:
; Adapted from SPHEREFIELD by David G. Grier, New York Unviersty 10/09/2008.
; Relevant MODIFICATION HISTORY follows
; Written by David G. Grier, New York University, 5/2007
; 06/20/2007: DGG Calculate \tau_n(\cos\theta) and \pi_n(\cos\theta)
;    according to recurrence relations in 
;    W. J. Wiscombe, Appl. Opt. 19, 1505-1509 (1980).
;    This is supposed to improve numerical accuracy.
; 04/03/2008: Bo Sun (Sephiroth), NYU: Calculate Lorenz-Mie a and b
;    coefficients using continued fractions rather than recursion.
;    Osman Akcakir from Arryx pointed out that the results are
;    more accurate in extreme cases.  Method described in
;    William J. Lentz, "Generating Bessel functions in Mie scattering
;    calculations using continued fractions," Appl. Opt. 15, 668-671
;    (1976).
; 04/11/2008: Sephiroth: Corrected small error in jump code for
;    repeated fractions in Mie coefficients.
; 01/16/2009: DGG: Silently correct for negative input parameters,
;    which are encountered sometimes in fits gone awry.
; 02/08/2009: DGG: Corrected rare infinite loop in jump code.
;    Explicitly cast constants to double.
;
; Copyright (c) 2007-2010 Bo Sun and David G. Grier.
;-

; Compute the Lorenz-Mie scattering coefficients:
; INPUTS:
;    x: n_m k a -- reduced radius of sphere in medium of index n_m
;    m: n_p/n_m -- relative index of particle

function sphere_coefficients, ap, np, nm, lambda

maxnc = 300                     ; maximum number of multipole terms
maxn = 2000.                    ; break out of aberrant loops

x = abs(2.d * !dpi * nm * ap / lambda) ; size parameter
m = np/nm                              ; relative refractive index

mx = dcomplex(x*m)              ; scaled size parameter.

; series convergence (Mishchenko, (5.237) page 159)
nc = floor(x + 4.05d * x^(1.d/3.d) + 2.d)

if nc gt maxnc then begin
   message, "excessively many multipole terms", /inf
   print, "ap:", ap
   print, "np:", np
   print, "nm:", nm
   print, "lambda:", lambda
   print, "truncating to "+strtrim(maxnc)
   nc = maxnc
endif

; Equation numbers refer to Lentz (1976)
; Mie coefficients
; start at order 0 for bookkeeping
a = dcomplexarr(1, nc+1)
b = dcomplexarr(1, nc+1)

; Deirmendjian's derivative term
D = dblarr(1, nc+1)

; Calculate D_n(mx)
for ii = 1.d, nc do begin       ; Eq. (8)
   u = ii + 0.5d
   term = 0

   an = double(2.d*u/mx)        ; coefficient for n = 1, Eq. (9)
   fnn = an                     ; numerator in Eq. (8)
   fdn = 1.d                    ; denominator in Eq. (8)
   pnn = an                     ; numerator continued fraction: [an...a1]
   pdn = 1d80                   ; denominator continued fraction: [an...a2]
   
                                ; continued fractions with very small
                                ; values can can affect accuracy.  See
                                ; Eqs. (10), (11) and associated discussion.
                                ; flag special cases
   jumpnn = 0                   ; numerator jump flag
   jumpdn = 0                   ; denominator jump flag
   
   n = 2.d                      ; build up confinued fractions from n = 2
   while (term eq 0) do begin 
      an = (-1.d)^(n+1.d) * 2.d * (u + n - 1.d) / mx ; Eq. (9)

      if (abs(pnn) lt 1d-6) then begin   ; small partial numerator
         xin = pnn * an + 1.d            ; \xi for numerator 
         amp2 = (-1.d)^n * 2.d * (u + n) / mx
         Zn = (amp2 * xin + pnn) / xin   ; Z for numerator, Eq. (11)
         jumpnn = 3                      ; flag special case
         pnn =  1.d
      endif

      if (abs(pdn) lt 1d-6) then begin   ; small partial denominator
         xid = pdn * an + 1.d            ; \xi for denominator
         amp2 = (-1.d)^n * 2.d * (u + n) / mx
         Zd = (amp2 * xid + pdn)/ xid    ; Z for denominator, Eq. (11)
         jumpdn = 2                      ; flag special case
         pdn = 1.d
      endif

      if (jumpnn eq 0) then begin ; no jump needed
         pnn = an + 1.d/pnn       ; partial numerator [an...a1]
         fnn *= pnn               ; numerator up to order n
      endif

      if (jumpnn eq 1) then begin ; jump required in last step
         pnn = Zn 
         fnn *= xin
      endif

      if (jumpdn eq 0) then begin ; no jump needed
         pdn = an + 1.d/pdn       ; partial numerator [an...a1]
         fdn *= pdn               ; denominator up to order n
      endif

      if (jumpdn eq 1) then begin ; jump required in last step
         pdn = Zd 
         fdn *= xid
      endif

      ; stop at required accuracy when no jumps remain to be resolved
      term = ((abs(pnn/pdn-1.d) le 1d-15) and ((jumpdn+jumpnn) eq 0))

      ; update flags
      if jumpnn gt 0 then jumpnn--
      if jumpdn gt 0 then jumpdn--
      n++
      if n gt maxn then begin
         ;message, "did not converge in first loop", /inf
         ;print, "order: ", ii
         ;print, "pnn:", pnn, jumpnn
         ;print, "pdn:", pdn, jumpdn
         term = 1
      endif
   endwhile

   D[ii] = fnn/fdn - ii/mx
endfor

; \psi(x), \xi(x)
ASpsi = dblarr(1, nc+1)
ASeta = dblarr(1, nc+1)
ASpsi[0] = sin(x)
ASeta[0] = - cos(x)

ASxi  = dcomplex(ASpsi, ASeta)

for ii = 1.d, nc do begin         ; Eq. (8)
   u = ii + 0.5d
   term = 0
   n = 2.d                      ; start iteration from n = 2
   an = double(2.d*u/x)         ; a_n for n = 1 
   fnn = an                     ; numerator in Eq. (8)
   fdn = 1.d                    ; denominator in Eq. (8)
   pnn = an                     ; continued fraction for numerator: [an...a1]
   pdn = 1d80                   ; continued fraction for denominator: [an...a2]
   jumpnn = 0                   ; numerator jump flag, Eq. (10), (11)
   jumpdn = 0                   ; denominator jump flag, Eq. (10), (11)

   while (term eq 0) do begin 
      an = (-1.d)^(n + 1.d) * 2.d * (u + n - 1.d) / x

      if (abs(pnn) lt 1d-6) then begin ; small partial numerator
         xin = pnn * an + 1.d          ; \xi for numerator 
         amp2 = (-1.d)^n * 2.d * (u + n) / x
         Zn = (amp2 * xin + pnn) / xin ; Z for numerator, Eq. (11)
         jumpnn = 3                    ; flag special case
         pnn =  1.d
      endif

      if (abs(pdn) lt 1d-6) then begin ; small partial denominator
         xid = pdn * an + 1.d          ; \xi for denominator
         amp2 = (-1.d)^n * 2.d * (u + n) / x
         Zd = (amp2 * xid + pdn) / xid ; Z for denominator, Eq. (11)
         jumpdn = 2                    ; flag special case
         pdn = 1.d
      endif

      if (jumpnn eq 0) then begin ; no jump needed
         pnn = an + 1.d/pnn       ; partial numerator [an...a1]
         fnn = fnn * pnn          ; numerator up to order n
      endif

      if (jumpnn eq 1) then begin ; jump required in last step
         pnn = Zn 
         fnn = fnn * xin
      endif

      if (jumpdn eq 0) then begin ; no jump required
         pdn = an + 1.d/pdn       ; partial denominator [an...a2]
         fdn = fdn * pdn          ; denominator to order n
      endif

      if (jumpdn eq 1) then begin ; jump required in last step
         pdn = Zd 
         fdn = fdn * xid
      endif

      ; stop at required accuracy with all jumps resolved
      term = ((abs(pnn/pdn-1.d) le 1d-15) and ((jumpdn+jumpnn) eq 0))

      ; update flags
      if jumpnn gt 0 then jumpnn--
      if jumpdn gt 0 then jumpdn--
      n++
      if n gt maxn then begin
         ;message, "did not converge in second loop", /inf
         term = 1
      endif
   endwhile

   ; Abramowitz and Stegun (10.1.31):
   ASpsi[ii] = ASpsi[ii-1] / fnn * fdn
   ASeta[ii] = ASeta[ii-1] * fdn / fnn - 1.d/ASpsi[ii-1] 
   ASxi[ii] = dcomplex(ASpsi[ii], ASeta[ii])
endfor

for n = 1.d, nc do begin
   fac  = D[n]/m + n/x
   a[n] = (fac * ASpsi[n] - ASpsi[n-1])/(fac * ASxi[n] - ASxi[n-1])
   fac  = m * D[n] + n/x
   b[n] = (fac * ASpsi[n] - ASpsi[n-1])/(fac * ASxi[n] - ASxi[n-1])
endfor

return, [a, b]
end
