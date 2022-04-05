;  returns tmp_harmonic, set_gap, k_value, set_gap+meas_gap_offset

function calc_p05_gap, energy,bestintensity=bestintensity
beamline=5
;HIGHEST_HARMONIC=30
meas_gap_offset= + 0. ; measured gap offset
minimum_gap = 9.8 - meas_gap_offset

if ( (beamline lt 0) or (beamline gt 21) ) then begin
  message, "Not beamline defined. Check Beamline property."
endif


; K(gap) = A*exp(B*(gap/LU)-C*(gap/LU)^2)
; gap = LU*( (B + sqrt(B^2 - 4*C*ln(K/A)))/(2*C) )

;PU01a, PU01b, PU02, PU03, PU04, PU05, PU06, PU07, PU08, PU09, PU10, PU11, PU21, PU22, PU23, PU24, PU25, PU61, PU62, PU63, PU64, PU65
coeff_A = [7.1567, 7.16677, 5.8368, 6.8358, 0, 6.7359, 6.9285, 6.8192, 7.2253, 7.2285, 6.4432, 7.1264, 0., 0., 0., 7.6315, 0., 0., 0., 0., 8.5921, 8.6096 ]
coeff_B = [-3.3815, -3.41482, -3.4729, -3.38148, 0, -3.3723, -3.3318, -3.3889, -3.5297, -3.5266, -3.3956, -3.4683, 0., 0., 0., -3.5277, 0., 0., 0., 0.,-3.5469, -3.5378]
coeff_C = [0.16065, 0.199182, 0.19491, 0.150098, 0, 0.15653, 0.12184, 0.17542, 0.23594, 0.36212, 0.17752, 0.32219, 0., 0., 0., 0.22369, 0., 0., 0., 0., 0.24743, 0.24427]
coeff_LU = [ 31.4, 31.4, 23., 29.0, 0, 29., 31.4 , 29., 29., 31.4, 29., 31.4, 0., 0.,  0., 0.,  0., 0.,  0., 0.,  0., 0. ] ; ACHTUNG NOCH FALSCH

; lambda = (lambda_undulator/(2*gamma^2*n))*(1+K^2/2)
; Energy = n*(CTE/lambda_undulator)*(1/(1+K^2/2))
;
; CTE = (h*c*2*gamma^2)
; h = 4.1356673*1e-15 eV*s
; c = 299792458*1e3 mm/s
; gamma = 1957*6.08 [GeV] ???
; lamda_undulator = 29 mm
; CTE = 351306.2997 eV*mm

CTE = 351063.0345219   ; eV*mm

;with 'n' the harmonic number
; Minimum energy for each harmonic is computed from the miminum gap 
; value, giving the maximum K. Minimum gap 9.8 (PL05).
  
tmp = coeff_B[beamline]*(minimum_gap/coeff_LU[beamline])+(-coeff_C[beamline]*(minimum_gap/coeff_LU[beamline])*(minimum_gap/coeff_LU[beamline]));
Kmax = coeff_A[beamline]*exp(tmp);

; Energy limits for harmonics:
; Minimum Energy (to be multiplied by Harmonic n):
  Emin=(CTE/coeff_LU[beamline])*(1/(1 + Kmax*Kmax/2));
; Emax corresponds to K = 0 => E = Harmonic n*CTE
  Emax= (CTE/coeff_LU[beamline]) 
    
  en=energy*1000.
  
  maxharmonic= fix(en / emin)
  minharmonic= fix(en / emax)
  
  if minharmonic lt 1 then minharmonic=1

  erg=fltarr(4,maxharmonic+1)
    
    
  ;  K(gap) = A*exp(B*(gap/LU)-C*(gap/LU)^2)
  ; gap = LU*( (B + sqrt(B^2 - 4*C*ln(K/A)))/(2*C) )
  
  for i=minharmonic, maxharmonic do begin
       tmp_harmonic=i
       if (2*(tmp_harmonic*(CTE/coeff_LU[beamline])/en) - 2 ge 0) then begin
          k_value = sqrt(2*(tmp_harmonic*(CTE/coeff_LU[beamline])/en) - 2);
          tmp_value = coeff_B[beamline]*coeff_B[beamline]-4*coeff_C[beamline]*alog(k_value/coeff_A[beamline]);
          tmp_value = sqrt(tmp_value);
          set_gap = coeff_LU[beamline] *( (coeff_B[beamline] + tmp_value) / (2*coeff_C[beamline]) );
          if finite(set_gap) then erg[*,i]= [ tmp_harmonic, set_gap, k_value, set_gap+meas_gap_offset ] 
        endif
  endfor
  w=where(erg[0,*] ne 0)
  if keyword_set(bestintensity) then begin
     maxw=max(w)
     print,erg[*,maxw]
     res=erg[3,maxw]
  endif else res=erg[*,w]
  return,res
  end
  