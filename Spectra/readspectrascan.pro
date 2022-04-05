function calck, gap
beamline=7

; K(gap) = A*exp(B*(gap/LU)-C*(gap/LU)^2)
; gap = LU*( (B + sqrt(B^2 - 4*C*ln(K/A)))/(2*C) )


; PU01a, PU01b, PU02, PU03, PU04, PU05, PU06, PU07, PU08, PU09, PU10,  PU11

coeff_A = [7.1567,   7.16677,  5.3769,   6.8358, 0,  6.7359,  7.11306359521,  6.8192,  6.6198,  7.04743,  6.6498,  7.1264]
coeff_B = [-3.3815, -3.41482,  -3.446, -3.38148, 0, -3.3723, -3.43948118481, -3.3889, -3.3723, -3.38194, -3.4854, -3.4683]
coeff_C = [0.16065, 0.199182, 0.26176, 0.150098, 0, 0.15653, 0.254105985615, 0.17542,  0.1442,  0.20397,  0.2868, 0.32219]
coeff_LU = [31.4,       31.4,     23.,     29.0, 0,     29.,           31.4,     29.,     29.,     31.4,     29.,    31.4]

; lambda = (lambda_undulator/(2*gamma^2*n))*(1+K^2/2)
; Energy = n*(CTE/lambda_undulator)*(1/(1+K^2/2))
;
; CTE = (h*c*2*gamma^2)
; h = 4.1356673*1e-15 eV*s
; c = 299792458*1e3 mm/s
; gamma = 1957*6.085 [GeV] ???
; lamda_undulator = 29 mm
; CTE = 351306.2997 eV*mm

CTE = 351640.678249   ; eV*mm

A=coeff_A[beamline]
B=coeff_B[beamline]
C=coeff_C[beamline]
LU=coeff_LU[beamline]

K = A*exp(B*(gap/LU)-C*(gap/LU)^2)

return,k
end

function readspectra, name, info=info, noerror=noerror
  print,name
  inf=file_info(name)
  if not inf.exists then begin
    print, 'Spectra is not existing: ',name
    if keyword_set(noerror) then return,0
    retall
  endif else if not inf.read then begin
    print, 'Spectra is not readable: ',name
    retall
  endif
  st=bytarr(inf.size)
  openr, lun, /get_lun, name
  readu,lun,st
  close,lun
  free_lun,lun
  stcon=strsplit(string(st),string(13b)+string(10b),/extract,count=co)
  s=''
  for i=0,4 do begin
    h0=strsplit(stcon[2*i],/extract,co=coh0)
    h1=strsplit(stcon[2*i+1],/extract,co=coh1)
    if coh0 ne coh1 then begin
      stop
    endif
    for j=0,coh1-1 do begin
      s+=(j eq 0 ? "" : ",")+h0[j]+"="+h1[j]
    endfor
  endfor
  info=s
  
  erg=fltarr(coh0,co-10)
  
  for i=0,co-11 do begin
    erg[*,i]=float(strsplit(stcon[i+10],/extract))
  endfor
  
  return,erg
end

function readspectrascan, name, energy=energy, kvalue=kvalue, gap=gap

inf=file_info(name+'.hdr')
if not inf.exists then begin
   print, 'Spectra-Scan is not existing: ',name
   retall
endif else if not inf.read then begin
  print, 'Spectra-Scan is not readable: ',name
  retall
endif

st=bytarr(inf.size)
openr, lun, /get_lun, name+'.hdr'
readu,lun,st
close,lun
free_lun,lun
stcon=strsplit(string(st),string(13b)+string(10b),/extract,count=co)
posconf=where(stcon eq '[SCANCONFIGS]')
posdesc=where(stcon eq '[DESCRIPTIONS]')
possuff=where(stcon eq '[SUFFICES]')
posfile=where(stcon eq '[FILES]')
confarr=strarr(2,posdesc-posconf-1)
for i=0,posdesc[0]-posconf[0]-2 do begin
  h=strsplit(stcon[i+1+posconf[0]],/extract,co=coh)
  confarr(0,i)=h[0]
  if coh eq 2 then  confarr(1,i)=h[1]
endfor

s='info=create_struct('
for i=0,posdesc[0]-posconf[0]-2 do begin
   s+=(i eq 0 ? "'" : ",'")+confarr(0,i)+"'"+','+"'"+confarr(1,i)+"'"
endfor
s+=')'
res=execute(s)

filearr=strarr(2,fix(info.points))

for i=0,fix(info.points)-1 do begin
  h=strsplit(stcon[i+posfile[0]+1],': ',/extract)
  filearr[0,i]=strtrim(h[0],2)
  filearr[1,i]=strtrim(h[1],2)
endfor

i=-1
repeat begin
  i++
  sp0=readspectra(name+'/'+filearr[1,i]+'.dc0',info=in,/noerror)
endrep until (n_elements(sp0) ne 1)

energy=reform(sp0[0,*])
kvalue=reform(float(filearr[0,*]))

gh=indgen(1000)/1000.*25+9.5
kh=calck(gh)
gap=interpol(gh,kh,kvalue)

erg=fltarr(fix(info.points),n_elements(energy))

for i=0,fix(info.points)-1 do begin
   sph=readspectra(name+'/'+filearr[1,i]+'.dc0',info=in,/noerror)
   erg[i,*]=(n_elements(sph) eq 1 ? 0 :sph[1,*])
endfor
return,erg
end
