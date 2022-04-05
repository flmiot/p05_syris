function untitled_1, input
res = input
return, res / 2
end

plot,en/1000,alog(erg[1,*])&oplot,[45,45],[0,40]

for i=0,200 do begin plot,en/1000,alog(erg[i,*]),title='k['+strtrim(i,2)+']='+strtrim(k[i],2)&oplot,[45,45],[0,40]&c=get_kbrd(0)&if c eq 's' then stop else if c eq 'm' then i=i-2
for i=0,200 do begin plot,en/1000,alog(erg[i,*]),title='k['+strtrim(i,2)+']='+strtrim(kv[i],2)&oplot,[45,45],[0,40]&c=get_kbrd(0)&if c eq 's' then stop else if c eq 'm' then i=i-2
for i=0,200 do begin plot,en/1000,alog(erg[i,*]),title='k['+strtrim(i,2)+']='+strtrim(kv[i],2)&oplot,[45,45],[0,40]&c=get_kbrd(1)&if c eq 's' then stop else if c eq 'm' then i=i-2
for i=0,200 do begin plot,en/1000,alog(erg[i,*]),title='k['+strtrim(i,2)+']='+strtrim(kv[i],2)&oplot,[45,45],[0,40]&c=get_kbrd(1)&if c eq 's' then stop else if c eq 'm' then i=i-2

for i=0,200 do begin plot,en/1000,alog(erg[i,*]),title='k['+strtrim(i,2)+']='+strtrim(kv[i],2)&oplot,en/1000,alog(erg1[i,*]),col=100&oplot,[45,45],[0,40]&c=get_kbrd(1)&if c eq 's' then stop else if c eq 'm' then i=i-2

erg1=readspectrascan('D:\data\spectra_p07\p07_k_2p2_1p2_filter_c_cu',en=en,kv=kv,ga=ga)

erg=readspectrascan('p07_hb_1p2_2p2_98m_1mm_1mm.hdr', en=en, kv=kv, ga=ga)
plot,en/1000,alog(erg[1,*])&oplot,[45,45],[0,40]
wshow

for i=0,200 do begin 
  ;plot,en/1000,alog(erg[i,*]),title='k['+strtrim(i,2)+']='+strtrim(kv[i],2)
  ;oplot,en/1000,alog(erg1[i,*]),col=100
  oplot,[45,45],[0,40]
  ;print,kv[i],ga[i]
  c=get_kbrd(1)
  if c eq 's' then stop else if c eq 'm' then i=i-2
endfor

; Plotting
plot,en/1000,alog(erg[1,*])&oplot,[45,45],[0,40]

; File handling
openr, lun, /get_lun, name+'.hdr'
readu,lun,st
close,lun
free_lun,lun