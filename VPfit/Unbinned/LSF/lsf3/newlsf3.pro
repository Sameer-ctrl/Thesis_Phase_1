PRO newlsf2, file=file, G=G

;file = 'fuv_G130M_1222_lsf.dat'
;G='G130M'



data=read_table(file)

sz=size(data)

cols=sz[1]
rows=sz[2]


wav_centres=intarr(cols)

for i=1,cols-1 DO BEGIN
wav_centres[i]=fix(data[i])
ENDFOR

index_val=findgen(rows-1)
index_val=index_val+1


for i=0,cols-1 DO BEGIN
print, i
str=string(wav_centres[i])
outfile='cos_' + G + '_lsf_' + strtrim(str,1) + '_ltpos_3.dat'

lsf_val=dblarr(rows-1)

for j=0, rows-2 DO BEGIN
lsf_val[j]=data[cols+i+j*cols]
endfor

writecol, outfile, index_val, lsf_val

ENDFOR

print, wav_centres
END
