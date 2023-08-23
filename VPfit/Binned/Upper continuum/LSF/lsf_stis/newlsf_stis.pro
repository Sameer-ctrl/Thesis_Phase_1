PRO newlsfstis, file=file

;file = 'STIS_LSF_E230H_2400.txt'

readcol, file, wav,lsf1,lsf2,lsf3,lsf4,stringskip=2 


sz=size(wav)

rows=sz[1]


index_val=findgen(rows-1)
index_val=index_val+1


outfile=file+'_1.dat'
writecol, outfile, index_val, lsf1

outfile=file+'_2.dat'
writecol, outfile, index_val, lsf2

outfile=file+'_3.dat'
writecol, outfile, index_val, lsf3


outfile=file+'_4.dat'
writecol, outfile, index_val, lsf4

END
