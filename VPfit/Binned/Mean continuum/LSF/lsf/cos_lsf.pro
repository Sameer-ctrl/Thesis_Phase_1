;; This code is to interpolate the cos lsf into the desired line center. 
;; Inputs: lam0     => Line center 
;;         ltpos    => Lifetime position of HST/COS (ltpos = 2 if data is taken after Nov. 2012) 
;;         grating  => The grating with which the data is obtained (G130M/G160M/G185M) 
;; 
;; Sowgat- 28/05/2013 
;; 	 
pro cos_lsf,lam0,ltpos=ltpos,grating=grating,plot=plot           
;; 
if (n_params() lt 1) then begin 
	print, 'CALLING SEQUENCE: '+ $ 
		'cos_lsf,lam0,ltpos=ltpos,grating=grating' 
return  
endif 
;; 
path1='~/IDLCODES/cos_lsf/' 
;; 
if ltpos eq 1 then begin 
	if grating eq 'G130M' then begin 
	arr=[1150.,1200.,1250.,1300.,1350.,1400.,1450.]
	format='(i,d,d,d,d,d,d,d)'  
	readcol,path1+'cos_lsf_G130M.dat',col0,col1,col2,col3,col4,col5,col6,col7,f=fmt 
	print, 'Reading data file: cos_lsf_G130M.dat'
	col = [1#col1,1#col2,1#col3,1#col4,1#col5,1#col6,1#col7]  
	endif 
	;; 
	if grating eq 'G160M' then begin 
	arr = [1450.,1500.,1550.,1600.,1650.,1700.,1750.]   
	format='(i,d,d,d,d,d,d,d)'   
	readcol,path1+'cos_lsf_G160M.dat',col0,col1,col2,col3,col4,col5,col6,col7,f=fmt  
	print, 'Reading data file: cos_lsf_G160M.dat'
	col = [1#col1,1#col2,1#col3,1#col4,1#col5,1#col6,1#col7]  
	endif  
	;;	
	if grating eq 'G185M' then begin 
	arr=[1700.,1800.,1900.,2000.,2100.,2200.,2300.,2400.,2500.,2600.,2700,$ 
	     2800.,2900.,3000.,3100.,3200]   
	format='(i,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d)'  
	readcol,path1+'cos_nuv_lsf.dat',col0,col1,col2,col3,col4,col5,col6,col7,col8,col9,$
	col10,col11,col12,col13,col14,col15,col16,f=fmt 
	print, 'Reading data file: cos_nuv_lsf.dat'
	col=[1#col1,1#col2,1#col3,1#col4,1#col5,1#col6,1#col7,1#col8,1#col9,$ 
	     1#col10,1#col11,1#col12,1#col13,1#col14,1#col15,1#col16] 
	endif 	 
	;; 
endif else begin 
	if grating eq 'G130M' then begin 
	arr=[1150.,1200.,1250.,1300.,1350.,1400.,1450.]
	format='(i,d,d,d,d,d,d,d)'  
	readcol,path1+'cos_lsf_ltpos2_G130M_1300.dat',col0,col1,col2,col3,col4,col5,col6,col7,f=fmt 
	print, 'Reading data file: cos_lsf_G130M.dat'
	col = [1#col1,1#col2,1#col3,1#col4,1#col5,1#col6,1#col7]  
	endif 
	;; 
	if grating eq 'G160M' then begin 
	arr = [1400.,1450.,1500.,1550.,1600.,1650.,1700.,1750.]   
	format='(i,d,d,d,d,d,d,d)'   
	readcol,path1+'cos_lsf_ltpos2_G160M_1600.dat',col0,col1,col2,col3,col4,col5,col6,col7,col8,f=fmt  
	print, 'Reading data file: cos_lsf_G160M.dat'
	col = [1#col1,1#col2,1#col3,1#col4,1#col5,1#col6,1#col7,1#col8]    
	endif  
	;;	
	if grating eq 'G185M' then begin 
	message,/info, '!!No data file available for G185M for ltpos = 2, so far!!'  
	return     
	endif  
endelse 
;; 
	ip_col=dblarr(n_elements(col0))
	for j=0, n_elements(col0)-1 do begin 	
		ip_col(j) = interpol(col(*,j),arr(*),lam0)
	endfor 
;; 
	result=[1#col0,1#ip_col]  

	title= 'COS/'+grating+'  LSF at'+string(lam0,format='(f8.2)')+ ' Angstrom'   

	if keyword_set(plot) then begin  
	window,10 
	plot, col0, ip_col,psym=10,thick=2,title=title  
	endif   	

	filenum='cos_'+grating+'_lsf_'+string(lam0,format='(f5.0)')+'dat' 
	writecol,filenum,col0,ip_col,fmt='(f20.5,e20.10)' 
	
	message,/info, ' Output is written as:   '+filenum    
;;  
END 	   
