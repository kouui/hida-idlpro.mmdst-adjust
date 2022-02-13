
; history
; 2022.02.13    k.u.    for display IQUV[R] 2d spectra

;; s4 (npix_x, npix_y, >=4), >=4 : I,Q,U,V[,R]
PRO showiquvr, s4, pmax=pmax, bin=bin, wid=wid, ialog=ialog, $
         labels=labels, x0s=x0s, right_bounds=right_bounds,$
         dd=dd, wexist=wexist

   if not keyword_set(wid) then wid=0
   if not keyword_set(pmax) then pmax=0.05
   if not keyword_set(bin) then bin=1

	ss = size(s4)
   if ss[0] ne 3 then throw_error, "ndim of s4 !=3"
   if ss[3] ne 4 and ss[3] ne 5 then throw_error, "not allowed size of the 3rd dimension of s4==", ss[3]
   
   if keyword_set(labels) and n_elements(labels) ne ss[3] $ 
   then throw_error, 'size of labels not match of size of 3rd dimension of s4'
   coms = ['I', 'Q', 'U', 'V'] & if ss[3] eq 5 then coms = [coms,['R']]

   nx0 = ss[1]
   ny0 = ss[2]
   ns0 = ss[3]
   nx = nx0 / bin
   ny = ny0 / bin
   
   dd = 2
   if keyword_set(wexist) then wset, wid  $
   else window,wid,xs=nx*ns0+dd*(ns0-1),ys=ny+dd*2
   
   ;; show I spectrum image
   x0 = 0
   x0s = indgen(ns0) & x0s[0] = x0
   ;x0s = [x0,0,0,0]
   y0 = 0
   i2d=congrid(s4[*,*,0],nx,ny)
   if keyword_set(ialog) then tvscls,alog(i2d),x0,y0+dd $
   else tvscls,i2d,x0,y0+dd,sgm=[-3,1]
   if keyword_set(islabel) then xyouts,x0+10,ny-30,coms[0],chars=3,charthick=2,/dev
   ;x0s = indgen(ns0) & x0s[0] = x0
   right_bounds = x0s & right_bounds[0] += nx
;   right_bounds = [x0 + nx,0,0,0]
   ;; show Q, U, V spectrum image
   for j=1,ns0-1 do begin
		x0 = x0 + nx + dd
		tvscls,congrid(s4[*,*,j],nx,ny),x0,y0+dd,vmin=-pmax,vmax=pmax
		if keyword_set(islabel) then xyouts,x0+10,ny-30,coms[j],chars=3,charthick=2,/dev
	   right_bounds[j] = x0 + nx
      x0s[j] = x0
   endfor
END