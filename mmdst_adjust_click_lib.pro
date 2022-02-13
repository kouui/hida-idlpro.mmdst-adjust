

; history
; 2022.02.13    k.u.    i


;; RETURN CURSOR CLICK COORDINATE
FUNCTION on_bclick, wid=wid, msg=msg, data=data
    if keyword_set(wid) then wset, wid
    if keyword_set(msg) then print, msg
    
    if keyword_set(data) then bcursor,xpos,ypos,/data,/down, button=button $
    else cursor,xpos,ypos,/dev,/down, button=button
    return, {click_struct, x:xpos, y:ypos, button:button}
END


;------------------------------------------------------------------------
;; SELECT y POS WITH CURSOR
;; s2d (ny, 4), 4 : I,Q,U,V
FUNCTION select_ypos_profile, s2d, wid=wid, ialog=ialog, msg_arr=msg_arr

   ss = size(s2d) & ny = ss[1]
   if ss[0] ne 2 then throw_error,'ndim of s2d != 2'
   if ss[2] ne 4 then throw_error,'size of 2rd dimension != 4 (iquv)'

   
   if not keyword_set(wid) then wid=0
   if not keyword_set(msg_arr) then msg_arr = ['click left button to select a x posiiton ...']
   nclick = n_elements(msg_arr)

   labels = ['I','Q','U','V']
   colors = ['ffffff'x, '66ff00'x, '9900ff'x, '00ffff'x] ; 'bbggrr'x
   ;colors = ['White', 'Red', 'Green', 'Yellow']

   window, wid, xs=1200, ys=600
   i = 0
   mx = 0.05
   yr = max(s2d[*,1:3])
   yr = max(mx,yr)
   yr = yr*[-1.,+1.]
   plot, s2d[*,i]*mx, color=colors[i], linestyle=0, yr = yr
   ;axis, yaxis=1, yr = [min(s2d[*,1:3]),max(s2d[*,1:3])]
   for i=1,3 do begin
      oplot, s2d[*,i], color=colors[i], linestyle=0
   endfor

   for i=0,3 do xyouts,100+50*i,100,labels[i],chars=2.5, color=colors[i], /device

   xpos_arr = []
   for kk=0,nclick-1 do begin
      ;; click and compute x position
      pos = click_event(msg=msg_arr[kk], /data)
      xpos = fix(pos[0])

      print, 'select xpos=', xpos
      
      xpos_arr = [xpos_arr,[xpos]]
      oplot,[xpos], [s2d[xpos,0]*mx], color=colors[0], psym=5, thick=4
   endfor



   if nclick eq 1 then return, xpos_arr[0]
   return, xpos_arr


END

;------------------------------------------------------------------------
;; SELECT X POS WITH CURSOR
;; s4 (npix_x, npix_y, 4), 4 : I,Q,U,V
FUNCTION select_xpos, s4, pmax=pmax, bin=bin, wid=wid, ialog=ialog, islabel=islabel, msg_arr=msg_arr

   if not keyword_set(wid) then wid=0
   if not keyword_set(pmax) then pmax=0.05
   if not keyword_set(bin) then bin=1
   if not keyword_set(msg_arr) then msg_arr = ['click left button to select a x posiiton ...']
   nclick = n_elements(msg_arr)
   coms = ['I', 'Q', 'U', 'V']

	ss = size(s4)
   if not ss[0] eq 3 then throw_error, "ndim of s4 !=3"
   if not ss[3] eq 4 then throw_error, "size of the 3rd dimension of s4 !=4"
   nx0 = ss[1]
   ny0 = ss[2]
   ns0 = ss[3]
   nx = nx0 / bin
   ny = ny0 / bin
   
   dd = 2
   window,wid,xs=nx*ns0+dd*(ns0-1),ys=ny+dd*2
   ;; show I spectrum image
   x0 = 0
   x0s = [x0,0,0,0]
   y0 = 0
   i2d=congrid(s4[*,*,0],nx,ny)
   if keyword_set(ialog) then tvscls,alog(i2d),x0,y0+dd $
   else tvscls,i2d,x0,y0+dd,sgm=[-3,1]
   if keyword_set(islabel) then xyouts,x0+10,ny-30,coms[0],chars=3,charthick=2,/dev
   right_bounds = [x0 + nx,0,0,0]
   ;; show Q, U, V spectrum image
   for j=1,ns0-1 do begin
		x0 = x0 + nx + dd
		tvscls,congrid(s4[*,*,j],nx,ny),x0,y0+dd,vmin=-pmax,vmax=pmax
		if keyword_set(islabel) then xyouts,x0+10,ny-30,coms[j],chars=3,charthick=2,/dev
	   right_bounds[j] = x0 + nx
      x0s[j] = x0
   endfor

   xpos_arr = []
   for kk=0,nclick-1 do begin
      ;; click and compute x position
      pos = click_event(msg=msg_arr[kk])
      xpos = pos[0]

      bias = 0
      for j=0,ns0-1 do begin
         if right_bounds[j] gt xpos then break
         bias = bias + (nx+dd)
      endfor

      xpos = xpos - bias
      
   ;;   print, j, xpos
      for j=0, ns0-1 do begin
         draw, (x0s[j]+xpos)*[1,1], [dd,dd+ny], color=150, thick=2
      endfor
      xpos = xpos*bin
      print, 'select xpos=', xpos
      
      xpos_arr = [xpos_arr,[xpos]]
   endfor

   if nclick eq 1 then return, xpos_arr[0]
   return, xpos_arr
END