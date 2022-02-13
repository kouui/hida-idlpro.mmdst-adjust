
; history
; 2022.02.13    k.u.    i

;------------------------------------------------------------------------
; SELECT XPOS AND DRAW LINE ONE IQUV

FUNCTION SELECTX_IQUV, s4, pmax=pmax, bin=bin, wid=wid, $
            labels=labels, ialog=ialog, msg_arr=msg_arr, $
            change=change, button=button, isclose=isclose, $
            isdraw=isdraw

    showiquvr, s4, pmax=pmax, bin=bin, wid=wid, ialog=ialog, labels=labels, x0s=x0s, right_bounds=right_bounds, dd=dd, wexist=change

    ss = size(s4)
    nx0 = ss[1]
    ny0 = ss[2]
    ns0 = ss[3]
    nx = nx0 / bin
    ny = ny0 / bin

    if not keyword_set(msg_arr) then msg_arr = ['click left button to select a x posiiton ...']
    nclick = n_elements(msg_arr)

    wset, wid


    if keyword_set(change) then begin

        cs = on_bclick(/change, button=button)
        xpos = fix(cs.x)
        
        ;; find device xpos in single image
        bias = 0
        for j=0,ns0-1 do begin
            if right_bounds[j] gt xpos then break
            bias = bias + (nx+dd)
        endfor
        xpos = xpos - bias

        ;; draw vertical line
        if keyword_set(isdraw) then $
        for j=0, ns0-1 do $
        draw, (x0s[j]+xpos)*[1,1], [dd,dd+ny], color=150, thick=2, line=2

        ;; correct xpos with pixel unit
        xpos = xpos*bin
        print, 'select xpos=', xpos

        return, [xpos]

    endif else begin

        xpos_arr = []
        for kk=0,nclick-1 do begin
            cs = on_bclick(msg=msg_arr[kk],button=button)
            xpos = fix(cs.x)
            
            ;; find device xpos in single image
            bias = 0
            for j=0,ns0-1 do begin
                if right_bounds[j] gt xpos then break
                bias = bias + (nx+dd)
            endfor
            xpos = xpos - bias

            ;; draw vertical line
            if keyword_set(isdraw) then $
            for j=0, ns0-1 do $
            draw, (x0s[j]+xpos)*[1,1], [dd,dd+ny], color=150, thick=2, line=2

            ;; correct xpos with pixel unit
            xpos = xpos*bin
            print, 'select xpos=', xpos
            xpos_arr = [xpos_arr,[xpos]]
        endfor

        if keyword_set(isclose) then wdelete, wid
        return, xpos_arr

    endelse
END


;;file = '/nwork/kouui/data-lvl1/dstpol/20220110.giant-prominence/spec/cal/sav.for-mmdst-adjust'
;;file = file + '/step5.s.ar.sav'
;;restore, file


;ret = SELECTX_IQUV(s[*,*,0:3], /isclose)
;ret = SELECTX_IQUV(s[*,*,0:3], /isdraw, button=button)

;showiquvr, s[*,*,0:3], pmax=pmax, bin=bin, wid=wid
;while 1b do begin
;    ret = SELECTX_IQUV(s[*,*,0:3], bin=bin, pmax=pmax, wid=wid, /isdraw, /change, button=button)
;    print, ret, button
;    if button eq 'right' then break
;endwhile

;END