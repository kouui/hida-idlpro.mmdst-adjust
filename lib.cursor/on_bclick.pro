
; history
; 2022.02.13    k.u.    return cursor information

;------------------------------------------------------------------------
;; RETURN CURSOR CLICK COORDINATE
FUNCTION on_bclick, wid=wid, msg=msg, data=data, change=change, button=button
    if keyword_set(wid) then wset, wid
    if keyword_set(msg) then print, msg
    
    if keyword_set(data) then bcursor,xpos,ypos,/data, button=button, change=change $
    else bcursor,xpos,ypos,/device, button=button, change=change
    
    return, {click_struct, x:xpos, y:ypos, button:button}
END