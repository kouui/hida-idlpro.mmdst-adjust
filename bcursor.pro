
; history
; 2022.02.13    k.u.    corsor + button type

PRO BCURSOR, x, y, change=change, down=down, nowait=nowait, up=up, wait=wait, data=data, device=device, normal=normal, button=button

    CURSOR, x, y, change=change, down=down, nowait=nowait, up=up, wait=wait, data=data, device=device, normal=normal

    case !mouse.button of
        0 : button='none' 
        1 : button='left'
        2 : button='middle'
        4 : button='right'
    endcase

END