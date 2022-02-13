
; history
; 2022.02.13    k.u.    created

FUNCTION MINVAL, v1, v2

    if v1 lt v2 then return, v1
    return, v2
END